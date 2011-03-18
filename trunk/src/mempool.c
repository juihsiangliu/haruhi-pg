#include "mempool.h"


static MempoolSet * _pool = NULL;
static MempoolSet ** _poolList = NULL;
static int _poolListSize = 0;


static const int _2M = 2*1024*1024;
static const int _4M = 4*1024*1024;
static const int _8M = 8*1024*1024;

/*
static int int_min(const int a,const int b)
{
	if(a<b) return a;
	else return b;
}


static int int_max(const int a,const int b)
{
	if(a<b) return b;
	else return a;
}
*/

//==================================================




static void dumpMempool(FILE *fp,const Mempool *pool)
{
	fprintf(fp,"===============\n");
	fprintf(fp,"pool name: %s\n",pool->poolName);
	fprintf(fp,"batchNumOfBlock: %d\n",pool->batchNumOfBlock);
	fprintf(fp,"blockSize: %d\n",pool->blockSize);

	fprintf(fp,"*** freeList ***\n");
	dumpDqueue(fp,pool->dqueue);
}



static Mempool *createMempool(const char *poolName,const int blockSize,const int numOfBlock)
{
	Mempool *pool = (Mempool *)malloc(sizeof(Mempool));
	strcpy(pool->poolName,poolName);
	pool->blockSize = blockSize;

//	const int maxHeapSize = 2*int_max(numOfBlock,1);
	const int maxHeapSize = 2*numOfBlock;

	pool->dqueue = createDqueue(maxHeapSize);
//	pool->batchNumOfBlock = int_max(numOfBlock,1);
	pool->batchNumOfBlock = numOfBlock;

	pool->count = 0;
	pool->malloc_count = 0;
	pthread_mutex_init(&(pool->mutex), NULL);


	return pool;
}




static void freeMempool(Mempool *pool)
{	
	Dqueue *dqueue = pool->dqueue;
	while(!isEmptyDqueue(dqueue))
	{
		void *mem = delTailDqueue(dqueue);
		free(mem);
		pool->malloc_count--;
	}
	freeDqueue(dqueue);
	pthread_mutex_destroy(&(pool->mutex));
	free(pool);
}




static void clearMempool(Mempool *pool)
{
	pthread_mutex_lock(&(pool->mutex));
	Dqueue *dqueue = pool->dqueue;
	while(!isEmptyDqueue(dqueue))
	{
		void *mem = delTailDqueue(dqueue);
		free(mem);
		pool->malloc_count--;
	}
	pthread_mutex_unlock(&(pool->mutex));
}




// return a mem block from pool to user 
static void *getMempool(Mempool *pool)
{
	pthread_mutex_lock(&(pool->mutex));

	if(isEmptyDqueue(pool->dqueue))
	{
		int i;
		for(i=0;i<pool->batchNumOfBlock;i++)
		{
			if(!isFullDqueue(pool->dqueue))
			{
				void *ptr = malloc(pool->blockSize);
				pool->malloc_count++;
				if(ptr == NULL)
				{
					fprintf(stderr,"error: not enought memory\n");
					exit(0);
				}
				insertTailDqueue(pool->dqueue,ptr);
			}
			else break;
		}
	}
	void *ret = delTailDqueue(pool->dqueue);
	pool->count++;

	pthread_mutex_unlock(&(pool->mutex));

	return ret;
}



static void retMempool(Mempool *pool, void *data)
{
	pthread_mutex_lock(&pool->mutex);
	if(isFullDqueue(pool->dqueue))
	{
		void *head = delHeadDqueue(pool->dqueue);
		free(head);
		pool->malloc_count--;
	}
	insertTailDqueue(pool->dqueue,data);
	pool->count--;
	pthread_mutex_unlock(&(pool->mutex));
}


// ======================================================



static MempoolSet *createMempoolSetKernelNew(const int minBlockSize,const int numOfMempool,const int *const numOfBlockList)
{
	MempoolSet *set = (MempoolSet *)malloc(sizeof(MempoolSet));
	set->poolSet = (Mempool **)malloc(sizeof(Mempool *)*numOfMempool);
	int i;
	set->minBlockSize = minBlockSize;
	set->maxBlockSize = minBlockSize;
//	set->hugeBlockSize = hugeBlockSize;
	set->numOfPool = numOfMempool;
	set->lgMinBlockSize = log2(set->minBlockSize);

	for(i=0;i<numOfMempool;i++)
	{
		char name[32];
		int n;
		n = sprintf(name,"normal:%d",set->maxBlockSize);
		set->poolSet[i] = createMempool(name,set->maxBlockSize,numOfBlockList[i]);
		set->maxBlockSize = set->maxBlockSize*2;
	}
	set->maxBlockSize = set->maxBlockSize/2;

	set->huge = createMempool("huge",-1,0);


	return set;
}





static void freeMempoolSetKernel(MempoolSet *ptr)
{
	int i;
	for(i=0;i<ptr->numOfPool;i++) freeMempool(ptr->poolSet[i]);
	free(ptr->poolSet);
	freeMempool(ptr->huge);
	free(ptr);
}



static void clearMempoolSetKernel(MempoolSet *ptr)
{
	int i;
	for(i=0;i<ptr->numOfPool;i++) clearMempool(ptr->poolSet[i]);
//	clearMempool(ptr->huge);
}





static void *getMempoolSetKernel(MempoolSet *ptr,const int blockSize)
{
	void *ret = NULL;
	if(blockSize > ptr->maxBlockSize)
	{
		pthread_mutex_lock(&ptr->huge->mutex);
		ret = malloc(blockSize);
		ptr->huge->count++;
		pthread_mutex_unlock(&ptr->huge->mutex);
		if(ret == NULL)
		{
			fprintf(stderr,"not enough memory ... alloc huge fail\n");
			exit(0);
		}
	}
	else
	{
		const double lgBlockSize = log2(blockSize);
		int index = ceil(lgBlockSize - ptr->lgMinBlockSize);
		if(index < 0) index = 0;
		ret = getMempool(ptr->poolSet[index]);
	}

//	memset(ret,0,sizeof(blockSize));
	return ret;
}





static void retMempoolSetKernel(MempoolSet *ptr,void *data, const int blockSize)
{
	if(blockSize > ptr->maxBlockSize)
	{
		pthread_mutex_lock(&ptr->huge->mutex);
		free(data);
		ptr->huge->count--;
		pthread_mutex_unlock(&ptr->huge->mutex);
	}
	else
	{
		const double lgBlockSize = log2(blockSize);
		int index = ceil(lgBlockSize - ptr->lgMinBlockSize);
		if(index < 0) index = 0;
		retMempool(ptr->poolSet[index],data);
	}
}




static void dumpMempoolSetKernel(FILE *fp, const MempoolSet *set)
{
	int i;
	for(i=0;i<set->numOfPool;i++) dumpMempool(fp,set->poolSet[i]);
	dumpMempool(fp,set->huge);
}



//=============================================================





void createMempoolSet(const int minBlockSize,const int numOfMempool,const int *const numOfBlockList,const int poolListSize)
{
	int i;
	_poolList = (MempoolSet **)malloc(sizeof(MempoolSet *)*poolListSize);
	for(i=0;i<poolListSize;i++)
	{
		_poolList[i] = createMempoolSetKernelNew(minBlockSize,numOfMempool,numOfBlockList);
	}
	_pool = _poolList[0];
	_poolListSize = poolListSize;
}




void freeMempoolSet(void)
{
	int i;
	for(i=0;i<_poolListSize;i++)
	{
		freeMempoolSetKernel(_poolList[i]);
	}
	free(_poolList);
}



void clearPidMempoolSet(const int pid)
{
	clearMempoolSetKernel(_poolList[pid%_poolListSize]);
}




void *getMempoolSet(const int blockSize)
{
	void *ptr;
	ptr =  getMempoolSetKernel(_pool,blockSize);
	return ptr;
}




void *getPidMempoolSet(const int blockSize,const int pid)
{
	void *ptr;
	ptr =  getMempoolSetKernel(_poolList[pid%_poolListSize],blockSize);
	return ptr;
}



void retMempoolSet(void *data, const int blockSize)
{
	retMempoolSetKernel(_pool,data,blockSize);
}




void retPidMempoolSet(void *data, const int blockSize,const int pid)
{
	retMempoolSetKernel(_poolList[pid%_poolListSize],data,blockSize);
}



void dumpMempoolSet(FILE *fp)
{
	dumpMempoolSetKernel(fp,_pool);
}


void usageMempoolSet(FILE *fp)
{
	int j;

//	fprintf(stderr,"pool list size:%d\n",_poolListSize);
	for(j=0;j<_poolListSize;j++)
	{
		int i;
		fprintf(fp,"poolListIndex: %d\n",j);
		for(i=0;i<_poolList[j]->numOfPool;i++)
		{
//			fprintf(fp,"*\tmalloc count:%2d\n",_poolList[j]->poolSet[i]->malloc_count);
			fprintf(fp,"\tpool:%2d count:%2d queueSize:%d\n",i,_poolList[j]->poolSet[i]->count,_poolList[j]->poolSet[i]->dqueue->size);
		}
		fprintf(fp,"\tpool:hg count:%2d\n",_poolList[j]->huge->count);
		fprintf(fp,"\n");
/*
		fprintf(fp,"\n================================\n");
		fprintf(fp,"memory usage stats\n");
		const int MB = 1024*1024;
		int i;
		int total = 0;
		for(i=0;i<_poolList[j]->numOfPool;i++)
		{
			const int val = _poolList[j]->poolSet[i]->batchNumOfBlock * ((float)_pool->poolSet[i]->blockSize/MB);
			total += val;
			fprintf(fp,"memPool %d byte: %dMB\n",_poolList[j]->poolSet[i]->blockSize,val);
		}

		const int hugeVal = _poolList[j]->huge->batchNumOfBlock * (_poolList[j]->huge->blockSize/MB);
		total += hugeVal;
		fprintf(fp,"memPool huge %d byte: %dMB\n",_poolList[j]->hugeBlockSize,hugeVal);
		fprintf(fp,"total memPool usage: %dMB\n",total);
		fprintf(fp,"================================\n");
*/
	}
}





