#include "parallel_common.h"


static pthread_mutex_t freeSinkList_mutex = PTHREAD_MUTEX_INITIALIZER;

static gdsl_element_t alloc_int(void *ptr)
{

	int *n = (int *) ptr;
	int *value = getMempoolSet(sizeof(int));

	// copy from n to value
	memcpy(value,n,sizeof(int));
	

	return (gdsl_element_t) value;
}


static void free_int(gdsl_element_t e)
{
	retMempoolSet(e,sizeof(int));
}


static int gdsl_queue_get_size_mutex(gdsl_queue_t queue)
{
	int n;
	n = gdsl_queue_get_size(queue);
	return n;
}



static void getInitSinkList(gdsl_queue_t freeSinkList, const ParallelETree *tree)
{
	int i;
	for(i=0;i<tree->size;i++)
//	for(i=tree->size-1;i>-1;i--)
	{
		if(tree->node[i]!=NULL)
		{
			if(tree->node[i]->type == lu)
			{
				gdsl_queue_insert(freeSinkList,&i);
			}
		}
	}
}




static void updateFreeSinkList(gdsl_queue_t freeSinkList, ParallelETree *tree, const int currentNodeIndex,const int treeInternalSize,const int *priorityTable)
{
	int rootIndex = currentNodeIndex/2;
	const int leftIndex = rootIndex*2;
	const int rightIndex = rootIndex*2 + 1;
	int *mark = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	memset(mark,0,sizeof(int)*(treeInternalSize+1));

	tree->node[currentNodeIndex]->visitLog = visit;

	if(tree->node[rootIndex]==NULL)
	{
		retMempoolSet(mark,sizeof(int)*(treeInternalSize+1));
		return;
	}


	// reorder the queue .... suck implement
	if( tree->node[leftIndex]->visitLog==visit && tree->node[rightIndex]->visitLog==visit )
	{
		int i,j;
		int nnz = 0;
		while(gdsl_queue_get_size(freeSinkList)!=0)
		{
			int *indexInQueue = gdsl_queue_get_head(freeSinkList);
			gdsl_queue_remove(freeSinkList);
			mark[*indexInQueue] = 1;
			nnz++;
			retMempoolSet(indexInQueue,sizeof(int));
		}
		mark[rootIndex] = 1;
		nnz++;
		for(i=0;i<nnz;i++)
		{
			int maxIndex = -1;
			int maxPri = -1;
			for(j=1;j<treeInternalSize+1;j++)
			{
				if( (mark[j]==1) && (priorityTable[j] > maxPri))
				{
					maxPri = priorityTable[j];
					maxIndex = j;
				}
			}
			mark[maxIndex] = 0;
			gdsl_queue_insert(freeSinkList,&maxIndex);
		}
	}


	if( tree->node[leftIndex]->visitLog==visit && tree->node[rightIndex]->visitLog==visit )
	{
		int temp1,temp2;
		temp1 = GSL_MIN(tree->node[leftIndex]->doneRowBegin,tree->node[rightIndex]->doneRowBegin);
		temp2 = GSL_MIN(tree->node[leftIndex]->rowBegin,tree->node[rightIndex]->rowBegin);
//
		tree->node[rootIndex]->doneRowBegin = GSL_MIN(temp1,temp2);
		temp1 = GSL_MAX(tree->node[leftIndex]->doneRowEnd,tree->node[rightIndex]->doneRowEnd);
		temp2 = GSL_MAX(tree->node[leftIndex]->rowEnd,tree->node[rightIndex]->rowEnd);
		tree->node[rootIndex]->doneRowEnd = GSL_MAX(temp1,temp2);

	}
	
	retMempoolSet(mark,sizeof(int)*(treeInternalSize+1));
}


/*
void getDoneRowInfo(ParallelETree *tree,const gdsl_queue_t freeSinkListSrc)
{
	int i;
	gdsl_queue_t freeSinkList = gdsl_queue_alloc("freeSinkListTemp",alloc_int,free_int);
	getInitSinkList(freeSinkList,tree);
	// set the doneRow information
	while(gdsl_queue_get_size(freeSinkList)!=0)
	{
		int *indexInQueue = gdsl_queue_get_head(freeSinkList);
		gdsl_queue_remove(freeSinkList);
		updateFreeSinkList(freeSinkList,tree,*indexInQueue);
		retMempoolSet(indexInQueue,sizeof(int));
	}

	// recover the setting of visiting the tree
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i] == NULL) continue;
		else tree->node[i]->visitLog = notvisit;
	}

	gdsl_queue_free(freeSinkList);
}
*/




void getDoneRowInfoNew(ParallelETree *tree,const int *postorder,const int *postorder_inv,const int treeInternalSize)
{	
	int i;

	int *priorityTable = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=0;i<treeInternalSize+1;i++) priorityTable[i] = treeInternalSize - postorder_inv[i];

	gdsl_queue_t freeSinkList = gdsl_queue_alloc("freeSinkListTemp",alloc_int,free_int);
	getInitSinkList(freeSinkList,tree);
	// set the doneRow information
	while(gdsl_queue_get_size(freeSinkList)!=0)
	{
		int *indexInQueue = gdsl_queue_get_head(freeSinkList);
		gdsl_queue_remove(freeSinkList);
		updateFreeSinkList(freeSinkList,tree,*indexInQueue,treeInternalSize,priorityTable);
		retMempoolSet(indexInQueue,sizeof(int));
	}

	// recover the setting of visiting the tree
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i] == NULL) continue;
		else tree->node[i]->visitLog = notvisit;
	}

	gdsl_queue_free(freeSinkList);
	retMempoolSet(priorityTable,sizeof(int)*(treeInternalSize+1));
}


// =======================================================



static pthread_mutex_t todolist_mutex = PTHREAD_MUTEX_INITIALIZER;
ToDoList *createToDoList(const int size,const int *postorder)
{
	ToDoList *ptr = getMempoolSet(sizeof(ToDoList));
	ptr->maxSize = size;
	ptr->currentSize = size;

	// the init list
	ptr->list = getMempoolSet(sizeof(int)*size);
/*
	int i;
	int j = 0;
	for(i= size/2 + 1;i<=size;i++) 
	{
		ptr->list[j] = i;
		j++;
	}
	for(i=1;i<= size/2;i++)
	{
		ptr->list[j] = i;
		j++;
	}
*/
	memcpy(ptr->list,postorder,sizeof(int)*size);
//	for(i=0;i<15;i++) printf("%d\n",ptr->list[i]);
	return ptr;
}



void freeToDoList(ToDoList *ptr)
{
	retMempoolSet(ptr->list,sizeof(int)*ptr->maxSize);
	retMempoolSet(ptr,sizeof(ToDoList));
}




// NOT thread safe
// N is zero based
static void pushNthToDoList(ToDoList *ptr,const int N, const int data)
{
	if(ptr->currentSize == ptr->maxSize)
	{
		printf("error!! todolist is full\n");
	}
	else
	{
		if(N == 0) // into first
		{
			memmove(&ptr->list[1],&ptr->list[0],sizeof(int)*(ptr->currentSize));
			ptr->list[0] = data;
			ptr->currentSize++;
		}
		else if(N == ptr->currentSize) // into end
		{
			ptr->list[ptr->currentSize] = data;
			ptr->currentSize++;
		}
		else
		{
			memmove(&ptr->list[N+1],&ptr->list[N],sizeof(int)*(ptr->currentSize-N));
			ptr->list[N] = data;
			ptr->currentSize++;
		}
	}
}





void pushBackToDoList(ToDoList *ptr, int data)
{
	pthread_mutex_lock(&todolist_mutex);
	if(ptr->currentSize == ptr->maxSize)
	{
		printf("error!! todolist is full\n");
	}
	else
	{
		ptr->list[ptr->currentSize] = data;
		ptr->currentSize++;
	}
	pthread_mutex_unlock(&todolist_mutex);
}



/*
static void pushByOrderToDoList(ToDoList *ptr, int data)
{
	pthread_mutex_lock(&todolist_mutex);
	if(ptr->currentSize == 15)
	{
		printf("error!! todolist is full\n");
	}
	else
	{
		if(ptr->currentSize == 0)
		{
			ptr->list[ptr->currentSize] = data;
			ptr->currentSize++;
		}
		else
		{
			int i = 0;
			int dataOrder = nodeOrder[data];
			for(i=ptr->currentSize-1;i>-1;i--)
			{
				int order =  nodeOrder[ptr->list[i]];
				if(dataOrder > order) break;
			}
			pushNthToDoList(ptr,i+1,data);
		}
	}
	pthread_mutex_unlock(&todolist_mutex);
}
*/





// NOT thread safe
static void pushFirstToDoList(ToDoList *ptr, int data)
{
	if(ptr->currentSize == ptr->maxSize)
	{
		printf("error!! todolist is full\n");
	}
	else
	{
		memmove(&ptr->list[1],&ptr->list[0],sizeof(int)*(ptr->currentSize));
		ptr->list[0] = data;
		ptr->currentSize++;
	}
}





// NOT thread safe
// N is zero based
static int getNthToDoList(ToDoList *ptr,const int N)
{
	int ret = -1;
	if(N == 0) // the fist one
	{
		if(ptr->currentSize == 0)
		{
	//		printf("error!! todolist is empty\n");
		}
		else
		{
			ret = ptr->list[0];
			ptr->currentSize--;
			memmove(&ptr->list[0],&ptr->list[1],sizeof(int)*ptr->currentSize);
		}	
	}
	else if(N == ptr->currentSize-1) // the last one
	{
		ret = ptr->list[N];
		ptr->currentSize--;
	}
	else
	{
		if(ptr->currentSize == 0)
		{
			printf("error get Nth!! todolist is empty\n");
		}
		else
		{
			ret = ptr->list[N];
			memmove(&ptr->list[N],&ptr->list[N+1],sizeof(int)*(ptr->currentSize-N+1));
			ptr->currentSize--;
		}
	}
	return ret;
}




int getFirstToDoList(ToDoList *ptr)
{
	pthread_mutex_lock(&todolist_mutex);
	int ret = getNthToDoList(ptr,0);
	pthread_mutex_unlock(&todolist_mutex);
	return ret;
}



/*
static int compareInt (const void * a, const void * b)
{
	int indexA = *(int *)a;
	int indexB = *(int *)b;
	return nodeOrder[indexA] - nodeOrder[indexB];
}




static void sortToDoList(ToDoList *ptr)
{
	pthread_mutex_lock(&todolist_mutex);
	if(ptr->currentSize > 0)
	{
		qsort(ptr->list,ptr->currentSize,sizeof(int),compareInt);
	}
	pthread_mutex_unlock(&todolist_mutex);
}
*/



// sort all parents to the head of the todolist
void sortParentsToHeadToDoList(ToDoList *ptr,const int key,const int *rootCurrentBeginList)
{
	pthread_mutex_lock(&todolist_mutex);
	int parent = key/2;
	int i;
	int ins = 0;
	while(parent != 0)
	{
		int index = -1;
		// find tht index in todolist;
		for(i=0;i<ptr->currentSize;i++)
		{
			if(ptr->list[i] == parent)
			{
				index = i;
				break;
			}
		}
		if(index!=-1)
		{
			int val = getNthToDoList(ptr,index);
//			pushNthToDoList(ptr,ins,val);
			pushFirstToDoList(ptr,val);
			ins++;
		}
		parent = parent / 2;
	}
	pthread_mutex_unlock(&todolist_mutex);
}





void dumpToDoList(FILE *fp,ToDoList *ptr)
{
	int i;
	pthread_mutex_lock(&todolist_mutex);
	fprintf(fp,"ToDoList: ");
	for(i=0;i<ptr->currentSize;i++)
	{
		fprintf(fp,"%d, ",ptr->list[i]);
	}
	fprintf(fp,"\n");
	pthread_mutex_unlock(&todolist_mutex);
}



// =================================

int _activeThread = 0;

static pthread_mutex_t activeThread_mutex = PTHREAD_MUTEX_INITIALIZER;

void initActiveThread()
{
	pthread_mutex_lock(&activeThread_mutex);
	_activeThread = 0;
	pthread_mutex_unlock(&activeThread_mutex);
}


void incActiveThread()
{
	pthread_mutex_lock(&activeThread_mutex);
	_activeThread = _activeThread+1;
	pthread_mutex_unlock(&activeThread_mutex);
}


void decActiveThread()
{
	pthread_mutex_lock(&activeThread_mutex);
	_activeThread = _activeThread-1;
	pthread_mutex_unlock(&activeThread_mutex);
}



int getActiveThread(void)
{
	return _activeThread;
}





void safeWaitFlagMatrix(int *flag,int index)
{
	while(flag[index]!=1)
	{
		usleep(10000);
	}
	flag[index] = 0;
	usleep(10000);
	return;
}



