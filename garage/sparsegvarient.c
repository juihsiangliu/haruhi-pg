#include "sparsegvarient.h"


/*
struct Queue_pid
{
	gdsl_queue_t queue;
	int pid;
};

typedef struct Queue_pid Queue_pid;
*/



static gdsl_element_t allocGControlInfoQueue(void *ptr)
{
	const GControlInfo *src = (GControlInfo *) ptr;
	GControlInfo *dest = createGControlInfo(src->vdsListSize,src->vgsListSize,src->gvNum);
	dest->sign = src->sign;
	dest->type = src->type;
	dest->gate = src->gate;
	dest->drain = src->drain;
	dest->source = src->source;
	dest->pid = src->pid;

	dest->vdsList = src->vdsList;
	dest->vgsList = src->vgsList;
	dest->partialIdsVxs = src->partialIdsVxs;
	return (gdsl_element_t) dest;
}




/*
static gdsl_element_t allocGControlInfoQueuePid(void *ptr)
{
	const GControlInfo *src = (GControlInfo *) ptr;
	const int pid = src->pid;
	GControlInfo *dest = createGControlInfoPid(src->vdsListSize,src->vgsListSize,src->gvNum,pid);
	dest->sign = src->sign;
	dest->type = src->type;
	dest->gate = src->gate;
	dest->drain = src->drain;
	dest->source = src->source;
	dest->pid = src->pid;

	dest->vdsList = src->vdsList;
	dest->vgsList = src->vgsList;
	dest->partialIdsVxs = src->partialIdsVxs;
	return (gdsl_element_t) dest;
}
*/





static void freeGControlInfoQueue(void *ptr)
{
	freeGControlInfo(ptr);
}



/*
static void freeGControlInfoQueuePid(void *ptr)
{
	GControlInfo *src = (GControlInfo *)ptr;
	const int pid = src->pid;
	freeGControlInfoPid(src,pid);
}
*/



static int dumpGControlInfoQueue(const gdsl_element_t E,gdsl_location_t LOCATION,void *USER_DATA)
{
	dumpGControlInfo(USER_DATA,E);
	return 0;
}


// E - one element in src queue
// user_data - the des_queue
static int copyGControlInfoQueue(const gdsl_element_t E, gdsl_location_t LOCATION, void *USER_DATA)
{
	gdsl_queue_t dest_queue = (gdsl_queue_t) USER_DATA;	
	gdsl_queue_insert(dest_queue,E);
}


/*
// E - one element in src queue
// user_data - queue_pid: des_queue + pid
static int copyPidGControlInfoQueue(const gdsl_element_t E, gdsl_location_t LOCATION, void *USER_DATA)
{
	Queue_pid *queue_pid = USER_DATA;
	gdsl_queue_t dest_queue = queue_pid->queue;
	GControlInfo *Ecopy = E;
	const int pid_backup = Ecopy->pid;
	Ecopy->pid = queue_pid->pid;
	gdsl_queue_insert(dest_queue,Ecopy);
	Ecopy->pid = pid_backup;
}
*/


//=======================================================


GControlInfo *createGControlInfo(const int vdsListSize, const int vgsListSize,const int gvNum)
{
	return createGControlInfoPid(vdsListSize,vgsListSize,gvNum,0);
}




GControlInfo *createGControlInfoPid(const int vdsListSize, const int vgsListSize,const int gvNum,const int pid)
{
	GControlInfo *dest = getPidMempoolSet(sizeof(GControlInfo),pid);
	dest->vdsListSize = vdsListSize;
	dest->vgsListSize = vgsListSize;
	dest->gvNum = gvNum;
	dest->sign = 0;
	dest->type = gm;
	dest->gate = -1;
	dest->drain = -1;
	dest->source = -1;
	dest->pid = pid;

	return dest;
}




void freeGControlInfo(GControlInfo *ptr)
{
	freeGControlInfoPid(ptr,0);
}






void freeGControlInfoPid(GControlInfo *ptr,const int pid)
{
	if(ptr!=NULL)
	{
		retPidMempoolSet(ptr,sizeof(GControlInfo),pid);
	}
}



void dumpGControlInfo(FILE *fp,const GControlInfo *ptr)
{
	fprintf(fp,"gate:%d,drain:%d,source:%d\n",ptr->gate,ptr->drain,ptr->source);
	fprintf(fp,"sign=%d\n",ptr->sign);
	if(ptr->type == ro_1) fprintf(fp,"type=ro_1\n");
	else fprintf(fp,"type=gm\n");
	int i;
	for(i=0;i<ptr->vdsListSize;i++) fprintf(fp,"%01.2f ",ptr->vdsList[i]);
	fprintf(fp,"\n");
	for(i=0;i<ptr->vgsListSize;i++) fprintf(fp,"%01.2f ",ptr->vgsList[i]);
	fprintf(fp,"\n");
}




/*
GControlInfo *copyGControlInfo(const GControlInfo *const src)
{
	GControlInfo *dest = getMempoolSet(sizeof(GControlInfo));
	memcpy(dest,src,sizeof(GControlInfo));
	return dest;
}
*/



SparseGControlInfo *createSparseGControlInfo(void)
{
	SparseGControlInfo *ptr = getMempoolSet(sizeof(SparseGControlInfo));
	ptr->row = -1;
	ptr->col = -1;
	ptr->queue = gdsl_queue_alloc("queue",allocGControlInfoQueue,freeGControlInfoQueue);
	ptr->rowLink = NULL;
	ptr->colLink = NULL;
	return ptr;
}





void freeSparseGControlInfo(SparseGControlInfo *ptr)
{
	gdsl_queue_free(ptr->queue);
	retMempoolSet(ptr,sizeof(SparseGControlInfo));
}




void dumpSparseGControlInfo(FILE *fp,const SparseGControlInfo *ptr)
{
	fprintf(fp,"row:%d,col:%d\n",ptr->row,ptr->col);
	gdsl_queue_map_forward(ptr->queue,dumpGControlInfoQueue,fp);
}




SparseGVarientTable *createSparseGVarientTable(const int row,const int col)
{
	SparseGVarientTable *ptr = (SparseGVarientTable *)getMempoolSet(sizeof(SparseGVarientTable));
	ptr->totalRow = row;
	ptr->totalCol = col;

	ptr->rowIndex = (SparseGControlInfo **)getMempoolSet(row*sizeof(SparseGControlInfo *));
	ptr->colIndex = (SparseGControlInfo **)getMempoolSet(col*sizeof(SparseGControlInfo *));

	int i;
	for(i=0;i<row;i++) ptr->rowIndex[i] = createSparseGControlInfo();
	for(i=0;i<col;i++) ptr->colIndex[i] = createSparseGControlInfo();
	return ptr;
}



void freeSparseGVarientTable(SparseGVarientTable *ptr)
{
	// free the memory row by row
	int i;
	SparseGControlInfo *currentEntry;
	SparseGControlInfo *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freeSparseGControlInfo(removeEntry);
		}
	}

	for(i=0;i<ptr->totalRow;i++) freeSparseGControlInfo(ptr->rowIndex[i]);
	for(i=0;i<ptr->totalCol;i++) freeSparseGControlInfo(ptr->colIndex[i]);
	
	retMempoolSet(ptr->rowIndex,ptr->totalRow*sizeof(SparseGControlInfo *));
	retMempoolSet(ptr->colIndex,ptr->totalCol*sizeof(SparseGControlInfo *));
	retMempoolSet(ptr,sizeof(SparseGVarientTable));
}




void insertSparseGVarientTable(SparseGVarientTable *table,GControlInfo *element,const int row,const int col)
{
	SparseGControlInfo *rowTarget = table->rowIndex[row];
	SparseGControlInfo *colTarget = table->colIndex[col];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= col) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= row) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink==NULL || rowTarget->rowLink->col!=col)
	{
		SparseGControlInfo *insert = createSparseGControlInfo();
		insert->row = row;
		insert->col = col;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;
		// insert element into insert->queue
		gdsl_queue_insert(insert->queue,element);
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else
	{
		gdsl_queue_insert(rowTarget->rowLink->queue,element);
	}
}




const SparseGControlInfo* getSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col)
{
	SparseGControlInfo *target = table->rowIndex[row]->rowLink;
	while(target != NULL)
	{
		if(target->col >= col) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=col) return NULL;
	else return target;
}




const gdsl_queue_t getQueueSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col)
{
	const SparseGControlInfo *target = getSparseGVarientTable(table,row,col);
	if(target == NULL) return NULL;
	return target->queue;
}



static pthread_mutex_t mutexCopy = PTHREAD_MUTEX_INITIALIZER;


gdsl_queue_t getCopyQueueSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col)
{
	const SparseGControlInfo *target = getSparseGVarientTable(table,row,col);
	if(target == NULL) return NULL;
	else
	{
		pthread_mutex_lock(&mutexCopy);
		gdsl_queue_t src_queue = target->queue;
		gdsl_queue_t des_queue =  gdsl_queue_alloc("queue",allocGControlInfoQueue,freeGControlInfoQueue);
		gdsl_queue_map_forward(src_queue,copyGControlInfoQueue,des_queue);
		pthread_mutex_unlock(&mutexCopy);
		return des_queue;
	}
}






/*
gdsl_queue_t getCopyPidQueueSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col,const int pid)
{
	const SparseGControlInfo *target = getSparseGVarientTable(table,row,col);
	if(target == NULL) return NULL;
	else
	{
//		pthread_mutex_lock(&mutexCopy);
		gdsl_queue_t src_queue = target->queue;
		gdsl_queue_t des_queue =  gdsl_queue_alloc("queue",allocGControlInfoQueuePid,freeGControlInfoQueuePid);
		Queue_pid *queue_pid = getPidMempoolSet(sizeof(Queue_pid),pid);
		queue_pid->queue = des_queue;
		queue_pid->pid = pid;
		gdsl_queue_map_forward(src_queue,copyPidGControlInfoQueue,queue_pid);
		retPidMempoolSet(queue_pid,sizeof(Queue_pid),pid);
//		pthread_mutex_unlock(&mutexCopy);
		return des_queue;
	}
}
*/





void freeCopyQueueSparseGVarientTable(gdsl_queue_t queue)
{
	pthread_mutex_lock(&mutexCopy);
	if(queue!=NULL) gdsl_queue_free(queue);
	pthread_mutex_unlock(&mutexCopy);
}





void dumpSparseGVarientTable(FILE *fp,const SparseGVarientTable *table)
{
	int i,j;
	for(i=0;i<table->totalRow;i++)
	{
		for(j=0;j<table->totalCol;j++)
		{
			const SparseGControlInfo *target = getSparseGVarientTable(table,i,j);
			if(target!=NULL) dumpSparseGControlInfo(fp,target);
		}
	}
}





