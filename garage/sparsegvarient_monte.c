#include "sparsegvarient_monte.h"


static gdsl_element_t allocGControlInfoMonteQueue(void *ptr)
{
	const GControlInfoMonte *src = (GControlInfoMonte *) ptr;
	GControlInfoMonte *dest = createGControlInfoMonte(src->vdsListSize,src->vgsListSize);
	dest->sign = src->sign;
	dest->gate = src->gate;
	dest->drain = src->drain;
	dest->source = src->source;

	dest->vdsList = src->vdsList;
	dest->vgsList = src->vgsList;
	dest->partialIdsVxs = src->partialIdsVxs;

	return (gdsl_element_t) dest;
}

static void freeGControlInfoMonteQueue(void *ptr)
{
	freeGControlInfoMonte(ptr);
}


static int dumpGControlInfoMonteQueue(const gdsl_element_t E,gdsl_location_t LPCATION,void *USER_DATA)
{
	dumpGControlInfoMonte(USER_DATA,E);
	return 0;
}


//=======================================================


GControlInfoMonte *createGControlInfoMonte(const int vdsListSize, const int vgsListSize)
{
	GControlInfoMonte *dest = getMempoolSet(sizeof(GControlInfoMonte));
	dest->vdsListSize = vdsListSize;
	dest->vgsListSize = vgsListSize;
	dest->sign = 0;
	dest->gate = -1;
	dest->drain = -1;
	dest->source = -1;
	return dest;
}


void freeGControlInfoMonte(GControlInfoMonte *ptr)
{
	if(ptr!=NULL)
	{
		retMempoolSet(ptr,sizeof(GControlInfoMonte));
	}
}



void dumpGControlInfoMonte(FILE *fp,const GControlInfoMonte *ptr)
{
	fprintf(fp,"gate:%d,drain:%d,source:%d\n",ptr->gate,ptr->drain,ptr->source);
	fprintf(fp,"sign=%d\n",ptr->sign);
	int i;
	for(i=0;i<ptr->vdsListSize;i++) fprintf(fp,"%01.2f ",ptr->vdsList[i]);
	fprintf(fp,"\n");
	for(i=0;i<ptr->vgsListSize;i++) fprintf(fp,"%01.2f ",ptr->vgsList[i]);
	fprintf(fp,"\n");
}





SparseGControlInfoMonte *createSparseGControlInfoMonte(void)
{
	SparseGControlInfoMonte *ptr = (SparseGControlInfoMonte *)malloc(sizeof(SparseGControlInfoMonte));
	ptr->row = -1;
	ptr->col = -1;
	ptr->queue = gdsl_queue_alloc("queue",allocGControlInfoMonteQueue,freeGControlInfoMonteQueue);
	ptr->rowLink = NULL;
	ptr->colLink = NULL;
	return ptr;
}


void freeSparseGControlInfoMonte(SparseGControlInfoMonte *ptr)
{
	gdsl_queue_free(ptr->queue);
	free(ptr);
}




void dumpSparseGControlInfoMonte(FILE *fp,const SparseGControlInfoMonte *ptr)
{
	fprintf(fp,"row:%d,col:%d\n",ptr->row,ptr->col);
	gdsl_queue_map_forward(ptr->queue,dumpGControlInfoMonteQueue,fp);
}




void setGControlInfoMonte(GControlInfoMonte *dest, const GControlInfo *src,const NonlinearMonteList *list)
{
	dest->vdsListSize = src->vdsListSize;
	dest->vgsListSize = src->vgsListSize;
	dest->sign = src->sign;
	dest->gate = src->gate;
	dest->drain = src->drain;
	dest->source = src->source;

	dest->vdsList = src->vdsList;
	dest->vgsList = src->vgsList;

	NonlinearInfoMonte *current = list->list;
	while(current!=NULL)
	{
		if(src->vdsList == current->vdsList)
			dest->vdsList = current->vdsList;
		if(src->vgsList == current->vgsList)
			dest->vgsList = current->vgsList;
		if(src->partialIdsVxs == current->quadPartialIdsVds)
			dest->partialIdsVxs = current->partialIdsVds;
		if(src->partialIdsVxs == current->quadPartialIdsVgs)
			dest->partialIdsVxs = current->partialIdsVgs;
		current = current->next;
	}
}



// ===============================================================


SparseGVarientTableMonte *createSparseGVarientTableMonte(const int row,const int col)
{
	SparseGVarientTableMonte *ptr = (SparseGVarientTableMonte *)malloc(sizeof(SparseGVarientTableMonte));
	ptr->totalRow = row;
	ptr->totalCol = col;

	ptr->rowIndex = (SparseGControlInfoMonte **)malloc(row*sizeof(SparseGControlInfoMonte *));
	ptr->colIndex = (SparseGControlInfoMonte **)malloc(col*sizeof(SparseGControlInfoMonte *));

	int i;
	for(i=0;i<row;i++) ptr->rowIndex[i] = createSparseGControlInfoMonte();
	for(i=0;i<col;i++) ptr->colIndex[i] = createSparseGControlInfoMonte();
	return ptr;
}



void freeSparseGVarientTableMonte(SparseGVarientTableMonte *ptr)
{
	// free the memory row by row
	int i;
	SparseGControlInfoMonte *currentEntry;
	SparseGControlInfoMonte *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freeSparseGControlInfoMonte(removeEntry);
		}
	}

	for(i=0;i<ptr->totalRow;i++) freeSparseGControlInfoMonte(ptr->rowIndex[i]);
	for(i=0;i<ptr->totalCol;i++) freeSparseGControlInfoMonte(ptr->colIndex[i]);
	
	free(ptr->rowIndex);
	free(ptr->colIndex);
	free(ptr);
}




void insertSparseGVarientTableMonte(SparseGVarientTableMonte *table,GControlInfoMonte *element,const int row,const int col)
{
	SparseGControlInfoMonte *rowTarget = table->rowIndex[row];
	SparseGControlInfoMonte *colTarget = table->colIndex[col];
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
		SparseGControlInfoMonte *insert = createSparseGControlInfoMonte();
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




const SparseGControlInfoMonte* getSparseGVarientTableMonte(const SparseGVarientTableMonte *table,const int row,const int col)
{
	SparseGControlInfoMonte *target = table->rowIndex[row]->rowLink;
	while(target != NULL)
	{
		if(target->col >= col) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=col) return NULL;
	else return target;
}




const gdsl_queue_t getQueueSparseGVarientTableMonte(const SparseGVarientTableMonte *table,const int row,const int col)
{
	const SparseGControlInfoMonte *target = getSparseGVarientTableMonte(table,row,col);
	if(target == NULL) return NULL;
	return target->queue;
}




void dumpSparseGVarientTableMonte(FILE *fp,const SparseGVarientTableMonte *table)
{
	int i,j;
	for(i=0;i<table->totalRow;i++)
	{
		for(j=0;j<table->totalCol;j++)
		{
			const SparseGControlInfoMonte *target = getSparseGVarientTableMonte(table,i,j);
			if(target!=NULL) dumpSparseGControlInfoMonte(fp,target);
		}
	}
}





// this function could only be called when the "nonlinear queue" is established
void setSparseGVarientTableMonte(SparseGVarientTableMonte *dest, const SparseGVarientTable *src,const NonlinearMonteList *list)
{
	int i,j;
	SparseGControlInfo *current;
	for(i=0;i<src->totalRow;i++)
	{
		current = src->rowIndex[i];
		while(current!=NULL)
		{
			const int queueSize = gdsl_queue_get_size(current->queue);
			for(j=0;j<queueSize;j++)
			{
				const GControlInfo *element = gdsl_queue_search_by_position(current->queue,j);
				GControlInfoMonte *monteElement = createGControlInfoMonte(element->vdsListSize,element->vgsListSize);
				setGControlInfoMonte(monteElement,element,list);
				insertSparseGVarientTableMonte(dest,monteElement,current->row,current->col);
			}
			current = current->rowLink;
		}
	}
}





