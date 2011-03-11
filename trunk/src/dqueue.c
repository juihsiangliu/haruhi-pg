#include "dqueue.h"


static int inc(Dqueue *ptr,int x)
{
	if(x == ptr->maxSize-1)
	{
		return 0;
	}
	else
	{
		return x+1;	
	}
}



static int dec(Dqueue *ptr,int x)
{
	if(x == 0)
	{
		return ptr->maxSize-1;
	}
	else
	{
		return x-1;
	}
}





Dqueue *createDqueue(const int maxSize)
{
	Dqueue *ptr = (Dqueue *)malloc(sizeof(Dqueue));
	ptr->maxSize = maxSize;
	ptr->head = 0;
	ptr->tail = 0;
	ptr->queue = malloc(sizeof(void *)*(maxSize));
	ptr->size = 0;
	memset(ptr->queue,0,sizeof(void *)*(maxSize));
	return ptr;
}



void freeDqueue(Dqueue *ptr)
{
	free(ptr->queue);
	free(ptr);
}




int isFullDqueue(Dqueue *ptr)
{
//	if(ptr->size == ptr->maxSize-1) return 1;
	if(ptr->size == ptr->maxSize) return 1;
	else return 0;
}




int isEmptyDqueue(Dqueue *ptr)
{
	if(ptr->size == 0) return 1;
	else return 0;
}




void dumpDqueue(FILE *fp,Dqueue *ptr)
{
	int current = 0;
	int index = ptr->head+1;
	fprintf(fp,"current size / max size = %d / %d\n",ptr->size,ptr->maxSize);
	fprintf(fp,"head:%d ,tail:%d\n",ptr->head,ptr->tail);
	while(current < ptr->size)
	{
		fprintf(fp,"%d, address:%p\n",current,ptr->queue[index%ptr->maxSize]);
		index++;
		current++;
	}
}





void insertTailDqueue(Dqueue *ptr, void *data)
{
	if(isFullDqueue(ptr))
	{
//		fprintf(stderr,"the queue is full\n");
//		exit(0);
		return;
	}
	else
	{
		int newTail = inc(ptr,ptr->tail);
		ptr->tail = newTail;
		ptr->queue[newTail] = data;
		ptr->size++;
	}
}




void *delHeadDqueue(Dqueue *ptr)
{
	if(isEmptyDqueue(ptr))
	{
//		fprintf(stderr,"the queue is empty\n");
//		exit(0);
		return NULL;
	}
	else
	{
		int newHead = inc(ptr,ptr->head);
		ptr->head = newHead;
		ptr->size--;
		void *ret = ptr->queue[newHead];
		ptr->queue[newHead] = NULL;
		return ret;
	}
}




void *delTailDqueue(Dqueue *ptr)
{
	if(isEmptyDqueue(ptr))
	{
//		fprintf(stderr,"the queue is empty\n");
//		exit(0);
		return NULL;
	}
	else
	{
		int newTail = dec(ptr,ptr->tail);
		void *ret = ptr->queue[ptr->tail];
		ptr->queue[ptr->tail] = NULL;
		ptr->tail = newTail;
		ptr->size--;
		return ret;
	}
}



/*

int main(void)
{
	Dqueue *ptr = createDqueue(4);
	int i = 0;
	dumpDqueue(stdout,ptr);
	for(i=0;i<4;i++)
	{
		int *mem = malloc(sizeof(16));
		insertTailDqueue(ptr,mem);
		dumpDqueue(stdout,ptr);
	}
	for(i=0;i<2;i++)
	{
		int *mem = delHeadDqueue(ptr);
		printf("obtained:%p\n",mem);
		dumpDqueue(stdout,ptr);
	}

	for(i=0;i<3;i++)
	{
		int *mem = malloc(sizeof(16));
		insertTailDqueue(ptr,mem);
		dumpDqueue(stdout,ptr);
	}
}


*/
