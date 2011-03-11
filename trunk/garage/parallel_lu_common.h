#ifndef PARALLEL_LU_COMMON_H
#define PARALLEL_LU_COMMON_H


#include <gdsl_types.h>
#include <gdsl_queue.h>
#include <gsl_math.h>
#include "mempool.h"
#include "partition_double.h"


gdsl_element_t alloc_int(void *ptr);
void free_int(gdsl_element_t e);
int gdsl_queue_get_size_mutex(gdsl_queue_t queue);
void getInitSinkList(gdsl_queue_t freeSinkList, const ParallelETree *tree);
void updateFreeSinkList(gdsl_queue_t freeSinkList, ParallelETree *tree, const int currentNodeIndex);
void getDoneRowInfo(ParallelETree *tree,const gdsl_queue_t freeSinkListSrc);



void safeWaitFlagMatrix(int *flag,int index);



void initActiveThread(void);
void incActiveThread(void);
void decActiveThread(void);
int getActiveThread(void);


struct ToDoList
{
	int currentSize;
	int list[15];
};

typedef struct ToDoList ToDoList;


ToDoList *createToDoList(void);
void freeToDoList(ToDoList *ptr);
int getFirstToDoList(ToDoList *ptr);
void dumpToDoList(FILE *fp,ToDoList *ptr);
void sortParentsToHeadToDoList(ToDoList *ptr,const int key);
void pushBackToDoList(ToDoList *ptr, int data);

#endif
