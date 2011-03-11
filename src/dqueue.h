#ifndef DQUEUE_H
#define DQUEUE_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Dqueue
{
	int maxSize;
	int size;
	int head;
	int tail;
	void **queue;
};

// the data is put in queue from head+1 to tail
// typically, the size can only insert to "maxSize -1"

typedef struct Dqueue Dqueue;


Dqueue *createDqueue(const int maxSize);
void freeDqueue(Dqueue *ptr);

void dumpDqueue(FILE *fp,Dqueue *ptr);

void insertTailDqueue(Dqueue *ptr, void *data);
void *delHeadDqueue(Dqueue *ptr);
void *delTailDqueue(Dqueue *ptr);

int isFullDqueue(Dqueue *ptr);
int isEmptyDqueue(Dqueue *ptr);


#endif
