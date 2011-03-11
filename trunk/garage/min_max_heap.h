#ifndef MIN_MAX_HEAP_H
#define MIN_MAX_HEAP_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct MinMaxHeap
{
	int maxSize;
	int size;
	void **heap;
};

typedef struct MinMaxHeap MinMaxHeap;


MinMaxHeap *createMinMaxHeap(const int maxSize);
void freeMinMaxHeap(MinMaxHeap *ptr);

void min_max_dump(FILE *fp,MinMaxHeap *ptr);

int min_max_is_empty(MinMaxHeap *ptr);
int min_max_is_full(MinMaxHeap *ptr);

void *min_max_get_min(MinMaxHeap *ptr);
void *min_max_get_max(MinMaxHeap *ptr);

void min_max_insert(MinMaxHeap *ptr, void *data);



#endif
