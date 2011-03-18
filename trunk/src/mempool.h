#ifndef MEMPOOL_H
#define MEMPOOL_H


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "dqueue.h"



struct Mempool
{
	char poolName[64];
	int batchNumOfBlock;
	int blockSize;

	// to the free block list  -- with heap data structure
	Dqueue *dqueue;
	int count; 
	int malloc_count;

	pthread_mutex_t mutex;
};

typedef struct Mempool Mempool;



struct MempoolSet
{
	int minBlockSize;
	int maxBlockSize;
//	int hugeBlockSize;
	int numOfPool;
	int lgMinBlockSize;

	Mempool **poolSet;
	Mempool *huge; // for the very big blockSize
};

typedef struct MempoolSet MempoolSet;



// =====================================================================

//void createMempoolSet(const int minBlockSize,const int numOfMempool,const int hugeBlockSize,const int *const numOfBlockList,const int poolListSize);
void createMempoolSet(const int minBlockSize,const int numOfMempool,const int *const numOfBlockList,const int poolListSize);
void freeMempoolSet(void);


void clearPidMempoolSet(const int pid);

void *getMempoolSet(const int blockSize);
void *getPidMempoolSet(const int blockSize,const int pid);

void retMempoolSet(void *data, const int blockSize);
void retPidMempoolSet(void *data, const int blockSize,const int pid);

void dumpMempoolSet(FILE *fp);
void usageMempoolSet(FILE *fp);


#endif
