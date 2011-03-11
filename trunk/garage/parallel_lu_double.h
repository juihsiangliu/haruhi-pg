#ifndef PARALLEL_LU_DOUBLE_H
#define PARALLEL_LU_DOUBLE_H

#include <string.h>
#include <pthread.h>
#include <assert.h>
#include <gdsl_types.h>
#include <gdsl_queue.h>
#include <gsl_math.h>

#include "sparsedoublematrix.h"
#include "partition_double.h"
#include "mempool.h"
#include "parallel_common.h"
#include "mymatrix.h"
#include "solvesparsedoublematrix.h"

struct ALUDouble
{
	// generate by createALU
	int aRow;
	int aCol;
	int lRow;
	int lCol;
	int uRow;
	int uCol;

	SparseDoubleMatrix *l;
	SparseDoubleMatrix *u;

	// used in parallel ic mode
	pthread_mutex_t mutex; // used to contronl the setting of csr_u & csr_u_count
	int csr_u_count;
	CSR_SparseDoubleMatrix *csr_u;
};

typedef struct ALUDouble ALUDouble;




struct ParallelDoneListDouble
{
	int treeInternalSize; // not include the "virtual" 0 node and all the dummy "null" leave 
	const int* orderList; // size is treeInternalSize , store in 0 base
	int* done; // size is treeInternalSize+1 , store in 1 base
	int* eachNodeCurrentDone; // size is treeInternalSize+1 , store in 1 base and init all elements to zero

	ALUDouble **alu; // the same as it in ParallelLUDoubleShareData
	const ParallelETree *tree; // the same as it in ParallelLUDoubleShareData
//	const SparseDoubleMatrix *a; // the same as it in ParallelLUDoubleShareData
	// if the ooc mode is used, the address of a will be changed during the process
};

typedef struct ParallelDoneListDouble ParallelDoneListDouble;



struct OOCInfo
{
	int node; // the index of aluList
	char type; // 'L' or 'U'
	char name[32];
	int rowBegin;
	int rowEnd; // from begin to end-1
	int count; // how many nodes to read this node
	int totalLength;
	// treeInternalSize = (totalLength - 1)/2

	pthread_mutex_t mutex;

	// only store at node 0
	int *postorder;
	int *baseRowL;
	int *baseColL;
	int *baseRowU;
	int *baseColU;

	SparseDoubleMatrix *a;
};

typedef struct OOCInfo OOCInfo;

struct OOCInfo *createOOCInfoList(ALUDouble **aluList, const int treeInternalSize, const int *postorder,enum OOCFlag oocFlag,const int *baseRowL, const int *baseColL, const int *baseRowU,const int *baseColU);
void freeOOCInfoList(struct OOCInfo *ptr);
void readByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList);
void pseudoWriteByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList);
void loadByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList);
void writeByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList);


struct ParallelLUDoubleShareData
{
	ParallelDoneListDouble *doneList;
	ToDoList *todolist;	

	// internal used mutex and cv
	pthread_mutex_t *mutex;
	pthread_cond_t *cond;
	
	const double *colNormA;
	double tol;
	int pid;
	int N;
	int rootCurrentBegin;
	int *rootCurrentBeginList; // one based ... as treeSize
	int currentEnd;
	int maxThreadNum;

	const int *baseRowL;
	const int *baseColL;
	const int *baseRowU;
	const int *baseColU;
	struct OOCInfo *oocInfoList;
	enum OOCFlag oocFlag;
};

typedef struct ParallelLUDoubleShareData ParallelLUDoubleShareData;


struct ThreadHandlerParDouble
{
	ParallelLUDoubleShareData **list;
	int threadNum;
};

typedef struct ThreadHandlerParDouble ThreadHandlerParDouble;


struct OOCInfo * parallelLUDouble(SparseDoubleMatrix *l,SparseDoubleMatrix *u, ParallelETree *tree, SparseDoubleMatrix *a, const SparseDoubleMatrix *p, const int threadNum,enum OOCFlag oocFlag);

// incomplete LU
struct OOCInfo * parallelILUDouble(SparseDoubleMatrix *l,SparseDoubleMatrix *u, ParallelETree *tree, SparseDoubleMatrix *a, const SparseDoubleMatrix *p, const int threadNum,const double tol,enum OOCFlag oocFlag);


void oocTriSolveSparseDoubleMatrix(double *x, struct OOCInfo *oocInfoList, const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const double *b);



struct PipeSolve
{
	int treeInternalSize;
	int *readFlag;
	int *shiftFlag;
	int *execFlag;
	int *freeFlag;

	const int *baseRow;
	const int *baseCol;
	const int *postorder;

	CSR_SparseDoubleMatrix **mtx_list;
	const OOCInfo *oocInfoList;
	const SparseDoubleMatrix *perm;
	int base;
};

typedef struct PipeSolve PipeSolve;



#endif
