#ifndef PARALLEL_LU_QUAD_H
#define PARALLEL_LU_QUAD_H

#include <gdsl_types.h>
#include <gdsl_queue.h>
#include <gsl_math.h>

#include "sparsequadmatrix.h"
#include "partition_double.h"
#include "parallel_lu_common.h"
#include "solvesparsequadmatrix.h"



struct ALUQuad
{
	SparseQuadMatrix *a;
	SparseQuadMatrix *l;
	SparseQuadMatrix *u;
};

typedef struct ALUQuad ALUQuad;



struct ParallelDoneList
{
	int orderList[15];
	int done[16]; // set in 1 base
	int gvNum;
	int eachNodeCurrentDone[16]; // set in 1 base and init all elements to zero

	ALUQuad **alu; // the same as it in ParallelLUQuadShareData
	const ParallelETree *tree; // the same as it in ParallelLUQuadShareData
	const SparseQuadMatrix *a; // the same as it in ParallelLUQuadShareData
};

typedef struct ParallelDoneList ParallelDoneList;




struct ParallelLUQuadShareDataNew
{
	ParallelDoneList *doneList;
	ToDoList *todolist;	

	// internal used mutex and cv
	pthread_mutex_t *mutex;
	pthread_cond_t *cond;

	int pid;
	int N;
	int rootCurrentBegin;
	int currentEnd;
};

typedef struct ParallelLUQuadShareDataNew ParallelLUQuadShareDataNew;



struct ThreadHandlerPar
{
	struct ParallelLUQuadShareDataNew **list;
	int threadNum;
};




void parallelLUQuad(SparseQuadMatrix *l,SparseQuadMatrix *u, ParallelETree *tree, const SparseQuadMatrix *a, const int threadNum);






#endif
