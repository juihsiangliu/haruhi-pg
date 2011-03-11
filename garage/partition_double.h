#ifndef PARTITION_DOUBLE_H
#define PARTITION_DOUBLE_H

#include <limits.h>
#include <gdsl_queue.h>
#include <gdsl_types.h>
#include <metis.h>

#include "sparsedoublematrix.h"
#include "mempool.h"


enum OOCFlag{ic,ooc,undef};


// ==================================================================
// 		start of the internal structures
// =================================================================

enum LockStat{lock,unlock};
enum VisitStat{visit,notvisit};

struct PartitionResult
{
	int partitionSize;
	int* partA; // element = 0 ~ size-1
	int* partB;

	int partASize;
	int partBSize;

	int *crossList;
	int crossListSize;

	// ===============
	enum VisitStat visitLog;
	// ===============  used for partitionMetis
	int partSerialBegin;
	int partSerialEnd;
};	

typedef struct PartitionResult PartitionResult;


// because it is assumed as complete binary tree, use array format to save this tree
struct EliminationTree
{
	int size; 
	// total node number of the elimination tree
	// the "empty head node" is also included
	PartitionResult **node;
	// the tree nodes
	
	// ============================

	int count;
	// used for internal insertation function
};

typedef struct EliminationTree EliminationTree;


// ==================================================================
// 		end of the internal structures
// ==================================================================


enum ParallelEtreeNodeType{undefine,lu,cross};

struct ParallelETreeNode
{
	int rowBegin;
	int rowEnd;
	int doneRowBegin;
	int doneRowEnd;
	enum VisitStat visitLog;
	enum ParallelEtreeNodeType type;
};

typedef struct ParallelETreeNode ParallelETreeNode;


struct ParallelETree
{
	int size;
	// total node number of the elimination tree
	// the "empty head node" is also included
	// also include the "null" dummy leaf
	ParallelETreeNode **node;
	// the tree nodes
};

typedef struct ParallelETree ParallelETree;

ParallelETree *createParallelETree(const int size);
void freeParallelETree(ParallelETree *ptr);



// =====================================================================================

// use amd directly to reorder the matrix a
void amdSparseDoubleMatrix(SparseDoubleMatrix *p,const SparseDoubleMatrix *a);



SparseDoubleMatrix * partitionSparseDoubleMatrix(SparseDoubleMatrix *p,SparseDoubleMatrix *pTrans,ParallelETree *tree,const SparseDoubleMatrix *a,const int goalPartition, enum OOCFlag oocFlag);

#endif
