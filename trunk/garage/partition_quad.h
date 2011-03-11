#ifndef PARTITION_QUAD_H
#define PARTITION_QUAD_H

#include <limits.h>
#include <gdsl_queue.h>
#include <gdsl_types.h>
#include <amd.h>

#include "partition_double.h"
#include "sparsedoublematrix.h"
#include "sparsequadmatrix.h"
#include "mempool.h"




// =====================================================================================

// amd reordering
void amdSparseQuadMatrix(SparseQuadMatrix *p,const SparseQuadMatrix *a);



void partitionSparseQuadMatrix(SparseQuadMatrix *p,SparseQuadMatrix *pTrans,ParallelETree *tree,SparseQuadMatrix *aRefine,const SparseQuadMatrix *a,const int goalPartition);

#endif
