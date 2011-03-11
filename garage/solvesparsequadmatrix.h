#ifndef SOLVESPARSEQUADMATRIX_H
#define SOLVESPARSEQUADMATRIX_H

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <math.h>

#include "quadelement.h"
#include "quadmatrix.h"
#include "mempool.h"
#include "sparsedoublematrix.h"
#include "sparsequadmatrix.h"
#include "partition_quad.h"
#include "parallel_lu_quad.h"

// a = lu
void luSparseQuadMatrix(SparseQuadMatrix *l, SparseQuadMatrix *u,const SparseQuadMatrix *a);
void luPidSparseQuadMatrix(SparseQuadMatrix *l, SparseQuadMatrix *u,const SparseQuadMatrix *a,const int pid);

//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveSparseQuadMatrix(QuadMatrix *x,const SparseQuadMatrix *p, const SparseQuadMatrix *pTrans, const SparseQuadMatrix *l,const SparseQuadMatrix *u, const QuadMatrix *b); 

// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b
void solveSparseQuadMatrix(QuadMatrix *x,const SparseQuadMatrix *A,const QuadMatrix *b,const int threadNum);

// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b with pre-defined permutation matrix - p and pTrans
void solveWithPermutationSparseQuadMatrix(QuadMatrix *x, const SparseQuadMatrix *p, const SparseQuadMatrix *pTrans, const SparseQuadMatrix *A,const QuadMatrix *b,const int threadNum,ParallelETree *tree);


// defined in parallel_lu_quad.h
extern void parallelLUQuad(SparseQuadMatrix *l,SparseQuadMatrix *u, ParallelETree *tree, const SparseQuadMatrix *a, const int threadNum);

#endif
