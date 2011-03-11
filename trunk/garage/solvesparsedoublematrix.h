#ifndef SOLVESPARSEDOUBLEMATRIX_H
#define SOLVESPARSEDOUBLEMATRIX_H


#include "mempool.h"
#include "sparsedoublematrix.h"
#include "partition_double.h"


#include <gsl/gsl_vector.h>


enum OrderMethod {orderUndef,orderMetis,orderAmd};


// a = lu
void luSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a);
void luPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid);

int symboic_nnz(const SparseDoubleMatrix *a);

// incomplete lu
void iluSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const double tol);
void iluPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid,const double tol);

//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *p, const SparseDoubleMatrix *pTrans, const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b); 

void triNoPSolveSparseDoubleMatrix(double *x, const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b); 

// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b
void solveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *A,const double *b);

// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b with pre-defined permutation matrix - p and pTrans
void solveWithPermutationSparseDoubleMatrix(double *x, const SparseDoubleMatrix *p, const SparseDoubleMatrix *pTrans, const SparseDoubleMatrix *A,const double *b);


// a = ll^t ... the fake version
void cholSparseDoubleMatrix(SparseDoubleMatrix *l, const SparseDoubleMatrix *a);

// used in "parallel metis lu solve"
int setGoalPartition(const SparseDoubleMatrix *mtx);

// =======================================================
//  CSR family
void triSolve_CSR_SparseDoubleMatrix(double *x,const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const CSR_SparseDoubleMatrix *l,const CSR_SparseDoubleMatrix *u,const double *b);
#endif
