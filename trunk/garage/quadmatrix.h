#ifndef QUADMATRIX_H
#define QUADMATRIX_H


#include "quadelement.h"
#include <pthread.h>
#include <gsl/gsl_math.h>

struct QuadMatrix
{
	int row;
	int col;
	int gvNum;
	QuadElement **data;	
};

typedef struct QuadMatrix QuadMatrix;


QuadMatrix *createQuadMatrix(const int row,const int col,const int gvNum);
QuadMatrix *createPidQuadMatrix(const int row,const int col,const int gvNum,const int pid);
void freeQuadMatrix(QuadMatrix *ptr);
void freePidQuadMatrix(QuadMatrix *ptr,const int pid);

// set ptr->data[rowIndex][colIndex] = element
// will allocate a zero entry to ptr->data[rowIndex][ColIndex] if element is empty
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null
void setQuadMatrix(QuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex);
void setPidQuadMatrix(QuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,const int pid);
// return ptr->data[rowIndex][colIndex]
// the modification of the return value will effect the inside element in QuadraticMatrix "ptr"
QuadElement *getPtrEntryQuadMatrix(const QuadMatrix *ptr,const int rowIndex,const int colIndex);
// return a copy of ptr->data[rowIndex][colIndex]
QuadElement *getCopyPtrEntryQuadMatrix(QuadElement *copy, const QuadMatrix *ptr,const int rowIndex,const int colIndex);
// swap ptr->data[row1][col1] and ptr->data[row2][col2]
// direct swap the pointer
void swapQuadMatrix(QuadMatrix *ptr,const int row1,const int col1,const int row2, const int col2);
void copyQuadMatrix(QuadMatrix *dest, const QuadMatrix *src);
void copyPidQuadMatrix(QuadMatrix *dest, const QuadMatrix *src,const int pid);


// c->data[x][y] will be NULL if it is empty ... maintain the sparse structure
// c = a + b
void addQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b);
void addPidQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b,const int pid);
// c = a - b
void subQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b);
void subPidQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b,const int pid);
// c = a * b
void mulQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b);
void mulPidQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b,const int pid);
// c = k * a // k is a scale
void scaleQuadMatrix(QuadMatrix *c, const double k, const QuadMatrix *a);
void scalePidQuadMatrix(QuadMatrix *c, const double k, const QuadMatrix *a,const int pid);
// c = q * a // q is a quadelement
void scaleqQuadMatrix(QuadMatrix *c, const QuadElement *q, const QuadMatrix *a);
// b = a'
void transposeQuadMatrix(QuadMatrix *b,const QuadMatrix *a);
void transposePidQuadMatrix(QuadMatrix *b,const QuadMatrix *a,const int pid);


// dump only 1 entry
void dumpEntryQuadMatrix(const QuadMatrix *ptr,const int rowIndex,const int colIndex);
// dump the whole matrix
void dumpQuadMatrix(const QuadMatrix *ptr);

// a = plu
void luQuadMatrix(QuadMatrix *p,QuadMatrix *l,QuadMatrix *u,const QuadMatrix *a);
void luPidQuadMatrix(QuadMatrix *l,QuadMatrix *u,const QuadMatrix *a,const int pid);
//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveQuadMatrix(QuadMatrix *x,const QuadMatrix *p,const QuadMatrix *l,const QuadMatrix *u,const QuadMatrix *b);
void triSolvePidQuadMatrix(QuadMatrix *x,const QuadMatrix *l,const QuadMatrix *u,const QuadMatrix *b,const int pid);
// directly use luQuadMatrix() and triSolveQuadMatrix() to solve Ax = b
void solveQuadMatrix(QuadMatrix *x,const QuadMatrix *A,const QuadMatrix *b);
// no permutation inside ...
void solvePidQuadMatrix(QuadMatrix *x,const QuadMatrix *A, const QuadMatrix *b,const int pid);


// return the ptr of the nth column in a
void getColCopyQuadMatrix(QuadMatrix *col,const int n,const QuadMatrix *a);
// return the ptr of the nth row in a
void getRowCopyQuadMatrix(QuadMatrix *row,const int n,const QuadMatrix *a);
//set the nth col in a as col
void setColQuadMatrix(QuadMatrix *a,const QuadMatrix *col,const int n);
//set the nth row in a as row
void setRowQuadMatrix(QuadMatrix *a,const QuadMatrix *row,const int n);

void meanQuadMatrix(gsl_matrix *meanMatrix,const QuadMatrix *a,const double *r);
void varQuadMatrix(gsl_matrix *varMatrix,const QuadMatrix *a,const double *r);

void resetQuadMatrix(QuadMatrix *ptr);
// ptr need be allocated before
// all entry will be set to zero
// the null entry will be also filled
void setZeroQuadMatrix(QuadMatrix *ptr);


// c = a'*b,  a and b are n x 1 vector
void innerQuadMatrix(QuadElement *c,const QuadMatrix *a, const QuadMatrix *b);



void delQuadMatrix(QuadMatrix *ptr,const int rowIndex, const int colIndex);
void delPidQuadMatrix(QuadMatrix *ptr,const int rowIndex, const int colIndex,const int pid);


#endif

