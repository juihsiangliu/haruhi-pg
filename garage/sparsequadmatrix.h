#ifndef SPARSEQUADMATRIX_H
#define SPARSEQUADMATRIX_H

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <math.h>

#include "quadelement.h"
#include "quadmatrix.h"
#include "mempool.h"
#include "sparsedoublematrix.h"

struct SparseQuadElement
{
	QuadElement *data;
	int row;
	int col;
	struct SparseQuadElement *rowLink;
	struct SparseQuadElement *colLink;
};

typedef struct SparseQuadElement SparseQuadElement;

SparseQuadElement * createSparseQuadElement(const int gvNum);
void freeSparseQuadElement(SparseQuadElement *ptr);
SparseQuadElement * createPidSparseQuadElement(const int gvNum,const int pid);
void freePidSparseQuadElement(SparseQuadElement *ptr,const int pid);

//=============================================

struct SparseQuadMatrix
{
	int totalRow;
	int totalCol;
	int gvNum;
	int nnz;
	SparseQuadElement **rowIndex;
	SparseQuadElement **colIndex;
};

typedef struct SparseQuadMatrix SparseQuadMatrix;

SparseQuadMatrix *createSparseQuadMatrix(const int row,const int col,const int gvNum);
void freeSparseQuadMatrix(SparseQuadMatrix *ptr);
SparseQuadMatrix *createPidSparseQuadMatrix(const int row,const int col,const int gvNum,const int pid);
void freePidSparseQuadMatrix(SparseQuadMatrix *ptr,const int pid);

// important: element can NOT be null
// set ptr->data[rowIndex][colIndex] = element
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null
void setSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex);
void setPidSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,const int pid);
void setFastSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,SparseQuadElement **baseRow,SparseQuadElement **baseCol);
void setFastPidSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,SparseQuadElement **baseRow,SparseQuadElement **baseCol,const int pid);

// will return null if the entry is empty
const QuadElement *getSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex);
const QuadElement *getFastRowSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex, SparseQuadElement ** baseRow);
SparseQuadElement *getPtrFastRowSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex,SparseQuadElement **baseRow);
const QuadElement *getFastColSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex, SparseQuadElement ** baseCol);

// will skip the entry if it is empty
void delSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex,const int colIndex);
void delPidSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex,const int colIndex,const int pid);
void delFastSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex, const int colIndex, SparseQuadElement **baseRow, SparseQuadElement **baseCol);
void delFastPidSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex, const int colIndex, SparseQuadElement **baseRow, SparseQuadElement **baseCol,const int pid);

// dump out row by row
void dumpSparseQuadMatrix(FILE *fp,const SparseQuadMatrix *ptr);
void dumpHeadSparseQuadMatrix(FILE *fp,const SparseQuadMatrix *ptr,const char *name);

// swap 2 rows // suck implement
void swapRowSparseQuadMatrix(SparseQuadMatrix *ptr, const int row1, const int row2);

void copySparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src);
void copyPidSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src,const int pid);

// clear the SparseMatrix to the initial state
void clearSparseQuadMatrix(SparseQuadMatrix *ptr);
void clearPidSparseQuadMatrix(SparseQuadMatrix *ptr,const int pid);

// c = a + b
void addSparseQuadMatrix(SparseQuadMatrix *c,const SparseQuadMatrix *a,const SparseQuadMatrix *b);
void addPidSparseQuadMatrix(SparseQuadMatrix *c,const SparseQuadMatrix *a,const SparseQuadMatrix *b,const int pid);

// c = a - b
void subSparseQuadMatrix(SparseQuadMatrix *c,const SparseQuadMatrix *a,const SparseQuadMatrix *b);
void subPidSparseQuadMatrix(SparseQuadMatrix *c,const SparseQuadMatrix *a,const SparseQuadMatrix *b,const int pid);

// c = a * b
void mulSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a, const SparseQuadMatrix *b);

// c = k * a , k is a scale
void scaleSparseQuadMatrix(SparseQuadMatrix *c,const double k,const SparseQuadMatrix *a);
void scalePidSparseQuadMatrix(SparseQuadMatrix *c,const double k,const SparseQuadMatrix *a,const int pid);

// c = trans(A)
void transSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a);

// c = a * b, (sparse multiply dense)
// a is a sparse matrix
// b is a n*1 QuadMatrix
// c is a n*1 QuadMatrix
void mulVecSparseQuadMatrix(QuadMatrix *c,const SparseQuadMatrix *a, const QuadMatrix *b);
void mulVecPidSparseQuadMatrix(QuadMatrix *c,const SparseQuadMatrix *a, const QuadMatrix *b,const int pid);

// set a as identity matrix
void identitySparseQuadMatrix(SparseQuadMatrix *a);
void identityPidSparseQuadMatrix(SparseQuadMatrix *a,const int pid);



// a(i,j) = a(i,j) + element
// will allocate the memory if necessary
void incSparseQuadMatrix(SparseQuadMatrix *a,const QuadElement *element,const int row, const int col);


// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
void decSparseQuadMatrix(SparseQuadMatrix *a,const QuadElement *element,const int row, const int col);

// dest and src should be allocated before calling this function
void quad2SparseQuadMatrix(SparseQuadMatrix *dest,const QuadMatrix *src);
void quad2PidSparseQuadMatrix(SparseQuadMatrix *dest,const QuadMatrix *src,const int pid);
void delQuad2PidSparseQuadMatrix(SparseQuadMatrix *dest, QuadMatrix *src,const int pid);

// dest and src should be allocated before calling this function
void toDenseSparseQuadMatrix(QuadMatrix *dest, const SparseQuadMatrix *src);
void toDensePidSparseQuadMatrix(QuadMatrix *dest, const SparseQuadMatrix *src,const int pid);

// dest(i,col) will be delete if isEmptyQuadElement(src(i,0)) is true
void setDenseCol2SparseQuadMatrix(SparseQuadMatrix *dest, const QuadMatrix *src,const int col);

// just set src(i,0) to dest(i,col) directly without delete step.
// ( so, the dest might be clear before calling this function)
void setDenseColQuick2SparseQuadMatrix(SparseQuadMatrix *dest, const QuadMatrix *src,const int col);
void setDenseColQuick2SparseQuadMatrixPid(SparseQuadMatrix *dest, const QuadMatrix *src,const int col,const int pid);


// the boundaries row1,col1,row2,col2 are included
void clearBlockSparseQuadMatrix(SparseQuadMatrix *dest, const int row1, const int col1, const int row2, const int col2);

// the boundaries row1,col1,row2,col2 are included
void appendBlockSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int row1, const int col1, const int row2, const int col2);

// get the m part of sparseQuadMatrix
void mSparseQuadMatrix(SparseDoubleMatrix *dest, const SparseQuadMatrix *src);

// set the m part of sparseQuadMatrix
void setMSparseQuadMatrix(SparseQuadMatrix *dest, const SparseDoubleMatrix *src);

// permutate 
void permutateSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *pRow, const SparseQuadMatrix *pCol, const SparseQuadMatrix *src);

// get the sub matrix
void getSubSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc);
void getSubPidSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc,const int pid);

// merge the sub matrix
void mergeSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest);
void mergePidSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest,const int pid);

#endif
