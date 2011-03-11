#ifndef SPARSEDOUBLEMATRIX_H
#define SPARSEDOUBLEMATRIX_H

#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <pthread.h>

#include "mymatrix.h"
#include "mempool.h"
#include "parallel.h"


struct SparseDoubleElement
{
	int row;
	int col;
	struct SparseDoubleElement *rowLink;
	struct SparseDoubleElement *colLink;
	double data;
};

typedef struct SparseDoubleElement SparseDoubleElement;

SparseDoubleElement * createSparseDoubleElement(const double data);
SparseDoubleElement * createPidSparseDoubleElement(const double data,const int pid);

void freeSparseDoubleElement(SparseDoubleElement *ptr);
void freePidSparseDoubleElement(SparseDoubleElement *ptr,const int pid);

//=============================================

struct SparseDoubleMatrix
{
	int totalRow;
	int totalCol;
	
	long long nnz;
	SparseDoubleElement **rowIndex;
	SparseDoubleElement **colIndex;

	SparseDoubleElement **rowCache;
	SparseDoubleElement **colCache;
};

typedef struct SparseDoubleMatrix SparseDoubleMatrix;

SparseDoubleMatrix *createSparseDoubleMatrix(const int row,const int col);
SparseDoubleMatrix *createPidSparseDoubleMatrix(const int row,const int col,const int pid);

void freeSparseDoubleMatrix(SparseDoubleMatrix *ptr);
void freePidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const int pid);

// important: element can NOT be null
// set ptr->data[rowIndex][colIndex] = element
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null
void setSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex);
void setPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,const int pid);


// will return 0 if the entry is empty
//double getSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex);
//double getFastRowSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex, SparseDoubleElement **baseRow);
//double getFastColSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex, SparseDoubleElement **baseCol);
//double getNewSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex,const char *type);
double getSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex,const char *type);


// will skip the entry if it is empty
void delSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex);
void delPidSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex,const int pid);

// dump out row by row
void dumpSparseDoubleMatrix(FILE *fp,const SparseDoubleMatrix *ptr);

void plotSparseDoubleMatrix(FILE *fp, const SparseDoubleMatrix *ptr);

// swap 2 rows // suck implement
void swapRowSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int row1, const int row2);

void copySparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src);
void copyPidSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src,const int pid);

// clear the SparseMatrix to the initial state
void clearSparseDoubleMatrix(SparseDoubleMatrix *ptr);
void clearPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const int pid);

// c = a + b
void addSparseDoubleMatrix(SparseDoubleMatrix *c,const SparseDoubleMatrix *a,const SparseDoubleMatrix *b);

// c = a - b
void subSparseDoubleMatrix(SparseDoubleMatrix *c,const SparseDoubleMatrix *a,const SparseDoubleMatrix *b);

// c = a * b
void mulSparseDoubleMatrix(SparseDoubleMatrix *c,const SparseDoubleMatrix *a,const SparseDoubleMatrix *b);

// c = k * a , k is a scale
void scaleSparseDoubleMatrix(SparseDoubleMatrix *c,const double k,const SparseDoubleMatrix *a);

// c = trans(A)
void transSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a);


// c = a * b, (sparse multiply dense)
// a is a sparse matrix
// b is a n*1 double Matrix
// c is a n*1 double Matrix
void mulVecSparseDoubleMatrix(double *c,const SparseDoubleMatrix *a, const double *b);


struct ParOfMul2
{
	int rowBegin;
	int rowEnd;;
	// process from rowBegin to rowEnd-1
	const SparseDoubleMatrix *a;
	const double *b;
	double *c;
};

void parallelMulVecSparseDoubleMatrix(double *c,const SparseDoubleMatrix *a, const double *b,const int threadNum);

// set a as identity matrix
void identitySparseDoubleMatrix(SparseDoubleMatrix *a);
void identityPidSparseDoubleMatrix(SparseDoubleMatrix *a,const int pid);


// a(i,j) = a(i,j) + element
// will allocate the memory if necessary
void incSparseDoubleMatrix(SparseDoubleMatrix *a,const double element,const int row, const int col);


// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
void decSparseDoubleMatrix(SparseDoubleMatrix *a,const double element,const int row, const int col);
SparseDoubleElement * decFastSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,SparseDoubleElement *baseRow, SparseDoubleElement *baseCol);

// dest and src should be allocated before calling this function
// the src is assumed to solve as myMatrix format
// the dimension of src and dest should be the same
void dense2SparseDoubleMatrix(SparseDoubleMatrix *dest,const double *src);
void dense2PidSparseDoubleMatrix(SparseDoubleMatrix *dest,const double *src,const int pid);

void sparse2DenseDoubleMatrix(double *dest,const SparseDoubleMatrix *src);


// dest(i,col) will be delete if isEmptyQuadElement(src(i,0)) is true
void setDenseCol2SparseDoubleMatrix(SparseDoubleMatrix *dest, const double *src,const int totalRow, const int targetCol);

// just set src(i,0) to dest(i,col) directly without delete step.
// ( so, the dest might be clear before calling this function)
void setDenseColQuick2SparseDoubleMatrix(SparseDoubleMatrix *dest, const double *src, const int totalRow, const int targetCol);


// the boundaries row1,col1,row2,col2 are included
void clearBlockSparseDoubleMatrix(SparseDoubleMatrix *dest, const int row1, const int col1, const int row2, const int col2);

// the boundaries row1,col1,row2,col2 are included
void appendBlockSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int row1, const int col1, const int row2, const int col2);


// merge 2 matrix dest and src
// append the whole src into dest, the "base point" of the left top corner is (ltRowDwst,ltColDest)
void mergeSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest);
void mergePidSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest,const int pid);

void getSubSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc);
void getPidSubSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc,const int pid);


// permutate 
void permutateSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *pRow, const SparseDoubleMatrix *pCol,const SparseDoubleMatrix *src);
// src will be deleted during permutation ... to save memory
void dpermutateSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *pRow, const SparseDoubleMatrix *pCol,SparseDoubleMatrix *src);

// inverse lowerTriangular
void invLTSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src);

// read the sparse matrix from file
// the file should store in binary, and csr format
// totalRow totalCol nnz
// rowptr...colind...val...
SparseDoubleMatrix *read_pid_SparseDoubleMatrix(const char *filename,const int pid);

// write the sparse matrix to file
// the file will store in binary, and csr format
// totalRow totalCol nnz
// rowptr...colind...val...
void write_SparseDoubleMatrix(const char *filename,const SparseDoubleMatrix *ptr);

// read the sparse matrx from file
// the file should store in txt, and ind format
// totalRow totalCol nnz
// row1 col1 val1
// row2 col2 val2
// ...
SparseDoubleMatrix *read_ind_SparseDoubleMatrix(const char *filename);
void write_ind_SparseDoubleMatrix(const char *filename,const SparseDoubleMatrix *ptr);

double colNormSparseDoubleMatrix(const SparseDoubleMatrix *ptr,const int col);



struct CSR_SparseDoubleMatrix
{
	void *buf; // mmap addr
	int pid;
	int totalRow;
	int totalCol;
	long long nnz;

	int *rowPtr; // 1 * totalrow+1
	int *col; // 1 * nnz
	double *val; // 1 * nnz
};

typedef struct CSR_SparseDoubleMatrix CSR_SparseDoubleMatrix;


CSR_SparseDoubleMatrix *create_CSR_SparseDoubleMatrix(int row, int col, long long nnz,int pid);
void free_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx);
CSR_SparseDoubleMatrix *read_to_CSR_SparseDoubleMatrix(const char *filename,const int pid);
CSR_SparseDoubleMatrix *linus_read_to_CSR_SparseDoubleMatrix(const char *filename,const int pid);
void linus_free_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx);
CSR_SparseDoubleMatrix *write_CSR_SparseDoubleMatrix(const char *filename,CSR_SparseDoubleMatrix *mtx);
void expand_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx,const long long extra);
void copy_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *dest, const CSR_SparseDoubleMatrix *src);
void dump_CSR_SparseDoubleMatrix(FILE *fp, const CSR_SparseDoubleMatrix *mtx);
void lu_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *l, CSR_SparseDoubleMatrix *u, const CSR_SparseDoubleMatrix *a);
CSR_SparseDoubleMatrix *sparse2CSR(const SparseDoubleMatrix *mtx, const int pid);
void clear_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx);
void mulVec_CSR_SparseDoubleMatrix(double *res,const CSR_SparseDoubleMatrix *mtx,const double *vec); // res and vec can not be the same
void parallel_mulVec_CSR_SparseDoubleMatrix(double *res,const CSR_SparseDoubleMatrix *mtx,const double *vec,const int threadNum); // res and vec can not be the same


struct CSC_SparseDoubleMatrix
{
	int pid;
	int totalRow;
	int totalCol;
	long long nnz;

	int *row;
	int *colPtr;
	double *val;
};

typedef struct CSC_SparseDoubleMatrix CSC_SparseDoubleMatrix;


CSC_SparseDoubleMatrix *create_CSC_SparseDoubleMatrix(int row, int col, long long nnz,int pid);
void free_CSC_SparseDoubleMatrix(CSC_SparseDoubleMatrix *mtx);
CSC_SparseDoubleMatrix *sparse2CSC(const SparseDoubleMatrix *mtx, const int pid);


struct Linus_read
{
	char *name;
	int size;
	void *buf;
};

typedef struct Linus_read Linus_read;




struct Parallel_mulVec
{
	int p;
	int threadNum;
	double *res;
	const CSR_SparseDoubleMatrix *mtx;
	const double *vec;
};

typedef struct Parallel_mulVec Parallel_mulVec;


#endif
