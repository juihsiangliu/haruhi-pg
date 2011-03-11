#ifndef QUADELEMENT_H
#define QUADELEMENT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <pthread.h>
#include "mempool.h"
#include "gsl_extern.h"
#include "mymatrix.h"


struct QuadElement
{
	int gvNum;
	double m;
	double alpha;
	double *beta;
	double *gamma;
//	gsl_matrix *beta;
//	gsl_matrix *gamma;
};
typedef struct QuadElement QuadElement;


QuadElement *createQuadElement(const int gvNum);
QuadElement *createPidQuadElement(const int gvNum,const int pid);
void freeQuadElement(QuadElement *ptr);
void freePidQuadElement(QuadElement *ptr,const int pid);

// c = a + b 
// if a & b are empty, reset c
void addQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b);
// c = a - b
// if a & b are empty, reset c
void subQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b);
// c = a * b
// if a || b is empty, reset c
void mulQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b);
void mulPidQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b,const int pid);
// c = a / b
// if a is empty, reset c
// if b is empty, do not handle this
void divQuadElement(QuadElement *c, const QuadElement *a, const QuadElement *b);
// c = s * a , s is a constant scalar
// if a is empty, it will reset c
void scaleQuadElement(QuadElement *c, const double s,const QuadElement *a);
// c = a + k, k is a constant
// a can be empty ~ will direct reset c and add k to c
void addConstantQuadElement(QuadElement *c,const double k,const QuadElement *a);

// dump out the QuadElement
void dumpQuadElement(const QuadElement *ptr);
// set the QuadElement, ptr must be allocated before
void setQuadElement(QuadElement *ptr,const double m,const double alpha,const double *beta,const double *gamma);
// copy, will reset dest if src is empty
void copyQuadElement(QuadElement *dest, const QuadElement *src);
void copyPidQuadElement(QuadElement *dest, const QuadElement *src,const int pid);
// reset
void resetQuadElement(QuadElement *ptr);
// test if it is empty (all entries are 0 or ptr is null)
inline int isEmptyQuadElement(const QuadElement *ptr);

double meanQuadElement(const QuadElement *a,const double *r);
double varQuadElement(const QuadElement *a,const double *r);
double meanPidQuadElement(const QuadElement *a,const double *r,const int pid);
double varPidQuadElement(const QuadElement *a,const double *r,const int pit);

// if ptr is empty ~ this function will allocate it
QuadElement* setZeroQuadElement(QuadElement *ptr,const int gvNum);
QuadElement* setZeroPidQuadElement(QuadElement *ptr,const int gvNum,const int pid);


// c = a * b
// if a || b is empty, reset c
void mulQuadElementQuick(QuadElement *c,const QuadElement *a,const QuadElement *b);

#endif
