#ifndef MY_MATRIX_H
#define MY_MATRIX_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>


void setMyMatrix(double *a,double val,const int rowA,const int colA,const int i, const int j);

double getMyMatrix(const double *a,const int rowA,const int colA,const int i, const int j);

int isNullMyMatrix(const double *a,const int rowA,const int colA);

void addMyMatrix(double *c,const double *a,const double *b, const int rowA,const int colA);

void subMyMatrix(double *c,const double *a,const double *b, const int rowA,const int colA);

void mulMyMatrix(double *c,const double *a,const double *b, const int rowA,const int colA, const int rowB, const int colB);

void scaleMyMatrix(double *c, const double k, const double *a,const int rowA,const int colA);

// the address of a and c can not be the same
void transMyMatrix(double *c, const double *a,const int rowA,const int colA);

double traceMyMatrix(double *a,const int rowA,const int colA);


void copyMyMatrix(double *dest, const double *src, const int row, const int col);


double absMaxMyMatrix(const double *src, const int row, const int col);

#endif
