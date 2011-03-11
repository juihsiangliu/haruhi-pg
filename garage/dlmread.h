#ifndef DLMREAD_H
#define DLMREAD_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_cblas.h"
#include "fgetl.h"

// dlmread(M,filename,delimiter,R,C)  reads numeric data from the ASCII delimited file filename, using the specified delimiter.
// The values R C specify the row and column where the upper-left corner of and the down-right corner the data lies in the file.
// R and C are zero based.
void dlmread(gsl_matrix *M,const char *filename,const char *delimiter,const int R1,const int C1,const int R2,const int C2);

#endif
