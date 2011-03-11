#ifndef PARALLEL_PCG_DOUBLE_H
#define PARALLEL_PCG_DOUBLE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>

#include "sparsedoublematrix.h"
#include "partition_double.h"
#include "mempool.h"
#include "parallel_common.h"
#include "mymatrix.h"
#include "solvesparsedoublematrix.h"
#include "parallel_lu_double.h"


void parallelPCG(const SparseDoubleMatrix *a,const double *b, double *sol,const int threadNum,enum OOCFlag oocFlag,enum OrderMethod orderMethod);


#endif






