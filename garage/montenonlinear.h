#ifndef MONTENONLINEAR_H
#define MONTENONLINEAR_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#include "transfernetlist.h"
#include "mymatrix.h"
#include "mempool.h"
#include "sparsegvarient_monte.h"
#include "sparsedoublematrix.h"
#include "solvesparsedoublematrix.h"

void monteNonlinearSimulation(const MonteNetlist *netlist,double *result,const int threadNum,const int dumpNodeIndex);


#endif
