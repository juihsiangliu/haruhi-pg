#ifndef MONTELINEAR_H
#define MONTELINEAR_H

#include "montesample.h"
#include "transfernetlist.h"
#include "sparsedoublematrix.h"
#include "solvesparsedoublematrix.h"
#include "mempool.h"
#include "parallel_lu_double.h"
#include "partition_double.h"

void monteLinearSimulation(const MonteNetlist *netlist, double *result,const int threadNum,const int dumpNodeIndex);

#endif
