#ifndef QUADLINEAR_H
#define QUADLINEAR_H

#include "parser.h"
#include "quadmatrix.h"
#include "sparsequadmatrix.h"
#include "partition_double.h"
#include "solvesparsequadmatrix.h"

// QuadMatrix *result should be nodeNum * stepNum
// result should be pre-allocated
//void quadLinearSimulation(const NetlistStampResultQuad *netlist,QuadMatrix *result,const int threadNum);


void sparseQuadLinearSimulation(const SparseNetlistQuad *netlist,QuadMatrix *result,const int threadNum,const int dumpNodeIndex);

#endif
