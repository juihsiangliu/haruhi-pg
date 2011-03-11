#ifndef QUADNONLINEAR_H
#define QUADNONLINEAR_H

#include "quadmatrix.h"
#include "parser.h"
#include "sparsequadmatrix.h"
#include "mempool.h"
#include "mymatrix.h"
#include "sparsegvarient.h"
#include "solvesparsequadmatrix.h"
#include "parallel.h"

// used for the nonlinear simulation
struct GVarientIJQuad
{
	QuadElement *gVarientElement; 
	// need be allocated in updateGVaroemtMatrixFinite()
	// will be updated in updateGVarientMatrixFinite() => map_fun_setGvarient()
	const QuadMatrix *v;
	const double* s; 
	int gvNum;
	int pid;
};

typedef struct GVarientIJQuad GVarientIJQuad;



// quadMatrix *result should be nodeNum * stepNum
// result should be pre-allocated
void quadNonlinearSimulation(const SparseNetlistQuad *netlist,QuadMatrix *result,const int threadNum,const int dumpNodeIndex);


#endif
