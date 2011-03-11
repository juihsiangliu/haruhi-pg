#ifndef PLOT_H
#define PLOT_H

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "quadmatrix.h"
#include "parser.h"
#include "transfernetlist.h"
#include "mymatrix.h"

// dump the result ( 2 x (nodeNum*stepSize) )
// mean node 1 (step 1~n), mean node 2 (step 1~n), ...
// var node 1 (step 1~n), var node 2 (step 1~n), ...

//void dumpQuadResult(FILE *fp,const QuadMatrix *result,const NetlistStampResultQuad *netlist);

void dumpSparseQuadResult(FILE *fp,const QuadMatrix *result,const SparseNetlistQuad *netlist);

// dump the result for a monteSample
// result ( 1 x (nodeNum*stepSize) )
// val node1 (step 1~n), val node 2 (step 1~n), ...
void dumpMonteResult(FILE *fp,double *result, const MonteNetlist *netlist);

#endif
