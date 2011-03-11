#ifndef TRANSFERNETLIST_H
#define TRANSFERNETLIST_H

#include "parser.h"
#include "montesample.h"
#include "sparsegvarient_monte.h"
#include "sparsedoublematrix.h"
#include "sparsequadmatrix.h"
#include "mempool.h"


struct MonteNetlist
{
	enum CircuitType type;
	SparseDoubleMatrix *a;
	SparseDoubleMatrix *b;
	SparseDoubleMatrix *c;

	// u is a nodeNum * stepNum 2d double array
	// each column should be allocated from "nodalVoltagePool"
	double **u; 

	int nodeNum;
	int gvNum;
	double endTime;
	double deltaT;
	int stepNum;
	SparseGVarientTableMonte* gVarientTable;
	NonlinearMonteList *nonlinearNodeList;
};

typedef struct MonteNetlist MonteNetlist;



MonteNetlist* createMonteNetlist(const SparseNetlistQuad *src, const MonteSample *sample);

void freeMonteNetlist(MonteNetlist *monteNetlist);


#endif
