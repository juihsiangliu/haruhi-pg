#ifndef PARSER_H
#define PARSER_H

#include "fgetl.h"
#include "sparsegvarient.h"
#include "sparsequadmatrix.h"
#include "quadmatrix.h"
#include "dlmread.h"
#include "quadelement.h"
#include "parser_util.h"
#include "mempool.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gdsl_2darray.h>
#include <gdsl_queue.h>
#include <gdsl_types.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// gVarientTable is a 2dArray, each entry in the array is a queue.
// And each element in queue is GControlInfo.
// It save the information that the mos after stamping.
// There is no create and free function for GControlInfo
// It is generated and allocated by calling parseQuadInput
// The free function is implemented in "freeNetlistStampResultQuad".



enum CircuitType{linear,nonlinear};


//===================================================================
//===================================================================
//===================================================================
//===================================================================
//        [NEW] FULL SPARSE FORMAT SUPPORT



struct SparseNetlistQuad
{
	enum CircuitType type;
	SparseQuadMatrix *a;
	SparseQuadMatrix *b;
	SparseQuadMatrix *c;
	QuadMatrix *u;
	double *s;
//	gsl_matrix *s;
	int nodeNum;
	int gvNum;
	double endTime;
	double deltaT;
	int stepNum;
	SparseGVarientTable *gVarientTable;

	int numOfLinear;
	int numOfNonLinear;


	NonlinearNodeList *nonlinearNodeList;
};

typedef struct SparseNetlistQuad SparseNetlistQuad;


// it will create SparseNetslistQuad
// only stamp type, s,nodeNum, gvNum, endTime, deltaT, stepNum
// used for allocate the memory pool
SparseNetlistQuad *symbolicParseSparseQuadInput(const char *filename);



// no create function is provided for api
// it will be allocated by the parse funtion
SparseNetlistQuad *parseSparseQuadInput(SparseNetlistQuad *ptr,const char *filename);
void freeSparseNetlistQuad(SparseNetlistQuad *ptr);



#endif
