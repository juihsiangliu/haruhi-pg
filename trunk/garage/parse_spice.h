#ifndef PARSE_SPICE_H
#define PARSE_SPICE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "sparsedoublematrix.h"
#include "mempool.h"
#include "union_find.h"

#define NAME_LEN (32)


enum ElementType  {vs,res,cs,comment}; // voltage source , resistor , current source , comment


struct SpiceLine
{
	enum ElementType type;
	int pos; // in the beginning ~ it is pos_hash, later it stores pos_ind
	int neg; // in the beginning ~ it is neg_hash, later it stores neg_ind
	double val;
};

typedef struct SpiceLine SpiceLine;


struct ShortEntry
{
	int pos_ind;
	int neg_ind;
};

typedef struct ShortEntry ShortEntry;



enum SolveMethod {direct,iterative};


struct SpiceMtx
{
	int nodeNum;
	SparseDoubleMatrix *gMtx;
	double *current;
	double tol; // used in ILU & PCG
	char fileNameList[NAME_LEN];
	char **compactNameTable;
	enum SolveMethod method;
	double *nodalVoltage; // the unknown

	// ------------------
	Set *set;
};

typedef struct SpiceMtx SpiceMtx;


SpiceMtx* parseSpice(const char *filename);
SpiceMtx* parseSpice_fp(FILE *fp);

void freeSpiceMtx(SpiceMtx *ptr);

void parse2LucyFormat(const char *filename);

#endif
