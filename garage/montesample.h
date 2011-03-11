#ifndef MONTESAMPLE_H
#define MONTESAMPLE_H

#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "parser.h"
#include "mempool.h"
#include "gsl_extern.h"
#include "quadelement.h"
#include "mymatrix.h"
#include "parser_util.h"
#include "sparsedoublematrix.h"




struct MonteSample
{
	double localVarSample;
	double *globalVarSample;
};

typedef struct MonteSample MonteSample;

MonteSample *createMonteSample(const int gvNum);
void freeMonteSample(MonteSample *ptr,const int gvNum);
void setMonteSample(MonteSample *sample,const gsl_matrix *weight,const gsl_rng *r);



//===========================================================

struct NonlinearInfoMonte
{
	const NonlinearInfo *quadAddress;

	int vdsListSize;
	int vgsListSize;
	
	// the following 3 will just save the address in "nonlinear list"
	double *vdsList; // y, row
	double *vgsList; // x, col
	double *ivTable;
	
	double cgs;
	double cgd;

	const QuadMatrix *quadPartialIdsVds;
	double *partialIdsVds;
	const QuadMatrix *quadPartialIdsVgs;
	double *partialIdsVgs;

	struct NonlinearInfoMonte *next;
};

typedef struct NonlinearInfoMonte NonlinearInfoMonte;

NonlinearInfoMonte *createAndSetNonlinearInfoMonte(const NonlinearInfo *src,const MonteSample *sample);
void freeNonlinearInfoMonte(NonlinearInfoMonte *src);

struct NonlinearMonteList
{
	NonlinearInfoMonte *list;
};

typedef struct NonlinearMonteList NonlinearMonteList;


NonlinearMonteList * createAndSetNonlinearMonteList(const NonlinearNodeList *src, const MonteSample *sample);
void freeNonlinearMonteList(NonlinearMonteList *ptr);


// --------------  misc -------------------



void getGVWeight(gsl_matrix *weight, const gsl_matrix *covariance);

gsl_rng * createRNG(void);

// given a QuadElement and a MonteSample, return the double val
double quadElement2Sample(const QuadElement *quad, const MonteSample *monte);






// ========== misc =========







#endif
