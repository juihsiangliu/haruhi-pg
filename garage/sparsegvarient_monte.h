#ifndef SPARSEGVARIENT_MONTE_H
#define SPARSEGVARIENT_MONTE_H


#include <gdsl_queue.h>
#include <gdsl_types.h>
#include <stdlib.h>
#include "quadmatrix.h"
#include "mempool.h"
#include "sparsegvarient.h"
#include "montesample.h"



struct GControlInfoMonte
{
	// these 2 init by create()
	int vdsListSize;
	int vgsListSize;
	//========================
	int sign;
	// save in 1 based
	int gate;
	int drain;
	int source;
	//================================================

	// The following three do not allocate memory in the create funtion.
	// The actual data is saved in the "nonlinear list", 
	// which is defined in monte_sample.h
	
	double *vdsList; // y, row
	double *vgsList; // x, col
	double *partialIdsVxs;
	// if type = gm, use paritalIdsVgs
	// if type = ro_1, use paritalIdsVds
};

typedef struct GControlInfoMonte GControlInfoMonte;

GControlInfoMonte *createGControlInfoMonte(const int vdsListSize, const int vgsListSize);
void freeGControlInfoMonte(GControlInfoMonte *ptr);
void dumpGControlInfoMonte(FILE *fp,const GControlInfoMonte *ptr);
void setGControlInfoMonte(GControlInfoMonte *dest, const GControlInfo *src,const NonlinearMonteList *list);

struct SparseGControlInfoMonte
{
	int row;
	int col;
	gdsl_queue_t queue;
	
	struct SparseGControlInfoMonte *rowLink;
	struct SparseGControlInfoMonte *colLink;
};

typedef struct SparseGControlInfoMonte SparseGControlInfoMonte;

SparseGControlInfoMonte *createSparseGControlInfoMonte(void);
void freeSparseGControlInfoMonte(SparseGControlInfoMonte *ptr);
void dumpSparseGControlInfoMonte(FILE *fp,const SparseGControlInfoMonte *ptr);

struct SparseGVarientTableMonte
{
	int totalRow;
	int totalCol;
	SparseGControlInfoMonte **rowIndex;
	SparseGControlInfoMonte **colIndex;
};

typedef struct SparseGVarientTableMonte SparseGVarientTableMonte;


SparseGVarientTableMonte *createSparseGVarientTableMonte(const int row,const int col);
void freeSparseGVarientTableMonte(SparseGVarientTableMonte *ptr);

// important: element can NOT be null
// set table->[rowIndex][colIndex] = element
// will copy a new element into the talbe->[rowIndex][colIndex]
// will automatically create the entry if table->[rowIndex][colIndex] is null && element!=null
void insertSparseGVarientTableMonte(SparseGVarientTableMonte *table,GControlInfoMonte *element,const int row,const int col);

const SparseGControlInfoMonte* getSparseGVarientTableMonte(const SparseGVarientTableMonte *table,const int row,const int col);
const gdsl_queue_t getQueueSparseGVarientTableMonte(const SparseGVarientTableMonte *table,const int row,const int col);

void dumpSparseGVarientTableMonte(FILE *fp,const SparseGVarientTableMonte *table);

// this function could only be called when the "nonlinear queue" is established
void setSparseGVarientTableMonte(SparseGVarientTableMonte *dest, const SparseGVarientTable *src,const NonlinearMonteList *list); 

#endif
