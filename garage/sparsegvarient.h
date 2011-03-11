#ifndef SPARSEGVARIENT_H
#define SPARSEGVARIENT_H


#include <gdsl_queue.h>
#include <gdsl_types.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <stdlib.h>
#include <pthread.h>
#include "quadmatrix.h"
#include "mempool.h"


enum MosRole{gm,ro_1}; // ro_1 = 1/ro

struct GControlInfo
{
	// these 3 init by create()
	int vdsListSize;
	int vgsListSize;
	int gvNum;
	int pid;
	//========================
	int sign;
	enum MosRole type;
	// save in 1 based
	int gate;
	int drain;
	int source;
	//================================================

	// The following three do not allocate memory in the create funtion.
	// The actual data is saved in the "nonlinear list", 
	// which is defined in parser_util.h

	double *vdsList; // y, row
	double *vgsList; // x, col
	QuadMatrix *partialIdsVxs; 
	// if type = gm, use paritalIdsVgs
	// if type = ro_1, use paritalIdsVds
};

typedef struct GControlInfo GControlInfo;

GControlInfo *createGControlInfo(const int vdsListSize, const int vgsListSize,const int gvNum);
GControlInfo *createGControlInfoPid(const int vdsListSize, const int vgsListSize,const int gvNum,const int pid);
void freeGControlInfo(GControlInfo *ptr);
void freeGControlInfoPid(GControlInfo *ptr,const int pid);
void dumpGControlInfo(FILE *fp,const GControlInfo *ptr);
//GControlInfo *copyGControlInfo(const GControlInfo *const src);


struct SparseGControlInfo
{
	int row;
	int col;
	gdsl_queue_t queue;
	
	struct SparseGControlInfo *rowLink;
	struct SparseGControlInfo *colLink;
};

typedef struct SparseGControlInfo SparseGControlInfo;

SparseGControlInfo *createSparseGControlInfo(void);
void freeSparseGControlInfo(SparseGControlInfo *ptr);
void dumpSparseGControlInfo(FILE *fp,const SparseGControlInfo *ptr);

struct SparseGVarientTable
{
	int totalRow;
	int totalCol;
	int gvNum;
	SparseGControlInfo **rowIndex;
	SparseGControlInfo **colIndex;
};

typedef struct SparseGVarientTable SparseGVarientTable;


SparseGVarientTable *createSparseGVarientTable(const int row,const int col);
void freeSparseGVarientTable(SparseGVarientTable *ptr);

// important: element can NOT be null
// set table->[rowIndex][colIndex] = element
// will copy a new element into the talbe->[rowIndex][colIndex]
// will automatically create the entry if table->[rowIndex][colIndex] is null && element!=null
void insertSparseGVarientTable(SparseGVarientTable *table,GControlInfo *element,const int row,const int col);

const SparseGControlInfo* getSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col);
const gdsl_queue_t getQueueSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col);
gdsl_queue_t getCopyQueueSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col);
//gdsl_queue_t getCopyPidQueueSparseGVarientTable(const SparseGVarientTable *table,const int row,const int col,const int pid);
void freeCopyQueueSparseGVarientTable(gdsl_queue_t queue);


void dumpSparseGVarientTable(FILE *fp,const SparseGVarientTable *table);

#endif
