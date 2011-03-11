#ifndef PARSER_UTIL_H
#define PARSER_UTIL_H

#include <stdio.h>
#include <string.h>
#include "quadelement.h"
#include "quadmatrix.h"




struct LinearNode
{
	double m;
	char alphaName[128];
	char betaName[128];
	char gammaName[128];
	QuadElement *data;

	struct LinearNode *next;
};

typedef struct LinearNode LinearNode;

LinearNode *createLinearNode(const double m, const char *alphaName, const char *betaName, const char *gammaName, const QuadElement *data);
void freeLinearNode(LinearNode *ptr);
void dumpLinearNode(FILE *fp, const LinearNode *ptr);


// implement by stack like data structure
struct LinearNodeList
{
	LinearNode *list;
};

typedef struct LinearNodeList LinearNodeList;

LinearNodeList *createLinearNodeList(void);
void freeLinearNodeList(LinearNodeList *ptr);
// return NULL if not found the match filenames
const QuadElement *getLinearNodeList(const LinearNodeList *ptr, const double m, const char *alphaName, const char *betaName, const char *gammaName);
void insertLinearNodeList(LinearNodeList *ptr, const double m, const char *alphaName, const char *betaName, const char *gammaName, const QuadElement *data);
void dumpLinearNodeList(FILE *fp, const LinearNodeList *ptr);



//====================================================================


struct NonlinearInfo
{
	int vdsListSize;
	int vgsListSize;
	double *vdsList; // y , row
	double *vgsList; // x , col

	QuadMatrix *ivTable;
	QuadElement *cgs;
	QuadElement *cgd;
	
	QuadMatrix *partialIdsVds;
	QuadMatrix *partialIdsVgs;
};

typedef struct NonlinearInfo NonlinearInfo;

NonlinearInfo *createNonlinearInfo(const int vdsListSize, const int vgsListSize, const double *vdsList, const double *vgsList, const QuadMatrix *ivTable, const QuadElement *cgs, const QuadElement *cgd);
void freeNonlinearInfo(NonlinearInfo *ptr);


struct NonlinearNode
{
	char mName[128];
	char alphaName[128];
	char betaName[128];
	char gammaName[128];
	char configName[128];
	NonlinearInfo *data;

	struct NonlinearNode *next;
};

typedef struct NonlinearNode NonlinearNode;

NonlinearNode *createNonlinearNode(const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const char *configName, const NonlinearInfo *data);
void freeNonlinearNode(NonlinearNode *ptr);
void dumpNonlinearNode(FILE *fp, const NonlinearNode *ptr);


struct NonlinearNodeList
{
	NonlinearNode *list;
};

typedef struct NonlinearNodeList NonlinearNodeList;

NonlinearNodeList *createNonlinearNodeList(void);
void freeNonlinearNodeList(NonlinearNodeList *ptr);
// return NULL if not found the match filenames
const NonlinearInfo *getNonlinearNodeList(const NonlinearNodeList *ptr, const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const char *configName); 
void insertNonlinearNodeList(NonlinearNodeList *ptr, const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const char *configName, const NonlinearInfo *data);
void dumpNonlinearNodeList(FILE *fp, const NonlinearNodeList *ptr);


//============================================================




// if type = gm, use paritalIdsVgs
// if type = ro_1, use paritalIdsVds
void setPartialIdsVds(QuadMatrix *partialIdsVds,const QuadMatrix *ivTable,const gsl_matrix* vdsList);
void setPartialIdsVgs(QuadMatrix *partialIdsVgs,const QuadMatrix *ivTable,const gsl_matrix* vgsList);

#endif
