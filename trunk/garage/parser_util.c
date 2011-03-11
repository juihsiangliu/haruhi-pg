#include "parser_util.h"


static int compareLinearNode(const LinearNode *node, const double m,const char *alphaName,const char *betaName, const char *gammaName)
{
	if(node == NULL)
	{
		return 0;
	}
	else
	{
		int ret = 1;
		if(node->m != m ) ret = 0;
		if(strcmp(alphaName,node->alphaName) != 0 ) ret = 0;
		if(strcmp(betaName,node->betaName) != 0 ) ret = 0;
		if(strcmp(gammaName,node->gammaName) != 0 ) ret = 0;

		return ret;
	}
}



static int compareNonlinearNode(const NonlinearNode *node, const char *mName,const char *alphaName,const char *betaName, const char *gammaName, const char *configName)
{
	if(node == NULL)
	{
		return 0;
	}
	else
	{
		int ret = 1;
		if(strcmp(mName,node->mName) != 0 ) ret = 0;
		if(strcmp(alphaName,node->alphaName) != 0 ) ret = 0;
		if(strcmp(betaName,node->betaName) != 0 ) ret = 0;
		if(strcmp(gammaName,node->gammaName) != 0 ) ret = 0;
		if(strcmp(configName,node->configName) != 0 ) ret = 0;

		return ret;
	}
}


//====================================================================


LinearNode *createLinearNode(const double m, const char *alphaName, const char *betaName, const char *gammaName, const QuadElement *data)
{
	LinearNode *ptr = (LinearNode *)malloc(sizeof(LinearNode));
	ptr->m = m;
	strcpy(ptr->alphaName,alphaName);
	strcpy(ptr->betaName,betaName);
	strcpy(ptr->gammaName,gammaName);
	ptr->data = createQuadElement(data->gvNum);
	copyQuadElement(ptr->data,data);
	ptr->next = NULL;
	return ptr;
}


void freeLinearNode(LinearNode *ptr)
{
	freeQuadElement(ptr->data);
	free(ptr);
}



void dumpLinearNode(FILE *fp, const LinearNode *ptr)
{
	fprintf(fp,"m = %g\n",ptr->m);
	fprintf(fp,"alphaName = %s\n",ptr->alphaName);
	fprintf(fp,"betaName = %s\n",ptr->betaName);
	fprintf(fp,"gammaName = %s\n",ptr->gammaName);
//	dumpQuadElement(ptr->data);
}



//=====================================================================


LinearNodeList *createLinearNodeList(void)
{
	LinearNodeList *ptr = (LinearNodeList *)malloc(sizeof(LinearNodeList));
	ptr->list = NULL;
	return ptr;
}


void freeLinearNodeList(LinearNodeList *ptr)
{
	LinearNode *currentLinearNode = ptr->list;
	while(currentLinearNode != NULL)
	{
		LinearNode *delLinearNode = currentLinearNode;
		currentLinearNode = currentLinearNode->next;
		freeLinearNode(delLinearNode);
	}
	free(ptr);
}




// return NULL if not found the match filenames
const QuadElement *getLinearNodeList(const LinearNodeList *ptr, const double m, const char *alphaName, const char *betaName, const char *gammaName)
{
	LinearNode *currentLinearNode = ptr->list;
	int ret;
	while(currentLinearNode != NULL)
	{
		ret = compareLinearNode(currentLinearNode,m,alphaName,betaName,gammaName);
		if(ret==1) return currentLinearNode->data;
		else currentLinearNode = currentLinearNode->next;	
	}
	return NULL;
}



void insertLinearNodeList(LinearNodeList *ptr, const double m, const char *alphaName, const char *betaName, const char *gammaName, const QuadElement *data)
{
	LinearNode *currentLinearNode = ptr->list;
	LinearNode *target = createLinearNode(m,alphaName,betaName,gammaName,data);
	target->next = ptr->list;
	ptr->list = target;
}



void dumpLinearNodeList(FILE *fp, const LinearNodeList *ptr)
{
	fprintf(fp,"=== LinearNodeList ===\n");
	LinearNode *currentLinearNode = ptr->list;
	while(currentLinearNode != NULL)
	{
		dumpLinearNode(fp,currentLinearNode);
		currentLinearNode = currentLinearNode->next;	
	}
}


//====================================================================



NonlinearInfo *createNonlinearInfo(const int vdsListSize, const int vgsListSize, const double *vdsList, const double *vgsList, const QuadMatrix *ivTable, const QuadElement *cgs, const QuadElement *cgd)
{
	NonlinearInfo *ptr = (NonlinearInfo *)malloc(sizeof(NonlinearInfo));
	ptr->vdsListSize = vdsListSize;
	ptr->vgsListSize = vgsListSize;

	ptr->vdsList = getMempoolSet(sizeof(double)*vdsListSize);
	ptr->vgsList = getMempoolSet(sizeof(double)*vgsListSize);

	memcpy(ptr->vdsList,vdsList,vdsListSize*sizeof(double));
	memcpy(ptr->vgsList,vgsList,vgsListSize*sizeof(double));
	ptr->ivTable = createQuadMatrix(vdsListSize,vgsListSize,ivTable->gvNum);
	copyQuadMatrix(ptr->ivTable,ivTable);
	ptr->cgs = createQuadElement(cgs->gvNum);
	ptr->cgd = createQuadElement(cgd->gvNum);
	copyQuadElement(ptr->cgs,cgs);
	copyQuadElement(ptr->cgd,cgd);

	ptr->partialIdsVds = createQuadMatrix(vdsListSize,vgsListSize,ivTable->gvNum);
	ptr->partialIdsVgs = createQuadMatrix(vdsListSize,vgsListSize,ivTable->gvNum);

	gsl_matrix *tempVdsList = gsl_matrix_calloc(1,ptr->vdsListSize);
	gsl_matrix *tempVgsList = gsl_matrix_calloc(1,ptr->vgsListSize);
	double_to_gsl_matrix(tempVdsList,ptr->vdsList);
	double_to_gsl_matrix(tempVgsList,ptr->vgsList);

	setPartialIdsVds(ptr->partialIdsVds,ptr->ivTable,tempVdsList);	
	setPartialIdsVgs(ptr->partialIdsVgs,ptr->ivTable,tempVgsList);	

	gsl_matrix_free(tempVdsList);
	gsl_matrix_free(tempVgsList);

	return ptr;
}


void freeNonlinearInfo(NonlinearInfo *ptr)
{
	retMempoolSet(ptr->vdsList,sizeof(double)*ptr->vdsListSize);
	retMempoolSet(ptr->vgsList,sizeof(double)*ptr->vgsListSize);
	freeQuadMatrix(ptr->ivTable);
	freeQuadMatrix(ptr->partialIdsVds);
	freeQuadMatrix(ptr->partialIdsVgs);
	freeQuadElement(ptr->cgs);
	freeQuadElement(ptr->cgd);
	free(ptr);
}






NonlinearNode *createNonlinearNode(const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const char *configName, const NonlinearInfo *data)
{
	NonlinearNode *ptr = (NonlinearNode *)malloc(sizeof(NonlinearNode));
	strcpy(ptr->mName,mName);
	strcpy(ptr->alphaName,alphaName);
	strcpy(ptr->betaName,betaName);
	strcpy(ptr->gammaName,gammaName);
	strcpy(ptr->configName,configName);
	ptr->data = createNonlinearInfo(data->vdsListSize,data->vgsListSize,data->vdsList,data->vgsList,data->ivTable,data->cgs,data->cgd);
	ptr->next = NULL;

	return ptr;
}



void freeNonlinearNode(NonlinearNode *ptr)
{
	freeNonlinearInfo(ptr->data);
	free(ptr);
}



void dumpNonlinearNode(FILE *fp, const NonlinearNode *ptr)
{
	fprintf(fp,"mName = %s\n",ptr->mName);
	fprintf(fp,"alphaName = %s\n",ptr->alphaName);
	fprintf(fp,"betaName = %s\n",ptr->betaName);
	fprintf(fp,"gammaName = %s\n",ptr->gammaName);
	fprintf(fp,"configName = %s\n",ptr->configName);
}




NonlinearNodeList *createNonlinearNodeList(void)
{
	NonlinearNodeList *ptr = (NonlinearNodeList *)malloc(sizeof(NonlinearNodeList));
	ptr->list = NULL;
	return ptr;
}


void freeNonlinearNodeList(NonlinearNodeList *ptr)
{
	NonlinearNode *currentNonlinearNode = ptr->list;
	while(currentNonlinearNode != NULL)
	{
		NonlinearNode *delNonlinearNode = currentNonlinearNode;
		currentNonlinearNode = currentNonlinearNode->next;
		freeNonlinearNode(delNonlinearNode);
	}
	free(ptr);
}



// return NULL if not found the match filenames
const NonlinearInfo *getNonlinearNodeList(const NonlinearNodeList *ptr, const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const char *configName)
{
	NonlinearNode *currentNonlinearNode = ptr->list;
	int ret;
	while(currentNonlinearNode != NULL)
	{
		ret = compareNonlinearNode(currentNonlinearNode,mName,alphaName,betaName,gammaName,configName);
		if(ret==1) return currentNonlinearNode->data;
		else currentNonlinearNode = currentNonlinearNode->next;
	}
	return NULL;
}



void insertNonlinearNodeList(NonlinearNodeList *ptr, const char *mName, const char *alphaName, const char *betaName, const char *gammaName, const char *configName, const NonlinearInfo *data)
{
//	NonlinearNode *currentNonlinearNode = ptr->list;
	NonlinearNode *target = createNonlinearNode(mName,alphaName,betaName,gammaName,configName,data);
	target->next = ptr->list;
	ptr->list = target;
}



void dumpNonlinearNodeList(FILE *fp, const NonlinearNodeList *ptr)
{
	fprintf(fp,"=== NonlinearNodeList ===\n");
	NonlinearNode *currentNonlinearNode = ptr->list;
	while(currentNonlinearNode != NULL)
	{
		dumpNonlinearNode(fp,currentNonlinearNode);
		currentNonlinearNode = currentNonlinearNode->next;
	}
}



//============================================================




void setPartialIdsVds(QuadMatrix *partialIdsVds,const QuadMatrix *ivTable,const gsl_matrix* vdsList)
{
	const int row = ivTable->row;
	const int col = ivTable->col;
	int i,j;
	for(i=0;i<row-1;i++)
	{
		const double unit = gsl_matrix_get(vdsList,0,i+1) - gsl_matrix_get(vdsList,0,i); 
		for(j=0;j<col-1;j++)
		{
			const QuadElement *y1 = getPtrEntryQuadMatrix(ivTable,i,j);
			const QuadElement *y2 = getPtrEntryQuadMatrix(ivTable,i,j+1);
			QuadElement *result = createQuadElement(y1->gvNum);
			subQuadElement(result,y2,y1);
			scaleQuadElement(result,(double)1/unit,result);
			setQuadMatrix(partialIdsVds,result,i,j);
			freeQuadElement(result);
		}
	}
}



void setPartialIdsVgs(QuadMatrix *partialIdsVgs,const QuadMatrix *ivTable,const gsl_matrix* vgsList)
{
	const int row = ivTable->row;
	const int col = ivTable->col;
	int i,j;
	for(i=0;i<row-1;i++)
	{
		const double unit = gsl_matrix_get(vgsList,0,i+1) - gsl_matrix_get(vgsList,0,i);
		for(j=0;j<col-1;j++)
		{
			const QuadElement *x1 = getPtrEntryQuadMatrix(ivTable,i,j);
			const QuadElement *x2 = getPtrEntryQuadMatrix(ivTable,i+1,j);
			QuadElement* result = createQuadElement(x1->gvNum);
			subQuadElement(result,x2,x1);
			scaleQuadElement(result,(double)1/unit,result);
			setQuadMatrix(partialIdsVgs,result,i,j);
			freeQuadElement(result);
		}
	}
}
