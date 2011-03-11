#include "montesample.h"





MonteSample *createMonteSample(const int gvNum)
{
	MonteSample *ptr = (MonteSample *)malloc(sizeof(MonteSample));
	ptr->localVarSample = 0.0;
	ptr->globalVarSample = getMempoolSet(sizeof(double)*gvNum);
	return ptr;
}


void freeMonteSample(MonteSample *ptr,const int gvNum)
{
	retMempoolSet(ptr->globalVarSample,sizeof(double)*gvNum);
	free(ptr);
}



void setMonteSample(MonteSample *sample,const gsl_matrix *weight,const gsl_rng* r)
{
	// init
	const int gvNum = gsl_matrix_row_num(weight);
	int i;

	// set local
	sample->localVarSample =  gsl_ran_gaussian(r,1);

	// set global
	gsl_matrix *uncorrelatedSample = gsl_matrix_alloc(gvNum,1);
	gsl_matrix *globalVarSample = gsl_matrix_alloc(gvNum,1);

	for(i=0;i<gvNum;i++) gsl_matrix_set(uncorrelatedSample,i,0,gsl_ran_gaussian(r,1));

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,weight,uncorrelatedSample,0.0,globalVarSample);

	gsl_matrix_to_double(sample->globalVarSample,globalVarSample);

	// free
	gsl_matrix_free(uncorrelatedSample);
	gsl_matrix_free(globalVarSample);
}


// =================================================================



NonlinearInfoMonte *createAndSetNonlinearInfoMonte(const NonlinearInfo *src, const MonteSample *sample)
{
	int i,j;
	NonlinearInfoMonte *ptr = (NonlinearInfoMonte *)malloc(sizeof(NonlinearInfoMonte));
	
	
	ptr->ivTable = getMempoolSet(sizeof(double)*src->vdsListSize*src->vgsListSize);
	ptr->partialIdsVds = getMempoolSet(sizeof(double)*src->vdsListSize*src->vgsListSize);
	ptr->partialIdsVgs = getMempoolSet(sizeof(double)*src->vdsListSize*src->vgsListSize);

	ptr->quadAddress = src;
	ptr->vdsListSize = src->vdsListSize;
	ptr->vgsListSize = src->vgsListSize;
	ptr->cgs = quadElement2Sample(src->cgs,sample); 
	ptr->cgd = quadElement2Sample(src->cgd,sample);
	ptr->quadPartialIdsVds = src->partialIdsVds;
	ptr->quadPartialIdsVgs = src->partialIdsVgs;
	ptr->next = NULL;

	ptr->vdsList = src->vdsList;
	ptr->vgsList = src->vgsList;

//	memcpy(ptr->vdsList,src->vdsList,sizeof(double)*ptr->vdsListSize);
//	memcpy(ptr->vgsList,src->vgsList,sizeof(double)*ptr->vgsListSize);
	
	for(i=0;i<ptr->vdsListSize;i++)
	{
		for(j=0;j<ptr->vgsListSize;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(src->ivTable,i,j);
			const double val = quadElement2Sample(element,sample);
			setMyMatrix(ptr->ivTable,val,ptr->vdsListSize,ptr->vgsListSize,i,j);
		}
	}

	for(i=0;i<ptr->vdsListSize;i++)
	{
		for(j=0;j<ptr->vgsListSize;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(src->partialIdsVds,i,j);
			const double val = quadElement2Sample(element,sample);
			setMyMatrix(ptr->partialIdsVds,val,ptr->vdsListSize,ptr->vgsListSize,i,j);
		}
	}

	for(i=0;i<ptr->vdsListSize;i++)
	{
		for(j=0;j<ptr->vgsListSize;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(src->partialIdsVgs,i,j);
			const double val = quadElement2Sample(element,sample);
			setMyMatrix(ptr->partialIdsVgs,val,ptr->vdsListSize,ptr->vgsListSize,i,j);
		}
	}
	
	return ptr;
}


void freeNonlinearInfoMonte(NonlinearInfoMonte *ptr)
{
	if(ptr!=NULL)
	{
		retMempoolSet(ptr->ivTable,sizeof(double)*ptr->vdsListSize*ptr->vgsListSize);
		retMempoolSet(ptr->partialIdsVds,sizeof(double)*ptr->vdsListSize*ptr->vgsListSize);
		retMempoolSet(ptr->partialIdsVgs,sizeof(double)*ptr->vdsListSize*ptr->vgsListSize);
		free(ptr);
	}
}






// =================================================================


NonlinearMonteList * createAndSetNonlinearMonteList(const NonlinearNodeList *src, const MonteSample *sample)
{
	NonlinearMonteList *ptr = (NonlinearMonteList *)malloc(sizeof(NonlinearMonteList));
	ptr->list = NULL;
	
	NonlinearInfoMonte *insert = NULL;
	const NonlinearNode *currentNonlinearNode = src->list;
	while(currentNonlinearNode!=NULL)
	{
		const NonlinearInfo *data = currentNonlinearNode->data;
		insert = createAndSetNonlinearInfoMonte(data,sample);
		insert->next = ptr->list;
		ptr->list = insert;
		currentNonlinearNode = currentNonlinearNode->next;
	}

	return ptr;
}


void freeNonlinearMonteList(NonlinearMonteList *ptr)
{
	if(ptr!=NULL)
	{
		NonlinearInfoMonte *next;
		NonlinearInfoMonte *del = ptr->list;
		while(del!=NULL)
		{
			next = del->next;
			freeNonlinearInfoMonte(del);
			del = next;
		}
		free(ptr);
	}
}



// ========================================================




// --------------  misc -------------------


void getGVWeight(gsl_matrix *weight, const gsl_matrix *covariance)
{
	int i,j;
	double val;
	const int row = gsl_matrix_row_num(covariance); 
	const int col = gsl_matrix_col_num(covariance); 

	gsl_matrix *s = gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(s,covariance);
	gsl_linalg_cholesky_decomp(s);
	gsl_matrix_set_zero(weight);
	for(i=0;i<row;i++)
	{
		for(j=0;j<=i;j++)
		{
			val = gsl_matrix_get(s,i,j);
			gsl_matrix_set(weight,i,j,val);
		}
	}
	gsl_matrix_free(s);
}



gsl_rng* createRNG(void)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r,time(0));
	return r;
}



// given a QuadElement and a MonteSample, return the double val
double quadElement2Sample(const QuadElement *quad, const MonteSample *monteSample)
{
	if(quad == NULL)
	{
		return 0.0;
	}
	else
	{
		int i,j;
		const int gvNum = quad->gvNum;

		double ret = 0.0;
	
		ret = quad->m;
	
		ret += quad->alpha * monteSample->localVarSample;
	
		for(i=0;i<gvNum;i++)
		{
			ret += getMyMatrix(quad->beta,gvNum,1,i,0) * getMyMatrix(monteSample->globalVarSample,gvNum,1,i,0);
		}
	
		double gv1,gv2,gamma;
		for(i=0;i<gvNum;i++)
		{
			gv1 = getMyMatrix(monteSample->globalVarSample,gvNum,1,i,0);
			for(j=0;j<gvNum;j++)
			{
				gv2 = getMyMatrix(monteSample->globalVarSample,gvNum,1,j,0);
				gamma = getMyMatrix(quad->gamma,gvNum,gvNum,i,j);
				ret += gv1*gv2*gamma;
			}
		}

		return ret;
	}
}







/*

static gdsl_element_t allocGControlInfo2D(void *ptr)
{
	MonteGControlInfo *src = (MonteGControlInfo *)ptr;
	MonteGControlInfo *dest = (MonteGControlInfo *)malloc(sizeof(MonteGControlInfo));
	dest->sign = src->sign;
	dest->type = src->type;
	dest->gate = src->gate;
	dest->drain = src->drain;
	dest->source = src->source;
	const int vdsListCol = gsl_matrix_col_num(src->vdsList);
	dest->vdsList = gsl_matrix_alloc(1,vdsListCol);
	gsl_matrix_memcpy(dest->vdsList,src->vdsList);
	const int vgsListCol = gsl_matrix_col_num(src->vgsList);
	dest->vgsList = gsl_matrix_alloc(1,vgsListCol);
	gsl_matrix_memcpy(dest->vgsList,src->vgsList);
	dest->ivTable = gsl_matrix_alloc(vdsListCol,vgsListCol);
	gsl_matrix_memcpy(dest->ivTable,src->ivTable);
	dest->partialIdsVds = gsl_matrix_alloc(vdsListCol,vgsListCol);
	gsl_matrix_memcpy(dest->partialIdsVds,src->partialIdsVds);
	dest->partialIdsVgs = gsl_matrix_alloc(vdsListCol,vgsListCol);
	gsl_matrix_memcpy(dest->partialIdsVgs,src->partialIdsVgs);
	return dest;
}


static void freeGControlInfo2D(void *ptr)
{
	if(ptr!=NULL)
	{
		MonteGControlInfo *src = (MonteGControlInfo *)ptr;
		gsl_matrix_free(src->vdsList);
		gsl_matrix_free(src->vgsList);
		gsl_matrix_free(src->ivTable);
		gsl_matrix_free(src->partialIdsVds);
		gsl_matrix_free(src->partialIdsVgs);
		free(src);
	}
}




static gdsl_element_t allocQueue2D(void *ptr)
{

	gdsl_queue_t src = (gdsl_queue_t) ptr;
	gdsl_queue_t dest = gdsl_queue_alloc("dest",allocGControlInfo2D,freeGControlInfo2D);
	gdsl_queue_t tmp = gdsl_queue_alloc("tmp",allocGControlInfo2D,freeGControlInfo2D);
	gdsl_element_t element;
	while(!gdsl_queue_is_empty(src))
	{
		element = gdsl_queue_remove(src);
		gdsl_queue_insert(dest,element);
		gdsl_queue_insert(tmp,element);
	}

	while(!gdsl_queue_is_empty(tmp))
	{
		element = gdsl_queue_remove(tmp);
		gdsl_queue_insert(src,element);
	}
	gdsl_queue_free(tmp);
	return dest;
}



static void freeQueue2D(void *ptr)
{
	if(ptr!=NULL)
	{
		gdsl_queue_t queue = (gdsl_queue_t) ptr;
		gdsl_queue_free(queue);
	}
}



static double quad2monte(const QuadElement *src,const double localVarSample, const gsl_matrix *globalVarSample)
{
	if(isEmptyQuadElement(src)) 
	{
		return 0;
	}
	else
	{
		double ret = 0.0;
		gsl_matrix *betaTerm = gsl_matrix_alloc(1,1);
		gsl_matrix *gammaTerm1 = gsl_matrix_alloc(1,src->gvNum);
		gsl_matrix *gammaTerm = gsl_matrix_alloc(1,1);

		ret += src->m;
	
		ret += src->alpha * localVarSample;
	
		gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0,src->beta,globalVarSample,0.0,betaTerm);
		ret += gsl_matrix_get(betaTerm,0,0);

		gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0,globalVarSample,src->gamma,0.0,gammaTerm1);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0,gammaTerm1,globalVarSample,0.0,gammaTerm);
		ret += gsl_matrix_get(gammaTerm,0,0);

		gsl_matrix_free(betaTerm);
		gsl_matrix_free(gammaTerm1);
		gsl_matrix_free(gammaTerm);
		return ret;
	}
}



static MonteGControlInfo * createMonteGControlInfo(const GControlInfo* src,const double localVarSample,const gsl_matrix *globalVarSample)
{
	int i,j;
	double val;
	const QuadElement *currentQuad;
	MonteGControlInfo *dest = (MonteGControlInfo *)malloc(sizeof(MonteGControlInfo));
	dest->sign = src->sign;
	dest->type = src->type;
	dest->gate = src->gate;
	dest->drain = src->drain;
	dest->source = src->source;
	const int vdsListCol = gsl_matrix_col_num(src->vdsList);
	dest->vdsList = gsl_matrix_alloc(1,vdsListCol);
	gsl_matrix_memcpy(dest->vdsList,src->vdsList);
	const int vgsListCol = gsl_matrix_col_num(src->vgsList);
	dest->vgsList = gsl_matrix_alloc(1,vgsListCol);
	gsl_matrix_memcpy(dest->vgsList,src->vgsList);
	dest->ivTable = gsl_matrix_alloc(vdsListCol,vgsListCol);
	dest->partialIdsVds = gsl_matrix_alloc(vdsListCol,vgsListCol);
	dest->partialIdsVgs = gsl_matrix_alloc(vdsListCol,vgsListCol);

	for(i=0;i<vdsListCol;i++)
	{
		for(j=0;j<vgsListCol;j++)
		{
			currentQuad = getPtrEntryQuadMatrix(src->ivTable,i,j);	
			val = quad2monte(currentQuad,localVarSample,globalVarSample);
			gsl_matrix_set(dest->ivTable,i,j,val);

			currentQuad = getPtrEntryQuadMatrix(src->partialIdsVds,i,j);
			val = quad2monte(currentQuad,localVarSample,globalVarSample);
			gsl_matrix_set(dest->partialIdsVds,i,j,val);
			
			currentQuad = getPtrEntryQuadMatrix(src->partialIdsVgs,i,j);
			val = quad2monte(currentQuad,localVarSample,globalVarSample);
			gsl_matrix_set(dest->partialIdsVgs,i,j,val);
		}
	}
	return dest;
}



static void freeMonteGControlInfo(MonteGControlInfo *ptr)
{
	if(ptr == NULL) return;

	gsl_matrix_free(ptr->vdsList);
	gsl_matrix_free(ptr->vgsList);
	gsl_matrix_free(ptr->ivTable);
	gsl_matrix_free(ptr->partialIdsVds);
	gsl_matrix_free(ptr->partialIdsVgs);
	free(ptr);
}




//===================================================================================
//===================================================================================
//===================================================================================



void getGVWeight(gsl_matrix *w, const NetlistStampResultQuad *netlist)
{
	int i,j;
	double val;
	const int gvNum = netlist->gvNum;
	gsl_matrix *s = gsl_matrix_alloc(gvNum,gvNum);
	gsl_matrix_memcpy(s,netlist->s);
	gsl_linalg_cholesky_decomp(s);
	gsl_matrix_set_zero(w);
	for(i=0;i<gvNum;i++)
	{
		for(j=0;j<=i;j++)
		{
			val = gsl_matrix_get(s,i,j);
			gsl_matrix_set(w,i,j,val);
		}
	}
	gsl_matrix_free(s);
}



gsl_rng* createRNG(void)
{
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r,time(0));
	return r;
}


void generateVarSample(double *localVarSample, gsl_matrix *globalVarSample,const gsl_matrix *w,const gsl_rng* r)
{
	// init
	const int gvNum = gsl_matrix_row_num(w);
	int i;

	// set local
	*localVarSample =  gsl_ran_gaussian(r,1);

	// set global
	gsl_matrix *uncorrelatedSample = gsl_matrix_alloc(gvNum,1);
	for(i=0;i<gvNum;i++) gsl_matrix_set(uncorrelatedSample,i,0,gsl_ran_gaussian(r,1));

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,w,uncorrelatedSample,0.0,globalVarSample);

	// free
	gsl_matrix_free(uncorrelatedSample);
}



MonteNetlist* createMonteNetlist(const NetlistStampResultQuad *netlist,const double localVarSample, const gsl_matrix *globalVarSample)
{
	// init & alloc
	int i,j,k;
	MonteNetlist *ptr = (MonteNetlist *)malloc(sizeof(MonteNetlist));
	ptr->a = gsl_matrix_alloc(netlist->nodeNum,netlist->nodeNum);
	ptr->b = gsl_matrix_alloc(netlist->nodeNum,netlist->nodeNum);
	ptr->c = gsl_matrix_alloc(netlist->nodeNum,netlist->nodeNum);
	ptr->u = gsl_matrix_alloc(netlist->nodeNum,netlist->stepNum);
	ptr->s = gsl_matrix_alloc(netlist->gvNum,netlist->gvNum);
	ptr->gVarientTable = gdsl_2darray_alloc("gVarientTable",netlist->nodeNum,netlist->nodeNum,allocQueue2D,freeQueue2D);
	for(i=0;i<netlist->nodeNum;i++)
	{
		for(j=0;j<netlist->nodeNum;j++)
		{
			gdsl_queue_t monteGControlInfoQueue = gdsl_queue_alloc("monteGControlInfoQueue",allocGControlInfo2D,freeGControlInfo2D);
			gdsl_2darray_set_content(ptr->gVarientTable,i,j,monteGControlInfoQueue);
		}
	}
	// assignment
	gsl_matrix_memcpy(ptr->s,netlist->s);
	ptr->nodeNum = netlist->nodeNum;
	ptr->gvNum = netlist->gvNum;
	ptr->endTime = netlist->endTime;
	ptr->deltaT = netlist->deltaT;
	ptr->stepNum = netlist->stepNum;
	for(i=0;i<netlist->nodeNum;i++)
	{
		for(j=0;j<netlist->nodeNum;j++)
		{
			// a, b, c, then gVarientTable
			const double a = quad2monte(getPtrEntryQuadMatrix(netlist->a,i,j),localVarSample,globalVarSample);
			gsl_matrix_set(ptr->a,i,j,a);	
			const double b = quad2monte(getPtrEntryQuadMatrix(netlist->b,i,j),localVarSample,globalVarSample);
			gsl_matrix_set(ptr->b,i,j,b);	
			const double c = quad2monte(getPtrEntryQuadMatrix(netlist->c,i,j),localVarSample,globalVarSample);
			gsl_matrix_set(ptr->c,i,j,c);	

			const gdsl_queue_t srcQueue = gdsl_2darray_get_content(netlist->gVarientTable,i,j);
			const int queueSize = gdsl_queue_get_size(srcQueue);
			gdsl_queue_t destQueue = gdsl_2darray_get_content(ptr->gVarientTable,i,j);			
			for(k=1;k<=queueSize;k++)
			{
				const GControlInfo *gControlInfo = gdsl_queue_search_by_position(srcQueue,k);
				MonteGControlInfo * monteGControlInfo = createMonteGControlInfo(gControlInfo,localVarSample,globalVarSample);
				gdsl_queue_insert(destQueue,monteGControlInfo);
				freeMonteGControlInfo(monteGControlInfo);
			}
		}

		for(j=0;j<netlist->stepNum;j++)
		{
			const double u = quad2monte(getPtrEntryQuadMatrix(netlist->u,i,j),localVarSample,globalVarSample);
			gsl_matrix_set(ptr->u,i,j,u);	
		}
	}

	return ptr;
}


void freeMonteNetlist(MonteNetlist *ptr)
{
	gsl_matrix_free(ptr->a);
	gsl_matrix_free(ptr->b);
	gsl_matrix_free(ptr->c);
	gsl_matrix_free(ptr->u);
	gsl_matrix_free(ptr->s);
	gdsl_2darray_free(ptr->gVarientTable);
	free(ptr);
}


*/
