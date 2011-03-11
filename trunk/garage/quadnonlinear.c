#include "quadnonlinear.h"



static pthread_mutex_t mutexSetJacobian = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mutexAll = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mutexG = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mutexSetGVar = PTHREAD_MUTEX_INITIALIZER;



inline static int toZeroBase(const int x)
{
	return x-1;
}


static int setGVarientForQueue(const gdsl_element_t E, gdsl_location_t LOCATION, void *USER_DATA) 
{
	//   =====   initial   =====
	GVarientIJQuad* currentGVarientIJ = (GVarientIJQuad *) USER_DATA;
	const int pid = currentGVarientIJ->pid;
	const GControlInfo *currentGControlInfo = (GControlInfo *) E;
	const int gateIndex = toZeroBase(currentGControlInfo->gate);
	const int drainIndex = toZeroBase(currentGControlInfo->drain);
	const int sourceIndex = toZeroBase(currentGControlInfo->source);
	QuadElement *vgs = createPidQuadElement(currentGVarientIJ->gvNum,pid);
	QuadElement *vds = createPidQuadElement(currentGVarientIJ->gvNum,pid);
	const QuadElement *vDrain = getPtrEntryQuadMatrix(currentGVarientIJ->v,drainIndex,0);
	const QuadElement *vGate = getPtrEntryQuadMatrix(currentGVarientIJ->v,gateIndex,0);
	if(sourceIndex !=-1)
	{
		const QuadElement *vSource = getPtrEntryQuadMatrix(currentGVarientIJ->v,sourceIndex,0);
		subQuadElement(vgs,vGate,vSource);
		subQuadElement(vds,vDrain,vSource);
	}
	else
	{
		copyQuadElement(vgs,vGate);
		copyQuadElement(vds,vDrain);
	}

//	dumpQuadMatrix(currentGVarientIJ->v);
	//   =====   prepare to curve fitting   =====
	// x = vgs, y = vds
	int upXIndex,downXIndex,upYIndex,downYIndex = 0;
	const double meanVgs = meanPidQuadElement(vgs,currentGVarientIJ->s,pid);
	const double meanVds = meanPidQuadElement(vds,currentGVarientIJ->s,pid);
	const double stdVgs = sqrt(varPidQuadElement(vgs,currentGVarientIJ->s,pid));
	const double stdVds = sqrt(varPidQuadElement(vds,currentGVarientIJ->s,pid));
	const double fixRatio = 3.0;  // use to fit +- 3 sigma
	const double deltaXAxis = getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,1) - getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,0);
	const double deltaYAxis = getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,1) - getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,0);

	// first, use max to find out the minimum required rowIndexUp
	// second, use min to guarantee the boundary
	// (the rowIndexDown may up-side-down)
	// set rowIndex
	const int rowIndexTemp = GSL_MAX_INT( (int) ((meanVds - getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,0)) / deltaYAxis) , 1);
	const int rowIndex = GSL_MIN_INT(rowIndexTemp,currentGControlInfo->vdsListSize-3);
	const int rowIndexUpTemp = GSL_MAX_INT(rowIndex+1,(int)((meanVds+fixRatio*stdVds-getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,0))/deltaYAxis)  );
	const int rowIndexUp = GSL_MIN_INT(rowIndexUpTemp,currentGControlInfo->vdsListSize-2); // the final entry of partial table may be null (lazy to compute ....)
	const int rowIndexDownTemp = GSL_MIN_INT(rowIndex-1,(int)((meanVds-fixRatio*stdVds-getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,0))/deltaYAxis)  );
	const int rowIndexDown = GSL_MAX_INT(rowIndexDownTemp,0);
	// set colIndex	
	const int colIndexTemp = GSL_MAX_INT((int) ((meanVgs - getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vdsListSize,0,0)) / deltaXAxis) , 1);
	const int colIndex = GSL_MIN_INT(colIndexTemp,currentGControlInfo->vgsListSize-3 );
	const int colIndexUpTemp = GSL_MAX_INT(colIndex+1,(int)((meanVgs+fixRatio*stdVgs-getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,0))/deltaXAxis)  );
	const int colIndexUp = GSL_MIN_INT(colIndexUpTemp,currentGControlInfo->vgsListSize-2); // the final entry of partial table may be null (lazy to compute ....)
	const int colIndexDownTemp = GSL_MIN_INT(colIndex-1,(int)((meanVgs-fixRatio*stdVgs-getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,0))/deltaXAxis)  );
	const int colIndexDown = GSL_MAX_INT(colIndexDownTemp,0);
//	fprintf(stderr,"before fit\n");
//	fprintf(stderr,"rowIndextTemp=%d meanVds=%lf\n",rowIndexTemp,meanVds);
	
	//   =====   actual fit   =====
	double Z[4][3];
	int i,j;
	for(i=0;i<4;i++) Z[i][0] = 1;
	Z[0][1] = getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndexDown);
	Z[1][1] = getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndexUp);
	Z[2][1] = Z[3][1] = getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndex);
	Z[0][2] = Z[1][2] = getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,colIndex);
	Z[2][2] = getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,colIndexDown);
	Z[3][2] = getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,colIndexUp);
	
//	fprintf(stderr,"after set z\n");
	QuadMatrix *ZQuad = createPidQuadMatrix(4,3,currentGVarientIJ->gvNum,pid);
	QuadElement *tempElement = createPidQuadElement(currentGVarientIJ->gvNum,pid);
	for(i=0;i<4;i++)
	{
		for(j=0;j<3;j++)
		{
			tempElement->m = Z[i][j];
			setPidQuadMatrix(ZQuad,tempElement,i,j,pid);
			//printf("%lf ",Z[i][j]);
		}
		//printf("\n");
	}
//	fprintf(stderr,"after z\n");
	
	freePidQuadElement(tempElement,pid);
// -------------------------------------------------------
	
	QuadMatrix *ZQuadTrans = createPidQuadMatrix(3,4,currentGVarientIJ->gvNum,pid);
	transposePidQuadMatrix(ZQuadTrans,ZQuad,pid);
	QuadMatrix *ZTransZQuad = createPidQuadMatrix(3,3,currentGVarientIJ->gvNum,pid);
	mulPidQuadMatrix(ZTransZQuad,ZQuadTrans,ZQuad,pid);
	QuadMatrix *yQuad = createPidQuadMatrix(4,1,currentGVarientIJ->gvNum,pid);
		
	const QuadElement *yDownIndex = getPtrEntryQuadMatrix(currentGControlInfo->partialIdsVxs,rowIndexDown,colIndex);
	const QuadElement *yUpIndex = getPtrEntryQuadMatrix(currentGControlInfo->partialIdsVxs,rowIndexUp,colIndex);
	const QuadElement *yIndexDown = getPtrEntryQuadMatrix(currentGControlInfo->partialIdsVxs,rowIndex,colIndexDown);
	const QuadElement *yIndexUp = getPtrEntryQuadMatrix(currentGControlInfo->partialIdsVxs,rowIndex,colIndexUp);
	setPidQuadMatrix(yQuad,yDownIndex,0,0,pid);
	setPidQuadMatrix(yQuad,yUpIndex,1,0,pid);
	setPidQuadMatrix(yQuad,yIndexDown,2,0,pid);
	setPidQuadMatrix(yQuad,yIndexUp,3,0,pid);

	QuadMatrix *ZTransyQuad = createPidQuadMatrix(3,1,currentGVarientIJ->gvNum,pid);
	mulPidQuadMatrix(ZTransyQuad,ZQuadTrans,yQuad,pid);

//================================================================================================
//	pthread_mutex_lock(&mutexSetGVar);
	QuadMatrix *bQuad = createPidQuadMatrix(3,1,currentGVarientIJ->gvNum,pid);
	solvePidQuadMatrix(bQuad,ZTransZQuad,ZTransyQuad,pid);
//	pthread_mutex_unlock(&mutexSetGVar);
//================================================================================================
	freePidQuadMatrix(ZQuad,pid);
	freePidQuadMatrix(ZQuadTrans,pid);
	freePidQuadMatrix(ZTransZQuad,pid);
	freePidQuadMatrix(yQuad,pid);
	freePidQuadMatrix(ZTransyQuad,pid);
// -------------------------------------------------------


	//   =====   update gVarientElement according to the fitting result   =====
	QuadElement *tempQuadElement = createPidQuadElement(currentGVarientIJ->gvNum,pid);
	QuadElement *tempQuadElement2 = createPidQuadElement(currentGVarientIJ->gvNum,pid);
	QuadElement *tempQuadElement3 = createPidQuadElement(currentGVarientIJ->gvNum,pid);
	
	addQuadElement(tempQuadElement,tempQuadElement,getPtrEntryQuadMatrix(bQuad,0,0));
	mulPidQuadElement(tempQuadElement2,getPtrEntryQuadMatrix(bQuad,1,0),vds,pid);
	mulPidQuadElement(tempQuadElement3,getPtrEntryQuadMatrix(bQuad,2,0),vgs,pid);
	addQuadElement(tempQuadElement,tempQuadElement,tempQuadElement2);
	addQuadElement(tempQuadElement,tempQuadElement,tempQuadElement3);
	scaleQuadElement(tempQuadElement,currentGControlInfo->sign,tempQuadElement);
	addQuadElement(currentGVarientIJ->gVarientElement,currentGVarientIJ->gVarientElement,tempQuadElement);

	freePidQuadElement(tempQuadElement,pid);
	freePidQuadElement(tempQuadElement2,pid);
	freePidQuadElement(tempQuadElement3,pid);

	freePidQuadElement(vgs,pid);
	freePidQuadElement(vds,pid);

	freePidQuadMatrix(bQuad,pid);
}



static void updateGVarientMatrixFinite(SparseQuadMatrix *gVarient,const QuadMatrix* v,const SparseNetlistQuad *netlist,const int pid)
{
	GVarientIJQuad *currentIJ = getPidMempoolSet(sizeof(GVarientIJQuad),pid);
	currentIJ->gVarientElement = createPidQuadElement(netlist->gvNum,pid);
	currentIJ->gvNum = netlist->gvNum;
	currentIJ->v = v;
	currentIJ->s = netlist->s;	
	currentIJ->pid = pid;

//	pthread_mutex_lock(&mutexSetGVar);
	clearPidSparseQuadMatrix(gVarient,pid);

	int i,j;
	for(i=0;i<netlist->nodeNum;i++)
	{
		for(j=0;j<netlist->nodeNum;j++)
		{
//			gdsl_queue_t gControlInfoQueue = getCopyQueueSparseGVarientTable(netlist->gVarientTable,i,j);
			gdsl_queue_t gControlInfoQueue = getQueueSparseGVarientTable(netlist->gVarientTable,i,j);
			if(gControlInfoQueue == NULL) continue;
			else
			{
				resetQuadElement(currentIJ->gVarientElement);
				// use map_forward to calculate
				gdsl_queue_map_forward(gControlInfoQueue,setGVarientForQueue,currentIJ);
				setPidSparseQuadMatrix(gVarient,currentIJ->gVarientElement,i,j,pid);
			}
//			freeCopyQueueSparseGVarientTable(gControlInfoQueue);
		}
	}
//	pthread_mutex_unlock(&mutexSetGVar);

	freePidQuadElement(currentIJ->gVarientElement,pid);
	retPidMempoolSet(currentIJ,sizeof(GVarientIJQuad),pid);
}




static void getJacobianFinite(SparseQuadMatrix *jacobian, const SparseNetlistQuad *netlist, const QuadMatrix *v,const int threadNum)
{
		const double DV = 0.006;
		SparseQuadMatrix *z = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
		copySparseQuadMatrix(z,netlist->c);
		scaleSparseQuadMatrix(z,1.0/netlist->deltaT,z);
		int i,j;
		// setup delta V
		QuadElement *deltaVElement = createQuadElement(netlist->gvNum);
		setQuadElement(deltaVElement,DV,0,NULL,0);
		QuadMatrix *deltaV = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
		// calculate by finite difference
		QuadMatrix *currentV = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
		SparseQuadMatrix *zSubACurrent = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
		QuadMatrix *fUp = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
		QuadMatrix *fDown = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
		QuadMatrix *deltaF = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
		QuadMatrix *result = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
		SparseQuadMatrix *aCurrent = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
		SparseQuadMatrix *gVarient = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
		// clear jacobian
		clearSparseQuadMatrix(jacobian);
		// calculate column by column
		for(j=0;j<netlist->nodeNum;j++)
		{
			resetQuadMatrix(deltaV);
			setQuadMatrix(deltaV,deltaVElement,j,0);
			
			addQuadMatrix(currentV,v,deltaV);
			updateGVarientMatrixFinite(gVarient,currentV,netlist,0);
			addSparseQuadMatrix(aCurrent,netlist->a,gVarient);
			subSparseQuadMatrix(zSubACurrent,z,aCurrent);
			mulVecSparseQuadMatrix(fUp,zSubACurrent,currentV);

			subQuadMatrix(currentV,v,deltaV);
			updateGVarientMatrixFinite(gVarient,currentV,netlist,0);
			addSparseQuadMatrix(aCurrent,netlist->a,gVarient);
			subSparseQuadMatrix(zSubACurrent,z,aCurrent);
			mulVecSparseQuadMatrix(fDown,zSubACurrent,currentV);

			subQuadMatrix(deltaF,fUp,fDown);
			scaleQuadMatrix(result,1.0/(2*DV),deltaF);
			setDenseColQuick2SparseQuadMatrix(jacobian,result,j);
		}
		
		freeSparseQuadMatrix(z);
		freeQuadElement(deltaVElement);
		freeQuadMatrix(deltaV);
		freeQuadMatrix(currentV);
		freeSparseQuadMatrix(zSubACurrent);
		freeQuadMatrix(fUp);
		freeQuadMatrix(fDown);
		freeQuadMatrix(deltaF);
		freeQuadMatrix(result);
		freeSparseQuadMatrix(aCurrent);
		freeSparseQuadMatrix(gVarient);
}



struct ParJacobian
{	
	// --------   const data ----------
	const SparseNetlistQuad *netlist;
	const QuadMatrix *v;
	int lowerBound;
	int upperBound;
	int pid;

	// ---------  used local variables ----------------
	// in mempool pid
	SparseQuadMatrix *jacobian;
	SparseQuadMatrix *z;
	QuadElement *deltaVElement;
	QuadMatrix *deltaV;
	QuadMatrix *currentV;
	SparseQuadMatrix *zSubACurrent;
	QuadMatrix *fUp;
	QuadMatrix *fDown;
	QuadMatrix *deltaF;
	QuadMatrix *result;
	SparseQuadMatrix *aCurrent;
	SparseQuadMatrix *gVarient; 
};



static struct ParJacobian *createParJacobian(const int nodeNum,const int gvNum,const int pid)
{
	struct ParJacobian *ptr = (struct ParJacobian *)malloc(sizeof(struct ParJacobian));
	ptr->pid = pid;

	ptr->jacobian = createPidSparseQuadMatrix(nodeNum,nodeNum,gvNum,pid);
	ptr->z = createPidSparseQuadMatrix(nodeNum,nodeNum,gvNum,pid);
	ptr->deltaVElement = createPidQuadElement(gvNum,pid);
	ptr->deltaV = createPidQuadMatrix(nodeNum,1,gvNum,pid);
	ptr->currentV = createPidQuadMatrix(nodeNum,1,gvNum,pid);
	ptr->zSubACurrent = createPidSparseQuadMatrix(nodeNum,nodeNum,gvNum,pid);
	ptr->fUp = createPidQuadMatrix(nodeNum,1,gvNum,pid);
	ptr->fDown = createPidQuadMatrix(nodeNum,1,gvNum,pid);
	ptr->deltaF = createPidQuadMatrix(nodeNum,1,gvNum,pid);
	ptr->result = createPidQuadMatrix(nodeNum,1,gvNum,pid);
	ptr->aCurrent = createPidSparseQuadMatrix(nodeNum,nodeNum,gvNum,pid);

	ptr->gVarient = createPidSparseQuadMatrix(nodeNum,nodeNum,gvNum,pid);
	
	return ptr;
}





static void freeParJacobian(struct ParJacobian *ptr,const int pid)
{
	freePidSparseQuadMatrix(ptr->jacobian,pid);
	freePidSparseQuadMatrix(ptr->z,pid);
	freePidQuadElement(ptr->deltaVElement,pid);
	freePidQuadMatrix(ptr->deltaV,pid);
	freePidQuadMatrix(ptr->currentV,pid);
	freePidSparseQuadMatrix(ptr->zSubACurrent,pid);
	freePidQuadMatrix(ptr->fUp,pid);
	freePidQuadMatrix(ptr->fDown,pid);
	freePidQuadMatrix(ptr->deltaF,pid);
	freePidQuadMatrix(ptr->result,pid);
	freePidSparseQuadMatrix(ptr->aCurrent,pid);
	
	freePidSparseQuadMatrix(ptr->gVarient,pid);

	free(ptr);
}


static void *getJacobianFiniteParallelKernel(void *par)
{
	struct ParJacobian *ptr = (struct ParJacobian *)par;
	const double DV = 0.006;
	int i,j;
	const int pid = ptr->pid;

//	clearSparseQuadMatrix(ptr->jacobian);

	//========================================
	scalePidSparseQuadMatrix(ptr->z,1.0/ptr->netlist->deltaT,ptr->netlist->c,pid);
	setQuadElement(ptr->deltaVElement,DV,0,NULL,0);
	//========================================

	// calculate column by column
	for(j=ptr->lowerBound;j<=ptr->upperBound;j++)
	{
//		fprintf(stderr,"j = %d\n",j);
		resetQuadMatrix(ptr->deltaV);
		setPidQuadMatrix(ptr->deltaV,ptr->deltaVElement,j,0,pid);
		
		addPidQuadMatrix(ptr->currentV,ptr->v,ptr->deltaV,pid);
		updateGVarientMatrixFinite(ptr->gVarient,ptr->currentV,ptr->netlist,pid);

		addPidSparseQuadMatrix(ptr->aCurrent,ptr->netlist->a,ptr->gVarient,pid);
		subPidSparseQuadMatrix(ptr->zSubACurrent,ptr->z,ptr->aCurrent,pid);
		mulVecPidSparseQuadMatrix(ptr->fUp,ptr->zSubACurrent,ptr->currentV,pid);

		subPidQuadMatrix(ptr->currentV,ptr->v,ptr->deltaV,pid);
		updateGVarientMatrixFinite(ptr->gVarient,ptr->currentV,ptr->netlist,pid);
		
		addPidSparseQuadMatrix(ptr->aCurrent,ptr->netlist->a,ptr->gVarient,pid);
		subPidSparseQuadMatrix(ptr->zSubACurrent,ptr->z,ptr->aCurrent,pid);
		mulVecPidSparseQuadMatrix(ptr->fDown,ptr->zSubACurrent,ptr->currentV,pid);

		subPidQuadMatrix(ptr->deltaF,ptr->fUp,ptr->fDown,pid);
		scalePidQuadMatrix(ptr->result,1.0/(2*DV),ptr->deltaF,pid);

		setDenseColQuick2SparseQuadMatrixPid(ptr->jacobian,ptr->result,j,pid);
	}
		

	pthread_exit(NULL);
}






static void getJacobianFiniteParallel(SparseQuadMatrix *jacobian, const SparseNetlistQuad *netlist, const QuadMatrix *v,const int threadNum)
{
	time_t t1,t2;
	
	time(&t1);
	clearSparseQuadMatrix(jacobian);
	int i;
	void *status;
	pthread_t threads[threadNum];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct ParJacobian ** par = (struct ParJacobian **)malloc(threadNum*sizeof(struct ParJacobian*));
	for(i=0;i<threadNum;i++)
	{
		par[i] = createParJacobian(netlist->nodeNum,netlist->gvNum,i+1);
		par[i]->netlist = netlist;
		par[i]->v = v;
		par[i]->lowerBound = block_low(i,threadNum,netlist->nodeNum);	
		par[i]->upperBound = block_high(i,threadNum,netlist->nodeNum);
		par[i]->pid = i+1;
	}
	time(&t2);
//	fprintf(stderr,"jacobian init time:%g\n",difftime(t2,t1));

	time(&t1);
	for(i=0;i<threadNum;i++)
	{
		pthread_create(&threads[i], &attr, getJacobianFiniteParallelKernel,par[i]);
	}

	pthread_attr_destroy(&attr);

	for(i=0;i<threadNum;i++)
	{
		pthread_join(threads[i],NULL);
	}
	time(&t2);
	fprintf(stderr,"jacobian compute time:%g\n",difftime(t2,t1));


	time(&t1);
	for(i=0;i<threadNum;i++)
	{
		int j,k;
		for(j=par[i]->lowerBound;j<=par[i]->upperBound;j++)
		{
			const SparseQuadElement *ptr = par[i]->jacobian->colIndex[j]->colLink;
			while(ptr!=NULL)
			{
				const QuadElement *val = ptr->data;
				int k = ptr->row;
				setSparseQuadMatrix(jacobian,val,k,j);
				ptr = ptr->colLink;
			}
		}
	}
	time(&t2);
//	fprintf(stderr,"jacobian set time:%g\n",difftime(t2,t1));

	for(i=0;i<threadNum;i++) freeParJacobian(par[i],i+1);
	free(par);

//	exit(0);
}




static void solveDCSparseQuad(QuadMatrix *vDC, const SparseQuadMatrix *aCurrent, const QuadMatrix *vCurrent, const int nextUIndex, const SparseNetlistQuad *netlist,const int threadNum)
{
	QuadMatrix *uNext = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	SparseQuadMatrix *z = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	QuadMatrix *zv = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *bu = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *zv_bu = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	SparseQuadMatrix *z_sub_a = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);

	getColCopyQuadMatrix(uNext,nextUIndex,netlist->u);
	scaleSparseQuadMatrix(z,1.0/netlist->deltaT,netlist->c);
	mulVecSparseQuadMatrix(zv,z,vCurrent);
	mulVecSparseQuadMatrix(bu,netlist->b,uNext);
	addQuadMatrix(zv_bu,zv,bu);
	subSparseQuadMatrix(z_sub_a,z,aCurrent);
	solveSparseQuadMatrix(vDC,z_sub_a,zv_bu,threadNum);
	
	freeQuadMatrix(uNext);
	freeSparseQuadMatrix(z);
	freeQuadMatrix(zv);
	freeQuadMatrix(bu);
	freeQuadMatrix(zv_bu);
	freeSparseQuadMatrix(z_sub_a);
}



// k = C/dt * vt + B *ut+1
static void getK(QuadMatrix *k, const SparseNetlistQuad *netlist, const QuadMatrix *vCurrent, const int currentIndex)
{
	const int nodeNum = netlist->nodeNum;
	const int gvNum = netlist->gvNum;
	QuadMatrix *uNext = createQuadMatrix(nodeNum,1,gvNum);
	SparseQuadMatrix *z = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	QuadMatrix *zv = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *bu = createQuadMatrix(nodeNum,1,gvNum);
	
	getColCopyQuadMatrix(uNext,currentIndex+1,netlist->u);
	scaleSparseQuadMatrix(z,1.0/netlist->deltaT,netlist->c);
	mulVecSparseQuadMatrix(zv,z,vCurrent);
	mulVecSparseQuadMatrix(bu,netlist->b,uNext);
	addQuadMatrix(k,zv,bu);

	freeQuadMatrix(uNext);
	freeSparseQuadMatrix(z);
	freeQuadMatrix(zv);
	freeQuadMatrix(bu);
}




//  right = quad_sub_quad(K, quad_x_quad(quad_sub_quad(scalar_x_quad(1/delta_t,C),ACurrent),v_current));     
static void getRightHandSide(QuadMatrix *right, const QuadMatrix *k, const SparseNetlistQuad *netlist, const SparseQuadMatrix *aCurrent, const QuadMatrix *vCurrent)
{
	const int nodeNum = netlist->nodeNum;
	const int gvNum = netlist->gvNum;
	SparseQuadMatrix *z = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *z_sub_a = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	QuadMatrix *z_sub_a_v = createQuadMatrix(nodeNum,1,gvNum);

	scaleSparseQuadMatrix(z,1.0/netlist->deltaT,netlist->c);
	subSparseQuadMatrix(z_sub_a,z,aCurrent);
	mulVecSparseQuadMatrix(z_sub_a_v,z_sub_a,vCurrent);
	subQuadMatrix(right,k,z_sub_a_v);

	freeSparseQuadMatrix(z);
	freeSparseQuadMatrix(z_sub_a);
	freeQuadMatrix(z_sub_a_v);
}





// ======================================================


// gsl_matrix *result should be nodeNum * stepNum
// result should be pre-allocated

void quadNonlinearSimulation(const SparseNetlistQuad *netlist,QuadMatrix *result,const int threadNum,const int dumpNodeIndex)
{
	SparseQuadMatrix *gVarient = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *aCurrent = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *jacobian = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	QuadMatrix *oneStepResult = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *vDC = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *vCurrent = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *k = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *right = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *deltaV = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	gsl_matrix* meanDeltaV = gsl_matrix_alloc(netlist->nodeNum,1);

	// set the initial conditional
	setZeroQuadMatrix(result);	
	// ... wait to put

	int i,j;
	const int maxIteration = 20;
	const double magicRatio = 0.1;
	const double error = 0.01;
	

	time_t t1,t2;
	time(&t1);

	// i is the current step
	for(i=0;i<netlist->stepNum-1;i++)
	{
		time_t t1,t2;
		time(&t1);
		copyQuadMatrix(oneStepResult,vCurrent);
		updateGVarientMatrixFinite(gVarient,oneStepResult,netlist,0);
		addSparseQuadMatrix(aCurrent,netlist->a,gVarient);
		getK(k,netlist,vCurrent,i);  // k = C/dt * vt + B *ut+1
		solveDCSparseQuad(vDC,aCurrent,oneStepResult,i+1,netlist,1);
		copyQuadMatrix(vCurrent,vDC);
		updateGVarientMatrixFinite(gVarient,vCurrent,netlist,0);
		addSparseQuadMatrix(aCurrent,netlist->a,gVarient);
		time(&t2);
		fprintf(stderr,"update time in for:%g\n",difftime(t2,t1));


		time(&t1);
		for(j=0;j<maxIteration;j++)
		{
			time_t t1,t2,t3,t4;
			time(&t1);
			if(threadNum > 1)
				getJacobianFiniteParallel(jacobian,netlist,vCurrent,threadNum);
			else
				getJacobianFinite(jacobian,netlist,vCurrent,1);
			time(&t2);
			fprintf(stderr,"getJacobian time:%g\n",difftime(t2,t1));
			
			time(&t1);
			getRightHandSide(right,k,netlist,aCurrent,vCurrent);
			time(&t3);
			solveSparseQuadMatrix(deltaV,jacobian,right,1);
			time(&t4);
			fprintf(stderr,"solve:%g\n",difftime(t4,t3));
			scaleQuadMatrix(deltaV,magicRatio,deltaV);
			addQuadMatrix(vCurrent,vCurrent,deltaV);

			meanQuadMatrix(meanDeltaV,deltaV,netlist->s);
			time(&t2);
			fprintf(stderr,"time block2:%g\n",difftime(t2,t1));
			if(gsl_matrix_abs_max(meanDeltaV) < error)
			{
				break;
			}
			else
			{
				updateGVarientMatrixFinite(gVarient,vCurrent,netlist,0);
				addSparseQuadMatrix(aCurrent,netlist->a,gVarient);
			}

			if(j==maxIteration)
			{
				fprintf(stderr,"can not converge\n");
				fprintf(stderr,"step: %d, iteration: %d\n",i,j);
			}
		}
		time(&t2);
		fprintf(stderr,"converge time:%g\n",difftime(t2,t1));
		setQuadMatrix(result,getPtrEntryQuadMatrix(vCurrent,dumpNodeIndex,0),0,i+1);
		fprintf(stderr,"step: %d\n",i);
	}
//	exit(0);
	time(&t2);
	fprintf(stderr,"actual sim time:%g\n",difftime(t2,t1));

	freeSparseQuadMatrix(gVarient);
	freeSparseQuadMatrix(aCurrent);
	freeSparseQuadMatrix(jacobian);
	freeQuadMatrix(oneStepResult);
	freeQuadMatrix(vDC);
	freeQuadMatrix(vCurrent);
	freeQuadMatrix(k);
	freeQuadMatrix(right);
	freeQuadMatrix(deltaV);
	gsl_matrix_free(meanDeltaV);
}


