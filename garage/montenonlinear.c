#include "montenonlinear.h"


inline static int toZeroBase(const int x)
{
	return x-1;
}

struct GVarientIJ
{
	int i;
	int j;
	double gVarientElement;
	// will be updated in updateGVarientMatrixFinite() => map_fun_setGvarient()
	const double *v; // use "myMatrix"
};

typedef struct GVarientIJ GVarientIJ;


//==================================================================================================

static map_func_setGVarient(const gdsl_element_t E, gdsl_location_t LOCATION, void *USER_DATA)
{
	// init & alloc
	GVarientIJ *currentGVarientIJ = (GVarientIJ *) USER_DATA;
	const GControlInfoMonte *currentGControlInfo = (GControlInfoMonte *) E;
	const int gateIndex = toZeroBase(currentGControlInfo->gate);
	const int drainIndex = toZeroBase(currentGControlInfo->drain);
	const int sourceIndex = toZeroBase(currentGControlInfo->source);
	double vgs,vds;
	
	if(sourceIndex!=-1)
	{
		vgs = currentGVarientIJ->v[gateIndex] - currentGVarientIJ->v[sourceIndex];
		vds = currentGVarientIJ->v[drainIndex] - currentGVarientIJ->v[sourceIndex];
	}
	else
	{
		vgs = currentGVarientIJ->v[gateIndex];
		vds = currentGVarientIJ->v[drainIndex];
	}

	const double deltaXAxis = getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,1) - getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,0);
	const double deltaYAxis = getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,1) - getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,0);

	// first, use max to find out the minimum required rowIndexUp
	// second, use min to guarantee the boundary
	// (the rowIndexDown may up-side-down)
	// set rowIndex
	const int rowIndexTemp = GSL_MAX_INT( (vds - getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,0)) / deltaYAxis,1);
	const int rowIndex = GSL_MIN_INT(rowIndexTemp,currentGControlInfo->vdsListSize-3);
	const int rowIndexUp = rowIndex + 1;
	const int rowIndexDown = rowIndex - 1;

	// set colIndex	
	const int colIndexTemp = GSL_MAX_INT( (vgs - getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,0)) / deltaXAxis,1);
	const int colIndex = GSL_MIN_INT(colIndexTemp,currentGControlInfo->vgsListSize-3 );
	const int colIndexUp = colIndex + 1;
	const int colIndexDown = colIndex - 1;

	gsl_matrix * Z = gsl_matrix_alloc(4,3);
	gsl_matrix * ZTransZ = gsl_matrix_alloc(3,3);
	gsl_vector * ZTransY = gsl_vector_alloc(3);
	gsl_vector *y = gsl_vector_alloc(4);
	gsl_vector *b = gsl_vector_alloc(3);

	// curve fitting
	int i;
	double val;
	for(i=0;i<4;i++)	gsl_matrix_set(Z,i,0,1);
	gsl_matrix_set(Z,0,1,getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndexDown));
	gsl_matrix_set(Z,1,1,getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndexUp));
	gsl_matrix_set(Z,2,1,getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndex));
	gsl_matrix_set(Z,3,1,getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,rowIndex));
	gsl_matrix_set(Z,0,2,getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,colIndex));
	gsl_matrix_set(Z,1,2,getMyMatrix(currentGControlInfo->vgsList,1,currentGControlInfo->vgsListSize,0,colIndex));
	gsl_matrix_set(Z,2,2,getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,colIndexDown));
	gsl_matrix_set(Z,3,2,getMyMatrix(currentGControlInfo->vdsList,1,currentGControlInfo->vdsListSize,0,colIndexUp));
	
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Z,Z,0.0,ZTransZ);
	
	const double yDownIndex = getMyMatrix(currentGControlInfo->partialIdsVxs,currentGControlInfo->vdsListSize,currentGControlInfo->vgsListSize,rowIndexDown,colIndex);	
	const double yUpIndex = getMyMatrix(currentGControlInfo->partialIdsVxs,currentGControlInfo->vdsListSize,currentGControlInfo->vgsListSize,rowIndexUp,colIndex);	
	const double yIndexDown = getMyMatrix(currentGControlInfo->partialIdsVxs,currentGControlInfo->vdsListSize,currentGControlInfo->vgsListSize,rowIndex,colIndexDown);	
	const double yIndexUp = getMyMatrix(currentGControlInfo->partialIdsVxs,currentGControlInfo->vdsListSize,currentGControlInfo->vgsListSize,rowIndex,colIndexUp);
	gsl_vector_set(y,0,yDownIndex);
	gsl_vector_set(y,1,yUpIndex);
	gsl_vector_set(y,2,yIndexDown);
	gsl_vector_set(y,3,yIndexUp);
	
	gsl_blas_dgemv(CblasTrans,1.0,Z,y,0.0,ZTransY);
	gsl_matrix_solve(b,ZTransZ,ZTransY);

	// update gVarientElement according to the fitting result
	double gVarientElement = gsl_vector_get(b,0) + gsl_vector_get(b,1)*vds + gsl_vector_get(b,2)*vgs;

	gVarientElement = currentGControlInfo->sign *  gVarientElement;
	currentGVarientIJ->gVarientElement += gVarientElement;

	// free
	gsl_matrix_free(Z);
	gsl_matrix_free(ZTransZ);
	gsl_vector_free(ZTransY);
	gsl_vector_free(y);
	gsl_vector_free(b);
}



static void updateGVarientMatrixFinite(SparseDoubleMatrix *gVarient,const double* v,const MonteNetlist *netlist)
{

	// init & alloc
	const int nodeNum = netlist->nodeNum;
	int i,j;
	GVarientIJ *currentIJ = (GVarientIJ *)malloc(sizeof(GVarientIJ));
	// pre - set
	currentIJ->v = v;
	// set

	for(i=0;i<nodeNum;i++)
	{
		for(j=0;j<nodeNum;j++)
		{
			const gdsl_queue_t currentQueue = getQueueSparseGVarientTableMonte(netlist->gVarientTable,i,j);

			if(currentQueue == NULL) continue;
			else
			{
				currentIJ->i = i;
				currentIJ->j = j;
				currentIJ->gVarientElement = 0.0;
				if(currentQueue!=NULL) gdsl_queue_map_forward(currentQueue,map_func_setGVarient,currentIJ);
				setSparseDoubleMatrix(gVarient,currentIJ->gVarientElement,i,j);
			}
		}
	}
	free(currentIJ);
}





static void getJacobianFinite(SparseDoubleMatrix *jacobian, const MonteNetlist *netlist, const double *v)
{
	// init & alloc
	const double DV = 0.06;
	const int nodeNum = netlist->nodeNum;
	int i,j;
	clearSparseDoubleMatrix(jacobian);

	SparseDoubleMatrix *z = createSparseDoubleMatrix(nodeNum,nodeNum);
	double *deltaV = getMempoolSet(sizeof(double)*nodeNum);
	double *currentV = getMempoolSet(sizeof(double)*nodeNum);
	SparseDoubleMatrix *zSubACurrent = createSparseDoubleMatrix(nodeNum,nodeNum);
	double *fUp = getMempoolSet(sizeof(double)*nodeNum);
	double *fDown = getMempoolSet(sizeof(double)*nodeNum);
	double *deltaF = getMempoolSet(sizeof(double)*nodeNum);
	double *result = getMempoolSet(sizeof(double)*nodeNum);
	SparseDoubleMatrix *aCurrent = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *gVarient = createSparseDoubleMatrix(nodeNum,nodeNum);
	// set
	scaleSparseDoubleMatrix(z,1.0/netlist->deltaT,netlist->c);
	for(j=0;j<nodeNum;j++)
	{
		memset(deltaV,0,nodeNum*sizeof(double));
		setMyMatrix(deltaV,DV,nodeNum,1,j,0);
		addMyMatrix(currentV,v,deltaV,nodeNum,1);
		updateGVarientMatrixFinite(gVarient,currentV,netlist);
		addSparseDoubleMatrix(aCurrent,netlist->a,gVarient);
		subSparseDoubleMatrix(zSubACurrent,z,aCurrent);
		mulVecSparseDoubleMatrix(fUp,zSubACurrent,currentV);

		subMyMatrix(currentV,v,deltaV,nodeNum,1);
		updateGVarientMatrixFinite(gVarient,currentV,netlist);
		addSparseDoubleMatrix(aCurrent,netlist->a,gVarient);
		subSparseDoubleMatrix(zSubACurrent,z,aCurrent);
		mulVecSparseDoubleMatrix(fDown,zSubACurrent,currentV);
		
		subMyMatrix(deltaF,fUp,fDown,nodeNum,1);
		scaleMyMatrix(result,1.0/(2*DV),deltaF,nodeNum,1);
		setDenseColQuick2SparseDoubleMatrix(jacobian,result,nodeNum,j);
	}
	
	freeSparseDoubleMatrix(z);
	retMempoolSet(deltaV,sizeof(double)*nodeNum);
	retMempoolSet(currentV,sizeof(double)*nodeNum);
	freeSparseDoubleMatrix(zSubACurrent);
	retMempoolSet(fUp,sizeof(double)*nodeNum);
	retMempoolSet(fDown,sizeof(double)*nodeNum);
	retMempoolSet(deltaF,sizeof(double)*nodeNum);
	retMempoolSet(result,sizeof(double)*nodeNum);
	freeSparseDoubleMatrix(aCurrent);
	freeSparseDoubleMatrix(gVarient);
}




struct ParJacobian
{	
	// --------   const data ----------
	const MonteNetlist *netlist;
	const double *v;
	int lowerBound;
	int upperBound;
	int pid;

	// ---------  used local variables ----------------
	// in mempool pid
	SparseDoubleMatrix *jacobian;
	SparseDoubleMatrix *z;
	double *deltaVElement;
	double *deltaV;
	double *currentV;
	SparseDoubleMatrix *zSubACurrent;
	double *fUp;
	double *fDown;
	double *deltaF;
	double *result;
	SparseDoubleMatrix *aCurrent;
	SparseDoubleMatrix *gVarient; 
};










static void solveDC(double *vDC,const SparseDoubleMatrix *aCurrent, const double *vCurrent, const int nextUIndex, const MonteNetlist *netlist)
{

	// init and alloc
	const int nodeNum = netlist->nodeNum;
	//double *uNext = getMempool(nodalVoltagePool);
	SparseDoubleMatrix *z = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *z_sub_a = createSparseDoubleMatrix(nodeNum,nodeNum);
	double *zv = getMempoolSet(sizeof(double)*nodeNum);
	double *bu = getMempoolSet(sizeof(double)*nodeNum);
	double *zv_bu = getMempoolSet(sizeof(double)*nodeNum);
	// set
	const double *uNext = netlist->u[nextUIndex];
	scaleSparseDoubleMatrix(z,1.0/netlist->deltaT,netlist->c);
	subSparseDoubleMatrix(z_sub_a,z,aCurrent);
	mulVecSparseDoubleMatrix(zv,z,vCurrent);
	mulVecSparseDoubleMatrix(bu,netlist->b,uNext);
	addMyMatrix(zv_bu,zv,bu,1,nodeNum);
	// solve
	solveSparseDoubleMatrix(vDC,z_sub_a,zv_bu);
	// free
	freeSparseDoubleMatrix(z);
	freeSparseDoubleMatrix(z_sub_a);
	retMempoolSet(zv,sizeof(double)*nodeNum);
	retMempoolSet(bu,sizeof(double)*nodeNum);
	retMempoolSet(zv_bu,sizeof(double)*nodeNum);
}





// k = C/dt * vt + B * ut+1
static void getK(double *k, const MonteNetlist *netlist,const double *vCurrent,const int currentIndex)
{
	// init & alloc
	const int nodeNum = netlist->nodeNum;
	const double *uNext = netlist->u[currentIndex+1];
	SparseDoubleMatrix *z = createSparseDoubleMatrix(nodeNum,nodeNum);
	double *zv = getMempoolSet(sizeof(double)*nodeNum);
	double *bu = getMempoolSet(sizeof(double)*nodeNum);
	// set
	scaleSparseDoubleMatrix(z,1.0/netlist->deltaT,netlist->c);	
	mulVecSparseDoubleMatrix(zv,z,vCurrent);
	mulVecSparseDoubleMatrix(bu,netlist->b,uNext);
	addMyMatrix(k,zv,bu,nodeNum,1);
	// free
	freeSparseDoubleMatrix(z);
	retMempoolSet(zv,sizeof(double)*nodeNum);
	retMempoolSet(bu,sizeof(double)*nodeNum);
}



//  right = quad_sub_quad(K, quad_x_quad(quad_sub_quad(scalar_x_quad(1/delta_t,C),ACurrent),v_current));     
static void getRightHandSide(double *right, const double *k, const MonteNetlist *netlist, const SparseDoubleMatrix *aCurrent, const double *vCurrent)
{
	// init & alloc
	const int nodeNum = netlist->nodeNum;
	SparseDoubleMatrix *z = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *z_sub_a = createSparseDoubleMatrix(nodeNum,nodeNum);
	double *z_sub_a_v = getMempoolSet(sizeof(double)*nodeNum);
	// set
	scaleSparseDoubleMatrix(z,1.0/netlist->deltaT,netlist->c);	
	subSparseDoubleMatrix(z_sub_a,z,aCurrent);
	mulVecSparseDoubleMatrix(z_sub_a_v,z_sub_a,vCurrent);
	subMyMatrix(right,k,z_sub_a_v,nodeNum,1);
	// free
	freeSparseDoubleMatrix(z);
	freeSparseDoubleMatrix(z_sub_a);
	retMempoolSet(z_sub_a_v,sizeof(double)*nodeNum);
}



//=====================================================================

void monteNonlinearSimulation(const MonteNetlist *netlist,double *result,const int threadNum,const int dumpNodeIndex)
{
	// init & alloc
	const int nodeNum = netlist->nodeNum;
	SparseDoubleMatrix *gVarient = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *aCurrent = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *jacobian = createSparseDoubleMatrix(nodeNum,nodeNum);
	double *oneStepResult = getMempoolSet(sizeof(double)*nodeNum);
	double *vDC = getMempoolSet(sizeof(double)*nodeNum);
	double *vCurrent = getMempoolSet(sizeof(double)*nodeNum);
	double *k = getMempoolSet(sizeof(double)*nodeNum);
	double *right = getMempoolSet(sizeof(double)*nodeNum);
	double *deltaV = getMempoolSet(sizeof(double)*nodeNum);

	SparseDoubleMatrix *p = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(nodeNum,nodeNum);

	int i,j;
	const int maxIteration = 20;
	const double magicRatio = 0.1;
	const double error = 0.01;
	const stepNum = netlist->stepNum;

	// set the initial condition
	memset(vCurrent,0,sizeof(double)*nodeNum);


	getJacobianFinite(jacobian,netlist,vCurrent);
	amdSparseDoubleMatrix(p,jacobian);
	transSparseDoubleMatrix(pTrans,p);


	// i is the current step
	for(i=0;i<stepNum-1;i++)
	{
//		memcpy(oneStepResult,vCurrent,sizeof(double)*netlist->nodeNum);

		updateGVarientMatrixFinite(gVarient,vCurrent,netlist);
		addSparseDoubleMatrix(aCurrent,netlist->a,gVarient);
		getK(k,netlist,vCurrent,i); // k = c/dt * vt + b*ut+1
		solveDC(vDC,aCurrent,vCurrent,i+1,netlist);
		copyMyMatrix(vCurrent,vDC,nodeNum,1);
		updateGVarientMatrixFinite(gVarient,vCurrent,netlist);
		addSparseDoubleMatrix(aCurrent,netlist->a,gVarient);

		for(j=0;j<maxIteration;j++)
		{
			getJacobianFinite(jacobian,netlist,vCurrent);
			getRightHandSide(right,k,netlist,aCurrent,vCurrent);
//			solveSparseDoubleMatrix(deltaV,jacobian,right,1);
			solveWithPermutationSparseDoubleMatrix(deltaV,p,pTrans,jacobian,right);
//			fprintf(stderr,"\t solve deltaV\n");
//			fprintf(stderr,"\t solve \n");
			scaleMyMatrix(deltaV,magicRatio,deltaV,nodeNum,1);
//			fprintf(stderr,"\t trivial 1\n");
			addMyMatrix(vCurrent,vCurrent,deltaV,nodeNum,1);
//			fprintf(stderr,"\t trivial 2\n");
			
//			fprintf(stderr,"before if\n");
			if(absMaxMyMatrix(deltaV,nodeNum,1) < error )
			{
				break;
			}
			else
			{
				updateGVarientMatrixFinite(gVarient,vCurrent,netlist);
				addSparseDoubleMatrix(aCurrent,netlist->a,gVarient);
			}

			if(j == maxIteration)
			{
				fprintf(stderr,"can not converge\n");
				fprintf(stderr,"step: %d, iteration:%d\n",i,j);
			}
		}
		result[i+1] = vCurrent[dumpNodeIndex];
		//copyMyMatrix(result[i+1],vCurrent,nodeNum,1);
	}
	
	// free
	freeSparseDoubleMatrix(gVarient);
	freeSparseDoubleMatrix(aCurrent);
	freeSparseDoubleMatrix(jacobian);
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	retMempoolSet(oneStepResult,sizeof(double)*nodeNum);
	retMempoolSet(vDC,sizeof(double)*nodeNum);
	retMempoolSet(vCurrent,sizeof(double)*nodeNum);
	retMempoolSet(k,sizeof(double)*nodeNum);
	retMempoolSet(right,sizeof(double)*nodeNum);
	retMempoolSet(deltaV,sizeof(double)*nodeNum);
}



