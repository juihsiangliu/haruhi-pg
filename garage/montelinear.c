#include "montelinear.h"



void monteLinearSimulation(const MonteNetlist *netlist, double *result,const int threadNum,const int dumpNodeIndex)
{
	// init & alloc
	const int nodeNum = netlist->nodeNum;
	const int stepNum = netlist->stepNum;
	int i;
	// -----
	SparseDoubleMatrix *y = createSparseDoubleMatrix(nodeNum,nodeNum);	
	SparseDoubleMatrix *yRefine = createSparseDoubleMatrix(nodeNum,nodeNum);	
	SparseDoubleMatrix *z = createSparseDoubleMatrix(nodeNum,nodeNum);	
	// used for lu
	SparseDoubleMatrix *p = createSparseDoubleMatrix(nodeNum,nodeNum);
	identitySparseDoubleMatrix(p);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(nodeNum,nodeNum);	
	identitySparseDoubleMatrix(pTrans);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(nodeNum,nodeNum);	
	SparseDoubleMatrix *u = createSparseDoubleMatrix(nodeNum,nodeNum);	
	// -----
	double *uCurrent = getMempoolSet(sizeof(double)*nodeNum);
	double *vCurrent = getMempoolSet(sizeof(double)*nodeNum);
	double *zv = getMempoolSet(sizeof(double)*nodeNum);
	double *bu = getMempoolSet(sizeof(double)*nodeNum);
	double *zvbu = getMempoolSet(sizeof(double)*nodeNum);
	double *vDC = getMempoolSet(sizeof(double)*nodeNum);

	// y = c/dt - a and chol(y)
	scaleSparseDoubleMatrix(z,1.0/netlist->deltaT,netlist->c); // z = c/dt
	subSparseDoubleMatrix(y,z,netlist->a); // y = z - a

	time_t reorderBegin,reorderEnd,luBegin,luEnd;	
	if(threadNum != 1)
	{
		const int goalPartition = 8;
//		const int goalPartition = 16;
		ParallelETree *tree2 = createParallelETree(goalPartition*4);
		time(&reorderBegin);
		partitionSparseDoubleMatrix(p,pTrans,tree2,yRefine,y,goalPartition);
		time(&reorderEnd);
		fprintf(stderr,"reordering time: %g\n",difftime(reorderEnd,reorderBegin));
		time(&luBegin);
		freeSparseDoubleMatrix(y);
		clearPidMempoolSet(0);
		parallelLUDouble(l,u,tree2,yRefine,p,threadNum,ooc);	
//		parallelILUDouble(l,u,tree2,yRefine,p,threadNum,0.0005);	
		time(&luEnd);
		fprintf(stderr,"lu time:%g\n",difftime(luEnd,luBegin));
		freeParallelETree(tree2);
		fprintf(stderr,"row: %d, nnz: %d, L->nnz: %d, U->nnz: %d\n",yRefine->totalRow,yRefine->nnz,l->nnz,u->nnz);
	}
	else
	{
		// directly partition
		time(&reorderBegin);
		amdSparseDoubleMatrix(p,y);
		transSparseDoubleMatrix(pTrans,p);
		permutateSparseDoubleMatrix(yRefine,p,pTrans,y);	
		time(&reorderEnd);
		fprintf(stderr,"reordering time: %g\n",difftime(reorderEnd,reorderBegin));
		time(&luBegin);
		freeSparseDoubleMatrix(y);
		clearPidMempoolSet(0);
	   	iluSparseDoubleMatrix(l,u,yRefine,0.0);
		time(&luEnd);
		fprintf(stderr,"lu time:%g\n",difftime(luEnd,luBegin));
		fprintf(stderr,"row: %d, nnz: %d, L->nnz: %d, U->nnz: %d\n",yRefine->totalRow,yRefine->nnz,l->nnz,u->nnz);
	}


	memset(vDC,0,sizeof(double)*netlist->nodeNum);

	for(i=0;i<netlist->stepNum-1;i++)
	{
		// z*vt
//		mulVecSparseDoubleMatrix(zv,z,vDC);
		parallelMulVecSparseDoubleMatrix(zv,z,vDC,threadNum);
		// b*ut
		const double *uCurrent = netlist->u[i];
//		mulVecSparseDoubleMatrix(bu,netlist->b,uCurrent);
		parallelMulVecSparseDoubleMatrix(bu,netlist->b,uCurrent,threadNum);
		// zv + bu
		addMyMatrix(zvbu,zv,bu,nodeNum,1);
		triSolveSparseDoubleMatrix(vDC,p,pTrans,l,u,zvbu);

//		parallelPCG(l,u,yRefine,vDC,zvbu,vDC,4);

		result[i+1] = getMyMatrix(vDC,nodeNum,1,dumpNodeIndex,0);
	}

	freeSparseDoubleMatrix(yRefine);
	freeSparseDoubleMatrix(z);
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);

	retMempoolSet(uCurrent,sizeof(double)*nodeNum);
	retMempoolSet(vCurrent,sizeof(double)*nodeNum);
	retMempoolSet(zv,sizeof(double)*nodeNum);
	retMempoolSet(bu,sizeof(double)*nodeNum);
	retMempoolSet(zvbu,sizeof(double)*nodeNum);
	retMempoolSet(vDC,sizeof(double)*nodeNum);
}	



/*
void monteLinearSimulation(const MonteNetlist *netlist,gsl_matrix *result,const int threadNum)
{
	// init & alloc
	const int nodeNum = netlist->nodeNum;
	const int gvNum = netlist->gvNum;
	const int stepNum = netlist->stepNum;
	int i;
	gsl_matrix *y = gsl_matrix_alloc(nodeNum,nodeNum);
	gsl_matrix *z = gsl_matrix_alloc(nodeNum,nodeNum);
	// -----
	gsl_vector *uCurrent = gsl_vector_alloc(nodeNum);
	gsl_vector *vCurrent = gsl_vector_alloc(nodeNum);
	gsl_vector *zv = gsl_vector_alloc(nodeNum);
	gsl_vector *bu = gsl_vector_alloc(nodeNum);
	gsl_vector *zvbu = gsl_vector_alloc(nodeNum);
	gsl_vector *vDC = gsl_vector_alloc(nodeNum);

	// y = c/dt - a and chol(y)
	gsl_matrix_memcpy(z,netlist->c);
	gsl_matrix_scale(z,1.0/netlist->deltaT);
	gsl_matrix_memcpy(y,z);
	gsl_matrix_sub(y,netlist->a);
	gsl_linalg_cholesky_decomp(y);

	// set the initial of v to result
	for(i=0;i<nodeNum;i++) gsl_matrix_set(result,i,0,0);

	// predict the nodal voltage of i+1 step
	for(i=0;i<stepNum-1;i++)
	{
		gsl_matrix_get_col(uCurrent,netlist->u,i);	
		gsl_matrix_get_col(vCurrent,result,i);	
		// z * vt
		gsl_blas_dgemv(CblasNoTrans,1.0,z,vCurrent,0.0,zv);
		// b * ut
		gsl_blas_dgemv(CblasNoTrans,1.0,netlist->b,uCurrent,0.0,bu);
		// z*vt + b*ut
		gsl_vector_set_zero(zvbu);
		gsl_vector_add(zvbu,zv);
		gsl_vector_add(zvbu,bu);
		gsl_linalg_cholesky_solve(y,zvbu,vDC);
		// set result
		gsl_matrix_set_col(result,i+1,vDC);
	}

	// free
	gsl_matrix_free(y);
	gsl_matrix_free(z);
	gsl_vector_free(uCurrent);
	gsl_vector_free(vCurrent);
	gsl_vector_free(zv);
	gsl_vector_free(bu);
	gsl_vector_free(zvbu);
	gsl_vector_free(vDC);
}
*/

