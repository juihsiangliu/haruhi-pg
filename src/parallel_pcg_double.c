#include "parallel_pcg_double.h"



static double get_abs_max(const double *vec,const int size)
{
	int i;
	double max = 0.0;
	for(i=0;i<size;i++)
	{
		const double val = fabs(vec[i]);
		if( val > max ) max = val;
	}
	return max;
}




static double get_abs_mean(const double *vec,const int size)
{
	int i;
	double sum = 0.0;
	for(i=0;i<size;i++)
	{
		const double val = fabs(vec[i]);
		sum += val;
	}
	return sum / size;
}





static int violate_max_error(const double *vec,const int size,const double tol)
{
	if(get_abs_max(vec,size) < tol) return 0;
	else return 1;
}





static void dumpVector(FILE *fp, const double *vec, const int size)
{
	int i;
	for(i=0;i<size;i++) fprintf(fp,"%lf ",vec[i]);
	fprintf(fp,"\n");
}




static void swap_before_after(double **a,double **b)
{
	double *tmp = *a;
	*a = *b;
	*b = tmp;
}




static int set_pcg_goal_partition(const SparseDoubleMatrix *mtx)
{
	int goalPartition;
	int totalRow = mtx->totalRow;

	if(totalRow < 2000000) goalPartition = 16;
	else if(totalRow < 4000000) goalPartition = 16;
	else if(totalRow < 8000000) goalPartition = 32;
	else if(totalRow < 16000000) goalPartition = 64;
	else goalPartition = 64;

	return goalPartition;
/*	int ret = setGoalPartition(a);
	ret /= 2;
	if(ret < 8) return 8;
	else return ret;
*/
}




// use ILU
void parallelPCG(const SparseDoubleMatrix *a, const double *b,double *sol,const int threadNum,enum OOCFlag oocFlag,enum OrderMethod orderMethod)
{
	FILE *fp1 = fopen("mean_error.txt","w");
	FILE *fp2 = fopen("max_error.txt","w");

	time_t t1,t2;
	const double tol = 1e-4;
//	const double max_err = 1e-5;
	const double max_err = 1e-4;

	CSR_SparseDoubleMatrix *a_csr = sparse2CSR(a,0);
	SparseDoubleMatrix *p = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(a->totalRow,a->totalCol);	
	identitySparseDoubleMatrix(p);
	identitySparseDoubleMatrix(pTrans);

	CSR_SparseDoubleMatrix *l_csr = NULL;
	CSR_SparseDoubleMatrix *u_csr = NULL;
	struct OOCInfo *oocInfoList = NULL;
	// ILU stage
	if( orderAmd == orderMethod) // it must be ic && threadNum = 1
	{
		time(&t1);
		SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(a->totalRow,a->totalCol);
		amdSparseDoubleMatrix(p,a);
		transSparseDoubleMatrix(pTrans,p);
		permutateSparseDoubleMatrix(aRefine,p,pTrans,a);
		time(&t2);
		fprintf(stderr,"reorder time: %g\n",difftime(t2,t1));
		time(&t1);
		SparseDoubleMatrix *l = createSparseDoubleMatrix(a->totalRow,a->totalRow);
		SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);
		iluSparseDoubleMatrix(l,u,aRefine,tol);
		freeSparseDoubleMatrix(aRefine);

		// un-comment following to perform CG
//		identitySparseDoubleMatrix(l);
//		identitySparseDoubleMatrix(u);
		
		l_csr = sparse2CSR(l,0);
		freeSparseDoubleMatrix(l);
		u_csr = sparse2CSR(u,0);
		freeSparseDoubleMatrix(u);
		time(&t2);
		fprintf(stderr,"ilu time: %g\n",difftime(t2,t1));
	}
	else if( orderMetis == orderMethod)
	{
		time(&t1);
		const int goalPartition = set_pcg_goal_partition(a);
		// construct elimination tree
		ParallelETree *tree = createParallelETree(goalPartition*4);
		SparseDoubleMatrix *aRefine = partitionSparseDoubleMatrix(p,pTrans,tree,a,goalPartition,oocFlag);
		time(&t2);
		fprintf(stderr,"reorder time: %g\n",difftime(t2,t1));
		time(&t1);
		SparseDoubleMatrix *l = createSparseDoubleMatrix(a->totalRow,a->totalRow);
		SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);
		oocInfoList = parallelILUDouble(l,u,tree,aRefine,p,threadNum,tol,oocFlag);	
		fprintf(stderr,"a->nnz:%d l->nnz:%d\n",aRefine->nnz,l->nnz);
		freeParallelETree(tree);
		if(ic == oocFlag)
		{
			freeSparseDoubleMatrix(aRefine);
			l_csr = sparse2CSR(l,0);
			u_csr = sparse2CSR(u,0);
		}
		freeSparseDoubleMatrix(l);
		freeSparseDoubleMatrix(u);
		time(&t2);
		fprintf(stderr,"ilu time: %g\n",difftime(t2,t1));
	}
	else
	{
		fprintf(stderr,"pcg_pre_method is not defined\n");
		exit(0);
	}


	int i;
	int k = 0;
	const int n = a->totalCol;

	double *x_init = getMempoolSet(sizeof(double)*n);
	memset(x_init,0,sizeof(double)*n);
	if(oocFlag == ic) triSolve_CSR_SparseDoubleMatrix(x_init,p,pTrans,l_csr,u_csr,b);
	else oocTriSolveSparseDoubleMatrix(x_init,oocInfoList,p,pTrans,b);

	double *rBefore = getMempoolSet(sizeof(double)*n);
	double *rAfter = getMempoolSet(sizeof(double)*n);

	double *xBefore = getMempoolSet(sizeof(double)*n);
	double *xAfter = getMempoolSet(sizeof(double)*n);
	
	double *zBefore = getMempoolSet(sizeof(double)*n);
	double *zAfter = getMempoolSet(sizeof(double)*n);
	
	double *pBefore = getMempoolSet(sizeof(double)*n);
	double *pAfter = getMempoolSet(sizeof(double)*n);

	// r0 = b-Ax0
	double *Ax0 = getMempoolSet(sizeof(double)*n);
	parallel_mulVec_CSR_SparseDoubleMatrix(Ax0,a_csr,x_init,threadNum);
	for(i=0;i<n;i++) rBefore[i] = b[i] - Ax0[i];
	retMempoolSet(Ax0,sizeof(double)*n);
	
	if(!violate_max_error(rBefore,n,max_err))
	{
		memcpy(sol,x_init,sizeof(double)*n);
		return;
	}

	// z0 = inv(M) * r0
	if(oocFlag == ic) triSolve_CSR_SparseDoubleMatrix(zBefore,p,pTrans,l_csr,u_csr,rBefore);
	else oocTriSolveSparseDoubleMatrix(zBefore,oocInfoList,p,pTrans,rBefore);
	
	// p0 = z0
	memcpy(pBefore,zBefore,sizeof(double)*n);

	memcpy(xBefore,x_init,sizeof(double)*n);
	double alpha;
	double beta;
	double *Ap = getMempoolSet(sizeof(double)*n);
	for(k=0;k<n;k++)
	{
		double tmp1;
		double tmp2;
		// Ap
		parallel_mulVec_CSR_SparseDoubleMatrix(Ap,a_csr,pBefore,threadNum);
		// alpha
		tmp1 = tmp2 = 0.0;
		for(i=0;i<n;i++) tmp1 += rBefore[i]*zBefore[i];
		for(i=0;i<n;i++) tmp2 += pBefore[i]*Ap[i];
		alpha = tmp1 / tmp2;
		// x_next = x_now + alpha * p
		for(i=0;i<n;i++) xAfter[i] = xBefore[i] + alpha*pBefore[i];
		// r_next = r_current - alpha * Ap
		for(i=0;i<n;i++) rAfter[i] = rBefore[i] - alpha*Ap[i];

	 	fprintf(fp1,"%g\n",get_abs_mean(rAfter,n));
	 	fprintf(fp2,"%g\n",get_abs_max(rAfter,n));
		// if r_next is small enough , break
		if(!violate_max_error(rAfter,n,max_err) || (k == n-1))
		{
			memcpy(sol,xAfter,sizeof(double)*n);
			break;
		}
		else
		{
//			fprintf(stderr,"iteration = %d\n",k);
		}
		// z_next = inv(M) * r_next
		if(ic == oocFlag) triSolve_CSR_SparseDoubleMatrix(zAfter,p,pTrans,l_csr,u_csr,rAfter);
		else oocTriSolveSparseDoubleMatrix(zAfter,oocInfoList,p,pTrans,rAfter);
		// beta
		tmp2 = 0.0;
		for(i=0;i<n;i++) tmp2 += rAfter[i]*zAfter[i];
		beta = tmp2 / tmp1;
		// p_next = z_next + beta * p_current
		for(i=0;i<n;i++) pAfter[i] = zAfter[i] + beta * pBefore[i];

		swap_before_after(&rBefore,&rAfter);
		swap_before_after(&xBefore,&xAfter);
		swap_before_after(&zBefore,&zAfter);
		swap_before_after(&pBefore,&pAfter);
	}
	fprintf(stderr,"max error = %g\n",get_abs_max(rAfter,n));
	
	retMempoolSet(Ap,sizeof(double)*n);
	retMempoolSet(rBefore,sizeof(double)*n);
	retMempoolSet(rAfter,sizeof(double)*n);
	retMempoolSet(xBefore,sizeof(double)*n);
	retMempoolSet(xAfter,sizeof(double)*n);
	retMempoolSet(zBefore,sizeof(double)*n);
	retMempoolSet(zAfter,sizeof(double)*n);
	retMempoolSet(pBefore,sizeof(double)*n);
	retMempoolSet(pAfter,sizeof(double)*n);
	retMempoolSet(x_init,sizeof(double)*n);
	free_CSR_SparseDoubleMatrix(a_csr);
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);

	if( orderMetis == orderMethod) freeOOCInfoList(oocInfoList);
	if( orderAmd == orderMethod || (orderMetis == orderMethod && ic == oocFlag) )
	{
		free_CSR_SparseDoubleMatrix(l_csr);
		free_CSR_SparseDoubleMatrix(u_csr);
	}

	if(k!=n-1) fprintf(stderr,"pcg converage in iteration: %d\n",k);
	else fprintf(stderr,"pcg can not converge until max iteration: %d\n",k);


	fclose(fp1);
	fclose(fp2);
}

