#include "solvesparsedoublematrix.h"



// ============================================




// a = plu
void luSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a)
{
	luPidSparseDoubleMatrix(l,u,a,0);
}




void luPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid)
{
	iluPidSparseDoubleMatrix(l,u,a,pid,0.0);
}




static double getColNorm(const SparseDoubleMatrix *a,const int col)
{
	double sum = 0;
	const SparseDoubleElement *ptr = a->colIndex[col]->colLink;
	while(ptr != NULL)
	{
		sum += ptr->data*ptr->data;
		ptr = ptr->colLink;
	}

	return sqrt(sum);
}




// return 1 if it can be drop, else 0
inline static int dropLij(const double lij,const double tol,const SparseDoubleMatrix *u,const int i,const int j,const double ujj,const SparseDoubleMatrix *a)
{
	if( i == j || tol == 0.0) return 0;

	if(tol != -1)
	{	
//		if( fabs(lij)*(ujj) > tol*colNorm[j] ) return 0;
		if( fabs(lij)*(ujj) > tol ) return 0;
		else return 1;
	}
	else // zero fill-in
	{
		const double aij = getSparseDoubleMatrix(a,i,j,"col");
		if(aij!=0) return 0;
		else return 1;
	}
}




// return 1 if it can be drop, else 0
inline static int dropUij(const double uij,const double tol,const SparseDoubleMatrix *u,const int i,const int j,const double uii,const SparseDoubleMatrix *a)
{
	if( i == j || tol == 0.0) return 0;

	if(tol != -1)
	{
//		if( fabs(uij) > tol*colNorm[j] ) return 0;
		if( fabs(uij) > tol*uii ) return 0;
		else return 1;
	}
	else // zero fill-in
	{
		const double aij = getSparseDoubleMatrix(a,i,j,"row");
		if(aij!=0) return 0;
		else return 1;
	}
}




void iluSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const double tol)
{
	iluPidSparseDoubleMatrix(l,u,a,0,tol);
}




void iluPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid,const double tol)
{
	// init & alloc
	int i,j;
	double element = 0.0;
		
	clearPidSparseDoubleMatrix(l,pid);
	clearPidSparseDoubleMatrix(u,pid);

//	preproceesing
	copyPidSparseDoubleMatrix(u,a,pid);
	identityPidSparseDoubleMatrix(l,pid);
	// perform lu
	for(i=0;i<u->totalRow-1;i++)
	{
		const SparseDoubleElement *uii_ptr = u->rowIndex[i]->rowLink;
		const SparseDoubleElement *eachRow = uii_ptr->colLink;
		const double drop_l_val = getSparseDoubleMatrix(u,i,i,"row");
		assert(uii_ptr->data!=0);
		while(eachRow != NULL)
		{
			double drop_u_val = getSparseDoubleMatrix(u,eachRow->row,eachRow->row,"row");
			const double scale = eachRow->data / uii_ptr->data;

			if(scale!=0)
			{	
				if(!dropLij(scale,tol,u,eachRow->row,i,drop_l_val,a))
					setPidSparseDoubleMatrix(l,scale,eachRow->row,i,pid);
				SparseDoubleElement *inRow = uii_ptr->rowLink;
				while(inRow != NULL)
				{
					element = scale * getSparseDoubleMatrix(u,i,inRow->col,"row");
					element = getSparseDoubleMatrix(u,eachRow->row,inRow->col,"row") - element;
					if(!dropUij(element,tol,u,eachRow->row,inRow->col,drop_u_val,a))
						setPidSparseDoubleMatrix(u,element,eachRow->row,inRow->col,pid);
					inRow = inRow->rowLink;
				}
			}
			const SparseDoubleElement *next = eachRow->colLink;
			delPidSparseDoubleMatrix(u,eachRow->row,i,pid);
			eachRow = next;
		}
	}
}

/*


static void copy_ai_to_ui(SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int i,const int pid)
{
	const SparseDoubleElement *ai_ptr = a->rowIndex[i]->rowLink;
	while(ai_ptr!=NULL)
	{
		setPidSparseDoubleMatrix(u,ai_ptr->data,i,ai_ptr->col,pid);
		ai_ptr = ai_ptr->rowLink;
	}
}




static void copy_ai_to_buf(double *buf,const int bufSize,const SparseDoubleMatrix *a,const int i)
{
	memset(buf,0,bufSize*sizeof(double));
	const SparseDoubleElement *ai_ptr = a->rowIndex[i]->rowLink;
	while(ai_ptr!=NULL)
	{
		buf[ai_ptr->col] = ai_ptr->data;
		ai_ptr = ai_ptr->rowLink;
	}
}



static void copy_buf_to_ui(SparseDoubleMatrix *u,const int i,const double *buf,const int bufSize,const int pid)
{
	int idx = 0;
	for(idx=0;idx<bufSize;idx++)
	{
		if(buf[idx]!=0) setPidSparseDoubleMatrix(u,buf[idx],i,idx,pid);
	}
}


static int first_nnz_in_buf(double *val,const double *buf,const int bufSize)
{
	int idx = 0;
	for(idx=0;idx<bufSize;idx++)
	{
		if(buf[idx]!=0) 
		{
			*val = buf[idx];
			return idx;
		}
	}
	*val = 0;
	return idx;
}



void iluPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid,const double tol)
{
	// init & alloc
	int i,j;
	double *buf = getPidMempoolSet(a->totalCol*sizeof(double),pid);
	clearPidSparseDoubleMatrix(l,pid);
	clearPidSparseDoubleMatrix(u,pid);
	for(i=0;i<a->totalRow;i++) // for each row to eliminate
	{
		setPidSparseDoubleMatrix(l,1.0,i,i,pid);
		copy_ai_to_buf(buf,a->totalCol,a,i);
		for(j=0;j<i;j++) // for each row before row_i
		{
			double ui;
			const SparseDoubleElement *uj_ptr = u->rowIndex[j]->rowLink;
			if(first_nnz_in_buf(&ui,buf,a->totalCol) == uj_ptr->col)
			{
				const double scale = ui / uj_ptr->data;
				setPidSparseDoubleMatrix(l,scale,i,j,pid);
				while(uj_ptr!=NULL)
				{
					buf[uj_ptr->col] -= scale*uj_ptr->data;
					uj_ptr = uj_ptr->rowLink;
				}
				buf[j] = 0.0;
			}
		}
		copy_buf_to_ui(u,i,buf,a->totalCol,pid);
	}
	retPidMempoolSet(buf,a->totalCol*sizeof(double),pid);
}

*/


int symboic_nnz(const SparseDoubleMatrix *a)
{
	// init & alloc
	
	SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);

//	preproceesing
	copySparseDoubleMatrix(u,a);
	int i,j;
	// perform lu
	for(i=0;i<u->totalRow-1;i++)
	{
		const SparseDoubleElement *uii_ptr = u->rowIndex[i]->rowLink;
		const SparseDoubleElement *eachRow = uii_ptr->colLink;
		while(eachRow != NULL)
		{
			const double scale = eachRow->data;

			if(scale!=0)
			{	
				SparseDoubleElement *inRow = uii_ptr->rowLink;
				while(inRow != NULL)
				{
					setSparseDoubleMatrix(u,1.0,eachRow->row,inRow->col);
					inRow = inRow->rowLink;
				}
			}
			const SparseDoubleElement *next = eachRow->colLink;
			delSparseDoubleMatrix(u,eachRow->row,i);
			eachRow = next;
		}
	}

	return u->nnz;
}






static double getSumInChol(const SparseDoubleMatrix *l, const int i,const int j)
{
	double sum = 0.0;
	int colI,colJ;
	const SparseDoubleElement *rowIPtr = l->rowIndex[i]->rowLink;
	const SparseDoubleElement *rowJPtr = l->rowIndex[j]->rowLink;
	while(rowIPtr!=NULL && rowJPtr!=NULL)
	{
		colI = rowIPtr->col;
		colJ = rowJPtr->col;
		if(colI>=j || colJ>=j)
		{
			break;
		}
		else
		{
			if(colI == colJ)
			{
				sum += rowIPtr->data * rowJPtr->data;
				rowIPtr = rowIPtr->rowLink;
				rowJPtr = rowJPtr->rowLink;
			}
			else if(colI > colJ)
			{
				rowJPtr = rowJPtr->rowLink;
			}
			else
			{
				rowIPtr = rowIPtr->rowLink;
			}
		}
	}
	return sum;
}


static void cdiv(SparseDoubleMatrix *l,const int j)
{
	SparseDoubleElement *ajjPtr = l->colIndex[j]->colLink;
	ajjPtr->data = sqrt(ajjPtr->data);
	const double ajj = ajjPtr->data;
	SparseDoubleElement *currentPtr = ajjPtr->colLink;
	while(currentPtr!=NULL)
	{
		int row = currentPtr->row;
		const double aij = currentPtr->data;
		currentPtr->data = aij/ajj;
		currentPtr = currentPtr->colLink;
	}
}




static void cmod(SparseDoubleMatrix *l,const int j,const int k,const double ajk)
{
	int i;
	/*
	for(i=j;i<l->totalRow;i++)
	{
		const double aij = getSparseDoubleMatrix(l,i,j);
		const double aik = getSparseDoubleMatrix(l,i,k);
		const double result = aij-aik*ajk;
		if(result!=0) setSparseDoubleMatrix(l,result,i,j);
		else delSparseDoubleMatrix(l,i,j);
	}
*/
	SparseDoubleElement *colKPtr = l->colIndex[k]->colLink;
	SparseDoubleElement *base = NULL;
	while(colKPtr!=NULL)
	{
		i = colKPtr->row;
		if(i < j)
		{
			colKPtr = colKPtr->colLink;
		}
		else
		{
			const double aik = colKPtr->data;
			const double result = aik*ajk;
			if(result!=0)
			{
//				decSparseDoubleMatrix(l,result,i,j);
				base = decFastSparseDoubleMatrix(l,result,i,j,colKPtr,base);
			}
			else
			{
				fprintf(stderr,"damn\n");
				delSparseDoubleMatrix(l,i,j);
			}
			colKPtr = colKPtr->colLink;
		}
	}
}




void cholSparseDoubleMatrix(SparseDoubleMatrix *l, const SparseDoubleMatrix *a)
{
	// copy the lower triangular from a to l
	clearSparseDoubleMatrix(l);
	int i,j,k;
	for(i=0;i<a->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = a->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			const int row = i;
			const int col = eachRow->col;
			if(col > row)
			{
				break;
			}
			else
			{
				setSparseDoubleMatrix(l,eachRow->data,row,col);		
				eachRow = eachRow->rowLink;
			}
		}
	}

/*
//	dumpSparseDoubleMatrix(stderr,l);
	// perfom sparse cholesky , left looking
	for(j=0;j<a->totalRow;j++) // for each row j
	{
		fprintf(stderr,"%d/%d\n",a->totalRow,j);
		SparseDoubleElement *currentPtr = l->rowIndex[j]->rowLink;
		while(currentPtr!=NULL)
		{
			k = currentPtr->col;
			const double ajk = currentPtr->data;
			if(k>=j) break;
			else
			{
				cmod(l,j,k,ajk);
				currentPtr = currentPtr->rowLink;
			}
		}
		cdiv(l,j);
	}
*/

	// cholesky ~ right looking
	for(k=0;k<a->totalCol;k++)
	{
		fprintf(stderr,"%d/%d\n",k,a->totalRow);
		cdiv(l,k);
		SparseDoubleElement *currentPtr = l->colIndex[k]->colLink;
		while(currentPtr!=NULL)
		{
			j = currentPtr->row;
			const double ajk = currentPtr->data;
			if(j>k)
			{
				cmod(l,j,k,ajk);
			}
			currentPtr = currentPtr->colLink;
		}
	}


//	dumpSparseDoubleMatrix(stderr,l);

	/* the dense implementation
	int i;
	int j;
	clearSparseDoubleMatrix(l);
	for(i=0;i<a->totalRow;i++)
	{
		for(j=0;j<=i;j++)
		{
			if(i==j)
			{
				const double aii = getSparseDoubleMatrix(a,i,i);
				const double sum = getSumInChol(l,i,i);
				const double result = sqrt(aii - sum);
				setSparseDoubleMatrix(l,result,i,i);
			}
			else
			{
				const double ljj = getSparseDoubleMatrix(l,j,j);
				const double aij = getSparseDoubleMatrix(a,i,j);
				const double sum = getSumInChol(l,i,j);
				const double result = (aij-sum)/ljj;
				if(result!=0) setSparseDoubleMatrix(l,result,i,j);
			}
		}
	}
	*/


	
}






//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b)
{
	CSR_SparseDoubleMatrix *l_csr = sparse2CSR(l,0);
	CSR_SparseDoubleMatrix *u_csr = sparse2CSR(u,0);
	triSolve_CSR_SparseDoubleMatrix(x,p,pTrans,l_csr,u_csr,b);
	free_CSR_SparseDoubleMatrix(u_csr);
	free_CSR_SparseDoubleMatrix(l_csr);
/*
	double *y = getMempoolSet(sizeof(double)*u->totalCol);
	double *pb = getMempoolSet(sizeof(double)*p->totalRow);

	double sum = 0.0;
	double element = 0.0;
	
	int i,k;
	SparseDoubleElement *eachRow;
	SparseDoubleElement *inRow;
	for(i=0;i<u->totalRow;i++)
	{
		int col = p->rowIndex[i]->rowLink->col;
		pb[i] = b[col];
	}
	
	// solve ly = pb
	for(i=0;i<l->totalRow;i++)
	{
		sum = 0.0;
		eachRow = l->rowIndex[i]->rowLink;
		inRow = eachRow;
		while(inRow->rowLink != NULL)
		{
			const double lij = inRow->data;
			const double yj = y[inRow->col];
			element = lij*yj;
			sum = sum+element;
			inRow = inRow->rowLink;
		}
		element = pb[i];
		element = (element - sum)/inRow->data;
		y[i] = element;
	}
	retMempoolSet(pb,sizeof(double)*p->totalRow);
	pb = NULL;

	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		sum = 0;
		eachRow = u->rowIndex[i]->rowLink;
		inRow = eachRow->rowLink;
		while(inRow != NULL)
		{
			const double uij = inRow->data;
			const double xj = x[inRow->col];
			element = uij*xj;
			sum = sum + element;
			inRow = inRow->rowLink;
		}
		element = y[i];
		element = (element-sum) / u->rowIndex[i]->rowLink->data;
		x[i] = element;
	}
	double *xTemp = getMempoolSet(sizeof(double)*u->totalRow);
	for(i=0;i<u->totalRow;i++)
	{
		int col = pTrans->rowIndex[i]->rowLink->col;
		xTemp[i] = x[col];
	}
	memcpy(x,xTemp,sizeof(double)*u->totalRow);
	retMempoolSet(xTemp,sizeof(double)*u->totalRow);

	retMempoolSet(y,sizeof(double)*u->totalCol);
*/
}





void triNoPSolveSparseDoubleMatrix(double *x, const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b)
{
	double *y = getMempoolSet(sizeof(double)*u->totalCol);

	double sum = 0.0;
	double element = 0.0;
	
	int i,k;
	SparseDoubleElement *eachRow;
	SparseDoubleElement *inRow;
	
	// solve ly = b
	for(i=0;i<l->totalRow;i++)
	{
		sum = 0.0;
		eachRow = l->rowIndex[i]->rowLink;
		inRow = eachRow;
		while(inRow->rowLink != NULL)
		{
			const double lij = inRow->data;
			const double yj = y[inRow->col];
			element = lij*yj;
			sum = sum+element;
			inRow = inRow->rowLink;
		}
		element = b[i];
		element = (element - sum)/inRow->data;
		y[i] = element;
	}

	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		sum = 0;
		eachRow = u->rowIndex[i]->rowLink;
		inRow = eachRow->rowLink;
		while(inRow != NULL)
		{
			const double uij = inRow->data;
			const double xj = x[inRow->col];
			element = uij*xj;
			sum = sum + element;
			inRow = inRow->rowLink;
		}
		element = y[i];
		element = (element-sum) / u->rowIndex[i]->rowLink->data;
		x[i] = element;
	}
	retMempoolSet(y,sizeof(double)*u->totalCol);
}



// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b
void solveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *a,const double *b)
{
	const int nodeNum = a->totalRow;
	SparseDoubleMatrix *p = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(nodeNum,nodeNum);

	amdSparseDoubleMatrix(p,a);
	transSparseDoubleMatrix(pTrans,p);
	mulSparseDoubleMatrix(aRefine,p,a);
	mulSparseDoubleMatrix(aRefine,aRefine,pTrans);

	luSparseDoubleMatrix(l,u,aRefine);
	triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	freeSparseDoubleMatrix(aRefine);
	freeSparseDoubleMatrix(pTrans);
}




void solveWithPermutationSparseDoubleMatrix(double *x,const SparseDoubleMatrix *p, const SparseDoubleMatrix *pTrans, const SparseDoubleMatrix *a, const double *b)
{
	const int nodeNum = a->totalRow;
	SparseDoubleMatrix *l = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(nodeNum,nodeNum);

	mulSparseDoubleMatrix(aRefine,p,a);
	mulSparseDoubleMatrix(aRefine,aRefine,pTrans);

	luSparseDoubleMatrix(l,u,aRefine);
	triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	freeSparseDoubleMatrix(aRefine);
}




static double myNorm2(const double *a,const int size)
{
	int i;
	double sum = 0.0;
	for(i=0;i<size;i++) sum += a[i]*a[i];

	return sqrt(sum);
}






int setGoalPartition(const SparseDoubleMatrix *mtx)
{
	int goalPartition;
	int totalRow = mtx->totalRow;

	if(totalRow < 800000) goalPartition = 8;
	else if(totalRow < 1600000) goalPartition = 16;
	else if(totalRow < 3200000) goalPartition = 32;
	else if(totalRow < 6400000) goalPartition = 64;
	else goalPartition = 64;

	return goalPartition;
}






// =======================================================
//  CSR family



void triSolve_CSR_SparseDoubleMatrix(double *x,const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const CSR_SparseDoubleMatrix *l,const CSR_SparseDoubleMatrix *u,const double *b)
{
	int i,j;
	double *pb = getMempoolSet(sizeof(double)*p->totalRow);
	memset(pb,0,sizeof(double)*p->totalRow);
	
	for(i=0;i<u->totalRow;i++)
	{
		const int col = p->rowIndex[i]->rowLink->col;
		pb[i] = b[col];
	}
	
	// solve ly = pb
	double *y = getMempoolSet(sizeof(double)*u->totalCol);
	memset(y,0,sizeof(double)*u->totalCol);
	for(i=0;i<l->totalRow;i++)
	{
		double sum = 0;
		for(j=l->rowPtr[i];j<(l->rowPtr[i+1]-1);j++)
		{
			const double lij = l->val[j];
			const double yj = y[l->col[j]];
			sum += lij*yj;
		}
		y[i] = (pb[i]-sum)/l->val[l->rowPtr[i+1]-1];
	}
	retMempoolSet(pb,sizeof(double)*p->totalRow);

	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		double sum = 0;
		double element = 0;
		for(j=u->rowPtr[i]+1;j<u->rowPtr[i+1];j++)
		{
			const double uij = u->val[j];
			const double xj = x[u->col[j]];
			sum += uij*xj;
		}
		x[i] = (y[i]-sum) / u->val[u->rowPtr[i]];
	}
	retMempoolSet(y,sizeof(double)*u->totalCol);
	
	// post ordering
	double *xTemp = getMempoolSet(sizeof(double)*u->totalRow);
	for(i=0;i<u->totalRow;i++)
	{
		const int col = pTrans->rowIndex[i]->rowLink->col;
		xTemp[i] = x[col];
	}
	memcpy(x,xTemp,sizeof(double)*u->totalRow);
	retMempoolSet(xTemp,sizeof(double)*u->totalRow);
}


