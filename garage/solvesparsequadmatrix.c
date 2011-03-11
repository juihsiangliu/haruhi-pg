#include "solvesparsequadmatrix.h"


// a = plu
void luSparseQuadMatrix(SparseQuadMatrix *l, SparseQuadMatrix *u,const SparseQuadMatrix *a)
{
	luPidSparseQuadMatrix(l,u,a,0);
}





/*
// a = plu
void luPidSparseQuadMatrix(SparseQuadMatrix *l, SparseQuadMatrix *u,const SparseQuadMatrix *a,const int pid)
{
		// init & alloc
//		SparseQuadMatrix *aTemp = createSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum);
		SparseQuadMatrix *aTemp = createPidSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum,pid);
		QuadElement *scale = createQuadElement(a->gvNum);
		QuadElement *element = createQuadElement(a->gvNum);
//		preproceesing
//		copySparseQuadMatrix(aTemp,a);
		copyPidSparseQuadMatrix(aTemp,a,pid);
		identitySparseQuadMatrix(l);
		int i,j;

		SparseQuadElement *lBaseRow = NULL;
		SparseQuadElement *lBaseCol = NULL;
		SparseQuadElement *aTempRow = NULL;
		SparseQuadElement *aTempCol = NULL;
		SparseQuadElement *aTempDelRow = NULL;
		SparseQuadElement *aTempDelCol = NULL;

		SparseQuadElement **rowCache = getPidMempoolSet(sizeof(SparseQuadElement *)*a->totalRow,pid);
		SparseQuadElement **colCache = getPidMempoolSet(sizeof(SparseQuadElement *)*a->totalCol,pid);
		for(i=0;i<a->totalRow;i++) rowCache[i] = NULL;
		for(i=0;i<a->totalCol;i++) colCache[i] = NULL;
		// perform lu
		for(i=0;i<a->totalRow-1;i++) 
		{
//			fprintf(stderr,"row %d\n",i);
			const SparseQuadElement *aii = aTemp->rowIndex[i]->rowLink;
			const SparseQuadElement *eachRow = aii->colLink;
			while(eachRow != NULL) // update E
			{
				divPidQuadElement(scale,eachRow->data,aii->data,pid);
				setFastSparseQuadMatrix(l,scale,eachRow->row,i,&lBaseRow,&lBaseCol);
				SparseQuadElement *inRow = aii->rowLink;
				while(inRow != NULL)
				{
					mulPidQuadElement(element,scale,getSparseQuadMatrix(aTemp,i,inRow->col),pid);
					subQuadElement(element,getSparseQuadMatrix(aTemp,eachRow->row,inRow->col),element);
					setFastPidSparseQuadMatrix(aTemp,element,eachRow->row,inRow->col,&rowCache[eachRow->row],&colCache[inRow->col],pid);
					inRow = inRow->rowLink;
				}
				const SparseQuadElement *next = eachRow->colLink;
				delFastPidSparseQuadMatrix(aTemp,eachRow->row,i,&aTempDelRow,&aTempDelCol,pid);
				eachRow = next;
			}
		}
		copySparseQuadMatrix(u,aTemp);

		freePidSparseQuadMatrix(aTemp,pid);
		freeQuadElement(scale);
		freeQuadElement(element);
}
*/


// a = plu
void luPidSparseQuadMatrix(SparseQuadMatrix *l, SparseQuadMatrix *u,const SparseQuadMatrix *a,const int pid)
{
	// init & alloc
	const int gvNum = a->gvNum;
	QuadElement* scale = createPidQuadElement(gvNum,pid);
	QuadElement* element = createPidQuadElement(gvNum,pid);

	clearPidSparseQuadMatrix(l,pid);
	clearPidSparseQuadMatrix(u,pid);

//	preproceesing
	copyPidSparseQuadMatrix(u,a,pid);
	identityPidSparseQuadMatrix(l,pid);
	int i,j;
	
	SparseQuadElement *uiCache = NULL;	
	SparseQuadElement *uEachRowCache = NULL;	
	
	SparseQuadElement **lRowCache = getPidMempoolSet(sizeof(SparseQuadElement *)*l->totalRow,pid);
	SparseQuadElement *lColCache = NULL;
	SparseQuadElement *uRowCache = NULL;
	SparseQuadElement **uColCache = getPidMempoolSet(sizeof(SparseQuadElement *)*u->totalCol,pid);
	SparseQuadElement **delRowCache = getPidMempoolSet(sizeof(SparseQuadElement *)*u->totalRow,pid);
	SparseQuadElement *delColCache = NULL;
	
	for(i=0;i<l->totalRow;i++) lRowCache[i] = NULL;
	for(i=0;i<u->totalCol;i++) uColCache[i] = NULL;
	for(i=0;i<u->totalRow;i++) delRowCache[i] = NULL;

	// perform lu
	for(i=0;i<u->totalRow-1;i++) // permutate
	{
		const SparseQuadElement *aii = u->rowIndex[i]->rowLink;
		const SparseQuadElement *eachRow = aii->colLink;
		while(eachRow != NULL) // update E
		{
			divPidQuadElement(scale,eachRow->data,aii->data,pid);
			if(!isEmptyQuadElement(scale))
			{
				setFastPidSparseQuadMatrix(l,scale,eachRow->row,i,&lRowCache[eachRow->row],&lColCache,pid);
				SparseQuadElement *inRow = aii->rowLink;
				while(inRow != NULL)
				{
					mulPidQuadElement(element,scale,getFastRowSparseQuadMatrix(u,i,inRow->col,&uiCache),pid);
					subQuadElement(element,getFastRowSparseQuadMatrix(u,eachRow->row,inRow->col,&uEachRowCache),element);
					setFastPidSparseQuadMatrix(u,element,eachRow->row,inRow->col,&uRowCache,&uColCache[inRow->col],pid);
					inRow = inRow->rowLink;
				}
			}
			const SparseQuadElement *next = eachRow->colLink;
			delFastPidSparseQuadMatrix(u,eachRow->row,i,&delRowCache[eachRow->row],&delColCache,pid);
			eachRow = next;
		}
	}
	
	retPidMempoolSet(lRowCache,sizeof(SparseQuadElement *)*l->totalRow,pid);
	retPidMempoolSet(uColCache,sizeof(SparseQuadElement *)*u->totalCol,pid);
	retPidMempoolSet(delRowCache,sizeof(SparseQuadElement *)*u->totalRow,pid);
}


//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveSparseQuadMatrix(QuadMatrix *x,const SparseQuadMatrix *p,const SparseQuadMatrix *pTrans,const SparseQuadMatrix *l,const SparseQuadMatrix *u, const QuadMatrix *b)
{
	QuadMatrix *y = createQuadMatrix(x->row,1,x->gvNum);
	QuadMatrix *pb = createQuadMatrix(x->row,1,x->gvNum);
	QuadElement *sum = createQuadElement(x->gvNum);
	QuadElement *element = createQuadElement(x->gvNum);
	
	int i,k;
	SparseQuadElement *eachRow;
	SparseQuadElement *inRow;
	setZeroQuadElement(sum,x->gvNum);
	mulVecSparseQuadMatrix(pb,p,b);
	
	// solve ly = pb
	for(i=0;i<l->totalRow;i++)
	{
		resetQuadElement(sum);
		eachRow = l->rowIndex[i]->rowLink;
		inRow = eachRow;
		while(inRow->rowLink != NULL)
		{
			const QuadElement *lij = inRow->data;
			const QuadElement *yj = getPtrEntryQuadMatrix(y,inRow->col,0);
			mulQuadElement(element,lij,yj);
			addQuadElement(sum,sum,element);
			inRow = inRow->rowLink;
		}
		copyQuadElement(element,getPtrEntryQuadMatrix(pb,i,0));
		subQuadElement(element,element,sum);
		divQuadElement(element,element,inRow->data);
		setQuadMatrix(y,element,i,0);
	}

	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		resetQuadElement(sum);
		eachRow = u->rowIndex[i]->rowLink;
		inRow = eachRow->rowLink;
		while(inRow != NULL)
		{
			const QuadElement *uij = inRow->data;
			const QuadElement *xj = getPtrEntryQuadMatrix(x,inRow->col,0);
			mulQuadElement(element,uij,xj);
			addQuadElement(sum,sum,element);
			inRow = inRow->rowLink;
		}
		copyQuadElement(element,getPtrEntryQuadMatrix(y,i,0));
		subQuadElement(element,element,sum);
		divQuadElement(element,element,u->rowIndex[i]->rowLink->data);
		setQuadMatrix(x,element,i,0);
	}

	mulVecSparseQuadMatrix(x,pTrans,x);

	freeQuadMatrix(y);
	freeQuadMatrix(pb);
	freeQuadElement(sum);
	freeQuadElement(element);

}



// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b
void solveSparseQuadMatrix(QuadMatrix *x,const SparseQuadMatrix *a,const QuadMatrix *b, const int threadNum)
{
	const int nodeNum = a->totalRow;
	const int gvNum = a->gvNum;
	SparseQuadMatrix *p = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *pTrans = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *l = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *u = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *aRefine = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	
	if(threadNum == 1)
	{
		amdSparseQuadMatrix(p,a);
		transSparseQuadMatrix(pTrans,p);
		permutateSparseQuadMatrix(aRefine,p,pTrans,a);
		luSparseQuadMatrix(l,u,aRefine);
	}
	else
	{
		const int goalPartition = 4;
		ParallelETree *tree2 = createParallelETree(goalPartition*2 + goalPartition+1);
		partitionSparseQuadMatrix(p,pTrans,tree2,aRefine,a,goalPartition);
	    parallelLUQuad(l,u,tree2,aRefine,threadNum);
		freeParallelETree(tree2);
	}
	triSolveSparseQuadMatrix(x,p,pTrans,l,u,b);
		
	freeSparseQuadMatrix(p);
	freeSparseQuadMatrix(l);
	freeSparseQuadMatrix(u);
	freeSparseQuadMatrix(aRefine);
	freeSparseQuadMatrix(pTrans);

}




void solveWithPermutationSparseQuadMatrix(QuadMatrix *x,const SparseQuadMatrix *p, const SparseQuadMatrix *pTrans, const SparseQuadMatrix *a, const QuadMatrix *b, const int threadNum, ParallelETree *tree)
{
	const int nodeNum = a->totalRow;
	const int gvNum = a->gvNum;
	SparseQuadMatrix *l = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *u = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);
	SparseQuadMatrix *aRefine = createSparseQuadMatrix(nodeNum,nodeNum,gvNum);

	if(threadNum == 1)
	{
		time_t t1,t2,t3;
		permutateSparseQuadMatrix(aRefine,p,pTrans,a);
		time(&t1);
		luSparseQuadMatrix(l,u,aRefine);
		time(&t2);
		triSolveSparseQuadMatrix(x,p,pTrans,l,u,b);
		time(&t3);
//		fprintf(stderr,"----- lu  time:%g\n",difftime(t2,t1));
//		fprintf(stderr,"----- tri time:%g\n",difftime(t3,t2));
	}
	else
	{
	/*
		permutateSparseQuadMatrix(aRefine,p,pTrans,a);
	    parallelLUQuad(l,u,tree,aRefine,threadNum);
		triSolveSparseQuadMatrix(x,p,pTrans,l,u,b);
	*/
		fprintf(stderr,"not support parallel solving here\n");
	}

	triSolveSparseQuadMatrix(x,p,pTrans,l,u,b);
	freeSparseQuadMatrix(l);
	freeSparseQuadMatrix(u);
	freeSparseQuadMatrix(aRefine);
}
