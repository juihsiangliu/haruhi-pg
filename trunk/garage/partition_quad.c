#include "partition_quad.h"


static void sparseQuadMatrix2AmdAi(int *ai,const SparseQuadMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	SparseQuadElement *current;
	for(i=0;i<totalRow;i++)
	{
		current = a->rowIndex[i]->rowLink;
		while(current!=NULL)
		{
			ai[index] = current->col;
			index++;
			current = current->rowLink;
		}
	}
}



static void sparseQuadMatrix2AmdAp(int *ap,const SparseQuadMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	int numInRow = 0;
	SparseQuadElement *current;
	ap[0] = 0;
	for(i=0;i<totalRow;i++)
	{
		current = a->rowIndex[i]->rowLink;
		numInRow = 0;
		while(current!=NULL)
		{
			numInRow++;
			current = current->rowLink;
		}
		ap[i+1] = ap[i] + numInRow;
	}	
}


//====================================================================

void amdSparseQuadMatrix(SparseQuadMatrix *p,const SparseQuadMatrix *a)
{
	int *pInt = (int *)malloc(sizeof(int)*a->totalRow);
	int *ap = (int *)malloc(sizeof(int)*a->totalRow+1);
	int *ai = (int *)malloc(sizeof(int)*a->nnz);
	int i;
	QuadElement *unit = createQuadElement(a->gvNum);
	unit->m = 1;

//	dumpSparseQuadMatrix(stdout,a);

	sparseQuadMatrix2AmdAi(ai,a);
	sparseQuadMatrix2AmdAp(ap,a);
	amd_order(a->totalRow,ap,ai,pInt,NULL,NULL);
/*
	for(i=0;i<a->totalCol;i++)
	{
		fprintf(stderr,"%d\n",pInt[i]);
	}
*/
	clearSparseQuadMatrix(p);

	for(i=0;i<a->totalCol;i++)
	{
		setSparseQuadMatrix(p,unit,i,pInt[i]);
	}

	freeQuadElement(unit);
	free(pInt);
	free(ap);
	free(ai);

}



void partitionSparseQuadMatrix(SparseQuadMatrix *p,SparseQuadMatrix *pTrans,ParallelETree *tree,SparseQuadMatrix *aRefine,const SparseQuadMatrix *a,const int goalPartition)
{
	SparseDoubleMatrix *pDouble = createSparseDoubleMatrix(p->totalRow,p->totalCol);
	SparseDoubleMatrix *pTransDouble = createSparseDoubleMatrix(pTrans->totalRow,pTrans->totalCol);
	SparseDoubleMatrix *aDouble = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *aDoubleRefine = createSparseDoubleMatrix(a->totalRow,a->totalCol);

	mSparseQuadMatrix(aDouble,a);
	partitionSparseDoubleMatrix(pDouble,pTransDouble,tree,aDoubleRefine,aDouble,goalPartition);
	setMSparseQuadMatrix(p,pDouble);
	setMSparseQuadMatrix(pTrans,pTransDouble);
	
//	mulSparseQuadMatrix(aRefine,p,a);
//	mulSparseQuadMatrix(aRefine,aRefine,pTrans);
	permutateSparseQuadMatrix(aRefine,p,pTrans,a);	

	freeSparseDoubleMatrix(pDouble);
	freeSparseDoubleMatrix(pTransDouble);
	freeSparseDoubleMatrix(aDouble);
	freeSparseDoubleMatrix(aDoubleRefine);
}

