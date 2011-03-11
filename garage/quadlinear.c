#include "quadlinear.h"


// result should be nodeNum * stepNum
// result should be pre-allocated
void sparseQuadLinearSimulation(const SparseNetlistQuad *netlist,QuadMatrix *result,const int threadNum,const int dumpNodeIndex)
{
	QuadMatrix *uCurrent = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *vCurrent = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *zv = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *bu= createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *zvbu = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	QuadMatrix *vDC = createQuadMatrix(netlist->nodeNum,1,netlist->gvNum);
	SparseQuadMatrix *y = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *yRefine = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *z = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	// used for lu decompostion
	SparseQuadMatrix *p = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *pTrans = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *l = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	SparseQuadMatrix *u = createSparseQuadMatrix(netlist->nodeNum,netlist->nodeNum,netlist->gvNum);
	// z = c/dt
	scaleSparseQuadMatrix(z,1.0/netlist->deltaT,netlist->c);
	// y = z - a	
	subSparseQuadMatrix(y,z,netlist->a);

	time_t reorderBegin,reorderEnd,luBegin,luEnd;

	if(threadNum == 0)
	{
		time(&reorderBegin);
		amdSparseQuadMatrix(p,y);
		transSparseQuadMatrix(pTrans,p);
		permutateSparseQuadMatrix(yRefine,p,pTrans,y);
		time(&reorderEnd);
		fprintf(stderr,"reordering time: %g\n",difftime(reorderEnd,reorderBegin));
		time(&luBegin);
		luSparseQuadMatrix(l,u,yRefine);
		time(&luEnd);
		fprintf(stderr,"lu time:%g\n",difftime(luEnd,luBegin));
/*
		const int goalPartition = 8;
		ParallelETree *tree2 = createParallelETree(goalPartition*2 + goalPartition+1);
		time(&reorderBegin);
		partitionSparseQuadMatrix(p,pTrans,tree2,yRefine,y,goalPartition);
		time(&reorderEnd);
		fprintf(stderr,"reordering time: %g\n",difftime(reorderEnd,reorderBegin));
		time(&luBegin);
	    luSparseQuadMatrix(l,u,yRefine);
		time(&luEnd);
		fprintf(stderr,"lu time:%g\n",difftime(luEnd,luBegin));
		freeParallelETree(tree2);
	//	dumpSparseQuadMatrix(stdout,yRefine);
*/
	}
	else
	{
//		const int goalPartition = threadNum;
		const int goalPartition = 8;
//		ParallelETree *tree2 = createParallelETree(goalPartition*2 + goalPartition+1);
		ParallelETree *tree2 = createParallelETree(goalPartition*4);
		time(&reorderBegin);
		partitionSparseQuadMatrix(p,pTrans,tree2,yRefine,y,goalPartition);
		time(&reorderEnd);
		fprintf(stderr,"reordering time: %g\n",difftime(reorderEnd,reorderBegin));
		time(&luBegin);
	    parallelLUQuad(l,u,tree2,yRefine,threadNum);
		time(&luEnd);
		fprintf(stderr,"lu time:%g\n",difftime(luEnd,luBegin));
		fprintf(stderr,"row: %d, nnz: %d, L->nnz: %d, U->nnz: %d\n",yRefine->totalRow,yRefine->nnz,l->nnz,u->nnz);
		freeParallelETree(tree2);

	}


	// set the initial of v to result (stamp result[0])
	setZeroQuadMatrix(vCurrent);	

	int i;
	// predict the nodal voltage of i+1 step
	
	i = netlist->stepNum;

	for(i=0;i<netlist->stepNum-1;i++)
	{
//		fprintf(stderr,"step: %d\n",i);
		// z * vt
//		getColCopyQuadMatrix(vCurrent,i,result);
		mulVecSparseQuadMatrix(zv,z,vCurrent);
		// b *ut
		getColCopyQuadMatrix(uCurrent,i,netlist->u);
		mulVecSparseQuadMatrix(bu,netlist->b,uCurrent);
		// z*vt + b*ut
		addQuadMatrix(zvbu,zv,bu);
		triSolveSparseQuadMatrix(vCurrent,p,pTrans,l,u,zvbu);
		setQuadMatrix(result,getPtrEntryQuadMatrix(vCurrent,dumpNodeIndex,0),0,i+1);
//		setColQuadMatrix(result,vDC,i+1);
//		dumpQuadMatrix(vCurrent);
	}

	freeSparseQuadMatrix(p);
	freeSparseQuadMatrix(pTrans);
	freeSparseQuadMatrix(l);
	freeSparseQuadMatrix(u);
	freeSparseQuadMatrix(y);
	freeSparseQuadMatrix(yRefine);
	freeSparseQuadMatrix(z);
	freeQuadMatrix(uCurrent);
	freeQuadMatrix(vCurrent);
	freeQuadMatrix(zv);
	freeQuadMatrix(bu);
	freeQuadMatrix(zvbu);
	freeQuadMatrix(vDC);
}
