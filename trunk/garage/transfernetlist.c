#include "transfernetlist.h"


MonteNetlist* createMonteNetlist(const SparseNetlistQuad *src, const MonteSample *sample)
{
	int i,j;
	const SparseQuadElement *sourceElement;
	const SparseGControlInfo *sourceGControlInfo;

	MonteNetlist *dest = (MonteNetlist *)malloc(sizeof(MonteNetlist));

	dest->nonlinearNodeList = createAndSetNonlinearMonteList(src->nonlinearNodeList,sample);

	dest->a = createSparseDoubleMatrix(src->nodeNum,src->nodeNum);
	dest->b = createSparseDoubleMatrix(src->nodeNum,src->nodeNum);
	dest->c = createSparseDoubleMatrix(src->nodeNum,src->nodeNum);

	dest->u = (double **)malloc(sizeof(double *)*src->stepNum);
	for(i=0;i<src->stepNum;i++) dest->u[i] = getMempoolSet(sizeof(double)*src->nodeNum);

	dest->gVarientTable = createSparseGVarientTableMonte(src->nodeNum,src->nodeNum);


	// set a
	for(i=0;i<src->nodeNum;i++)
	{
		sourceElement = src->a->rowIndex[i]->rowLink;
		while(sourceElement != NULL)
		{
			const double val = quadElement2Sample(sourceElement->data,sample);
			setSparseDoubleMatrix(dest->a,val,sourceElement->row,sourceElement->col);
			sourceElement = sourceElement->rowLink;
		}
	}

	// set b
	for(i=0;i<src->nodeNum;i++)
	{
		sourceElement = src->b->rowIndex[i]->rowLink;
		while(sourceElement != NULL)
		{
			const double val = quadElement2Sample(sourceElement->data,sample);
			setSparseDoubleMatrix(dest->b,val,sourceElement->row,sourceElement->col);
			sourceElement = sourceElement->rowLink;
		}
	}


	// set c
	for(i=0;i<src->nodeNum;i++)
	{
		sourceElement = src->c->rowIndex[i]->rowLink;
		while(sourceElement != NULL)
		{
			const double val = quadElement2Sample(sourceElement->data,sample);
			setSparseDoubleMatrix(dest->c,val,sourceElement->row,sourceElement->col);
			sourceElement = sourceElement->rowLink;
		}
	}

	for(i=0;i<src->stepNum;i++)
	{
		for(j=0;j<src->nodeNum;j++)
		{
			const QuadElement *quadElement = getPtrEntryQuadMatrix(src->u,j,i);
			const double val = quadElement2Sample(quadElement,sample);
			setMyMatrix(dest->u[i],val,src->nodeNum,1,j,0);
		}
	}

	for(i=0;i<src->nodeNum;i++)
	{
		sourceGControlInfo = src->gVarientTable->rowIndex[i]->rowLink;
		while(sourceGControlInfo!=NULL)
		{
			int j;
			const gdsl_queue_t queue = sourceGControlInfo->queue;
			const int queueSize = gdsl_queue_get_size(queue);
			// operate the data in queue
			for(j=1;j<=queueSize;j++)
			{
				const GControlInfo *gControlInfo = gdsl_queue_search_by_position(queue,j);
				GControlInfoMonte *gControlInfoMonte = createGControlInfoMonte(gControlInfo->vdsListSize,gControlInfo->vgsListSize);
				setGControlInfoMonte(gControlInfoMonte,gControlInfo,dest->nonlinearNodeList);
				insertSparseGVarientTableMonte(dest->gVarientTable,gControlInfoMonte,sourceGControlInfo->row,sourceGControlInfo->col);
				freeGControlInfoMonte(gControlInfoMonte); // add on 2010/05/15 , not sure the accuracy
			}
			sourceGControlInfo = sourceGControlInfo->rowLink;
		}
	}

	dest->type = src->type;
	dest->nodeNum = src->nodeNum;
	dest->gvNum = src->gvNum;
	dest->endTime = src->endTime;
	dest->deltaT = src->deltaT;
	dest->stepNum = src->stepNum;

	return dest;
}



void freeMonteNetlist(MonteNetlist *monteNetlist)
{

	int i;
	freeSparseDoubleMatrix(monteNetlist->a);	
	freeSparseDoubleMatrix(monteNetlist->b);	
	freeSparseDoubleMatrix(monteNetlist->c);

	for(i=0;i<monteNetlist->stepNum;i++) retMempoolSet(monteNetlist->u[i],sizeof(double)*monteNetlist->nodeNum);
	free(monteNetlist->u);

	freeSparseGVarientTableMonte(monteNetlist->gVarientTable);

	freeNonlinearMonteList(monteNetlist->nonlinearNodeList);

	free(monteNetlist);
}
