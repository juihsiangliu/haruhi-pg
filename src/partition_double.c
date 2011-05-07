#include "partition_double.h"


static void sparseDoubleMatrix2AmdAi(int *ai,const SparseDoubleMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	SparseDoubleElement *current;
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



static void sparseDoubleMatrix2AmdAp(int *ap,const SparseDoubleMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	int numInRow = 0;
	SparseDoubleElement *current;
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



static void sparseDoubleMatrix2MetisADJNCY(idxtype *ai,const SparseDoubleMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	const SparseDoubleElement *current;
	for(i=0;i<totalRow;i++)
	{
		current = a->rowIndex[i]->rowLink;
		while(current!=NULL)
		{
			if(current->col!=i)
			{
				ai[index] = current->col;
				index++;
			}
			current = current->rowLink;
		}
	}
}


/*
static void sparseDoubleMatrix2MetisXADJ(idxtype *ap,const SparseDoubleMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	int numInRow = 0;
	SparseDoubleElement *current;
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
	
	for(i=1;i<=totalRow;i++) ap[i] = ap[i] - i;
}
*/



static void sparseDoubleMatrix2MetisXADJ(idxtype *ap,const SparseDoubleMatrix *a)
{
	const int totalRow = a->totalRow;
	int i;
	int index = 0;
	int numInRow = 0;
	SparseDoubleElement *current;
	ap[0] = 0;
	for(i=0;i<totalRow;i++)
	{
		current = a->rowIndex[i]->rowLink;
		numInRow = 0;
		while(current!=NULL)
		{
			if(current->col!=i) numInRow++;
			current = current->rowLink;
		}
		ap[i+1] = ap[i] + numInRow;
	}	
	
}





static void postorder2Permutation(SparseDoubleMatrix *p, const int *orderList)
{
	int i;	
	clearSparseDoubleMatrix(p);
	// the permutation matrix
	for(i=0;i<p->totalRow;i++)
	{
		setSparseDoubleMatrix(p,1,i,orderList[i]);
	}
}




static PartitionResult* createPartitionResult(const int partitionSize)
{
	PartitionResult *ptr = getMempoolSet(sizeof(PartitionResult));
	ptr->partitionSize = partitionSize;
	ptr->partA = getMempoolSet(sizeof(int)*partitionSize);
	ptr->partB = getMempoolSet(sizeof(int)*partitionSize);

	ptr->partASize = 0;
	ptr->partBSize = 0;

	ptr->crossList = getMempoolSet(sizeof(int)*2*partitionSize);
	ptr->crossListSize = 0;

	ptr-> visitLog = notvisit;

	ptr->partSerialBegin = ptr->partSerialEnd = -1;
	return ptr;
}


static void freePartitionResult(PartitionResult *ptr)
{
	retMempoolSet(ptr->partA,sizeof(int)*ptr->partitionSize);
	retMempoolSet(ptr->partB,sizeof(int)*ptr->partitionSize);
	retMempoolSet(ptr->crossList,sizeof(int)*2*ptr->partitionSize);
	retMempoolSet(ptr,sizeof(PartitionResult));
}




static void insertToPartitionResult(PartitionResult *result, const int val, const char flag) // flag = 'A' or flag = 'B'
{
	switch(flag)
	{
		case 'A':
			if(result->partASize == result->partitionSize)
			{
				fprintf(stderr,"error in insertToPartitionResult\n");
				exit(0);
			}
			else
			{
				result->partA[result->partASize] = val;
				result->partASize++;
			}
			break;
		case 'B':
			if(result->partBSize == result->partitionSize)
			{
				fprintf(stderr,"error in insertToPartitionResult\n");
				exit(0);
			}
			else
			{
				result->partB[result->partBSize] = val;
				result->partBSize++;
			}
			break;
		case 'C':  // add to cross list
			if(result->crossListSize == 2*result->partitionSize)
			{
				fprintf(stderr,"error in insertToPartitionResult\n");
				exit(0);
			}
			else
			{
				result->crossList[result->crossListSize] = val;
				result->crossListSize++;
			}
			break;
		default:
			fprintf(stderr,"error in insertToPartitionResult, wrong flag\n");
			exit(0);
	}
}



static gdsl_element_t alloc_partitionResult(void *ptr)
{
	PartitionResult *n = (PartitionResult *) ptr;
	PartitionResult *value = createPartitionResult(n->partitionSize);

	// copy from n to value
	value->partASize = n->partASize;
	value->partBSize = n->partBSize;
	value->crossListSize = n->crossListSize;

	memcpy(value->partA,n->partA,n->partASize*sizeof(int));
	memcpy(value->partB,n->partB,n->partBSize*sizeof(int));
	memcpy(value->crossList,n->crossList,n->crossListSize*sizeof(int));

	value->partSerialBegin = n->partSerialBegin;
	value->partSerialEnd = n->partSerialEnd;

	return (gdsl_element_t) value;
}


static void free_partitionResult(gdsl_element_t e)
{
	freePartitionResult(e);
}



static PartitionResult *allocAndCopyPartitionResult(const  PartitionResult *src,const char flag)
{
	PartitionResult *dest = createPartitionResult(src->partitionSize);
	switch(flag)
	{
		case 'A':
			dest->partASize = src->partASize;
			memcpy(dest->partA,src->partA,src->partASize*sizeof(int));
			break;
		case 'B':
			dest->partBSize = src->partBSize;
			memcpy(dest->partB,src->partB,src->partBSize*sizeof(int));
			break;
		case 'C':
			dest->crossListSize = src->crossListSize;
			memcpy(dest->crossList,src->crossList,src->crossListSize*sizeof(int));
			break;
		default:
			dest->partASize = src->partASize;
			memcpy(dest->partA,src->partA,src->partASize*sizeof(int));
			dest->partBSize = src->partBSize;
			memcpy(dest->partB,src->partB,src->partBSize*sizeof(int));
			dest->crossListSize = src->crossListSize;
			memcpy(dest->crossList,src->crossList,src->crossListSize*sizeof(int));
	}

	return dest;
}



/*

// ==================================================  the utility in KL


// assuming ai and bj are exchanged, ai is original in partA and bj is original in partB
static void updataD
(int *D,const SparseDoubleMatrix *g,const int ai,const int bj,const PartitionResult* partition,const enum LockStat* lockTable)
{
	int i;

	// updatd Dx'
	for(i=0;i<partition->partASize;i++)
	{
		const int x =partition->partA[i];
		if(lockTable[x]!=lock)
		{
			const int cxai = getSparseDoubleMatrix(g,x,ai);  // Cx,ai
			const int cxbj = getSparseDoubleMatrix(g,x,bj);  // Cx,bj
			D[x] = D[x] + 2*cxai - 2*cxbj;
		}
	}

	// updatd Dy'
	for(i=0;i<partition->partBSize;i++)
	{
		const int y =partition->partB[i];
		if(lockTable[y]!=lock)
		{
			const int cybj = getSparseDoubleMatrix(g,y,bj);  // Cy,bj
			const int cyai = getSparseDoubleMatrix(g,y,ai);  // Cy,ai
			D[y] = D[y] + 2*cybj - 2*cyai;
		}
	}
}


// ==================================================

static int getCutSize(const PartitionResult *partition,const SparseDoubleMatrix *g)
{
	int cutSize = 0;
	int i,j;

	
	for(i=0;i<partition->partASize;i++)
	{
		const int indexA = partition->partA[i];
		for(j=0;j<partition->partBSize;j++)
		{
			const int indexB = partition->partB[j];
			const double connectivity = getSparseDoubleMatrix(g,indexA,indexB);
			if(connectivity != 0) cutSize++;
		}
	}

	return cutSize;
}
*/


// generate the crossList
// A' = A - crossList
// B' = B - crossList
static void refinePartitionResult(PartitionResult *result,const SparseDoubleMatrix *g)
{
	// init
	const int origASize = result->partASize;
	const int origBSize = result->partBSize;
	int *crossLogA = getMempoolSet(sizeof(int)*origASize);	
	int *crossLogB = getMempoolSet(sizeof(int)*origBSize);
	memset(crossLogA,0,sizeof(int)*origASize);
	memset(crossLogB,0,sizeof(int)*origBSize);
	int i,j;

	// calculate the cross terms
	//==============================================
	// A: 1 B: 2
	int *mark = getMempoolSet(sizeof(int)*g->totalRow);
	memset(mark,0,sizeof(int)*g->totalRow);
	int *indexLogA = getMempoolSet(sizeof(int)*g->totalRow);
	int *indexLogB = getMempoolSet(sizeof(int)*g->totalRow);
	for(i=0;i<g->totalRow;i++) indexLogA[i] = indexLogB[i] = -1;
	
	for(i=0;i<origASize;i++)
	{
		mark[result->partA[i]] = 1;
		indexLogA[result->partA[i]] = i;
	}
	for(i=0;i<origBSize;i++)
	{
		mark[result->partB[i]] = 2;
		indexLogB[result->partB[i]] = i;
	}


	for(i=0;i<g->totalRow;i++)
	{
		SparseDoubleElement *ptr = g->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
			if(mark[i] == 1 && mark[col] == 2)
			{
				crossLogA[indexLogA[i]] = 1;
				crossLogB[indexLogB[col]] = 1;
			}
			else if(mark[i] == 2 && mark[col] == 1)
			{
				crossLogA[indexLogA[col]] = 1;
				crossLogB[indexLogB[i]] = 1;
			}
			ptr = ptr->rowLink;
		}
	}
	retMempoolSet(mark,sizeof(int)*g->totalRow);
	retMempoolSet(indexLogA,sizeof(int)*g->totalRow);
	retMempoolSet(indexLogB,sizeof(int)*g->totalRow);

	//==============================================

	// output the refined result to temp
	PartitionResult *temp = createPartitionResult(result->partitionSize);
	for(i=0;i<origASize;i++)
	{
		const int indexA = result->partA[i];
		if(crossLogA[i] == 0) insertToPartitionResult(temp,indexA,'A'); 
		else insertToPartitionResult(temp,indexA,'C'); 
	}
	for(i=0;i<origBSize;i++)
	{
		const int indexB = result->partB[i];
		if(crossLogB[i] == 0) insertToPartitionResult(temp,indexB,'B'); 
		else insertToPartitionResult(temp,indexB,'C'); 
	}


	srand(0);
	// random shuffle
	for(i=0;i<4*temp->crossListSize;i++)
	{
		const int index1 = rand()%temp->crossListSize;
		const int index2 = rand()%temp->crossListSize;
		int tmp = temp->crossList[index1];
		temp->crossList[index1] = temp->crossList[index2];
		temp->crossList[index2] = tmp;
	}


	// temp stats
	// A = 1, B = 2, C = 3
	int *tempStat = getMempoolSet(sizeof(int)*g->totalCol);
	memset(tempStat,0,sizeof(int)*g->totalCol);
	for(i=0;i<temp->partASize;i++) tempStat[temp->partA[i]] = 1;
	for(i=0;i<temp->partBSize;i++) tempStat[temp->partB[i]] = 2;
	for(i=0;i<temp->crossListSize;i++) tempStat[temp->crossList[i]] = 3;

	time_t t1,t2;
	time(&t1);
//	printf("cross Size = %d\n",temp->crossListSize);
	for(i=0;i<temp->crossListSize;i++)
	{
		const int indexC = temp->crossList[i];
		int flag = 0;
		/*
		for(j=0;j<origASize;j++)
		{
			if(indexC == result->partA[j])
			{
				flag = 1;
				break;
			}
		}
		for(j=0;j<origBSize;j++)
		{
			if(indexC == result->partB[j])
			{
				flag = 2;
				break;
			}
		}
		*/
		flag = tempStat[indexC];

		assert(flag!=0);

		if(flag == 1) // cross term is original in partA
		{
			// check does it have connection with "tempPartB"
			// add to temp partA if it is ok
			// set temp crossList[i] = -1;
			int k;
			int connectFlag = 0;
			const SparseDoubleElement *ptr = g->rowIndex[indexC]->rowLink;
			while(ptr!=NULL)
			{
				const int col = ptr->col;
				if(tempStat[col] == 2)
				{
					connectFlag = 1;
					break;
				}
				ptr = ptr->rowLink;
			}
			
			if(connectFlag == 0)
			{
				insertToPartitionResult(temp,indexC,'A');
				tempStat[indexC] = 1;
				temp->crossList[i] = -1;
			}
		}
		else // flag == 2, cross term is original in partB
		{
			// check does it have connection with "tempPartA"
			// add to temp partB if it is ok
			// set temp crossList[i] = -1;
			int k;
			int connectFlag = 0;
			const SparseDoubleElement *ptr = g->rowIndex[indexC]->rowLink;
			while(ptr!=NULL)
			{
				const int col = ptr->col;
				if(tempStat[col] == 1)
				{
					connectFlag = 1;
					break;
				}
				ptr = ptr->rowLink;
			}
			
			if(connectFlag == 0)
			{
				insertToPartitionResult(temp,indexC,'B');
				tempStat[indexC] = 2;
				temp->crossList[i] = -1;
			}
		}
	}
	retMempoolSet(tempStat,sizeof(int)*g->totalCol);
	time(&t2);
//	fprintf(stderr,"test time:%g\n",difftime(t2,t1));



	// copy from temp to result
	result->partASize = temp->partASize;
	result->partBSize = temp->partBSize;
	result->crossListSize = temp->crossListSize;

	memcpy(result->partA,temp->partA,temp->partASize*sizeof(int));
	memcpy(result->partB,temp->partB,temp->partBSize*sizeof(int));
	memcpy(result->crossList,temp->crossList,temp->crossListSize*sizeof(int));

	
	// cross refinement
	j = 0;
	for(i=0;i<temp->crossListSize;i++)
	{
		if(temp->crossList[i] == -1)
		{
			result->crossListSize--;
		}
		else
		{
			result->crossList[j] = temp->crossList[i];
			j++;
		}
	}

	if(result->crossListSize == 0)
	{
		insertToPartitionResult(result,result->partA[result->partASize-1],'C');
		result->partASize--;
		insertToPartitionResult(result,result->partB[result->partBSize-1],'C');
		result->partBSize--;
	}



	// free
	freePartitionResult(temp);
	retMempoolSet(crossLogA,sizeof(int)*origASize);	
	retMempoolSet(crossLogB,sizeof(int)*origBSize);
}





static void refinePartitionResultBackup(PartitionResult *result,const SparseDoubleMatrix *g)
{
	// init
	const int origASize = result->partASize;
	const int origBSize = result->partBSize;
	int *crossLogA = getMempoolSet(sizeof(int)*origASize);	
	int *crossLogB = getMempoolSet(sizeof(int)*origBSize);
	memset(crossLogA,0,sizeof(int)*origASize);
	memset(crossLogB,0,sizeof(int)*origBSize);
	int i,j;

	// calculate the cross terms
	for(i=0;i<origASize;i++)
	{
		const int indexA = result->partA[i];
		for(j=0;j<origBSize;j++)
		{
			const int indexB = result->partB[j];
			const double connectivity = getSparseDoubleMatrix(g,indexA,indexB,"row");
			if(connectivity != 0)
			{
				crossLogA[i] = 1;
				crossLogB[j] = 1;
			}
		}
	}

	// output the refined result to temp
	PartitionResult *temp = createPartitionResult(result->partitionSize);
	for(i=0;i<origASize;i++)
	{
		const int indexA = result->partA[i];
		if(crossLogA[i] == 0) insertToPartitionResult(temp,indexA,'A'); 
		else insertToPartitionResult(temp,indexA,'C'); 
	}
	for(i=0;i<origBSize;i++)
	{
		const int indexB = result->partB[i];
		if(crossLogB[i] == 0) insertToPartitionResult(temp,indexB,'B'); 
		else insertToPartitionResult(temp,indexB,'C'); 
	}

	// copy from temp to result
	result->partASize = temp->partASize;
	result->partBSize = temp->partBSize;
	result->crossListSize = temp->crossListSize;

	memcpy(result->partA,temp->partA,temp->partASize*sizeof(int));
	memcpy(result->partB,temp->partB,temp->partBSize*sizeof(int));
	memcpy(result->crossList,temp->crossList,temp->crossListSize*sizeof(int));

	// free
	freePartitionResult(temp);
	retMempoolSet(crossLogA,sizeof(int)*origASize);	
	retMempoolSet(crossLogB,sizeof(int)*origBSize);
}






// ==================================================

/*
static void fastKL(PartitionResult *result, const SparseDoubleMatrix *g,const PartitionResult *init)
{
	
	// ===============   below is the KL body    ===================
	int i,j,k;
//	const int gSize = 2 * init->partitionSize;
	const int gSize = g->totalRow;

	int *I = getMempoolSet(sizeof(int)*gSize);
	int *E = getMempoolSet(sizeof(int)*gSize);
	int *D = getMempoolSet(sizeof(int)*gSize);
	memset(I,0,sizeof(int)*gSize);
	memset(E,0,sizeof(int)*gSize);
	memset(D,0,sizeof(int)*gSize);

	enum LockStat *lockTable = getMempoolSet(gSize*sizeof(enum LockStat));
	for(i=0;i<gSize;i++) lockTable[i] = unlock;	

	// calculate I
	for(i=0;i<init->partASize;i++)
	{
		const int indexSrc = init->partA[i];
		for(j=0;j<init->partASize;j++)
		{
			if(i==j) continue;
			else
			{
				const int indexDest = init->partA[j];
				const double connectivity = getSparseDoubleMatrix(g,indexSrc,indexDest);
				if(connectivity != 0) I[indexSrc]++;
			}
		}
	}
	for(i=0;i<init->partBSize;i++)
	{
		const int indexSrc = init->partB[i];
		for(j=0;j<init->partBSize;j++)
		{
			if(i==j) continue;
			else
			{
				const int indexDest = init->partB[j];
				const double connectivity = getSparseDoubleMatrix(g,indexSrc,indexDest);
				if(connectivity != 0) I[indexSrc]++;
			}
		}
	}

	// calculate E
	for(i=0;i<init->partASize;i++)
	{
		const int indexA = init->partA[i];
		for(j=0;j<init->partBSize;j++)
		{
			const int indexB = init->partB[j];
			const double connectivity = getSparseDoubleMatrix(g,indexA,indexB);
			if(connectivity != 0) E[indexA]++;
		}
	}
	for(i=0;i<init->partBSize;i++)
	{
		const int indexB = init->partB[i];
		for(j=0;j<init->partASize;j++)
		{
			const int indexA = init->partA[j];
			const double connectivity = getSparseDoubleMatrix(g,indexA,indexB);
			if(connectivity != 0) E[indexB]++;
		}
	}

	// calculate D
	for(i=0;i<gSize;i++) D[i] = E[i] - I[i];

	// actual swap
	while(1)
	{
		int gain = INT_MIN;
		int candidateA = -1;
		int candidateB = -1;
		for(j=0;j<init->partASize;j++)
		{
			const int valInPartA = init->partA[j];
			if(lockTable[valInPartA] == lock) continue;
			for(k=0;k<init->partBSize;k++)
			{
				const int valInPartB = init->partB[k];
				if(lockTable[valInPartB] == lock) continue;
				else
				{
					const int currentGain = D[valInPartA] + D[valInPartB] - 2*getSparseDoubleMatrix(g,valInPartA,valInPartB);
					if(currentGain>gain && currentGain>0)
					{
						gain = currentGain;
						candidateA = valInPartA;
						candidateB = valInPartB;
					}
				}
			}
		}
		if(gain > 0)
		{
			lockTable[candidateA] = lock;
			lockTable[candidateB] = lock;
			updataD(D,g,candidateA,candidateB,init,lockTable);
			insertToPartitionResult(result,candidateB,'A');
			insertToPartitionResult(result,candidateA,'B');
//			fprintf(stderr,"part a: %d\n",candidateB);
//			fprintf(stderr,"part b: %d\n",candidateA);
		}
		else
		{
			// insert the remaining unlock element to result
			for(i=0;i<init->partASize;i++)
			{
				const int valInPartA = init->partA[i];
				if(lockTable[valInPartA] == lock) continue;
				else
				{
//					fprintf(stderr,"part a: %d\n",valInPartA);
					insertToPartitionResult(result,valInPartA,'A');
				}
			}
			for(i=0;i<init->partBSize;i++)
			{
				const int valInPartB = init->partB[i];
				if(lockTable[valInPartB] == lock) continue;
				else
				{
//					fprintf(stderr,"part b: %d\n",valInPartB);
					insertToPartitionResult(result,valInPartB,'B');
				}
			}

			break;
		}
	}
	refinePartitionResult(result,g);
	// ===============   above is the KL body    ===================

//	for(i=0;i<result->crossListSize;i++) fprintf(stderr,"cross: %d\n",result->crossList[i]);
	
	retMempoolSet(I,sizeof(int)*gSize);
	retMempoolSet(E,sizeof(int)*gSize);
	retMempoolSet(D,sizeof(int)*gSize);
	retMempoolSet(lockTable,gSize*sizeof(enum LockStat));
}
*/

// ==================================================





static EliminationTree *createEliminationTree(const int size)
{
	EliminationTree *ptr = getMempoolSet(sizeof(EliminationTree));
	ptr->size = size;
	ptr->node = getMempoolSet(size*sizeof(PartitionResult *));
	memset(ptr->node,0,sizeof(PartitionResult *)*size);

	ptr->count = 1;

	return ptr;
}



static void freeEliminationTree(EliminationTree *ptr)
{
	int i;
	for(i=0;i<ptr->size;i++)
	{
		if(ptr->node[i]!=NULL) freePartitionResult(ptr->node[i]);
	}
	retMempoolSet(ptr->node,ptr->size*sizeof(PartitionResult *));
	retMempoolSet(ptr,sizeof(EliminationTree));
}



static void insertEliminationTree(EliminationTree *tree,const PartitionResult *node)
{
	if(tree->count == tree->size)
	{
		fprintf(stderr,"Elimination tree construction error\n");
		exit(0);
	}
	else
	{
		tree->node[tree->count] = createPartitionResult(node->partitionSize);
		tree->node[tree->count]->partASize = node->partASize;
		tree->node[tree->count]->partBSize = node->partBSize;
		tree->node[tree->count]->crossListSize = node->crossListSize;

		memcpy(tree->node[tree->count]->partA,node->partA,node->partASize*sizeof(int));
		memcpy(tree->node[tree->count]->partB,node->partB,node->partBSize*sizeof(int));
		memcpy(tree->node[tree->count]->crossList,node->crossList,node->crossListSize*sizeof(int));
		
		tree->count++;
	}
}



static void postOrder(int *orderList,const EliminationTree *tree)
{
	int orderListCounter = 0;
	int currentIndex = 1;
	PartitionResult* currentNode;
	while(1)
	{
//		fprintf(stderr,"currentIndex=%d\n",currentIndex);
		currentNode = tree->node[currentIndex];

		if(currentNode == NULL) break;
		
		const PartitionResult *currentLeft = tree->node[2*currentIndex];
		const PartitionResult *currentRight = tree->node[2*currentIndex+1];
		if(currentLeft!=NULL && currentLeft->visitLog != visit)
		{
			currentIndex = currentIndex*2;
		}
		else if(currentRight!=NULL && currentRight->visitLog != visit)
		{
			currentIndex = currentIndex*2 + 1;
		}
		else if(currentNode->visitLog!=visit) // print out the "N"
		{
			currentNode->visitLog = visit;
			if(currentNode->partASize!=0)
			{
				int j;
				for(j=0;j<currentNode->partASize;j++)
				{
					orderList[orderListCounter] = currentNode->partA[j];
					orderListCounter++;
//					fprintf(stderr,"%d ",currentNode->partA[j]);
				}
//				fprintf(stderr,"\n");
			}
			else if(currentNode->crossListSize!=0)
			{
				int j;
				for(j=0;j<currentNode->crossListSize;j++)
				{
					orderList[orderListCounter] = currentNode->crossList[j];
					orderListCounter++;
//					fprintf(stderr,"%d ",currentNode->crossList[j]);
				}
//				fprintf(stderr,"\n");
			}
			else
			{
//				fprintf(stderr,"error in tree\n");
//				exit(0);
			}
		}
		else
		{
			currentIndex = currentIndex/2;
		}
	}
}



/*
static void partition(EliminationTree *tree,const SparseDoubleMatrix *g,const int goalPartition)
{
	int i;
	SparseDoubleMatrix *gTrans = createSparseDoubleMatrix(g->totalRow,g->totalCol);
	SparseDoubleMatrix *gMap = createSparseDoubleMatrix(g->totalRow,g->totalCol);
	transSparseDoubleMatrix(gTrans,g);
	addSparseDoubleMatrix(gMap,g,gTrans);

	for(i=0;i<gMap->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = gMap->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			const int setRow = i;
			const int setCol = eachRow->col;
			const double val = 1.0;
			setSparseDoubleMatrix(gMap,val,setRow,setCol);
			eachRow = eachRow->rowLink;
		}
	}

	// initial partitinon
	PartitionResult *init = createPartitionResult((g->totalRow+1)/2);

	for(i=0;i<gMap->totalRow/2;i++)
	{
		insertToPartitionResult(init,i,'A');
	}
	for(i=gMap->totalRow/2;i<gMap->totalRow;i++)
	{
		insertToPartitionResult(init,i,'B');
	}

	PartitionResult *result = createPartitionResult((g->totalRow+1)/2 );	
	fprintf(stderr,"\tbefore fast KL\n");
	fastKL(result,gMap,init);
	fprintf(stderr,"\tend fast KL\n");
	
	PartitionResult *nodeToTree;
	nodeToTree = allocAndCopyPartitionResult(result,'C');
	insertEliminationTree(tree,nodeToTree);
	freePartitionResult(nodeToTree);

	gdsl_queue_t partitionResultQueue = gdsl_queue_alloc("partitionList",alloc_partitionResult,free_partitionResult); 

	PartitionResult *resultToQueueA = createPartitionResult(result->partASize);
	for(i=0;i<result->partASize;i++)
	{
		const int val = result->partA[i];
		insertToPartitionResult(resultToQueueA,val,'A');
	}
	gdsl_queue_insert(partitionResultQueue,resultToQueueA);
	freePartitionResult(resultToQueueA);

	PartitionResult *resultToQueueB = createPartitionResult(result->partBSize);
	for(i=0;i<result->partBSize;i++)
	{
		const int val = result->partB[i];
		insertToPartitionResult(resultToQueueB,val,'A');
	}
	gdsl_queue_insert(partitionResultQueue,resultToQueueB);
	freePartitionResult(resultToQueueB);

	int currentPartition;
	for(currentPartition=2;currentPartition<goalPartition && gdsl_queue_get_size(partitionResultQueue)>0 ;currentPartition++)
	{
		PartitionResult *currentNode = gdsl_queue_get_head(partitionResultQueue);
		PartitionResult *currentInit = createPartitionResult((currentNode->partASize+1)/2);
		PartitionResult *currentResult = createPartitionResult((currentNode->partASize+1)/2);

		for(i=0;i<currentNode->partASize/2;i++)
		{
			const int val = currentNode->partA[i];
			insertToPartitionResult(currentInit,val,'A');
		}
		for(i=currentNode->partASize/2;i<currentNode->partASize;i++)
		{
			const int val = currentNode->partA[i];
			insertToPartitionResult(currentInit,val,'B');
		}

		fastKL(currentResult,gMap,currentInit);
		
		// insert to elimination tree
		nodeToTree = allocAndCopyPartitionResult(currentResult,'C');
		insertEliminationTree(tree,nodeToTree);
		freePartitionResult(nodeToTree);

		// insert to partition queue
		PartitionResult *currentResultA = createPartitionResult(currentResult->partASize);
		for(i=0;i<currentResult->partASize;i++)
		{
			const int val = currentResult->partA[i];
			insertToPartitionResult(currentResultA,val,'A');
		}
		gdsl_queue_insert(partitionResultQueue,currentResultA);
		freePartitionResult(currentResultA);

		PartitionResult *currentResultB = createPartitionResult(currentResult->partBSize);
		for(i=0;i<currentResult->partBSize;i++)
		{
			const int val = currentResult->partB[i];
			insertToPartitionResult(currentResultB,val,'A');
		}
		gdsl_queue_insert(partitionResultQueue,currentResultB);
		freePartitionResult(currentResultB);
		
		gdsl_queue_remove(partitionResultQueue);

		freePartitionResult(currentResult);
		freePartitionResult(currentInit);
		freePartitionResult(currentNode);
	}	


	while(gdsl_queue_get_size(partitionResultQueue) > 0)
	{
		PartitionResult *currentNode = gdsl_queue_get_head(partitionResultQueue);
		
		// insert to elimination tree
		nodeToTree = allocAndCopyPartitionResult(currentNode,'A');
		insertEliminationTree(tree,nodeToTree);
		freePartitionResult(nodeToTree);

		gdsl_queue_remove(partitionResultQueue);
		freePartitionResult(currentNode);
	}


	gdsl_queue_free(partitionResultQueue);
	freePartitionResult(init);
	freeSparseDoubleMatrix(gTrans);
	freeSparseDoubleMatrix(gMap);
	freePartitionResult(result);
}
*/


static void paritionMetis(EliminationTree *tree, const SparseDoubleMatrix *a,const int goalPartition)
{
//	SparseDoubleMatrix *aTrans = createSparseDoubleMatrix(a->totalRow,a->totalCol);
//	SparseDoubleMatrix *aMap = createSparseDoubleMatrix(a->totalRow,a->totalCol);
//	transSparseDoubleMatrix(aTrans,a);
//	addSparseDoubleMatrix(aMap,a,aTrans);
//	freeSparseDoubleMatrix(aTrans);

	const SparseDoubleMatrix *aMap = a;

	int i,j;
	int n = a->totalRow;
	int wgtflag = 2;
	int numflag = 0;
	int nparts = goalPartition;
	int options[5] = {0,0,0,0,0};
	int edgecut = 0;
	int *part = getMempoolSet(sizeof(int)*n);

	time_t t1,t2;
	time(&t1);
	idxtype *xadj = getMempoolSet(sizeof(idxtype)*(a->totalRow+1));
	idxtype *adjncy = getMempoolSet(sizeof(idxtype)*(a->nnz-a->totalRow));
	idxtype *vwgt = getMempoolSet(sizeof(idxtype)*(a->totalRow));
	sparseDoubleMatrix2MetisXADJ(xadj,aMap);
	sparseDoubleMatrix2MetisADJNCY(adjncy,aMap);
	memset(vwgt,0,sizeof(idxtype)*a->totalRow);
	for(i=0;i<a->totalRow;i++) 
	{
		if(a->rowIndex[i]->rowLink->rowLink!=NULL) vwgt[i] = 1;
	}

//	fprintf(stderr,"metis pre done\n");
//	METIS_PartGraphRecursive(&n,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,&nparts,options,&edgecut,part);
	METIS_PartGraphKway(&n,xadj,adjncy,vwgt,NULL,&wgtflag,&numflag,&nparts,options,&edgecut,part);
//	fprintf(stderr,"metis kernel done\n");
	retMempoolSet(xadj,sizeof(idxtype)*(a->totalRow+1));
	retMempoolSet(adjncy,sizeof(idxtype)*(a->nnz-a->totalRow));
	retMempoolSet(vwgt,sizeof(idxtype)*(a->totalRow));
	time(&t2);
	fprintf(stderr,"kernel metis time:%g\n",difftime(t2,t1));

	gdsl_queue_t partitionResultQueue = gdsl_queue_alloc("partitionList",alloc_partitionResult,free_partitionResult); 

	PartitionResult* initPartitionResult = createPartitionResult(n);
	for(i=0;i<n;i++) insertToPartitionResult(initPartitionResult,i,'A');
	initPartitionResult->partSerialBegin = 0;
	initPartitionResult->partSerialEnd = goalPartition-1;
	gdsl_queue_insert(partitionResultQueue,initPartitionResult);
	freePartitionResult(initPartitionResult);


	for(i=0;i<goalPartition-1;i++)
	{
		PartitionResult* currentPartitionResult = gdsl_queue_get_head(partitionResultQueue);
		PartitionResult* currentToRefine = createPartitionResult(currentPartitionResult->partASize);
		const int partSerialBegin = currentPartitionResult->partSerialBegin;
		const int partSerialEnd = currentPartitionResult->partSerialEnd;
		const int separate = partSerialBegin + (partSerialEnd - partSerialBegin +1)/2;
/*
		fprintf(stderr,"begin:%d ,end:%d, mid:%d\n",partSerialBegin,partSerialEnd,separate);
		fprintf(stderr,"partASize:%d\n",currentPartitionResult->partASize);
		for(j=0;j<n;j++) fprintf(stderr,"%d ",part[j]);
		fprintf(stderr,"\n");
*/
		for(j=0;j<currentPartitionResult->partASize;j++)
		{
			const int data = currentPartitionResult->partA[j];
			if(part[data]<separate) // for part A
			{
				insertToPartitionResult(currentToRefine,data,'A'); // flag = 'A' or flag = 'B'
			}
			else if(part[data]>=separate) // for partB
			{
				insertToPartitionResult(currentToRefine,data,'B'); // flag = 'A' or flag = 'B'
			}
		}
		time_t t1,t2,t3;
		
		time(&t1);
		refinePartitionResult(currentToRefine,aMap);
		time(&t2);

		PartitionResult *nodeToTree;
		nodeToTree = allocAndCopyPartitionResult(currentToRefine,'C');
		insertEliminationTree(tree,nodeToTree);
		freePartitionResult(nodeToTree);
		time(&t3);
		
//		fprintf(stderr,"t2-t1:%g\tt3-t2:%g\n",difftime(t2,t1),difftime(t3,t2));

		PartitionResult *currentToRefineA = createPartitionResult(currentToRefine->partASize);
		for(j=0;j<currentToRefine->partASize;j++)
		{
			const int val = currentToRefine->partA[j];
			insertToPartitionResult(currentToRefineA,val,'A');
		}
		currentToRefineA->partSerialBegin = partSerialBegin;
		currentToRefineA->partSerialEnd = separate-1;
		gdsl_queue_insert(partitionResultQueue,currentToRefineA);
		freePartitionResult(currentToRefineA);

		PartitionResult *currentToRefineB = createPartitionResult(currentToRefine->partBSize);
		for(j=0;j<currentToRefine->partBSize;j++)
		{
			const int val = currentToRefine->partB[j];
			insertToPartitionResult(currentToRefineB,val,'A');
		}
		currentToRefineB->partSerialBegin = separate;
		currentToRefineB->partSerialEnd = partSerialEnd;
		gdsl_queue_insert(partitionResultQueue,currentToRefineB);
		freePartitionResult(currentToRefineB);

		gdsl_queue_remove(partitionResultQueue);
		freePartitionResult(currentPartitionResult);
		freePartitionResult(currentToRefine);
	}
//	freeSparseDoubleMatrix(aMap);
//	fprintf(stderr,".....................\n");

	for(i=0;i<goalPartition;i++)
	{
		PartitionResult* currentPartitionResult = gdsl_queue_get_head(partitionResultQueue);
		insertEliminationTree(tree,currentPartitionResult);
		gdsl_queue_remove(partitionResultQueue);
		freePartitionResult(currentPartitionResult);
	}
	gdsl_queue_free(partitionResultQueue);
		
	retMempoolSet(part,sizeof(int)*n);
}



static void getPermutation(SparseDoubleMatrix *p,const EliminationTree *tree)
{
	int i;
	// the order list
	int *orderList = getMempoolSet(p->totalRow*sizeof(int));
	postOrder(orderList,tree);
	postorder2Permutation(p,orderList);
	retMempoolSet(orderList,p->totalRow*sizeof(int));
}




static void renameingETree(EliminationTree *tree, const SparseDoubleMatrix *p)
{
	int i;
	// the tree renaming
	for(i=1;i<tree->size;i++)	
	{
		PartitionResult *currentNode = tree->node[i];
		if(currentNode==NULL) continue;
		if(currentNode->partASize!=0)
		{
			int j;
			for(j=0;j<currentNode->partASize;j++)
			{
				int oldVal = currentNode->partA[j];
				int newVal = p->colIndex[oldVal]->colLink->row;
				currentNode->partA[j] = newVal;
//				fprintf(stderr,"%d ",currentNode->partA[j]);
			}
//			fprintf(stderr,"\n");
		}
		else if(currentNode->crossListSize!=0)
		{
			int j;
			for(j=0;j<currentNode->crossListSize;j++)
			{
				int oldVal = currentNode->crossList[j];
				int newVal = p->colIndex[oldVal]->colLink->row;
				currentNode->crossList[j] = newVal;
//				fprintf(stderr,"%d ",currentNode->crossList[j]);
			}
//			fprintf(stderr,"\n");
		}
	}
}


static ParallelETreeNode *createParallelETreeNode(void)
{
	ParallelETreeNode *ptr = getMempoolSet(sizeof(ParallelETreeNode));
	ptr->rowBegin = ptr->rowEnd = -1;
	ptr->doneRowBegin = INT_MAX;
	ptr->doneRowEnd = INT_MIN;
	ptr->visitLog = notvisit; 
	ptr->type = undefine;
	return ptr;
}


static void freeParallelETreeNode(ParallelETreeNode *ptr)
{
	if(ptr!=NULL)
	{
		retMempoolSet(ptr,sizeof(ParallelETreeNode));
	}
}


ParallelETree *createParallelETree(const int size)
{
	ParallelETree *ptr = getMempoolSet(sizeof(ParallelETree));
	ptr->size = size;
	ptr->node = getMempoolSet(size*sizeof(ParallelETreeNode *));
	memset(ptr->node,0,sizeof(ParallelETreeNode *)*size);

	return ptr;
}


void freeParallelETree(ParallelETree *ptr)
{
	int i;
	for(i=0;i<ptr->size;i++)
	{
		freeParallelETreeNode(ptr->node[i]);
	}
	retMempoolSet(ptr->node,ptr->size*sizeof(ParallelETreeNode *));
	retMempoolSet(ptr,sizeof(ParallelETree));
}


static void copyFromEliminationTreeToParallelETree(ParallelETree *dest, const EliminationTree *src)
{
	if(dest->size != src->size)
	{
		fprintf(stderr,"convert from E-tree to Parall-E-tree error\n");
		exit(0);
	}
	int i,j;
	for(i=0;i<src->size;i++)
	{
		const PartitionResult *srcPtr = src->node[i];
		if(srcPtr==NULL) continue;
		else
		{
			dest->node[i] = createParallelETreeNode();
			if(srcPtr->partASize!=0)
			{
				dest->node[i]->rowBegin = srcPtr->partA[0];
				dest->node[i]->rowEnd = srcPtr->partA[srcPtr->partASize-1];
				dest->node[i]->type = lu;
			}
			else if(srcPtr->crossListSize!=0)
			{
				dest->node[i]->rowBegin = srcPtr->crossList[0];
				dest->node[i]->rowEnd = srcPtr->crossList[srcPtr->crossListSize-1];
				dest->node[i]->type = cross;
			}
		}
	}
}




static void roughPartition(ParallelETree *tree, SparseDoubleMatrix *p, const SparseDoubleMatrix *g, const int goalPartition)
{
	EliminationTree *treeTemp= createEliminationTree(goalPartition*4);

	time_t t1,t2;

//	fprintf(stderr,"rough partition - doing metis\n");
	time(&t1);
	paritionMetis(treeTemp,g,goalPartition);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));

//	fprintf(stderr,"partition phase 2\n");
	time(&t1);
	getPermutation(p,treeTemp);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));

//	fprintf(stderr,"partition phase 3\n");
	time(&t1);
	renameingETree(treeTemp,p);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));

//	fprintf(stderr,"partition phase 4\n");
	time(&t1);
	copyFromEliminationTreeToParallelETree(tree,treeTemp);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));
//	fprintf(stderr,"partition phase 5\n");
	time(&t1);
	freeEliminationTree(treeTemp);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));
//	fprintf(stderr,"rough partition end\n");
}


// assuming g is already permutated by "p" which is generated by roughPartition
static void refinePartition(SparseDoubleMatrix *pRefine,const ParallelETree *tree, const SparseDoubleMatrix *g)
{
	int i;

	const SparseDoubleMatrix *gMap = g;

	clearSparseDoubleMatrix(pRefine);

	// for all non-empty nodes in tree
	for(i=1;i<tree->size;i++)
	{
		if(tree->node[i]!=NULL)
		{
			const int rowBegin = tree->node[i]->rowBegin;
			const int rowEnd = tree->node[i]->rowEnd;
			const int num = rowEnd - rowBegin + 1;
			if(rowBegin == -1  && rowEnd == -1)
			{
				// no need to refine
			}
			else
			{
				SparseDoubleMatrix *subGMap = createSparseDoubleMatrix(num,num);
				getSubSparseDoubleMatrix(subGMap,gMap,rowBegin,rowBegin,rowEnd,rowEnd);
				int *pInt = getMempoolSet(sizeof(int)*subGMap->totalRow);
				int *ap = getMempoolSet(sizeof(int)*subGMap->totalRow+1);
				int *ai = getMempoolSet(sizeof(int)*subGMap->nnz);
				sparseDoubleMatrix2AmdAi(ai,subGMap);
				sparseDoubleMatrix2AmdAp(ap,subGMap);
				time_t t1,t2;
				time(&t1);
				amd_order(num,ap,ai,pInt,NULL,NULL);
				time(&t2);
//				fprintf(stderr,"kernel amd time:%g\n",difftime(t2,t1));

				int j;
				SparseDoubleMatrix *pCurrentNode = createSparseDoubleMatrix(num,num);
				postorder2Permutation(pCurrentNode,pInt);
				mergeSparseDoubleMatrix(pRefine,pCurrentNode,pRefine->totalRow,pRefine->totalCol,rowBegin,rowBegin);
				freeSparseDoubleMatrix(pCurrentNode);
	
				freeSparseDoubleMatrix(subGMap);
				retMempoolSet(pInt,sizeof(int)*subGMap->totalRow);
				retMempoolSet(ap,sizeof(int)*subGMap->totalRow+1);
				retMempoolSet(ai,sizeof(int)*subGMap->nnz);
			}
		}
	}
}



void amdSparseDoubleMatrix(SparseDoubleMatrix *p,const SparseDoubleMatrix *a)
{
	int *pInt = getMempoolSet(sizeof(int)*a->totalRow);
	int *ap = getMempoolSet(sizeof(int)*a->totalRow+1);
	int *ai = getMempoolSet(sizeof(int)*a->nnz);
	int i;
	const double unit = 1.0;

	sparseDoubleMatrix2AmdAi(ai,a);
	sparseDoubleMatrix2AmdAp(ap,a);
	const int aTotalRow = a->totalRow;

	amd_order(a->totalRow,ap,ai,pInt,NULL,NULL); 
	clearSparseDoubleMatrix(p);

	for(i=0;i<a->totalCol;i++)
	{
		setSparseDoubleMatrix(p,unit,i,pInt[i]);
	}

	retMempoolSet(pInt,sizeof(int)*a->totalRow);
	retMempoolSet(ap,sizeof(int)*a->totalRow+1);
	retMempoolSet(ai,sizeof(int)*a->nnz);
}



void static updatePermutationMatrix(SparseDoubleMatrix *pRefineP,const SparseDoubleMatrix *pRefine, const SparseDoubleMatrix *p)
{
	SparseDoubleMatrix *destTemp = createSparseDoubleMatrix(pRefineP->totalRow,pRefineP->totalCol);
	int i;
	// permutate row
	for(i=0;i<pRefine->totalRow;i++)
	{
		const int colForRowI = pRefine->rowIndex[i]->rowLink->col;
		const SparseDoubleElement *ptr = p->rowIndex[colForRowI]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
			setSparseDoubleMatrix(destTemp,ptr->data,i,col);			
			ptr = ptr->rowLink;
		}
	}
	
	copySparseDoubleMatrix(pRefineP,destTemp);
	freeSparseDoubleMatrix(destTemp);
}





SparseDoubleMatrix * partitionSparseDoubleMatrix(SparseDoubleMatrix *p,SparseDoubleMatrix *pTrans,ParallelETree *tree,const SparseDoubleMatrix *a,const int goalPartition,enum OOCFlag oocFlag)
{
	time_t t1,t2;
	const int pid = getpid();
	const int nodeNum = a->totalRow;

//	fprintf(stderr,"phase 1\n");
	time(&t1);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(nodeNum,nodeNum);
	roughPartition(tree,p,a,goalPartition);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));

//	fprintf(stderr,"phase 1.5\n");
	time(&t1);
	transSparseDoubleMatrix(pTrans,p);
	permutateSparseDoubleMatrix(aRefine,p,pTrans,a);	
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));
	
//	fprintf(stderr,"phase 2\n");
	time(&t1);
	SparseDoubleMatrix *pRefine = createSparseDoubleMatrix(nodeNum,nodeNum);
	refinePartition(pRefine,tree,aRefine);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));
	
//	fprintf(stderr,"phase 3\n");
	time(&t1);
	updatePermutationMatrix(p,pRefine,p);
	freeSparseDoubleMatrix(pRefine);
	transSparseDoubleMatrix(pTrans,p);
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));
	
//	fprintf(stderr,"phase 4\n");
	time(&t1);
	permutateSparseDoubleMatrix(aRefine,p,pTrans,a);	
	time(&t2);
//	fprintf(stderr,"time: %g\n",difftime(t2,t1));

	return aRefine;

/*
	int i;
	fprintf(stderr,"===========================\n");
	fprintf(stderr,"rowBegin, rowEnd setting\n");
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i]!=NULL)	fprintf(stderr,"i:%d,rowBegin:%d, rowEnd:%d\n",i,tree->node[i]->rowBegin,tree->node[i]->rowEnd); 	
	}
	fprintf(stderr,"===========================\n");
*/
}






