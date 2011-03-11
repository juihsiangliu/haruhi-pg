#include "parallel_lu_double.h"


//static int flagEachNode[16]; // one based
static int *flagEachNode; // one based
// the node[i] is ready to be processed




static pthread_mutex_t donelist_mutex = PTHREAD_MUTEX_INITIALIZER;

static ParallelDoneListDouble *createDoneList(ALUDouble **alu, const ParallelETree *tree, const SparseDoubleMatrix *a, const int treeInternalSize,const int *postorder)
{
	ParallelDoneListDouble *ptr = getMempoolSet(sizeof(ParallelDoneListDouble));
	ptr->treeInternalSize = treeInternalSize;

	// used to check if all the necessary nodes are done in setNode function
	ptr->orderList = postorder;

	// done
	ptr->done = getMempoolSet(sizeof(int)*(ptr->treeInternalSize+1));
	memset(ptr->done,0,sizeof(int)*(ptr->treeInternalSize+1));

	// eachnodeCurrentDone
	ptr->eachNodeCurrentDone = getMempoolSet(sizeof(int)*(ptr->treeInternalSize+1));
	memset(ptr->eachNodeCurrentDone,0,sizeof(int)*(ptr->treeInternalSize+1));

	ptr->alu = alu;
	ptr->tree = tree;
	return ptr;
}




static void freeDoneList(ParallelDoneListDouble *ptr)
{
	retMempoolSet(ptr->eachNodeCurrentDone,sizeof(int)*(ptr->treeInternalSize+1));
	retMempoolSet(ptr->done,sizeof(int)*(ptr->treeInternalSize+1));
	retMempoolSet(ptr,sizeof(ParallelDoneListDouble));
}




static int checkNextParallelDoneList(ParallelDoneListDouble *ptr,const int current)
{
	int nextIndex = ptr->orderList[current+1];
	return ptr->done[nextIndex];
}




static int checkNextAndDec(ParallelDoneListDouble *ptr, const int current)
{
	int ret = checkNextParallelDoneList(ptr,current);
	if(ret !=1 ) decActiveThread();
	return ret;
}






static void set_1_sort_ParallelDoneList(ParallelDoneListDouble *ptr,ToDoList *todoList,const int index,const int *rootCurrentBeginList,const int maxThreadNum)
{
	pthread_mutex_lock(&donelist_mutex);
	ptr->done[index] = 1;

	if(maxThreadNum != 1)
	{
		if(index%2 == 1) // index is in the right
		{
			if(ptr->done[index-1] == 1) // left is complete
			{
				sortParentsToHeadToDoList(todoList,index,rootCurrentBeginList);
			}
		}
		else
		{
			sortParentsToHeadToDoList(todoList,index,rootCurrentBeginList);
		}
	}

	pthread_mutex_unlock(&donelist_mutex);
}





static int isCompleteParallelDoneList(ParallelDoneListDouble *doneList)
{
	pthread_mutex_lock(&donelist_mutex);
	int ret = 1;
	int i;
	for(i=1;i<doneList->treeInternalSize;i++)
	{
		if(doneList->done[i] == 0)
		{
			ret = 0;
			break;
		}
	}
	pthread_mutex_unlock(&donelist_mutex);
	return ret;
}



// =================================================
// N is tree node
// for each tree node N, they have their own baseRow and baseCol
// rowFix = *baseRow + row 
// colFox = *baseCol + col
static void getBaseOfPartialU(const SparseDoubleMatrix *u, const int srcOrder, const ParallelDoneListDouble *ptr, int *baseRow, int *baseCol, const int N)
{
	const int uRow = u->totalRow;
	const int uCol = u->totalCol;

	ALUDouble **aluList = ptr->alu;
	*baseRow = 0;
	*baseCol = 0;
	int i = 0;
	while(ptr->orderList[i]!=N)
	{
		const int order = ptr->orderList[i];

		*baseRow = *baseRow + aluList[order]->uRow;
		i++;
	}
	*baseCol = uCol - aluList[srcOrder]->uCol;
}






// N is tree node
// for each tree node N, they have their own baseRow and baseCol
// rowFix = *baseRow + row 
// colFox = *baseCol + col
static void getBaseOfPartialL(const SparseDoubleMatrix *l, const ParallelDoneListDouble *ptr, int *baseRow, int *baseCol, const int N)
{
	const int lRow = l->totalRow;
	const int lCol = l->totalCol;
	const int treeInternalSize = ptr->treeInternalSize;
	int i;
	ALUDouble **aluList = ptr->alu;
	
	int *eachCol = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=1;i<=treeInternalSize;i++) 
	{
		eachCol[i] = aluList[i]->lCol;
	}

	*baseRow = 0;
	*baseCol = 0;
	i = 0;
	while(ptr->orderList[i]!=N)
	{
		const int order = ptr->orderList[i];
		if(aluList[order]->lRow != aluList[order]->lCol) // the cross node
		{
			*baseCol = *baseCol - eachCol[order*2] - eachCol[order*2+1];
		}
		*baseRow = *baseRow + aluList[order]->lRow;
		*baseCol = *baseCol + aluList[order]->lCol;
		i++;
	}
	const int order = ptr->orderList[i];
	
	const SparseDoubleMatrix *partialL = aluList[order]->l;
	if(aluList[order]->lRow != aluList[order]->lCol) // the cross node
	{
		*baseCol = *baseCol - eachCol[order*2] - eachCol[order*2+1];
	}

	retMempoolSet(eachCol,sizeof(int)*(treeInternalSize+1));
}






// =================================================

static ParallelLUDoubleShareData *createParallelLUDoubleShareData(ParallelDoneListDouble *doneList, ToDoList *todolist, const int pid, const int N,const int currentEnd,pthread_mutex_t *mutex, pthread_cond_t *cond, const double tol,const int *baseRowL,const int *baseColL, const int *baseRowU, const int *baseColU,struct OOCInfo *oocInfoList,enum OOCFlag oocFlag,int *rootCurrentBeginList,int threadNum)
{
	ParallelLUDoubleShareData *ptr = getMempoolSet(sizeof(ParallelLUDoubleShareData));

	ptr->doneList = doneList;
	ptr->todolist = todolist;
	ptr->pid = pid;
	ptr->N = N;
	ptr->rootCurrentBegin = rootCurrentBeginList[N];
	ptr->currentEnd = currentEnd;
	ptr->mutex = mutex;
	ptr->cond = cond;
	ptr->tol = tol;
	ptr->maxThreadNum = threadNum;

	ptr->rootCurrentBeginList = rootCurrentBeginList;
	ptr->baseRowL = baseRowL;
	ptr->baseColL = baseColL;
	ptr->baseRowU = baseRowU;
	ptr->baseColU = baseColU;
	ptr->oocInfoList = oocInfoList;
	ptr->oocFlag = oocFlag;

	return ptr;
}




static void freeParallelLUDoubleShareData(ParallelLUDoubleShareData *ptr)
{
	retMempoolSet(ptr,sizeof(ParallelLUDoubleShareData));
}




// =================================================



static ALUDouble **createALUList(const ParallelETree *tree,const int totalRow,const int threadNum,enum OOCFlag oocFlag)
{
	int i;
	ALUDouble **aluList = getMempoolSet(sizeof(ALUDouble *)*tree->size);
	memset(aluList,0,sizeof(ALUDouble *)*tree->size);
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i] == NULL)
		{
			continue;
		}
		else
		{
			aluList[i] = getPidMempoolSet(sizeof(ALUDouble),i);
			int aRow = tree->node[i]->rowEnd - tree->node[i]->rowBegin + 1; // the row of cross node
			int aCol;
			if(tree->node[i]->type == cross)
			{
				const int leftIndex = 2*i;
				const int rightIndex = 2*i + 1;
				const int aRowNew = aRow + (tree->node[i]->doneRowEnd - tree->node[i]->doneRowBegin +1);
				if(i%2==0) // left node
				{
					int parentIndex = GSL_MAX(1,i/2);
					if( parentIndex == 1) // it is already root
					{
						aCol = totalRow;
					}
					else
					{
						aCol = aluList[parentIndex]->uCol;
					}
				}
				else // right node
				{
					aCol = totalRow - tree->node[i]->doneRowBegin;
				}

				aluList[i]->aRow = aRowNew;
				aluList[i]->aCol = aCol;
				aluList[i]->lRow = aRow;
				aluList[i]->lCol = aRowNew;
				aluList[i]->uRow = aRow;
				aluList[i]->uCol = aCol;
				aluList[i]->csr_u_count = 0;
				pthread_mutex_init(&aluList[i]->mutex,NULL);
			}
			else if(tree->node[i]->type == lu)
			{
				if(i%2 == 0) // left node
				{
					const int rootIndex = i/2;
					aCol = totalRow - tree->node[i]->rowBegin;
				}
				else
				{
					aRow = tree->node[i]->rowEnd - tree->node[i]->rowBegin + 1;
					const int rootIndex = i/2;
					const int leftIndex = i-1;
					const int leftRow = tree->node[leftIndex]->rowEnd - tree->node[leftIndex]->rowBegin + 1;
					const int rootCol = tree->node[rootIndex]->rowEnd - tree->node[rootIndex]->doneRowBegin + 1;
					aCol = totalRow - tree->node[i]->rowBegin;
				}
				
				aluList[i]->aRow = aRow;
				aluList[i]->aCol = aCol;
				aluList[i]->lRow = aRow;
				aluList[i]->lCol = aRow;
				aluList[i]->uRow = aRow;
				aluList[i]->uCol = aCol;
				aluList[i]->csr_u_count = 0;
				pthread_mutex_init(&aluList[i]->mutex,NULL);
			}
			else
			{
				fprintf(stderr,"in ALU list (%d): undefine\n",i);
				exit(0);
			}
//			fprintf(stderr,"i:%d, row:%d, col:%d\n",i,aluList[i]->l->totalRow,aluList[i]->l->totalCol);
//			fprintf(stderr,"i:%d, rowBegin:%d, rowEnd:%d, aRow:%d, aCol:%d\n",i,tree->node[i]->rowBegin,tree->node[i]->rowEnd,aRow,aCol);
		}
		clearPidMempoolSet(i);
	}
	return aluList;
}





static void freePidALU(ALUDouble *alu,const int pid)
{
	if(alu!=NULL)
	{
		pthread_mutex_destroy(&alu->mutex);
		if(alu->l!=NULL) freePidSparseDoubleMatrix(alu->l,pid);
		if(alu->u!=NULL) freePidSparseDoubleMatrix(alu->u,pid);
		retPidMempoolSet(alu,sizeof(ALUDouble),pid);
		clearPidMempoolSet(pid);
	}
}


// =================================================


static void storeLU_dec(ParallelLUDoubleShareData *ptr,const int treeNodeID)
{
	ALUDouble **alu = ptr->doneList->alu;
	if(ptr->maxThreadNum == 1)
	{
		writeByOOCInfo(ptr->oocInfoList,treeNodeID,'L',alu);
		writeByOOCInfo(ptr->oocInfoList,treeNodeID,'U',alu);
		decActiveThread();
	}
	else
	{
		decActiveThread();
		writeByOOCInfo(ptr->oocInfoList,treeNodeID,'L',alu);
		writeByOOCInfo(ptr->oocInfoList,treeNodeID,'U',alu);
	}
}


// =================================================


static void *threadLU(void *par)
{
	time_t t1,t2;
	ParallelLUDoubleShareData *ptr = ((ParallelLUDoubleShareData *)par);
	const int treeNodeID = ptr->N;	
	ALUDouble **alu = ptr->doneList->alu;
	const int pid = ptr->pid;

	flagEachNode[treeNodeID] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);

	alu[treeNodeID]->l = createPidSparseDoubleMatrix(alu[treeNodeID]->lRow,alu[treeNodeID]->lCol,pid);
	alu[treeNodeID]->u = createPidSparseDoubleMatrix(alu[treeNodeID]->uRow,alu[treeNodeID]->uCol,pid);

	if(ptr->oocFlag == ooc)
	{
		readByOOCInfo(ptr->oocInfoList,0,0,alu); // a
//		ptr->doneList->a = ptr->oocInfoList->a;
	}
	
//	fprintf(stderr,"lu begin: %d\n",treeNodeID);
	time(&t1);
	const int ltRowSrc = ptr->doneList->tree->node[treeNodeID]->rowBegin;
	const int ltColSrc = ltRowSrc;
	const int rbRowSrc = ptr->doneList->tree->node[treeNodeID]->rowEnd;
	const int rbColSrc = ltColSrc + alu[treeNodeID]->u->totalCol -1;
	
	SparseDoubleMatrix *a = createPidSparseDoubleMatrix(alu[treeNodeID]->u->totalRow,alu[treeNodeID]->u->totalCol,pid);
	getPidSubSparseDoubleMatrix(a,ptr->oocInfoList->a,ltRowSrc,ltColSrc,rbRowSrc,rbColSrc,pid);
	
	if(ptr->oocFlag == ooc)  pseudoWriteByOOCInfo(ptr->oocInfoList,0,0,alu); // a;

	iluPidSparseDoubleMatrix(alu[treeNodeID]->l,alu[treeNodeID]->u,a,pid,ptr->tol);

	freePidSparseDoubleMatrix(a,pid);
	time(&t2);
//	fprintf(stderr,"lu %d time:%g\n",treeNodeID,difftime(t2,t1));

	clearPidMempoolSet(pid);

	if(ptr->oocFlag == ooc)
	{
		storeLU_dec(ptr,treeNodeID);
	}
	else
	{
		decActiveThread();
	}

	set_1_sort_ParallelDoneList(ptr->doneList,ptr->todolist,treeNodeID,ptr->rootCurrentBeginList,ptr->maxThreadNum);
}



//====================================================================



// i and j is the original index of U in setNode (not offset fixed)
static inline int dropPartialUij(const double uij,const double tol,const SparseDoubleMatrix *partialU,const int i,const int j,const int *baseRowU, const int *baseColU,const int current,const double uii)
{
	if(tol == 0) return 0;
	const int fixI = baseRowU[current] + i;
	const int fixJ = baseColU[current] + j;
	if(fixI == fixJ) return 0;

//	const double uii = getSparseDoubleMatrix(partialU,i,i,"row");
//	if(fabs(uii) >= 1e6) return 0;

//	printf("current = %d, %d %d\n",current,fixI,fixJ);
	if( fabs(uij)/fabs(uii) >= tol) return 0;
	else return 1;
}



static inline int dropPartialUij2(const double uij,const double tol,const SparseDoubleMatrix *partialU,const int i,const int j,const int *baseRowU, const int *baseColU,const int current,const double uii)
{
//	return dropPartialUij(uij,tol,partialU,i,j,baseRowU,baseColU,current,uii);
	if(tol == 0) return 0;
	const int fixI = baseRowU[current] + i;
	const int fixJ = baseColU[current] + j;

//	if(fixI == fixJ) return 0;
	
	if( fabs(uij)/fabs(uii) >= tol) return 0;
	else return 1;

}



// i and j is the original index of L in setNode (not offset fixed)
// fixJ = baseColL[N] + j
//
// fixURowInd = baseRowU[current] + x
// fixUColInd = baseColU[current] + y
//
// to get U(fixJ,fixJ) <=> to solve x and y then search partialU(x,y)
// fixJ = fixURowInd = fixUColInd
static int dropPartialLij(const double lij,const double tol,const SparseDoubleMatrix *partialU,const int i,const int j,const int *baseRowL, const int *baseColL, const int *baseRowU, const int *baseColU, const int current,const int N,const double u_fixJ_fixJ)
{
	if(tol == 0) return 0;

	const int fixI = baseRowL[N] + i;
	const int fixJ = baseColL[N] + j;
	if(fixI == fixJ) return 0;

	const int x = fixJ - baseRowU[current];
	const int y = fixJ - baseColU[current];
//	const double u_fixJ_fixJ = getSparseDoubleMatrix(partialU,x,y,"row");
//	const double u_fixJ_fixJ = partialU->rowIndex[x]->rowLink->data;

//	if( fabs(u_fixJ_fixJ) >= 1e6) return 0;
	if( fabs(lij) / fabs(u_fixJ_fixJ) >= tol) return 0;
	else return 1;

}


static int dropPartialLij2(const double lij,const double tol,const SparseDoubleMatrix *partialU,const int i,const int j,const int *baseRowL, const int *baseColL, const int *baseRowU, const int *baseColU,const int N,const double u_fixJ_fixJ)
{
	return dropPartialLij(lij,tol,partialU,i,j,baseRowL,baseColL,baseRowU,baseColU,N,N,u_fixJ_fixJ);
}



//====================================================================

static void setNode_phase_1(const CSR_SparseDoubleMatrix *partial_u,SparseDoubleMatrix *u_link,SparseDoubleMatrix *l_link,int *base,const int offset,const int col_offset,const int N,const int pid,const int current,ParallelDoneListDouble *doneList,ParallelLUDoubleShareData *ptr,const double tol)
{
	const int *baseRowL = ptr->baseRowL;
	const int *baseColL = ptr->baseColL;
	const int *baseRowU = ptr->baseRowU;
	const int *baseColU = ptr->baseColU;

	int i,j;
	double scale = 0.0;
	double newLinkU = 0.0;
	double scaledPartialUData = 0.0;
	for(i=0;i<partial_u->totalRow;i++) // sweep partial_u row by row
	{
		for(j=0;j<u_link->totalRow;j++) // sweep u_link row by row
		{
			int partial_u_ptr = partial_u->rowPtr[i];
			scale = getSparseDoubleMatrix(u_link,j,*base+i,"col")/partial_u->val[partial_u_ptr];
			if(scale == 0.0)
			{
				delPidSparseDoubleMatrix(u_link,j,*base+i,pid);
				continue;
			}
//			if( fabs(scale/partial_u->val[partial_u_ptr]) >= tol)
				setPidSparseDoubleMatrix(l_link,scale,j,*base+partial_u->col[partial_u_ptr]+col_offset,pid);
			
//			const double u_drop_val = getSparseDoubleMatrix(u_link,row,baseRowU[N]+row-baseColU[N],"row");
			int k;
			for(k=partial_u->rowPtr[i];k<partial_u->rowPtr[i+1];k++)
			{
				const int row = j;
				const int col = *base + partial_u->col[k] + col_offset;
				const double data = partial_u->val[k];
				double linkU = getSparseDoubleMatrix(u_link,row,col,"row");
				scaledPartialUData = scale * data;
				newLinkU = linkU - scaledPartialUData;
//				if(!dropPartialUij(newLinkU,tol,u_link,row,col,baseRowU,baseColU,N)) 
					setPidSparseDoubleMatrix(u_link,newLinkU,row,col,pid);
			}
			delPidSparseDoubleMatrix(u_link,j,*base+i,pid);
		}
	}
	*base = *base + offset;
	doneList->eachNodeCurrentDone[N] = current;
}






static void setNode_phase_2(SparseDoubleMatrix *l_link,SparseDoubleMatrix *u_link,ParallelLUDoubleShareData *ptr,const int N,const int pid,const int base,const double tol)
{
	const int *baseRowL = ptr->baseRowL;
	const int *baseColL = ptr->baseColL;
	const int *baseRowU = ptr->baseRowU;
	const int *baseColU = ptr->baseColU;
	int i,j;
	
	double scale = 0.0;
	double newLinkU = 0.0;
	double scaledPartialUData = 0.0;
	for(i=0;i<u_link->totalRow;i++) // for each row in link (to be pivot row)
	{
		for(j=i+1;j<u_link->totalRow;j++) // for each row under the pivot row
		{
			const SparseDoubleElement *pivotRowPtr = u_link->rowIndex[i]->rowLink;
			scale = getSparseDoubleMatrix(u_link,j,base+i,"col") / pivotRowPtr->data;
			if(scale == 0.0)
			{
				delPidSparseDoubleMatrix(u_link,j,base + i,pid);
				continue;
			}
			if( fabs(scale/pivotRowPtr->data) >= tol)
				setPidSparseDoubleMatrix(l_link,scale,j,pivotRowPtr->col,pid);

			const int row = j;
			const double u_drop_val = getSparseDoubleMatrix(u_link,row,baseRowU[N]+row-baseColU[N],"row");
			while(pivotRowPtr!=NULL) // for each col in pivot row
			{
				const int col = pivotRowPtr->col;
				const double data = pivotRowPtr->data;
				double linkU = getSparseDoubleMatrix(u_link,row,col,"row");
				scaledPartialUData = scale * data;
				newLinkU = linkU - scaledPartialUData;
//				if( !dropPartialUij2(newLinkU,tol,u_link,row,col,baseRowU,baseColU,N,u_drop_val) )
				if( fabs(newLinkU/u_drop_val) >= tol ) 
					setPidSparseDoubleMatrix(u_link,newLinkU,row,col,pid);

				pivotRowPtr = pivotRowPtr->rowLink;
			}
			delPidSparseDoubleMatrix(u_link,j,base + i,pid);
		}
	}
	
}




static void set_L_diag(ParallelDoneListDouble *doneList,const int N,const int pid)
{
	int i;
	int j = 1;
	for(i=doneList->alu[N]->lRow-1;i>=0;i--)
	{
		const double element = 1.0;
		const int row = i;
		const int col = doneList->alu[N]->lCol - j;
		setPidSparseDoubleMatrix(doneList->alu[N]->l,element,row,col,pid);
		j++;
	}
}




static SparseDoubleMatrix* set_l_link(ParallelLUDoubleShareData *ptr)
{
	SparseDoubleMatrix *l_link = ptr->doneList->alu[ptr->N]->l;
	return l_link;
}




static SparseDoubleMatrix* set_u_link(ParallelLUDoubleShareData *ptr)
{
	ParallelDoneListDouble *doneList = ptr->doneList;

	const int rowCrossBegin = doneList->tree->node[ptr->N]->rowBegin;
	const int rowCrossEnd = doneList->tree->node[ptr->N]->rowEnd;
	const int currentDoneRowBegin = doneList->tree->node[ptr->N]->doneRowBegin;
	
	if(ptr->oocFlag == ooc)
	{
		readByOOCInfo(ptr->oocInfoList,0,0,doneList->alu); // a
//		doneList->a = ptr->oocInfoList->a;
	}

	SparseDoubleMatrix *u_link = doneList->alu[ptr->N]->u;
	getPidSubSparseDoubleMatrix(u_link,ptr->oocInfoList->a,rowCrossBegin,currentDoneRowBegin,rowCrossEnd,ptr->oocInfoList->a->totalCol-1,ptr->pid);

	if(ptr->oocFlag == ooc)  pseudoWriteByOOCInfo(ptr->oocInfoList,0,0,doneList->alu); // a;

	return u_link;
}





static void get_partialU_and_offset(ParallelLUDoubleShareData *ptr,int *offset,int *col_offset,const int current)
{
	*offset = 0; // how many rows of the cross term or the lu term
	*col_offset = 0;

	ParallelDoneListDouble *doneList = ptr->doneList;
	// get alu[current]->csr_u, offset, col_fooset
	if(current > doneList->treeInternalSize/2 )  // lu node
	{
		if(ptr->oocFlag == ooc) readByOOCInfo(ptr->oocInfoList,current,'U',doneList->alu);			
		else
		{
			pthread_mutex_lock(&doneList->alu[current]->mutex);
			if(doneList->alu[current]->csr_u_count == 0)
			{
				doneList->alu[current]->csr_u = sparse2CSR(doneList->alu[current]->u,ptr->pid);
				doneList->alu[current]->csr_u_count++;
			}
			else doneList->alu[current]->csr_u_count++;
			pthread_mutex_unlock(&doneList->alu[current]->mutex);
		}
		*offset = doneList->alu[current]->uRow;
	}
	else  // internal node
	{
		const int left = current*2;
		const int right = left+1;
		const int leftRow = doneList->alu[left]->aRow;
		const int rightRow = doneList->alu[right]->aRow;
		if(ptr->oocFlag == ooc) readByOOCInfo(ptr->oocInfoList,current,'U',doneList->alu);
		else
		{
			pthread_mutex_lock(&doneList->alu[current]->mutex);
			if(doneList->alu[current]->csr_u_count == 0)
			{
				doneList->alu[current]->csr_u = sparse2CSR(doneList->alu[current]->u,ptr->pid);
				doneList->alu[current]->csr_u_count++;
			}
			else doneList->alu[current]->csr_u_count++;
			pthread_mutex_unlock(&doneList->alu[current]->mutex);
		}
					
		*offset = doneList->alu[current]->uRow;
		*col_offset = -1 * (leftRow + rightRow);
	}
}




static void free_partialU(ParallelLUDoubleShareData *ptr,const int current)
{
	if(ptr->oocFlag == ooc) pseudoWriteByOOCInfo(ptr->oocInfoList,current,'U',ptr->doneList->alu);
	else
	{
		pthread_mutex_lock(&ptr->doneList->alu[current]->mutex);
		if(ptr->doneList->alu[current]->csr_u_count == 1)
		{
			free_CSR_SparseDoubleMatrix(ptr->doneList->alu[current]->csr_u);
			ptr->doneList->alu[current]->csr_u_count--;
		}
		else ptr->doneList->alu[current]->csr_u_count--;
		pthread_mutex_unlock(&ptr->doneList->alu[current]->mutex);
	}
}




inline static void wait_done(ParallelLUDoubleShareData *ptr,const int idle_counter,const int rootCurrent)
{
	if(idle_counter > 5) sleep(idle_counter);
	do
	{
		flagEachNode[ptr->N] = 1;
		pushBackToDoList(ptr->todolist,ptr->N);
		pthread_cond_wait(ptr->cond,ptr->mutex);
	}
	while(checkNextAndDec(ptr->doneList,rootCurrent)!=1);
}





static void *setNode(void *par)
{
	int i,j;
	time_t t1,t2,t3;
	time(&t1);
	ParallelLUDoubleShareData *ptr = (ParallelLUDoubleShareData *)par;
	ParallelDoneListDouble *doneList = ptr->doneList;
	flagEachNode[ptr->N] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);
	
//	fprintf(stderr,"setNode begin: %d\n",ptr->N);
	const int N = ptr->N;
	const int pid = ptr->pid;
	const double tol = ptr->tol;
	const double *colNormA = ptr->colNormA;
	
	doneList->alu[N]->l = createPidSparseDoubleMatrix(doneList->alu[N]->lRow,doneList->alu[N]->lCol,ptr->pid);
	doneList->alu[N]->u = createPidSparseDoubleMatrix(doneList->alu[N]->uRow,doneList->alu[N]->uCol,ptr->pid);

	set_L_diag(doneList,N,pid);
	
	SparseDoubleMatrix *l_link = set_l_link(ptr);
	SparseDoubleMatrix *u_link = set_u_link(ptr);

	int idle_counter = 0;
	int *rootCurrent = &(ptr->rootCurrentBeginList[N]);
	int current = -1;
	int base = 0;
	while(current != ptr->currentEnd)
	{
		if(checkNextParallelDoneList(doneList,*rootCurrent) == 1)
		{
			time(&t1);
			idle_counter = 0;
			*rootCurrent = *rootCurrent + 1;
			current = doneList->orderList[*rootCurrent];
//			fprintf(stderr,"setNode %d process node %d\n",N,current);
			int offset; // how many rows of the cross term or the lu term
			int col_offset;
			get_partialU_and_offset(ptr,&offset,&col_offset,current); // partial_u is stored in alu[current]->csr_u
			setNode_phase_1(doneList->alu[current]->csr_u,u_link,l_link,&base,offset,col_offset,N,pid,current,doneList,ptr,tol);
			free_partialU(ptr,current);
			time(&t2);
//			fprintf(stderr,"setnode %d process %d time:%g\n",N,current,difftime(t2,t1));
		}
		else
		{
			idle_counter++;
			if(ptr->oocFlag == ooc)
			{
//				fprintf(stderr,"node %d is going to sleep\n",N);
				storeLU_dec(ptr,ptr->N);
				wait_done(ptr,idle_counter,*rootCurrent);
				
				loadByOOCInfo(ptr->oocInfoList,N,'L',doneList->alu);
				loadByOOCInfo(ptr->oocInfoList,N,'U',doneList->alu);
				l_link = doneList->alu[N]->l;
				u_link = doneList->alu[N]->u;

//				fprintf(stderr,"node %d is continue to work\n",N);
				time(&t1);
			}
			else // ic
			{
//				fprintf(stderr,"node %d is going to sleep\n",N);
				decActiveThread();
				clearPidMempoolSet(pid);
				wait_done(ptr,idle_counter,*rootCurrent);
				
//				fprintf(stderr,"node %d is continue to work\n",N);
				time(&t1);
			}
		}
	}
	time(&t2);
//	fprintf(stderr,"set node %d skew time:%g\n",N,difftime(t2,t1));
//	fprintf(stderr,"node %d enter the final stage\n",N);

	setNode_phase_2(l_link,u_link,ptr,N,pid,base,tol);
	// ====================================================================
	time(&t3);
//	fprintf(stderr,"node %d final time: %g\n",N,difftime(t3,t2));

	if(ptr->oocFlag == ooc) storeLU_dec(ptr,N);
	else decActiveThread();

	clearPidMempoolSet(pid);
	doneList->eachNodeCurrentDone[ptr->N] = ptr->N;
	set_1_sort_ParallelDoneList(doneList,ptr->todolist,ptr->N,ptr->rootCurrentBeginList,ptr->maxThreadNum);
	pthread_exit(0);
}





// =================================================


static void threadHandlerNew(void *par)
{
	struct ThreadHandlerParDouble *thread_handle_ptr = par;

	ParallelLUDoubleShareData **list = thread_handle_ptr->list;
	ParallelDoneListDouble *doneList = list[0]->doneList;
	ToDoList *todolist = list[0]->todolist; 
	const int threadNum = thread_handle_ptr->threadNum;

//	dumpToDoList(stdout,todolist);
	initActiveThread();
	while(!isCompleteParallelDoneList(doneList))
	{
		while(getActiveThread() < threadNum)
		{
			const int indexInToDoList = getFirstToDoList(todolist);
			if(indexInToDoList == -1) // already empty
			{
				if(isCompleteParallelDoneList(doneList))
				{
					break;
				}
				else
				{
					sleep(1);
				}
			}
			else
			{
				safeWaitFlagMatrix(flagEachNode,indexInToDoList);
//				fprintf(stderr,"wake up node %d\n",indexInToDoList);	
				incActiveThread();
//				dumpToDoList(stderr,todolist);
				pthread_cond_signal(list[indexInToDoList-1]->cond);
			}
		}
		usleep(10000);
	}
}





//=====================================================================



static void assembleU(SparseDoubleMatrix *u, ParallelDoneListDouble *ptr, const int orderListSize, const int *baseRow, const int *baseCol,struct OOCInfo* oocInfoList,enum OOCFlag oocFlag)
{
	clearSparseDoubleMatrix(u);
	const int uRow = u->totalRow;
	const int uCol = u->totalCol;
	ALUDouble **aluList = ptr->alu;
	time_t t1,t2,t3,t4;
	
	int i;	
	for(i=0;i<orderListSize;i++)
	{
		if(oocFlag == ic)
		{
			const int order = ptr->orderList[i];
			const SparseDoubleMatrix *partialU = aluList[order]->u;
			mergeSparseDoubleMatrix(u,partialU,uRow,uCol,baseRow[order],baseCol[order]);
			freePidSparseDoubleMatrix(aluList[order]->u,order);
			aluList[order]->u = NULL;
			clearPidMempoolSet(order);
		}
		else
		{
			fprintf(stderr,"error at assemble U\n");
			exit(0);
		}
	}
}







static void assembleL(SparseDoubleMatrix *l, ParallelDoneListDouble *ptr, const int orderListSize, const int *baseRow, const int *baseCol,struct OOCInfo * oocInfoList,enum OOCFlag oocFlag)
{
	identitySparseDoubleMatrix(l);
	const int lRow = l->totalRow;
	const int lCol = l->totalCol;
	ALUDouble **aluList = ptr->alu;
	time_t t1,t2,t3,t4;

	int i;
	for(i=orderListSize-1;i>=0;i--)
	{
		if(oocFlag == ic)
		{
			const int order = ptr->orderList[i];
			const SparseDoubleMatrix *partialL = aluList[order]->l;
			mergeSparseDoubleMatrix(l,partialL,lRow,lCol,baseRow[order],baseCol[order]);	
			freePidSparseDoubleMatrix(aluList[order]->l,order);
			aluList[order]->l = NULL;
			clearPidMempoolSet(order);
		}
		else
		{
			fprintf(stderr,"error at assemble U\n");
			exit(0);
		}
	}
}










// =================================================


static int findBeforeBegin(const int *postorder_inv,int key,const int treeInternalSize)
{
	if(key > treeInternalSize/2)
	{
		return 0;
	}
	else
	{
		int i = key;
		while(key < treeInternalSize)
		{
			// the left child
			key = key * 2;
		}
		key = key / 2;
		return postorder_inv[key]-1;
	}
}


// =================================================


struct OOCInfo * parallelLUDouble(SparseDoubleMatrix *l,SparseDoubleMatrix *u, ParallelETree *tree, SparseDoubleMatrix *a, const SparseDoubleMatrix *p, const int threadNum,enum OOCFlag oocFlag)
{
	return parallelILUDouble(l,u,tree,a,p,threadNum,0.0,oocFlag);
}




struct OOCInfo * parallelILUDouble(SparseDoubleMatrix *l,SparseDoubleMatrix *u, ParallelETree *tree, SparseDoubleMatrix *a, const SparseDoubleMatrix *p, const int threadNum,const double tol,enum OOCFlag oocFlag)
{
	const int aRow = a->totalRow;
	const int aCol = a->totalCol;

	const int treeInternalSize = (tree->size/2)-1;
	// pre-processing
	int i;
	int status;
	
	// the super nodal tree & its postorder traversal
	int *postorder = getMempoolSet(sizeof(int)*treeInternalSize);
	int *arrayTree = getMempoolSet(sizeof(int)*tree->size);
	memset(arrayTree,0,sizeof(int)*tree->size);
	for(i=1;i<=treeInternalSize;i++) arrayTree[i] = i;
	arrayToPostOrder(postorder,arrayTree,treeInternalSize);
	retMempoolSet(arrayTree,sizeof(int)*tree->size);

	// postorder inverse table
	int *postorder_inv = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=0;i<treeInternalSize;i++)
	{	
		const int index = postorder[i];
		postorder_inv[index] = i;
	}


	// used to communcate the state of each threads
	flagEachNode = getMempoolSet(sizeof(int)*treeInternalSize+1);
	memset(flagEachNode,0,sizeof(int)*treeInternalSize+1);

	getDoneRowInfoNew(tree,postorder,postorder_inv,treeInternalSize);
	
	ALUDouble **aluList = createALUList(tree,aRow,threadNum,oocFlag);
	ParallelDoneListDouble *doneList = createDoneList(aluList,tree,a,treeInternalSize,postorder);

	ToDoList *todolist = createToDoList(treeInternalSize,postorder);

	// set the par entries ... used for the parallel lu
	// mutex and cv list for node
	pthread_mutex_t mutex[treeInternalSize];
	pthread_cond_t cond[treeInternalSize];
	for(i=0;i<treeInternalSize;i++)
	{
		pthread_mutex_init(&mutex[i],NULL);
		pthread_cond_init(&cond[i],NULL);
	}


	int *baseRowL = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	int *baseColL = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	int *baseRowU = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	int *baseColU = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=0;i<treeInternalSize;i++)
	{
		const int order = postorder[i];
		getBaseOfPartialL(l,doneList,&baseRowL[order],&baseColL[order],order);
		getBaseOfPartialU(u,order,doneList,&baseRowU[order],&baseColU[order],order);
	}

	// used in the ooc mode
	struct OOCInfo *oocInfoList = createOOCInfoList(aluList,treeInternalSize,postorder,oocFlag,baseRowL,baseColL,baseRowU,baseColU);
	oocInfoList->a = a;
	if(oocFlag == ooc)
	{
		write_SparseDoubleMatrix(oocInfoList[0].name,a); // a
		freeSparseDoubleMatrix(a);
		clearPidMempoolSet(0);
	}


	ParallelLUDoubleShareData **parList = getMempoolSet(sizeof(ParallelLUDoubleShareData *)*treeInternalSize);
	int *rootCurrentBeginList = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=1;i<=treeInternalSize;i++)
	{
		rootCurrentBeginList[i] = findBeforeBegin(postorder_inv,i,treeInternalSize);
	}
		
	for(i=0;i<treeInternalSize;i++)
	{
		if(i<treeInternalSize/2)
			parList[i] = createParallelLUDoubleShareData(doneList,todolist,i+1,i+1,2*(i+1)+1,&mutex[i],&cond[i],
			                        tol,baseRowL,baseColL,baseRowU,baseColU,oocInfoList,oocFlag,rootCurrentBeginList,threadNum);
		else parList[i] = createParallelLUDoubleShareData(doneList,todolist,i+1,i+1,0,&mutex[i],&cond[i],
		                              tol,baseRowL,baseColL,baseRowU,baseColU,oocInfoList,oocFlag,rootCurrentBeginList,threadNum);
	}

	// send out to each pthread
	pthread_t pid[treeInternalSize];
	for(i=0;i<treeInternalSize;i++)
	{
		if(i<treeInternalSize/2) pthread_create(&pid[i],NULL,setNode,parList[i]);
		else pthread_create(&pid[i],NULL,threadLU,parList[i]);
	}
	sleep(1);

	struct ThreadHandlerParDouble thread_handle;
	thread_handle.list = parList;
	thread_handle.threadNum = threadNum;
	threadHandlerNew(&thread_handle);

	for(i=0;i<treeInternalSize;i++)
	{
		pthread_join(pid[i],NULL);
	}
	
	if(oocFlag == ic)
	{
		// copy to the result
//		fprintf(stderr,"copy the result\n");
//		fprintf(stderr,"berfore the assemble\n");
		assembleL(l,doneList,treeInternalSize,baseRowL,baseColL,oocInfoList,oocFlag);
		assembleU(u,doneList,treeInternalSize,baseRowU,baseColU,oocInfoList,oocFlag);
//		fprintf(stderr,"after the assemble\n");
	}
	
	retMempoolSet(rootCurrentBeginList,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseRowL,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseColL,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseRowU,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseColU,sizeof(int)*(treeInternalSize+1));

	retMempoolSet(postorder,sizeof(int)*(treeInternalSize));
	retMempoolSet(postorder_inv,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(flagEachNode,sizeof(int)*(treeInternalSize+1));
	for(i=0;i<tree->size;i++) freePidALU(aluList[i],i);
	for(i=0;i<treeInternalSize;i++) freeParallelLUDoubleShareData(parList[i]);
	retMempoolSet(parList,treeInternalSize*sizeof(ParallelLUDoubleShareData*));
	retMempoolSet(aluList,sizeof(ALUDouble *)*tree->size);
	freeDoneList(doneList);
	freeToDoList(todolist);

	return oocInfoList;
}



//=================================

struct OOCInfo *createOOCInfoList(ALUDouble **aluList, const int treeInternalSize, const int *postorder,enum OOCFlag oocFlag,const int *baseRowL,const int *baseColL, const int *baseRowU, const int *baseColU)
{
	const int totalLength = 1 + 2*treeInternalSize;
	struct OOCInfo *ptr = getMempoolSet(sizeof(struct OOCInfo)*(totalLength));
	memset(ptr,0,sizeof(struct OOCInfo)*(1+2*treeInternalSize));

	const int pid = getpid();

	ptr[0].totalLength = totalLength;
	ptr[0].count = 0;
	ptr[0].postorder = getMempoolSet(sizeof(int)*treeInternalSize);
	memcpy(ptr[0].postorder,postorder,sizeof(int)*treeInternalSize);
	sprintf(ptr[0].name,"a.mtx.%d",pid);
	pthread_mutex_init(&ptr[0].mutex, NULL);
	ptr[0].baseRowL = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	memcpy(ptr[0].baseRowL,baseRowL,sizeof(int)*(treeInternalSize+1));
	ptr[0].baseColL = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	memcpy(ptr[0].baseColL,baseColL,sizeof(int)*(treeInternalSize+1));
	ptr[0].baseRowU = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	memcpy(ptr[0].baseRowU,baseRowU,sizeof(int)*(treeInternalSize+1));
	ptr[0].baseColU = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	memcpy(ptr[0].baseColU,baseColU,sizeof(int)*(treeInternalSize+1));
	int i;
	
	int rowBegin = 0;
	for(i=0;i<treeInternalSize;i++)
	{
		const int index = postorder[i];
		ptr[index].node = postorder[i];
		ptr[index].type = 'L';
		sprintf(ptr[index].name,"lnode%d.%d",postorder[i],pid);
		ptr[index].rowBegin = rowBegin;
		rowBegin += aluList[index]->lRow;
		ptr[index].rowEnd = rowBegin;
		ptr[index].count = 0;
		ptr[index].totalLength = totalLength;
		pthread_mutex_init(&ptr[index].mutex, NULL);
	}

	rowBegin = 0;
	for(i=0;i<treeInternalSize;i++)
	{
		const int index = treeInternalSize + postorder[i];
		ptr[index].node = postorder[i];
		ptr[index].type = 'U';
		sprintf(ptr[index].name,"unode%d.%d",postorder[i],pid);
		ptr[index].rowBegin = rowBegin;
		rowBegin += aluList[postorder[i]]->uRow;
		ptr[index].rowEnd = rowBegin;
		ptr[index].count = 0;
		ptr[index].totalLength = totalLength;
		pthread_mutex_init(&ptr[index].mutex, NULL);
	}

	return ptr;
}


void freeOOCInfoList(struct OOCInfo *ptr)
{
	const int totalLength = ptr[0].totalLength;
	int i;
	for(i=0;i<totalLength;i++) pthread_mutex_destroy(&ptr[i].mutex); // 0 is dummy
	const int treeInternalSize = (totalLength-1)/2;
	retMempoolSet(ptr[0].baseRowL,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(ptr[0].baseColL,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(ptr[0].baseRowU,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(ptr[0].baseColU,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(ptr[0].postorder,sizeof(int)*treeInternalSize);
	retMempoolSet(ptr,sizeof(struct OOCInfo)*totalLength);
}



void readByOOCInfo(struct OOCInfo *ptr,const int node,const char type, ALUDouble **aluList)
{	
	time_t t1,t2;
	time(&t1);
	int index;
	if(node == 0) index = 0;
	else if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;

	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 0)
	{
		if(index == 0)
		{
			ptr->a = read_pid_SparseDoubleMatrix(ptr[index].name,node);
		}
		else
		{
			if(type == 'L') aluList[node]->l = read_pid_SparseDoubleMatrix(ptr[index].name,node);
			else // U
			{
				aluList[node]->csr_u = linus_read_to_CSR_SparseDoubleMatrix(ptr[index].name,node); // read as csr format
			}
		}
//		fprintf(stderr,"log: read %d from file successfully.\n",node);
	}
	else if(ptr[index].count > 0)
	{
//		fprintf(stderr,"log: read %d ... already in mem.\n",node);
	}
	else
	{
		fprintf(stderr,"error: read %d from file failed\n",node);
		exit(0);
	}
	ptr[index].count++;
	pthread_mutex_unlock(&ptr[index].mutex);
	time(&t2);
//	fprintf(stderr,"node %d:read time %g\n",node,difftime(t2,t1));
}




// logically, this function is used when the matrix is not modified since it is read to memory
void pseudoWriteByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList)
{
	int index;
	if(node == 0) index = 0;
	else if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;
	
	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 1) // no others occupy this file
	{
		if(index == 0)
		{
			freePidSparseDoubleMatrix(ptr->a,node);
			ptr->a = NULL;
		}
		else
		{
			if(type == 'L')
			{
				freePidSparseDoubleMatrix(aluList[node]->l,node);
				aluList[node]->l = NULL;
			}
			else // U
			{
//				free_CSR_SparseDoubleMatrix(aluList[node]->csr_u);
				linus_free_CSR_SparseDoubleMatrix(aluList[node]->csr_u);
			}
		}
		clearPidMempoolSet(node);
//		fprintf(stderr,"log: pseudo write %d to file successfully.\n",node);
	}
	else if(ptr[index].count > 1)
	{
		// nothing ... just let it go ...
//		fprintf(stderr,"log: pseudo write counter = %d, skip write %d to file.\n",ptr[index].count,node);
	}
	else
	{
		fprintf(stderr,"error: pseudo write %d ... counter illegal....\n",node);
		exit(0);
	}
	ptr[index].count--;
	pthread_mutex_unlock(&ptr[index].mutex);
	
}




void loadByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList)
{	
	int index;
	if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;

	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 0)
	{
		if(type == 'L') aluList[node]->l = read_pid_SparseDoubleMatrix(ptr[index].name,node);
		else  aluList[node]->u = read_pid_SparseDoubleMatrix(ptr[index].name,node);
//		fprintf(stderr,"log: load %d from file successfully.\n",node);
	}
	else
	{
		fprintf(stderr,"error: load %d ... can not load the same file twice....\n",node);
		exit(0);
	}
	pthread_mutex_unlock(&ptr[index].mutex);
}






// the file will be overwritten by this function.
// typicallly, load <=> write , read <=> pseudoWrite
void writeByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList)
{
	time_t t1,t2;
	time(&t1);
	int index;
	if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;

	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 0) // no others occupy this file
	{
		if(type == 'L')
		{
			write_SparseDoubleMatrix(ptr[index].name,aluList[node]->l);
			freePidSparseDoubleMatrix(aluList[node]->l,node);
			aluList[node]->l = NULL;
		}
		else
		{
			write_SparseDoubleMatrix(ptr[index].name,aluList[node]->u);
			freePidSparseDoubleMatrix(aluList[node]->u,node);
			aluList[node]->u = NULL;
		}
		clearPidMempoolSet(node);
//		fprintf(stderr,"log: write %d to file successfully.\n",node);
	}
/*
	else if(ptr[index].count > 0)
	{
		fprintf(stderr,"log: when writing the matrix %d to file. alread read by other nodes.\n",node);
	}
*/	else // count != 0
	{
		fprintf(stderr,"error: write node %d failed\n",node);
	}
	pthread_mutex_unlock(&ptr[index].mutex);
	time(&t2);
//	fprintf(stderr,"node %d:write time %g\n",node,difftime(t2,t1));
}







static void shiftBase(CSR_SparseDoubleMatrix *mtx,const int newRow, const int newCol, const int baseRow,const int baseCol)
{
	// backup the old rowPtr
	int *rowPtr_backup = getPidMempoolSet(sizeof(int)*(mtx->totalRow+1),mtx->pid);
	memcpy(rowPtr_backup,mtx->rowPtr,sizeof(int)*(mtx->totalRow+1));

	// set new row
	retPidMempoolSet(mtx->rowPtr,sizeof(int)*(mtx->totalRow+1),mtx->pid);
	mtx->rowPtr = getPidMempoolSet(sizeof(int)*(newRow+1),mtx->pid);
	int i;
	for(i=0;i<baseRow;i++) mtx->rowPtr[i] = 0;
	memcpy(&mtx->rowPtr[i],rowPtr_backup,sizeof(int)*(mtx->totalRow+1));
	retPidMempoolSet(rowPtr_backup,sizeof(int)*(mtx->totalRow+1),mtx->pid);

	// set new col
	long long j;
	for(j=0;j<mtx->nnz;j++) mtx->col[j] += baseCol;

	// set the header
	mtx->totalRow = newRow;
	mtx->totalCol = newCol;
}







static void *free_mtx(void *par)
{
	PipeSolve *ptr = (PipeSolve *) par;
	int i,pid;
	*(ptr->freeFlag) = -1;
	time_t t1,t2;

	if(ptr->base == 0)
	{
		for(i=0;i<ptr->treeInternalSize;i++)
		{
			while(*(ptr->execFlag) <= *(ptr->freeFlag)) usleep(1000);
			time(&t1);
			pid = ptr->mtx_list[i]->pid;
			linus_free_CSR_SparseDoubleMatrix(ptr->mtx_list[i]);
//			free_CSR_SparseDoubleMatrix(ptr->mtx_list[i]);
			*(ptr->freeFlag) = i;
			clearPidMempoolSet(pid);
			time(&t2);
//			fprintf(stderr,"free %d time: %g\n",i,difftime(t2,t1));
		}
	}
	else
	{
		for(i=ptr->treeInternalSize-1;i>=0;i--)
		{
			while(*(ptr->execFlag) <= *(ptr->freeFlag)) usleep(1000);
			time(&t1);
			pid = ptr->mtx_list[i]->pid;
			linus_free_CSR_SparseDoubleMatrix(ptr->mtx_list[i]);
//			free_CSR_SparseDoubleMatrix(ptr->mtx_list[i]);
			*(ptr->freeFlag) = *(ptr->freeFlag) + 1;
			clearPidMempoolSet(pid);
			time(&t2);
//			fprintf(stderr,"free %d time: %g\n",i,difftime(t2,t1));
		}
	}

	pthread_exit(0);
}





static void *read_mtx(void *par)
{
	long long nnz = 0;
	PipeSolve *ptr = (PipeSolve *) par;
	int i;
	*(ptr->readFlag) = -1;
	time_t t1,t2;

	if(ptr->base == 0)
	{
		int postorderCounter = 0;
		for(i=0;i<ptr->treeInternalSize;i++)
		{
			while(*(ptr->readFlag) - *(ptr->freeFlag) >= 1) usleep(1000);
			time(&t1);
			ptr->mtx_list[i] = linus_read_to_CSR_SparseDoubleMatrix(ptr->oocInfoList[ptr->base+ptr->postorder[postorderCounter]].name,ptr->postorder[postorderCounter]);
//			ptr->mtx_list[i] = read_to_CSR_SparseDoubleMatrix(ptr->oocInfoList[ptr->base+ptr->postorder[postorderCounter]].name,ptr->postorder[postorderCounter]);
			postorderCounter++;
			*(ptr->readFlag) = i;
			time(&t2);
//			fprintf(stderr,"read %d time:%g\n",i,difftime(t2,t1));
			nnz += ptr->mtx_list[i]->nnz;
		}
//		fprintf(stderr,"nnz of factor = %lld\n",nnz);
	}
	else
	{
		int postorderCounter = ptr->treeInternalSize-1;
		for(i=ptr->treeInternalSize-1;i>=0;i--)
		{
			while(*(ptr->readFlag) - *(ptr->freeFlag) >= 1) usleep(1000);
			time(&t1);
			ptr->mtx_list[i] = linus_read_to_CSR_SparseDoubleMatrix(ptr->oocInfoList[ptr->base+ptr->postorder[postorderCounter]].name,ptr->postorder[postorderCounter]);
//			ptr->mtx_list[i] = read_to_CSR_SparseDoubleMatrix(ptr->oocInfoList[ptr->base+ptr->postorder[postorderCounter]].name,ptr->postorder[postorderCounter]);
			postorderCounter--;
			*(ptr->readFlag) = *(ptr->readFlag)+1;
			time(&t2);
//			fprintf(stderr,"read %d time:%g\n",i,difftime(t2,t1));
		}
	}

	pthread_exit(0);
}










static PipeSolve* createPipeSolve(int treeInternalSize,int *readFlag,int *shiftFlag,int *execFlag,int *freeFlag,const int *baseRow,const int*baseCol,const int *postorder,CSR_SparseDoubleMatrix **mtx_list, const OOCInfo *oocInfoList, const SparseDoubleMatrix *perm,const int base)
{
	PipeSolve *ptr = getMempoolSet(sizeof(PipeSolve));

	ptr->treeInternalSize = treeInternalSize;
	ptr->readFlag = readFlag;
	ptr->shiftFlag = shiftFlag;
	ptr->execFlag = execFlag;
	ptr->freeFlag = freeFlag;
	ptr->baseRow = baseRow;
	ptr->baseCol = baseCol;
	ptr->postorder = postorder;
	ptr->mtx_list = mtx_list;
	ptr->oocInfoList = oocInfoList;
	ptr->perm = perm;
	ptr->base = base;
	
	return ptr;
}





static void solve_Ly_b(double *y,struct OOCInfo *oocInfoList, const double *b,const SparseDoubleMatrix *p,const int *postorder,const int *baseRowL, const int *baseColL)
{
	time_t t1,t2;
	memset(y,0,sizeof(double)*p->totalRow);
	int i,j,k=0;
	double sum = 0.0;
	int postorderCounter = 0;
	const int treeInternalSize = (oocInfoList[0].totalLength-1)/2;
	int readFlag = -1;
	int shiftFlag = -1;
	int execFlag = -1;
	int freeFlag = -1;

	double *pb = getMempoolSet(sizeof(double)*p->totalRow);
	for(i=0;i<p->totalRow;i++)
	{
		const int col = p->rowIndex[i]->rowLink->col;
		pb[i] = b[col];
	}

	CSR_SparseDoubleMatrix **l_list = getMempoolSet(sizeof(CSR_SparseDoubleMatrix *)*treeInternalSize);
	// create thread
	PipeSolve* pipePar = createPipeSolve(treeInternalSize,&readFlag,&shiftFlag,&execFlag,&freeFlag,baseRowL,baseColL,postorder,l_list,oocInfoList,p,0);
	pthread_t pid[2];
	pthread_create(&pid[0],NULL,read_mtx,pipePar);
	pthread_create(&pid[1],NULL,free_mtx,pipePar);
	CSR_SparseDoubleMatrix *l;
	// solve ly = pb
	for(i=0;i<p->totalRow;i++)
	{
		if( oocInfoList[postorder[postorderCounter]].rowBegin == i)
		{
			while( readFlag <= execFlag) usleep(1000);
			time(&t1);
			l = l_list[k++];
		}
		sum = 0.0;
		const int forEnd = l->rowPtr[i+1-baseRowL[postorder[postorderCounter]]]-1;
		for( j=l->rowPtr[i-baseRowL[postorder[postorderCounter]]]; j<forEnd ; j++)
		{
			const double lij = l->val[j];
			const double yj = y[l->col[j]+baseColL[postorder[postorderCounter]]];
			sum += lij*yj;
		}
		y[i] = (pb[i]-sum)/l->val[j];
		if( oocInfoList[postorder[postorderCounter]].rowEnd == i+1)
		{
			postorderCounter++;	
			execFlag++;
			time(&t2);
//			fprintf(stderr,"solve %d time %g\n",postorderCounter-1,difftime(t2,t1));
		}
	}
	retMempoolSet(pipePar,sizeof(PipeSolve));
	retMempoolSet(l_list,sizeof(CSR_SparseDoubleMatrix *)*treeInternalSize);
	retMempoolSet(pb,sizeof(double)*p->totalRow);
	pthread_join(pid[0],NULL);
	pthread_join(pid[1],NULL);
}






static void solve_Ux_y(double *x,struct OOCInfo *oocInfoList, const double *y,const SparseDoubleMatrix *pTrans,const int *postorder,const int *baseRowU,const int *baseColU)
{
	time_t t1,t2;
	const int base = oocInfoList[0].totalLength/2; // used for u matrix
	const int treeInternalSize = (oocInfoList[0].totalLength-1)/2;
	memset(x,0,sizeof(double)*pTrans->totalRow);
	double sum = 0.0;
	int i,j,k=treeInternalSize-1;
	int postorderCounter = treeInternalSize - 1;
	int readFlag = -1;
	int shiftFlag = -1;
	int execFlag = -1;
	int freeFlag = -1;
	
	CSR_SparseDoubleMatrix **u_list = getMempoolSet(sizeof(CSR_SparseDoubleMatrix *)*treeInternalSize);
	// create thread
	PipeSolve* pipePar = createPipeSolve(treeInternalSize,&readFlag,&shiftFlag,&execFlag,&freeFlag,baseRowU,baseColU,postorder,u_list,oocInfoList,pTrans,base);
	pthread_t pid[2];
	pthread_create(&pid[0],NULL,read_mtx,pipePar);
	pthread_create(&pid[1],NULL,free_mtx,pipePar);
	CSR_SparseDoubleMatrix *u;
	double *xTemp = getMempoolSet(sizeof(double)*pTrans->totalRow);
	// solve ux = y
	for(i=pTrans->totalRow-1;i>=0;i--)
	{
		if( oocInfoList[base + postorder[postorderCounter]].rowEnd-1 == i)
		{
			while( readFlag <= execFlag) usleep(1000);
			time(&t1);
			u = u_list[k--];
		}
		sum = 0.0;
		const int forEnd = u->rowPtr[i+1-baseRowU[postorder[postorderCounter]]];
		for(j=u->rowPtr[i-baseRowU[postorder[postorderCounter]]]+1 ; j<forEnd ; j++)
		{
			const double uij = u->val[j];
			const double xj = xTemp[u->col[j]+baseColU[postorder[postorderCounter]]];
			sum += uij * xj;
		}
		xTemp[i] = (y[i]-sum) / u->val[u->rowPtr[i-baseRowU[postorder[postorderCounter]]]];
		if( oocInfoList[base + postorder[postorderCounter]].rowBegin == i)
		{
			postorderCounter--;
			execFlag++;
			time(&t2);
//			fprintf(stderr,"solve %d time %g\n",postorderCounter+1,difftime(t2,t1));
		}
	}

	for(i=0;i<pTrans->totalRow;i++)
	{
		int col = pTrans->rowIndex[i]->rowLink->col;
		x[i] = xTemp[col];
	}
	retMempoolSet(xTemp,sizeof(double)*pTrans->totalRow);
	retMempoolSet(u_list,sizeof(CSR_SparseDoubleMatrix *)*treeInternalSize);
	retMempoolSet(pipePar,sizeof(PipeSolve));
	pthread_join(pid[0],NULL);
	pthread_join(pid[1],NULL);
}









//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void oocTriSolveSparseDoubleMatrix(double *x, struct OOCInfo *oocInfoList, const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const double *b)
{
	const int *postorder = oocInfoList[0].postorder;
	const int *baseRowL = oocInfoList[0].baseRowL;
	const int *baseColL = oocInfoList[0].baseColL;
	const int *baseRowU = oocInfoList[0].baseRowU;
	const int *baseColU = oocInfoList[0].baseColU;

	double *y = getMempoolSet(sizeof(double)*p->totalCol);
	solve_Ly_b(y,oocInfoList,b,p,postorder,baseRowL,baseColL);
	solve_Ux_y(x,oocInfoList,y,pTrans,postorder,baseRowU,baseColU);
	retMempoolSet(y,sizeof(double)*p->totalCol);
}











