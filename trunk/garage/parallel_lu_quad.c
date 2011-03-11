#include "parallel_lu_quad.h"


static int flagEachNode[17]; // one based
// the node[i] is ready to be processed


// =================================================


static pthread_mutex_t donelist_mutex = PTHREAD_MUTEX_INITIALIZER;

static ParallelDoneList *createDoneList(ALUQuad **alu, const ParallelETree *tree, const SparseQuadMatrix *a,const int gvNum) // assume the # of lead of tree node is 8
{
	ParallelDoneList *ptr = getMempoolSet(sizeof(ParallelDoneList));

	// the order used for checking in setNode()
	ptr->orderList[0] = 8;
	ptr->orderList[1] = 9;
	ptr->orderList[2] = 4;
	ptr->orderList[3] = 10;
	ptr->orderList[4] = 11;
	ptr->orderList[5] = 5;
	ptr->orderList[6] = 2;
	ptr->orderList[7] = 12;
	ptr->orderList[8] = 13;
	ptr->orderList[9] = 6;
	ptr->orderList[10] = 14;
	ptr->orderList[11] = 15;
	ptr->orderList[12] = 7;
	ptr->orderList[13] = 3;
	ptr->orderList[14] = 1;

	memset(ptr->done,0,sizeof(int)*16);
	memset(ptr->eachNodeCurrentDone,0,sizeof(int)*16);

	ptr->gvNum = gvNum;
	ptr->alu = alu;
	ptr->tree = tree;
	ptr->a = a;
	return ptr;
}




static void freeParallelDoneList(ParallelDoneList *ptr)
{
	retMempoolSet(ptr,sizeof(ParallelDoneList));
}




static int checkNextParallelDoneList(ParallelDoneList *ptr,const int current)
{
	int nextIndex = ptr->orderList[current+1];
	return ptr->done[nextIndex];
}



static void set_1_sort_ParallelDoneList(ParallelDoneList *ptr,ToDoList *todoList,const int index)
{
	pthread_mutex_lock(&donelist_mutex);
	ptr->done[index] = 1;

	if(index%2 == 1) // index is in the right
	{
		if(ptr->done[index-1] == 1) // left is complete
		{
			sortParentsToHeadToDoList(todoList,index);
		}
	}
	else
	{
		sortParentsToHeadToDoList(todoList,index);
	}
	pthread_mutex_unlock(&donelist_mutex);
}




static int isCompleteParallelDoneList(ParallelDoneList *doneList)
{
	pthread_mutex_lock(&donelist_mutex);
	int ret = 1;
	int i;
	for(i=1;i<16;i++)
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



static ParallelLUQuadShareDataNew *createParallelLUQuadShareDataNew(ParallelDoneList *doneList, ToDoList *todolist, const int pid, const int N,const int rootCurrentBegin,const int currentEnd,pthread_mutex_t *mutex, pthread_cond_t *cond)
{
	ParallelLUQuadShareDataNew *ptr = getMempoolSet(sizeof(ParallelLUQuadShareDataNew));

	ptr->doneList = doneList;
	ptr->todolist = todolist;
	ptr->pid = pid;
	ptr->N = N;
	ptr->rootCurrentBegin = rootCurrentBegin;
	ptr->currentEnd = currentEnd;
	ptr->mutex = mutex;
	ptr->cond = cond;

	return ptr;
}



static void freeParallelLUQuadShareDataNew(ParallelLUQuadShareDataNew *ptr)
{
	retMempoolSet(ptr,sizeof(ParallelLUQuadShareDataNew));
}





// =================================================




static ALUQuad **createALUList(const ParallelETree *tree,const int totalRow,const int gvNum)
{
	int i;
	ALUQuad **aluList = getMempoolSet(sizeof(ALUQuad *)*tree->size);
	memset(aluList,0,sizeof(ALUQuad *)*tree->size);
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i] == NULL)
		{
			continue;
		}
		else
		{
			aluList[i] = getPidMempoolSet(sizeof(ALUQuad),i);
			int aRow = tree->node[i]->rowEnd - tree->node[i]->rowBegin + 1;
			int aCol;
			if(tree->node[i]->type == cross)
			{
				const int leftIndex = 2*i;
				const int rightIndex = 2*i + 1;
				aRow = aRow + (tree->node[i]->doneRowEnd - tree->node[i]->doneRowBegin +1);
				if(i%2==0) // left node
				{
					int parentIndex = GSL_MAX(1,i/2);
					if( parentIndex == 1) // it is already root
					{
						aCol = totalRow;
					}
					else
					{
						aCol = aluList[parentIndex]->a->totalCol;
					}
				}
				else // right node
				{
					aCol = totalRow - tree->node[i]->doneRowBegin;
				}

				aluList[i]->a = createPidSparseQuadMatrix(aRow,aCol,gvNum,i);
				aluList[i]->l = createPidSparseQuadMatrix(aRow,aRow,gvNum,i);
				aluList[i]->u = createPidSparseQuadMatrix(aRow,aCol,gvNum,i);
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
				
				aluList[i]->a = createPidSparseQuadMatrix(aRow,aCol,gvNum,i);
				aluList[i]->l = createPidSparseQuadMatrix(aRow,aRow,gvNum,i);
				aluList[i]->u = createPidSparseQuadMatrix(aRow,aCol,gvNum,i);
			}
			else
			{
				fprintf(stderr,"in ALU list: undefine\n");
				exit(0);
			}

//			fprintf(stderr,"i:%d, rowBegin:%d, rowEnd:%d, aRow:%d, aCol:%d\n",i,tree->node[i]->rowBegin,tree->node[i]->rowEnd,aRow,aCol);
			fprintf(stderr,"i:%d, uRow:%d, uCol:%d\n",i,aluList[i]->u->totalRow,aluList[i]->u->totalCol);
		}
	}
	return aluList;
}



static void freePidALU(ALUQuad *alu,const int pid)
{
	if(alu!=NULL)
	{
		freePidSparseQuadMatrix(alu->a,pid);
		freePidSparseQuadMatrix(alu->l,pid);
		freePidSparseQuadMatrix(alu->u,pid);
		retPidMempoolSet(alu,sizeof(ALUQuad),pid);
	}
}


// =================================================


static void *threadLU(void *par)
{
	time_t t1,t2;
	ParallelLUQuadShareDataNew *ptr = (ParallelLUQuadShareDataNew *)par;
	const int treeNodeID = ptr->N;	
	ALUQuad **alu = ptr->doneList->alu;

	flagEachNode[treeNodeID] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);

	fprintf(stderr,"lu begin: %d\n",treeNodeID);
	time(&t1);
	const int pid = ptr->pid;

	luPidSparseQuadMatrix(alu[treeNodeID]->l,alu[treeNodeID]->u,alu[treeNodeID]->a,pid);
	time(&t2);
	fprintf(stderr,"lu %d time:%g\n",treeNodeID,difftime(t2,t1));

	set_1_sort_ParallelDoneList(ptr->doneList,ptr->todolist,treeNodeID);
	decActiveThread();
}




static pthread_mutex_t setnode_merge_mutex = PTHREAD_MUTEX_INITIALIZER;


static void *setNode(void *par)
{
	int i,j;
	time_t t1,t2,t3;
	time(&t1);
	ParallelLUQuadShareDataNew *ptr = (ParallelLUQuadShareDataNew *)par;
	ParallelDoneList *doneList = ptr->doneList;
	ToDoList *todolist = ptr->todolist;

	flagEachNode[ptr->N] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);
//	fprintf(stderr,"setNode begin: %d\n",ptr->N);

	const int N = ptr->N;
	const int L = 2*N;
	const int R = L+1;
	const int rootCurrentBegin = ptr->rootCurrentBegin;
	const int currentEnd = ptr->currentEnd;;
	const int pid = ptr->pid;

	const int gvNum = doneList->gvNum;
	const int rowLink = doneList->tree->node[N]->rowEnd - doneList->tree->node[N]->rowBegin +1;
	const int colLink = doneList->alu[N]->u->totalCol;
	const int rowCrossBegin = doneList->tree->node[N]->rowBegin;
	const int rowCrossEnd = doneList->tree->node[N]->rowEnd;
	const int currentDoneRowBegin = doneList->tree->node[N]->doneRowBegin;
	
	SparseQuadMatrix *l_link = createPidSparseQuadMatrix(rowLink,doneList->alu[N]->l->totalCol,gvNum,pid);
	SparseQuadMatrix *u_link = createPidSparseQuadMatrix(rowLink,colLink,gvNum,pid);
	getSubPidSparseQuadMatrix(u_link,doneList->a,rowCrossBegin,currentDoneRowBegin,rowCrossEnd,doneList->a->totalCol-1,pid);

	SparseQuadElement **l_link_row_cache = getPidMempoolSet(sizeof(SparseQuadElement *)*l_link->totalRow,pid);
	SparseQuadElement **l_link_col_cache = getPidMempoolSet(sizeof(SparseQuadElement *)*l_link->totalCol,pid);
	SparseQuadElement **u_link_row_cache = getPidMempoolSet(sizeof(SparseQuadElement *)*u_link->totalRow,pid);
	SparseQuadElement **u_link_col_cache = getPidMempoolSet(sizeof(SparseQuadElement *)*u_link->totalCol,pid);
	for(i=0;i<l_link->totalRow;i++) l_link_row_cache[i] = NULL;
	for(i=0;i<l_link->totalCol;i++) l_link_col_cache[i] = NULL;
	for(i=0;i<u_link->totalRow;i++) u_link_row_cache[i] = NULL;
	for(i=0;i<u_link->totalCol;i++) u_link_col_cache[i] = NULL;

	int rootCurrent = rootCurrentBegin;
	int current = -1;
	int base = 0;
	while(current != currentEnd)
	{
		if(checkNextParallelDoneList(doneList,rootCurrent) == 1)
		{
			rootCurrent++;
			current = doneList->orderList[rootCurrent];
//			printf("setNode %d process node %d\n",N,current);
			const int left = current*2;
			const int right = left+1;
			int offset; // how many rows of the cross term or the lu term
			SparseQuadMatrix *partial_u;
			if(current > 7)
			{
				offset = doneList->alu[current]->u->totalRow;
				partial_u = doneList->alu[current]->u;
			}
			else
			{
				const int leftRow = doneList->alu[left]->u->totalRow;
				const int rightRow = doneList->alu[right]->u->totalRow;
				offset = doneList->alu[current]->u->totalRow - leftRow - rightRow;
				const int aTotal = doneList->a->totalCol;
				const SparseQuadMatrix *const uTemp = doneList->alu[current]->u;
				partial_u = createPidSparseQuadMatrix(uTemp->totalRow-leftRow-rightRow,doneList->a->totalCol-base,gvNum,pid);
				fprintf(stderr,"size of partial U - row:%d, col:%d, nnz:%d\n",partial_u->totalRow,partial_u->totalCol,partial_u->nnz);
				getSubPidSparseQuadMatrix(partial_u,uTemp,leftRow+rightRow, leftRow+rightRow,uTemp->totalRow-1,uTemp->totalCol-1,pid);
			}
//			printf("current = %d, offset =%d, row=%d\n",current,offset,doneList->alu[current]->u->totalRow);

			int i,j;
			QuadElement *scale = createPidQuadElement(gvNum,pid);
			QuadElement *newLinkU = createPidQuadElement(gvNum,pid);
			QuadElement *scaledPartialUData = createPidQuadElement(gvNum,pid);
			for(i=0;i<partial_u->totalRow;i++) // sweep partial_u row by row
			{
				for(j=0;j<u_link->totalRow;j++) // sweep u_link row by row
				{
					const SparseQuadElement *partial_u_ptr = partial_u->rowIndex[i]->rowLink;
					if(partial_u_ptr == NULL)
					{
						printf(" ------- WTF!!!!!! -----------\n");
					}
					divPidQuadElement(scale,getSparseQuadMatrix(u_link,j,base+i),partial_u_ptr->data,pid);
					if(isEmptyQuadElement(scale))
					{
						delPidSparseQuadMatrix(u_link,j,base + i,pid);
						continue;
					}
//					setPidSparseQuadMatrix(l_link,scale,j,base + partial_u_ptr->col,pid);
					setFastPidSparseQuadMatrix(l_link,scale,j,base + partial_u_ptr->col,&l_link_row_cache[j],&l_link_col_cache[partial_u_ptr->col],pid);
					
					SparseQuadElement *linkUPtr = NULL;
					while(partial_u_ptr!=NULL) // sweep each element in partial row i
					{
						const int row = j;
						const int col = base + partial_u_ptr->col;
						const QuadElement *const data = partial_u_ptr->data;
						const QuadElement *linkU = getFastRowSparseQuadMatrix(u_link,row,col,&linkUPtr);

						mulPidQuadElement(scaledPartialUData,scale,data,pid);
						subQuadElement(newLinkU,linkU,scaledPartialUData);
						setFastPidSparseQuadMatrix(u_link,newLinkU,row,col,&u_link_row_cache[row],&u_link_col_cache[col],pid);
						partial_u_ptr = partial_u_ptr->rowLink;
					}
					delPidSparseQuadMatrix(u_link,j,base + i,pid);
				}
			}
			freePidQuadElement(scaledPartialUData,pid);
			freePidQuadElement(newLinkU,pid);
			freePidQuadElement(scale,pid);
			
			base = base + offset;
			if(current <=7)
			{
				freePidSparseQuadMatrix(partial_u,pid);
			}
			doneList->eachNodeCurrentDone[N] = current;
//			printf("%d - currentDone:%d\n",N,ptr->currentDone);
		}
		else
		{
//			printf("node %d is going to sleep\n",N);
			decActiveThread();
			flagEachNode[N] = 1;
			pushBackToDoList(todolist,N);
//			pushByOrderToDoList(todolist,N);

			pthread_cond_wait(ptr->cond,ptr->mutex);
//			printf("node %d is continue to work\n",N);
			time(&t1);
		}
	}


	if(N==1)
	{
		for(i=4;i<doneList->tree->size;i++)
		{
			freePidALU(doneList->alu[i],i);
		}
	}

	time(&t2);
//	printf("set node %d skew time:%g\n",N,difftime(t2,t1));
//	printf("node %d enter the final stage\n",N);
//	printf("base = %d, row = %d, col =%d\n",base,u_link->totalRow,u_link->totalCol);
	
	QuadElement *scale = createPidQuadElement(gvNum,pid);
	QuadElement *newLinkU = createPidQuadElement(gvNum,pid);
	QuadElement *scaledPartialUData = createPidQuadElement(gvNum,pid);
	for(i=0;i<u_link->totalRow;i++) // for each row in link (to be pivot row)
	{
		for(j=i+1;j<u_link->totalRow;j++) // for each row under the pivot row
		{
			const SparseQuadElement *pivotRowPtr = u_link->rowIndex[i]->rowLink;
			divPidQuadElement(scale,getSparseQuadMatrix(u_link,j,base+i),pivotRowPtr->data,pid);
			if(isEmptyQuadElement(scale))
			{
				delPidSparseQuadMatrix(u_link,j,base + i,pid);
				continue;
			}
			setFastPidSparseQuadMatrix(l_link,scale,j,pivotRowPtr->col,&l_link_row_cache[j],&l_link_col_cache[pivotRowPtr->col],pid);

			SparseQuadElement *linkUPtr = NULL;
			while(pivotRowPtr!=NULL) // for each col in pivot row
			{
				const int row = j;
				const int col = pivotRowPtr->col;
				const QuadElement * const data = pivotRowPtr->data;
				const QuadElement *linkU = getFastRowSparseQuadMatrix(u_link,row,col,&linkUPtr);

				mulPidQuadElement(scaledPartialUData,scale,data,pid);
				subQuadElement(newLinkU,linkU,scaledPartialUData);
				setFastPidSparseQuadMatrix(u_link,newLinkU,row,col,&u_link_row_cache[row],&u_link_col_cache[col],pid);
				pivotRowPtr = pivotRowPtr->rowLink;
			}
			delPidSparseQuadMatrix(u_link,j,base + i,pid);
		}
	}
	freePidQuadElement(scaledPartialUData,pid);
	freePidQuadElement(newLinkU,pid);
	freePidQuadElement(scale,pid);
	
	retPidMempoolSet(l_link_row_cache,sizeof(SparseQuadElement *)*l_link->totalRow,pid);
	retPidMempoolSet(l_link_col_cache,sizeof(SparseQuadElement *)*l_link->totalCol,pid);
	retPidMempoolSet(u_link_row_cache,sizeof(SparseQuadElement *)*u_link->totalRow,pid);
	retPidMempoolSet(u_link_col_cache,sizeof(SparseQuadElement *)*u_link->totalCol,pid);


	// ====================================================================
	QuadElement *element = createPidQuadElement(gvNum,pid);
	setQuadElement(element,1.0,0,NULL,NULL);
	for(i=0;i<l_link->totalRow;i++)
	{
		const int row = i;
		const int col = doneList->alu[L]->l->totalRow + doneList->alu[R]->l->totalRow + i;
		setPidSparseQuadMatrix(l_link,element,row,col,pid);
	}
	freePidQuadElement(element,pid);

//	pthread_mutex_lock(&setnode_merge_mutex);
	mergePidSparseQuadMatrix(doneList->alu[N]->l,doneList->alu[L]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,0,0,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->l,doneList->alu[R]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow,doneList->alu[L]->l->totalRow,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->l,l_link,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow+doneList->alu[R]->l->totalRow,0,N);
//	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseQuadMatrix(l_link,pid);

//	pthread_mutex_lock(&setnode_merge_mutex);
	mergePidSparseQuadMatrix(doneList->alu[N]->u,doneList->alu[L]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,0,0,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->u,doneList->alu[R]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow,doneList->alu[L]->u->totalRow,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->u,u_link,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow+doneList->alu[R]->u->totalRow,0,N);
//	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseQuadMatrix(u_link,pid);
	
	if(N==1)
	{
		freePidALU(doneList->alu[2],2);
		freePidALU(doneList->alu[3],3);
	}
	
	// ====================================================================

	time(&t3);
//	printf("node %d final time: %g\n",N,difftime(t3,t2));

			
	set_1_sort_ParallelDoneList(doneList,ptr->todolist,N);
	decActiveThread();
	clearPidMempoolSet(pid);
	doneList->eachNodeCurrentDone[N] = N;
	pthread_exit(0);
}



// ======================================================================================



static void *setNodeMix(void *par)
{
	int i,j,k;
	time_t t1,t2,t3;
	time(&t1);
	ParallelLUQuadShareDataNew *ptr = (ParallelLUQuadShareDataNew *)par;
	ParallelDoneList *doneList = ptr->doneList;
	ToDoList *todolist = ptr->todolist;

	flagEachNode[ptr->N] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);
	fprintf(stderr,"setNode begin: %d\n",ptr->N);

	const int N = ptr->N;
	const int L = 2*N;
	const int R = L+1;
	const int rootCurrentBegin = ptr->rootCurrentBegin;
	const int currentEnd = ptr->currentEnd;;
	const int pid = ptr->pid;

	const int gvNum = doneList->gvNum;
	const int rowLink = doneList->tree->node[N]->rowEnd - doneList->tree->node[N]->rowBegin +1;
	const int colLink = doneList->alu[N]->u->totalCol;
	const int rowCrossBegin = doneList->tree->node[N]->rowBegin;
	const int rowCrossEnd = doneList->tree->node[N]->rowEnd;
	const int currentDoneRowBegin = doneList->tree->node[N]->doneRowBegin;

	SparseQuadMatrix *link = createSparseQuadMatrix(rowLink,colLink,gvNum);
	getSubSparseQuadMatrix(link,doneList->a,rowCrossBegin,currentDoneRowBegin,rowCrossEnd,doneList->a->totalCol-1);
	freeSparseQuadMatrix(link);
	
	SparseQuadMatrix *l_link = createPidSparseQuadMatrix(rowLink,doneList->alu[N]->l->totalCol,gvNum,pid);
	SparseQuadMatrix *u_link = createPidSparseQuadMatrix(rowLink,colLink,gvNum,pid);
	copyPidSparseQuadMatrix(u_link,link,pid);

	int rootCurrent = rootCurrentBegin;
	int current = -1;
	int base = 0;
	while(current != currentEnd)
	{
		if(checkNextParallelDoneList(doneList,rootCurrent) == 1)
		{
			rootCurrent++;
			current = doneList->orderList[rootCurrent];
			printf("setNode %d process node %d\n",N,current);
			const int left = current*2;
			const int right = left+1;
			int offset; // how many rows of the cross term or the lu term
			SparseQuadMatrix *partial_u;
			if(current > 7)
			{
				offset = doneList->alu[current]->u->totalRow;
				partial_u = doneList->alu[current]->u;
			}
			else
			{
				const int leftRow = doneList->alu[left]->u->totalRow;
				const int rightRow = doneList->alu[right]->u->totalRow;
				offset = doneList->alu[current]->u->totalRow - leftRow - rightRow;
				const int aTotal = doneList->a->totalCol;
				const SparseQuadMatrix *const uTemp = doneList->alu[current]->u;
				partial_u = createSparseQuadMatrix(uTemp->totalRow-leftRow-rightRow,doneList->a->totalCol-base,gvNum);
				getSubSparseQuadMatrix(partial_u,uTemp,leftRow+rightRow, leftRow+rightRow,uTemp->totalRow-1,uTemp->totalCol-1);
			}
//			printf("current = %d, offset =%d, row=%d\n",current,offset,doneList->alu[current]->u->totalRow);

			int i,j;
			QuadElement *scale = createPidQuadElement(gvNum,pid);
			QuadElement *newLinkU = createPidQuadElement(gvNum,pid);
			QuadElement *scaledPartialUData = createPidQuadElement(gvNum,pid);

			for(i=0;i<partial_u->totalRow;i++) // sweep partial_u row by row
			{
				for(j=0;j<u_link->totalRow;j++) // sweep u_link row by row
				{
					const SparseQuadElement *partial_u_ptr = partial_u->rowIndex[i]->rowLink;
					if(partial_u_ptr == NULL)
					{
						printf(" ------- WTF!!!!!! -----------\n");
					}
					divPidQuadElement(scale,getSparseQuadMatrix(u_link,j,base+i),partial_u_ptr->data,pid);
					if(isEmptyQuadElement(scale))
					{
						delPidSparseQuadMatrix(u_link,j,base + i,pid);
						continue;
					}
					setPidSparseQuadMatrix(l_link,scale,j,base + partial_u_ptr->col,pid);

					SparseQuadElement *linkUPtr = NULL;
					while(partial_u_ptr!=NULL) // sweep each element in partial row i
					{
						const int row = j;
						const int col = base + partial_u_ptr->col;
						const QuadElement *const data = partial_u_ptr->data;
//						const QuadElement *const linkU = getSparseQuadMatrix(u_link,row,col);
//
						linkUPtr = getPtrFastRowSparseQuadMatrix(u_link,row,col,&linkUPtr);
						QuadElement *linkU = NULL;
						if(linkUPtr!=NULL) linkU = linkUPtr->data;

						mulPidQuadElement(scaledPartialUData,scale,data,pid);
						subQuadElement(newLinkU,linkU,scaledPartialUData);
//						setPidSparseQuadMatrix(u_link,newLinkU,row,col,pid);
						setFastPidSparseQuadMatrix(u_link,newLinkU,row,col,&linkUPtr,&linkUPtr,pid);
						partial_u_ptr = partial_u_ptr->rowLink;
					}
					delPidSparseQuadMatrix(u_link,j,base + i,pid);
				}
			}
			freePidQuadElement(scaledPartialUData,pid);
			freePidQuadElement(newLinkU,pid);
			freePidQuadElement(scale,pid);
			
			base = base + offset;
			if(current <=7)
			{
				freeSparseQuadMatrix(partial_u);
			}
		}
		else
		{
			printf("node %d is going to sleep\n",N);
			decActiveThread();
			flagEachNode[N] = 1;
//			pushBackToDoList(todolist,N);
			pushByOrderToDoList(todolist,N);

			pthread_cond_wait(ptr->cond,ptr->mutex);
			printf("node %d is continue to work\n",N);
			time(&t1);
		}
	}

	time(&t2);
	printf("set node %d skew time:%g\n",N,difftime(t2,t1));
	printf("node %d enter the final stage\n",N);

	QuadMatrix *l_link_dense = createPidQuadMatrix(rowLink,doneList->alu[N]->l->totalCol,gvNum,pid);
	QuadMatrix *u_link_dense = createPidQuadMatrix(rowLink,colLink,gvNum,pid);
	toDensePidSparseQuadMatrix(l_link_dense,l_link,pid);
	toDensePidSparseQuadMatrix(u_link_dense,u_link,pid);

	QuadElement *scale = createPidQuadElement(gvNum,pid);
	QuadElement *newLinkU = createPidQuadElement(gvNum,pid);
	QuadElement *scaledPartialUData = createPidQuadElement(gvNum,pid);
	for(i=0;i<u_link_dense->row;i++) // for each row in link (to be pivot row)
	{
		for(j=i+1;j<u_link_dense->row;j++) // for each row under the pivot row
		{
			const QuadElement *const pivot = getPtrEntryQuadMatrix(u_link_dense,i,base+i);	
			const QuadElement *const eachRowHeadData = getPtrEntryQuadMatrix(u_link_dense,j,base+i);
			divPidQuadElement(scale,eachRowHeadData,pivot,pid);
			if(isEmptyQuadElement(scale))
			{
				delPidQuadMatrix(u_link_dense,j,base + i,pid);
				continue;
			}
			setPidQuadMatrix(l_link_dense,scale,j,base+i,pid);
			for(k=base;k<u_link_dense->col;k++)
			{
				const QuadElement * const data1 = getPtrEntryQuadMatrix(u_link_dense,i,k);
				const QuadElement * const data2 = getPtrEntryQuadMatrix(u_link_dense,j,k);
				mulPidQuadElement(scaledPartialUData,scale,data1,pid);
				subQuadElement(newLinkU,data2,scaledPartialUData);
				setPidQuadMatrix(u_link_dense,newLinkU,j,k,pid);
			}
			delPidQuadMatrix(u_link_dense,j,base+i,pid);
		}
	}
	freePidQuadElement(scaledPartialUData,pid);
	freePidQuadElement(newLinkU,pid);
	freePidQuadElement(scale,pid);

	quad2PidSparseQuadMatrix(l_link,l_link_dense,pid);
	quad2PidSparseQuadMatrix(u_link,u_link_dense,pid);
	freePidQuadMatrix(l_link_dense,pid);
	freePidQuadMatrix(u_link_dense,pid);

	// ====================================================================
	pthread_mutex_lock(&setnode_merge_mutex);
	QuadElement *element = createPidQuadElement(gvNum,pid);
	setQuadElement(element,1.0,0,NULL,NULL);
	for(i=0;i<l_link->totalRow;i++)
	{
		const int row = i;
		const int col = doneList->alu[L]->l->totalRow + doneList->alu[R]->l->totalRow + i;
		setPidSparseQuadMatrix(l_link,element,row,col,pid);
	}
	freePidQuadElement(element,pid);

	mergeSparseQuadMatrix(doneList->alu[N]->l,doneList->alu[L]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,0,0);
	mergeSparseQuadMatrix(doneList->alu[N]->l,doneList->alu[R]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow,doneList->alu[L]->l->totalRow);
	mergeSparseQuadMatrix(doneList->alu[N]->l,l_link,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow+doneList->alu[R]->l->totalRow,0);
	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseQuadMatrix(l_link,pid);

	pthread_mutex_lock(&setnode_merge_mutex);
	mergeSparseQuadMatrix(doneList->alu[N]->u,doneList->alu[L]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,0,0);
	mergeSparseQuadMatrix(doneList->alu[N]->u,doneList->alu[R]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow,doneList->alu[L]->u->totalRow);
	mergeSparseQuadMatrix(doneList->alu[N]->u,u_link,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow+doneList->alu[R]->u->totalRow,0);
	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseQuadMatrix(u_link,pid);
	
	if(N==1)
	{
		for(i=2;i<doneList->tree->size;i++)
		{
			freeALU(doneList->alu[i]);
		}
	}
	
	// ====================================================================

	time(&t3);
	printf("node %d final time: %g\n",N,difftime(t3,t2));

			
	set_1_sort_ParallelDoneList(doneList,ptr->todolist,N);
	decActiveThread();

	pthread_exit(0);
}


// ======================================================================================


static int isAllAncestorsDone(ParallelDoneList *doneList, const int currentFinishNode)
{
	 int ret = 1;
	 int current = currentFinishNode;
	 while(current>1)
	 {
	 	if(doneList->done[current]!=1)
		{
			ret = 0;
			break;
		}
	 	current = current/2;
	 }
	 return ret;
}



static int indexInOrderList(ParallelDoneList *doneList, const int key)
{
	int i;
	int index = -1;
	for(i=0;i<15;i++)
	{
		if(key == doneList->orderList[i])
		{
			index = i;
		}
	}
	return index;
}




static int isAllAncestorsPartialDone(ParallelDoneList *doneList,const int currentFinishNode)
{
	int ret = 1;
	int current = currentFinishNode;
	int baseIndex = indexInOrderList(doneList,currentFinishNode);
	while(current>0)
	{
		int targetIndex = indexInOrderList(doneList,doneList->eachNodeCurrentDone[current]);
//		printf("current:%d , currentFinishNode:%d\n",current,currentFinishNode);
		if(targetIndex < baseIndex)
		{
			ret = 0;
			break;
		}
		current = current/2;
	}
//	printf("ret=%d\n",ret);
	return ret;
}




static int freeALU_flag[16] = {0};

static void freeUnNecessaryALU(ParallelDoneList *doneList, const int currentFinishNode)
{
//	if(currentFinishNode > 7) return;

	int i = 0;
	int current = doneList->orderList[i];
	while(1)
	{
		int parent = current/2;
		if(freeALU_flag[current]!=1 && doneList->done[parent]==1 && isAllAncestorsPartialDone(doneList,parent)) // not free yet
		{
			freePidALU(doneList->alu[current],current);
			clearPidMempoolSet(current);
			freeALU_flag[current] = 1;
			printf("free alu current:%d\n",current);
		}

		if(current == currentFinishNode)
		{
			break;
		}
		else
		{
			i++;
			current = doneList->orderList[i];	
		}
	}
	
}




static void freeExcept1_2and3(ParallelDoneList *doneList)
{
	int i;
	for(i=4;i<doneList->tree->size;i++)
	{
		if(freeALU_flag[i]!=1)
		{
			freePidALU(doneList->alu[i],i);
			clearPidMempoolSet(i);
			freeALU_flag[i] = 1;
		}
	}
}





static void *setNodeDense(void *par)
{
	int i,j,k;
	time_t t1,t2,t3;
	time(&t1);
	ParallelLUQuadShareDataNew *ptr = (ParallelLUQuadShareDataNew *)par;
	ParallelDoneList *doneList = ptr->doneList;
	ToDoList *todolist = ptr->todolist;

	flagEachNode[ptr->N] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);
	fprintf(stderr,"setNode begin: %d\n",ptr->N);

	const int N = ptr->N;
	const int L = 2*N;
	const int R = L+1;
	const int rootCurrentBegin = ptr->rootCurrentBegin;
	const int currentEnd = ptr->currentEnd;;
	const int pid = ptr->pid;

	const int gvNum = doneList->gvNum;
	const int rowLink = doneList->tree->node[N]->rowEnd - doneList->tree->node[N]->rowBegin +1;
	const int colLink = doneList->alu[N]->u->totalCol;
	const int rowCrossBegin = doneList->tree->node[N]->rowBegin;
	const int rowCrossEnd = doneList->tree->node[N]->rowEnd;
	const int currentDoneRowBegin = doneList->tree->node[N]->doneRowBegin;

	SparseQuadMatrix *l_link = createPidSparseQuadMatrix(rowLink,doneList->alu[N]->l->totalCol,gvNum,pid);
	SparseQuadMatrix *u_link = createPidSparseQuadMatrix(rowLink,colLink,gvNum,pid);
	getSubPidSparseQuadMatrix(u_link,doneList->a,rowCrossBegin,currentDoneRowBegin,rowCrossEnd,doneList->a->totalCol-1,pid);

	QuadMatrix *l_link_dense = createPidQuadMatrix(rowLink,doneList->alu[N]->l->totalCol,gvNum,pid);
	QuadMatrix *u_link_dense = createPidQuadMatrix(rowLink,colLink,gvNum,pid);
	toDensePidSparseQuadMatrix(l_link_dense,l_link,pid);
	clearPidSparseQuadMatrix(l_link,pid);
	toDensePidSparseQuadMatrix(u_link_dense,u_link,pid);
	clearPidSparseQuadMatrix(u_link,pid);
	
	int rootCurrent = rootCurrentBegin;
	int current = -1;
	int base = 0;
	while(current != currentEnd)
	{
		if(checkNextParallelDoneList(doneList,rootCurrent) == 1)
		{
			rootCurrent++;
			current = doneList->orderList[rootCurrent];
			printf("setNode %d process node %d\n",N,current);
			const int left = current*2;
			const int right = left+1;
			int offset; // how many rows of the cross term or the lu term
			SparseQuadMatrix *partial_u;
			if(current > 7)
			{
				offset = doneList->alu[current]->u->totalRow;
				partial_u = doneList->alu[current]->u;
			}
			else
			{
				const int leftRow = doneList->alu[left]->u->totalRow;
				const int rightRow = doneList->alu[right]->u->totalRow;
				offset = doneList->alu[current]->u->totalRow - leftRow - rightRow;
				const int aTotal = doneList->a->totalCol;
				const SparseQuadMatrix *const uTemp = doneList->alu[current]->u;
				partial_u = createPidSparseQuadMatrix(uTemp->totalRow-leftRow-rightRow,doneList->a->totalCol-base,gvNum,pid);
				getSubPidSparseQuadMatrix(partial_u,uTemp,leftRow+rightRow, leftRow+rightRow,uTemp->totalRow-1,uTemp->totalCol-1,pid);
			}

			int i,j;
			QuadElement *scale = createPidQuadElement(gvNum,pid);
			QuadElement *newLinkU = createPidQuadElement(gvNum,pid);
			QuadElement *scaledPartialUData = createPidQuadElement(gvNum,pid);
			for(i=0;i<partial_u->totalRow;i++) // sweep partial_u row by row
			{
				for(j=0;j<u_link_dense->row;j++) // sweep u_link_dense row by row
				{
					const SparseQuadElement *partial_u_ptr = partial_u->rowIndex[i]->rowLink;
					if(partial_u_ptr == NULL)
					{
						printf(" ------- WTF!!!!!! -----------\n");
					}
					divPidQuadElement(scale,getPtrEntryQuadMatrix(u_link_dense,j,base+i),partial_u_ptr->data,pid);
					if(isEmptyQuadElement(scale))
					{
						delPidQuadMatrix(u_link_dense,j,base + i,pid);
						continue;
					}
					setPidQuadMatrix(l_link_dense,scale,j,base + partial_u_ptr->col,pid);
					while(partial_u_ptr!=NULL) // sweep each element in partial row i
					{
						const int row = j;
						const int col = base + partial_u_ptr->col;
						const QuadElement *const data = partial_u_ptr->data;
						
						QuadElement *linkU = NULL;
						linkU = getPtrEntryQuadMatrix(u_link_dense,row,col);
						mulPidQuadElement(scaledPartialUData,scale,data,pid);
						subQuadElement(newLinkU,linkU,scaledPartialUData);
						setPidQuadMatrix(u_link_dense,newLinkU,row,col,pid);
						partial_u_ptr = partial_u_ptr->rowLink;
					}
					delPidQuadMatrix(u_link_dense,j,base + i,pid);
				}
			}
			freePidQuadElement(scaledPartialUData,pid);
			freePidQuadElement(newLinkU,pid);
			freePidQuadElement(scale,pid);
			
			base = base + offset;
			if(current <=7)
			{
				freePidSparseQuadMatrix(partial_u,pid);
			}
			doneList->eachNodeCurrentDone[N] = current;
			if(N==1)
			{
				freeUnNecessaryALU(doneList,current);
			}
		}
		else
		{
			printf("node %d is going to sleep\n",N);
			decActiveThread();
			flagEachNode[N] = 1;
			pushBackToDoList(todolist,N);
//			pushByOrderToDoList(todolist,N);

			pthread_cond_wait(ptr->cond,ptr->mutex);
			printf("node %d is continue to work\n",N);
			time(&t1);
		}
	}

	if(N==1)
	{
		freeExcept1_2and3(doneList);
	}

	time(&t2);
	printf("set node %d skew time:%g\n",N,difftime(t2,t1));
	printf("node %d enter the final stage\n",N);

	QuadElement *scale = createPidQuadElement(gvNum,pid);
	QuadElement *newLinkU = createPidQuadElement(gvNum,pid);
	QuadElement *scaledPartialUData = createPidQuadElement(gvNum,pid);
	for(i=0;i<u_link_dense->row;i++) // for each row in link (to be pivot row)
	{
		for(j=i+1;j<u_link_dense->row;j++) // for each row under the pivot row
		{
			const QuadElement *const pivot = getPtrEntryQuadMatrix(u_link_dense,i,base+i);	
			const QuadElement *const eachRowHeadData = getPtrEntryQuadMatrix(u_link_dense,j,base+i);
			divPidQuadElement(scale,eachRowHeadData,pivot,pid);
			if(isEmptyQuadElement(scale))
			{
				delPidQuadMatrix(u_link_dense,j,base + i,pid);
				continue;
			}
			setPidQuadMatrix(l_link_dense,scale,j,base+i,pid);
			for(k=base;k<u_link_dense->col;k++)
			{
				const QuadElement * const data1 = getPtrEntryQuadMatrix(u_link_dense,i,k);
				const QuadElement * const data2 = getPtrEntryQuadMatrix(u_link_dense,j,k);
				mulPidQuadElement(scaledPartialUData,scale,data1,pid);
				subQuadElement(newLinkU,data2,scaledPartialUData);
				setPidQuadMatrix(u_link_dense,newLinkU,j,k,pid);
			}
			delPidQuadMatrix(u_link_dense,j,base+i,pid);
		}
	}
	freePidQuadElement(scaledPartialUData,pid);
	freePidQuadElement(newLinkU,pid);
	freePidQuadElement(scale,pid);

	delQuad2PidSparseQuadMatrix(l_link,l_link_dense,pid);
	delQuad2PidSparseQuadMatrix(u_link,u_link_dense,pid);
	freePidQuadMatrix(l_link_dense,pid);
	freePidQuadMatrix(u_link_dense,pid);

	// ====================================================================
	QuadElement *element = createPidQuadElement(gvNum,pid);
	setQuadElement(element,1.0,0,NULL,NULL);
	for(i=0;i<l_link->totalRow;i++)
	{
		const int row = i;
		const int col = doneList->alu[L]->l->totalRow + doneList->alu[R]->l->totalRow + i;
		setPidSparseQuadMatrix(l_link,element,row,col,pid);
	}
	freePidQuadElement(element,pid);

	time_t begin,end;

	time(&begin);
//	pthread_mutex_lock(&setnode_merge_mutex);
	mergePidSparseQuadMatrix(doneList->alu[N]->l,doneList->alu[L]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,0,0,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->l,doneList->alu[R]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow,doneList->alu[L]->l->totalRow,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->l,l_link,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow+doneList->alu[R]->l->totalRow,0,N);
//	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseQuadMatrix(l_link,pid);

//	pthread_mutex_lock(&setnode_merge_mutex);
	mergePidSparseQuadMatrix(doneList->alu[N]->u,doneList->alu[L]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,0,0,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->u,doneList->alu[R]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow,doneList->alu[L]->u->totalRow,N);
	mergePidSparseQuadMatrix(doneList->alu[N]->u,u_link,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow+doneList->alu[R]->u->totalRow,0,N);
//	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseQuadMatrix(u_link,pid);
	time(&end);
	printf("merge time:%g\n",difftime(end,begin));
	
	if(N==1)
	{
		freePidALU(doneList->alu[2],2);
		freePidALU(doneList->alu[3],3);
		clearPidMempoolSet(2);
		clearPidMempoolSet(3);
	}
	// ====================================================================

	time(&t3);
	printf("node %d final time: %g\n",N,difftime(t3,t2));

	set_1_sort_ParallelDoneList(doneList,ptr->todolist,N);
	decActiveThread();
	clearPidMempoolSet(pid);
	doneList->eachNodeCurrentDone[N] = N;
	pthread_exit(0);
}



// =================================================

// =================================================


static void threadHandlerNew(void *par)
{
	struct ThreadHandlerPar *thread_handle_ptr = par;

	ParallelLUQuadShareDataNew **list = thread_handle_ptr->list;
	ParallelDoneList *doneList = list[0]->doneList;
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
//				printf("wake up node %d\n",indexInToDoList);	
				incActiveThread();
//				dumpToDoList(stdout,todolist);
				pthread_cond_signal(list[indexInToDoList-1]->cond);
			}
		}
		usleep(500000);
//		usleep(100000);
	}
}







// =================================================


void parallelLUQuad(SparseQuadMatrix *l,SparseQuadMatrix *u, ParallelETree *tree, const SparseQuadMatrix *a, const int threadNum)
{
	const int treeInternalSize = (tree->size/2)-1;
	// pre-processing
	int i;
	int status;
	memset(flagEachNode,0,sizeof(int)*17);


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


	getDoneRowInfoNew(tree,postorder,postorder_inv,treeInternalSize);
	
	ALUQuad **aluList = createALUList(tree,a->totalRow,a->gvNum);
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i]!=NULL)
		{
			if(tree->node[i]->type == lu)
			{
				const int ltRowSrc = tree->node[i]->rowBegin;
				const int ltColSrc = ltRowSrc;
				const int rbRowSrc = tree->node[i]->rowEnd;
				const int rbColSrc = ltColSrc + aluList[i]->u->totalCol -1 ;
				getSubPidSparseQuadMatrix(aluList[i]->a,a,ltRowSrc,ltColSrc,rbRowSrc,rbColSrc,i);
			}
		}
	}
	ParallelDoneList *doneList = createDoneList(aluList,tree,a,a->gvNum);
	ToDoList *todolist = createToDoList(15,postorder);

	// set the par entries ... used for the parallel lu
	const int treeSize = 15;
	// mutex and cv list for node
	pthread_mutex_t mutex[treeSize];
	pthread_cond_t cond[treeSize];
	for(i=0;i<treeSize;i++)
	{
		pthread_mutex_init(&mutex[i],NULL);
		pthread_cond_init(&cond[i],NULL);
	}

	ParallelLUQuadShareDataNew **parList = getMempoolSet(sizeof(ParallelLUQuadShareDataNew *)*treeSize);
	if(threadNum == 1)
	{
		//used for cross nodes
		parList[0] = createParallelLUQuadShareDataNew(doneList,todolist,1,1,-1, 3,&mutex[0],&cond[0]);
		parList[1] = createParallelLUQuadShareDataNew(doneList,todolist,1,2,-1, 5,&mutex[1],&cond[1]);
		parList[2] = createParallelLUQuadShareDataNew(doneList,todolist,1,3, 6, 7,&mutex[2],&cond[2]);
		parList[3] = createParallelLUQuadShareDataNew(doneList,todolist,1,4,-1, 9,&mutex[3],&cond[3]);
		parList[4] = createParallelLUQuadShareDataNew(doneList,todolist,1,5, 2,11,&mutex[4],&cond[4]);
		parList[5] = createParallelLUQuadShareDataNew(doneList,todolist,1,6, 6,13,&mutex[5],&cond[5]);
		parList[6] = createParallelLUQuadShareDataNew(doneList,todolist,1,7, 9,15,&mutex[6],&cond[6]);
		// used for "LU node"
		for(i=7;i<15;i++) parList[i] = createParallelLUQuadShareDataNew(doneList,todolist,1,i+1,0,0,&mutex[i],&cond[i]);
	}
	else
	{
		//used for cross nodes
		parList[0] = createParallelLUQuadShareDataNew(doneList,todolist,1,1,-1, 3,&mutex[0],&cond[0]);
		parList[1] = createParallelLUQuadShareDataNew(doneList,todolist,2,2,-1, 5,&mutex[1],&cond[1]);
		parList[2] = createParallelLUQuadShareDataNew(doneList,todolist,3,3, 6, 7,&mutex[2],&cond[2]);
		parList[3] = createParallelLUQuadShareDataNew(doneList,todolist,4,4,-1, 9,&mutex[3],&cond[3]);
		parList[4] = createParallelLUQuadShareDataNew(doneList,todolist,5,5, 2,11,&mutex[4],&cond[4]);
		parList[5] = createParallelLUQuadShareDataNew(doneList,todolist,6,6, 6,13,&mutex[5],&cond[5]);
		parList[6] = createParallelLUQuadShareDataNew(doneList,todolist,7,7, 9,15,&mutex[6],&cond[6]);
		// used for "LU node"
		for(i=7;i<15;i++) parList[i] = createParallelLUQuadShareDataNew(doneList,todolist,i+1,i+1,0,0,&mutex[i],&cond[i]);
	}

	pthread_t pid[15];
	for(i=0;i<15;i++)
	{
		if(i==0) pthread_create(&pid[i],NULL,setNode,parList[i]);
		else if(i<7) pthread_create(&pid[i],NULL,setNode,parList[i]);
		else pthread_create(&pid[i],NULL,threadLU,parList[i]);
	}
	sleep(1);

	struct ThreadHandlerPar thread_handle;
	thread_handle.list = parList;
	thread_handle.threadNum = threadNum;
	threadHandlerNew(&thread_handle);

	// copy to the result
	copySparseQuadMatrix(l,parList[0]->doneList->alu[1]->l);
	copySparseQuadMatrix(u,parList[0]->doneList->alu[1]->u);

	retMempoolSet(postorder,sizeof(int)*(treeInternalSize));
	retMempoolSet(postorder_inv,sizeof(int)*(treeInternalSize+1));
//	for(i=0;i<tree->size;i++) freeALU(aluList[i]);
	for(i=0;i<2;i++) freePidALU(aluList[i],i);
	clearPidMempoolSet(1);
	for(i=0;i<treeSize;i++) freeParallelLUQuadShareDataNew(parList[i]);
	retMempoolSet(parList,treeSize*sizeof(ParallelLUQuadShareDataNew*));
	retMempoolSet(aluList,sizeof(ALUQuad *)*tree->size);
	freeParallelDoneList(doneList);
	freeToDoList(todolist);
}

