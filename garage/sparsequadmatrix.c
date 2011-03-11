#include "sparsequadmatrix.h"


// find out the pivot from currentRow and update the pivot "p" and the goal matrix "a"
static void updateSparseQuadMatrix(SparseQuadMatrix *a,SparseQuadMatrix *p, const int currentRow)
{
	SparseQuadElement *currentNode = a->colIndex[currentRow]->colLink;;
	// move until currentNode->row >= currentRow
	while(currentNode!=NULL && currentNode->row < currentRow) currentNode = currentNode->colLink;

	if(currentNode == NULL)
	{
		fprintf(stderr,"error: pivot in LU\n");
		exit(0);
	}
	else if(currentNode->row == currentRow)
	{
		// no need to pivot
		return;	
	}
	else
	{
		// update the permutation matrix
		swapRowSparseQuadMatrix(p,currentRow,currentNode->row);
		// swap currentRow and the pivot row
		swapRowSparseQuadMatrix(a,currentRow,currentNode->row);
	}
}


static int nnzInColForOpSparseQuadMatrix(const SparseQuadMatrix *a,const int col)
{
	int nnz = 0;

	SparseQuadElement *element = a->colIndex[col]->colLink;
	while(element!=NULL)
	{
		if(element->row <= col) element = element->colLink;
		else
		{
			nnz++;
			element = element->colLink;
		}
	}

	return nnz;
}


static int indexOfIthElementInColSparseQuadMatrix(const SparseQuadMatrix *a, const int col, const int i)
{
	int j;
	SparseQuadElement *element = a->colIndex[col]->colLink;
	for(j=0;j<i;j++)
	{
		if(element == NULL) return -1;
		else element = element->colLink;
	}
	return element->row;
}


// include Aii
static int nnzUntilAiiSparseQuadMatrix(const SparseQuadMatrix *a, const int row)
{
	int nnz = 0;

	SparseQuadElement *element = a->colIndex[row]->colLink;
	while(element!=NULL)
	{
		if(element->row > row) break;
		else
		{
			nnz++;
			element = element->colLink;
		}
	}
	return nnz;
}







//========================================================


SparseQuadElement * createSparseQuadElement(const int gvNum)
{
	SparseQuadElement *ptr = getMempoolSet(sizeof(SparseQuadElement));

	ptr->row = ptr->col = -1;
	ptr->rowLink = ptr->colLink = NULL;
	ptr->data = createQuadElement(gvNum);
	
	return ptr;
}

void freeSparseQuadElement(SparseQuadElement *ptr)
{
	freeQuadElement(ptr->data);
	retMempoolSet(ptr,sizeof(SparseQuadElement));
}




SparseQuadElement * createPidSparseQuadElement(const int gvNum,const int pid)
{
	SparseQuadElement *ptr = getPidMempoolSet(sizeof(SparseQuadElement),pid);

	ptr->row = ptr->col = -1;
	ptr->rowLink = ptr->colLink = NULL;
	ptr->data = createPidQuadElement(gvNum,pid);
	
	return ptr;
}

void freePidSparseQuadElement(SparseQuadElement *ptr,const int pid)
{
	freePidQuadElement(ptr->data,pid);
	retPidMempoolSet(ptr,sizeof(SparseQuadElement),pid);
}


//========================================================



SparseQuadMatrix *createSparseQuadMatrix(const int row,const int col,const int gvNum)
{
	SparseQuadMatrix *ptr = (SparseQuadMatrix *)getMempoolSet(sizeof(SparseQuadMatrix));
	ptr->totalRow = row;
	ptr->totalCol = col;
	ptr->gvNum = gvNum;
	ptr->nnz = 0;

	ptr->rowIndex = (SparseQuadElement **)getMempoolSet(row*sizeof(SparseQuadElement *));
	ptr->colIndex = (SparseQuadElement **)getMempoolSet(col*sizeof(SparseQuadElement *));
	
	int i;
	for(i=0;i<row;i++)	ptr->rowIndex[i] = createSparseQuadElement(gvNum);
	for(i=0;i<col;i++)	ptr->colIndex[i] = createSparseQuadElement(gvNum);
	return ptr;
}




SparseQuadMatrix *createPidSparseQuadMatrix(const int row,const int col,const int gvNum,const int pid)
{
	SparseQuadMatrix *ptr = (SparseQuadMatrix *)getPidMempoolSet(sizeof(SparseQuadMatrix),pid);
	ptr->totalRow = row;
	ptr->totalCol = col;
	ptr->gvNum = gvNum;
	ptr->nnz = 0;

	ptr->rowIndex = (SparseQuadElement **)getPidMempoolSet(row*sizeof(SparseQuadElement *),pid);
	ptr->colIndex = (SparseQuadElement **)getPidMempoolSet(col*sizeof(SparseQuadElement *),pid);
	
	int i;
	for(i=0;i<row;i++)	ptr->rowIndex[i] = createPidSparseQuadElement(gvNum,pid);
	for(i=0;i<col;i++)	ptr->colIndex[i] = createPidSparseQuadElement(gvNum,pid);
	return ptr;
}


void freeSparseQuadMatrix(SparseQuadMatrix *ptr)
{
	// free the memory row by row
	int i;
	SparseQuadElement *currentEntry;
	SparseQuadElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry != NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freeSparseQuadElement(removeEntry);
		}
	}

	for(i=0;i<ptr->totalRow;i++) freeSparseQuadElement(ptr->rowIndex[i]);
	for(i=0;i<ptr->totalCol;i++) freeSparseQuadElement(ptr->colIndex[i]);
	
	retMempoolSet(ptr->rowIndex,ptr->totalRow*sizeof(SparseQuadElement *));
	retMempoolSet(ptr->colIndex,ptr->totalCol*sizeof(SparseQuadElement *));
	retMempoolSet(ptr,sizeof(SparseQuadMatrix));
}





void freePidSparseQuadMatrix(SparseQuadMatrix *ptr,const int pid)
{
	// free the memory row by row
	int i;
	SparseQuadElement *currentEntry;
	SparseQuadElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry != NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freePidSparseQuadElement(removeEntry,pid);
		}
	}

	for(i=0;i<ptr->totalRow;i++) freePidSparseQuadElement(ptr->rowIndex[i],pid);
	for(i=0;i<ptr->totalCol;i++) freePidSparseQuadElement(ptr->colIndex[i],pid);
	
	retPidMempoolSet(ptr->rowIndex,ptr->totalRow*sizeof(SparseQuadElement *),pid);
	retPidMempoolSet(ptr->colIndex,ptr->totalCol*sizeof(SparseQuadElement *),pid);
	retPidMempoolSet(ptr,sizeof(SparseQuadMatrix),pid);
}





// important: element can NOT be null
// set ptr->data[rowIndex][colIndex] = element
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null && element!=null
void setSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex)
{
	SparseQuadElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseQuadElement *colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseQuadElement *insert = createSparseQuadElement(element->gvNum);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	
		copyQuadElement(insert->data,element);
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else copyQuadElement(rowTarget->rowLink->data,element);
}







// important: element can NOT be null
// set ptr->data[rowIndex][colIndex] = element
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null && element!=null
void setPidSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,const int pid)
{
	SparseQuadElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseQuadElement *colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseQuadElement *insert = createPidSparseQuadElement(element->gvNum,pid);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	
		copyQuadElement(insert->data,element);
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else copyQuadElement(rowTarget->rowLink->data,element);
}




void setFastSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,SparseQuadElement **baseRow,SparseQuadElement **baseCol)
{
	SparseQuadElement *rowTarget = NULL;
	SparseQuadElement *colTarget = NULL;
	
	// check *baseRow is available and legal
	if((*baseRow)!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col < colIndex)
		{
			rowTarget = *baseRow;
		}
	}
	// check *baseCol is available and legal
	if((*baseCol)!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row < rowIndex)
		{
			colTarget = *baseCol;
		}
	}
	// not set by rowBase 
	if(rowTarget == NULL)  rowTarget = ptr->rowIndex[rowIndex];
//	else printf("rowTarget is useful\n");

	// not set by colBase 
	if(colTarget == NULL)  colTarget = ptr->colIndex[colIndex];
//	else printf("colTarget is useful\n");

	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseQuadElement *insert = createSparseQuadElement(element->gvNum);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	
		copyQuadElement(insert->data,element);
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else copyQuadElement(rowTarget->rowLink->data,element);

	*baseRow = rowTarget;
	*baseCol = colTarget;
}






void setFastPidSparseQuadMatrix(SparseQuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,SparseQuadElement **baseRow,SparseQuadElement **baseCol,const int pid)
{
	SparseQuadElement *rowTarget = NULL;
	SparseQuadElement *colTarget = NULL;
	
	// check *baseRow is available and legal
	if((*baseRow)!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col < colIndex)
		{
			rowTarget = *baseRow;
		}
	}
	// check *baseCol is available and legal
	if((*baseCol)!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row < rowIndex)
		{
			colTarget = *baseCol;
		}
	}
	// not set by rowBase 
	if(rowTarget == NULL)  rowTarget = ptr->rowIndex[rowIndex];
	// not set by colBase 
	if(colTarget == NULL)  colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseQuadElement *insert = createPidSparseQuadElement(element->gvNum,pid);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	
		copyQuadElement(insert->data,element);
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else copyQuadElement(rowTarget->rowLink->data,element);

	*baseRow = rowTarget;
	*baseCol = colTarget;
}






static const SparseQuadElement *getNewSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex)
{
	SparseQuadElement *target = ptr->rowIndex[rowIndex]->rowLink;
	// set target
	while(target != NULL)
	{
		if(target->col >= colIndex ) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=colIndex) return NULL;
	else return target;
}





const QuadElement *getSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex)
{
	const SparseQuadElement *target = getNewSparseQuadMatrix(ptr,rowIndex,colIndex);
	if(target == NULL) return NULL;
	else return target->data;
}





const QuadElement *getFastRowSparseQuadMatrix(const SparseQuadMatrix *ptr, const int rowIndex, const int colIndex,SparseQuadElement **baseRow)
{
	SparseQuadElement *target = ptr->rowIndex[rowIndex]->rowLink;

	if(*baseRow!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col > target->col && (*baseRow)->col <= colIndex)
		{
			target = *baseRow;
		}
	}

	// set target
	while(target != NULL)
	{
		if(target->col >= colIndex ) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=colIndex)
	{
		return NULL;
	}
	else
	{
		*baseRow = target;
		return target->data;
	}
}








const QuadElement *getFastColSparseQuadMatrix(const SparseQuadMatrix *ptr,const int rowIndex, const int colIndex,SparseQuadElement **baseCol)
{
	SparseQuadElement *target = ptr->colIndex[colIndex]->colLink;

	if(baseCol!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row > target->row && (*baseCol)->col <= rowIndex)
		{
			target = *baseCol;
		}
	}

	// set target
	while(target != NULL)
	{
		if(target->row >= rowIndex) break;
		else target = target->colLink;
	}

	if(target == NULL || target->row != rowIndex)
	{
		return NULL;
	}
	else 
	{
		*baseCol = target;
		return target->data;
	}
}




void delSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex,const int colIndex)
{
	SparseQuadElement *rowPrev = ptr->rowIndex[rowIndex];
	SparseQuadElement *colPrev = ptr->colIndex[colIndex];
	SparseQuadElement *del;
	// set rowTarget
	while(rowPrev->rowLink != NULL)
	{
		if(rowPrev->rowLink->col >= colIndex ) break;
		else rowPrev = rowPrev->rowLink;
	}
	// set colTarget
	while(colPrev->colLink != NULL)
	{
		if(colPrev->colLink->row >= rowIndex ) break;
		else colPrev = colPrev->colLink;
	}
	// return directly if the entry is null
	if(rowPrev->rowLink ==NULL || rowPrev->rowLink->col!=colIndex) return;
	else 
	{
		del = rowPrev->rowLink;
		// update the information of prev
		rowPrev->rowLink = del->rowLink;
		colPrev->colLink = del->colLink;
		freeSparseQuadElement(del);
		ptr->nnz--;
	}
}




void delPidSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex,const int colIndex,const int pid)
{
	SparseQuadElement *rowPrev = ptr->rowIndex[rowIndex];
	SparseQuadElement *colPrev = ptr->colIndex[colIndex];
	SparseQuadElement *del;
	// set rowTarget
	while(rowPrev->rowLink != NULL)
	{
		if(rowPrev->rowLink->col >= colIndex ) break;
		else rowPrev = rowPrev->rowLink;
	}
	// set colTarget
	while(colPrev->colLink != NULL)
	{
		if(colPrev->colLink->row >= rowIndex ) break;
		else colPrev = colPrev->colLink;
	}
	// return directly if the entry is null
	if(rowPrev->rowLink ==NULL || rowPrev->rowLink->col!=colIndex) return;
	else 
	{
		del = rowPrev->rowLink;
		// update the information of prev
		rowPrev->rowLink = del->rowLink;
		colPrev->colLink = del->colLink;
		freePidSparseQuadElement(del,pid);
		ptr->nnz--;
	}
}




void delFastSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex, const int colIndex, SparseQuadElement **baseRow, SparseQuadElement **baseCol)
{
	SparseQuadElement *rowTarget = NULL;
	SparseQuadElement *colTarget = NULL;
	SparseQuadElement *del;
	
	// check *baseRow is available and legal
	if((*baseRow)!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col < colIndex)
		{
			rowTarget = *baseRow;
		}
	}
	// check *baseCol is available and legal
	if((*baseCol)!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row < rowIndex)
		{
			colTarget = *baseCol;
		}	
	}
	// not set by rowBase 
	if(rowTarget == NULL)  rowTarget = ptr->rowIndex[rowIndex];
	// not set by colBase 
	if(colTarget == NULL)  colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// return directly if the entry is null
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex) 
	{
		// ...	
	}
	else 
	{
		del = rowTarget->rowLink;
		// update the information of prev
		rowTarget->rowLink = del->rowLink;
		colTarget->colLink = del->colLink;
		freeSparseQuadElement(del);
		ptr->nnz--;
	}
	*baseRow = rowTarget;
	*baseCol = colTarget;
}






void delFastPidSparseQuadMatrix(SparseQuadMatrix *ptr, const int rowIndex, const int colIndex, SparseQuadElement **baseRow, SparseQuadElement **baseCol,const int pid)
{
	SparseQuadElement *rowTarget = NULL;
	SparseQuadElement *colTarget = NULL;
	SparseQuadElement *del;
	
	// check *baseRow is available and legal
	if((*baseRow)!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col < colIndex)
		{
			rowTarget = *baseRow;
		}
	}
	// check *baseCol is available and legal
	if((*baseCol)!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row < rowIndex)
		{
			colTarget = *baseCol;
		}
	}
	// not set by rowBase 
	if(rowTarget == NULL)  rowTarget = ptr->rowIndex[rowIndex];
	// not set by colBase 
	if(colTarget == NULL)  colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// return directly if the entry is null
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex) 
	{
		// ...	
	}
	else 
	{
		del = rowTarget->rowLink;
		// update the information of prev
		rowTarget->rowLink = del->rowLink;
		colTarget->colLink = del->colLink;
		freePidSparseQuadElement(del,pid);
		ptr->nnz--;
	}
	*baseRow = rowTarget;
	*baseCol = colTarget;
}




void dumpHeadSparseQuadMatrix(FILE *fp,const SparseQuadMatrix *ptr,const char *name)
{
	int i;
	fprintf(fp,"%s\n",name);
	fprintf(fp,"row:%d,col:%d\n",ptr->totalRow,ptr->totalCol);
	fprintf(fp,"nnz= %d\n",ptr->nnz);
}





void dumpSparseQuadMatrix(FILE *fp,const SparseQuadMatrix *ptr)
{
	int i;
	fprintf(fp,"nnz= %d\n",ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseQuadElement *target = ptr->rowIndex[i]->rowLink;
		while(target!=NULL)
		{
			fprintf(fp,"row: %d, col:%d\n",target->row,target->col);
			fprintf(fp,"node: %p, rowLink: %p, colLink: %p\n",target,target->rowLink,target->colLink);
			dumpQuadElement(target->data);
			target = target->rowLink;
		}
	}
}




void swapRowSparseQuadMatrix(SparseQuadMatrix *ptr, const int row1,const int row2)
{
	SparseQuadMatrix *tmpMatrix = createSparseQuadMatrix(ptr->totalRow,ptr->totalCol,ptr->gvNum);

	while(ptr->rowIndex[row1]->rowLink != NULL)
	{
		const int rowNew = row2;
		const int col = ptr->rowIndex[row1]->rowLink->col;
		const QuadElement *tmp = ptr->rowIndex[row1]->rowLink->data;
		setSparseQuadMatrix(tmpMatrix,tmp,rowNew,col);
		delSparseQuadMatrix(ptr,row1,col);
	}
	while(ptr->rowIndex[row2]->rowLink != NULL)
	{
		const int rowNew = row1;
		const int col = ptr->rowIndex[row2]->rowLink->col;
		const QuadElement *tmp = ptr->rowIndex[row2]->rowLink->data;
		setSparseQuadMatrix(tmpMatrix,tmp,rowNew,col);
		delSparseQuadMatrix(ptr,row2,col);
	}
	// dump tempMatrix->row1 to ptr->row1
	while(tmpMatrix->rowIndex[row1]->rowLink != NULL)
	{
		const int col = tmpMatrix->rowIndex[row1]->rowLink->col;
		const QuadElement *tmp = tmpMatrix->rowIndex[row1]->rowLink->data;
		setSparseQuadMatrix(ptr,tmp,row1,col);
		delSparseQuadMatrix(tmpMatrix,row1,col);
	}
	// dump tempMatrix->row2 to ptr->row2
	while(tmpMatrix->rowIndex[row2]->rowLink != NULL)
	{
		const int col = tmpMatrix->rowIndex[row2]->rowLink->col;
		const QuadElement *tmp = tmpMatrix->rowIndex[row2]->rowLink->data;
		setSparseQuadMatrix(ptr,tmp,row2,col);
		delSparseQuadMatrix(tmpMatrix,row2,col);
	}
	freeSparseQuadMatrix(tmpMatrix);
}


void clearSparseQuadMatrix(SparseQuadMatrix *ptr)
{
	clearPidSparseQuadMatrix(ptr,0);
}





void clearPidSparseQuadMatrix(SparseQuadMatrix *ptr,const int pid)
{
	// free the memory row by row
	int i;
	SparseQuadElement *currentEntry;
	SparseQuadElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			removeEntry = currentEntry;
			currentEntry = removeEntry->rowLink;
			freePidSparseQuadElement(removeEntry,pid);
		}
	}

	for(i=0;i<ptr->totalRow;i++) ptr->rowIndex[i]->rowLink = NULL;
	for(i=0;i<ptr->totalCol;i++) ptr->colIndex[i]->colLink = NULL;
	ptr->nnz = 0;
}




void copySparseQuadMatrix(SparseQuadMatrix *dest,const SparseQuadMatrix *src)
{
	copyPidSparseQuadMatrix(dest,src,0);
}





void copyPidSparseQuadMatrix(SparseQuadMatrix *dest,const SparseQuadMatrix *src,const int pid)
{
	if(dest == src) return;

	clearPidSparseQuadMatrix(dest,pid);
	int i;	
	SparseQuadElement **rowCache = getMempoolSet(sizeof(SparseQuadElement *)*src->totalRow);
	for(i=0;i<src->totalRow;i++) rowCache[i] = NULL;
	SparseQuadElement **colCache = getMempoolSet(sizeof(SparseQuadElement *)*src->totalCol);
	for(i=0;i<src->totalCol;i++) colCache[i] = NULL;
	SparseQuadElement *currentEntry;
	for(i=0;i<src->totalRow;i++)
	{
		currentEntry = src->rowIndex[i];
		while(currentEntry->rowLink!=NULL)
		{
			const int insCol = currentEntry->rowLink->col;
			setFastPidSparseQuadMatrix(dest,currentEntry->rowLink->data,i,insCol,&rowCache[i],&colCache[insCol],pid);
			currentEntry = currentEntry->rowLink;
		}
	}
	retMempoolSet(rowCache,sizeof(SparseQuadElement *)*src->totalRow);
	retMempoolSet(colCache,sizeof(SparseQuadElement *)*src->totalCol);
}






void addSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a,const SparseQuadMatrix *b)
{
	addPidSparseQuadMatrix(c,a,b,0);
}





void addPidSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a,const SparseQuadMatrix *b,const int pid)
{
	SparseQuadMatrix *cResult = createPidSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum,pid);
	QuadElement *entryC = createPidQuadElement(a->gvNum,pid);
	SparseQuadElement *entryA;
	SparseQuadElement *entryB;
	int i;
	for(i=0;i<a->totalRow;i++)
	{
		// add row by row ,
		// increase min(col1,col2) until row a or row b is null
		// then insert the remaining of that row
		entryA = a->rowIndex[i]->rowLink;
		entryB = b->rowIndex[i]->rowLink;
		while(entryA!=NULL && entryB!=NULL)
		{
			if(entryA->col < entryB->col)
			{
				const QuadElement *targetA = entryA->data;
				setPidSparseQuadMatrix(cResult,targetA,i,entryA->col,pid);
				entryA = entryA->rowLink;
			}
			else if(entryA->col > entryB->col)
			{
				const QuadElement *targetB = entryB->data;
				setPidSparseQuadMatrix(cResult,targetB,i,entryB->col,pid);
				entryB = entryB->rowLink;
			}
			else
			{
				const QuadElement *targetA = entryA->data;
				const QuadElement *targetB = entryB->data;
				addQuadElement(entryC,targetA,targetB);
				setPidSparseQuadMatrix(cResult,entryC,i,entryA->col,pid);
				entryA = entryA->rowLink;
				entryB = entryB->rowLink;
			}
		}
		// insert the remaining of a
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			setPidSparseQuadMatrix(cResult,targetA,i,entryA->col,pid);
			entryA = entryA->rowLink;
		}
		// insert the remaining of b
		while(entryB!=NULL)
		{
			const QuadElement *targetB = entryB->data;
			setPidSparseQuadMatrix(cResult,targetB,i,entryB->col,pid);
			entryB = entryB->rowLink;
		}
	}
	copyPidSparseQuadMatrix(c,cResult,pid);
	
	freePidSparseQuadMatrix(cResult,pid);
	freePidQuadElement(entryC,pid);
}




// c = a - b
void subSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a,const SparseQuadMatrix *b)
{
	subPidSparseQuadMatrix(c,a,b,0);
}





void subPidSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a,const SparseQuadMatrix *b,const int pid)
{
	SparseQuadMatrix *cResult = createPidSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum,pid);
	QuadElement *entryC = createPidQuadElement(a->gvNum,pid);
	SparseQuadElement *entryA;
	SparseQuadElement *entryB;
	int i;
	for(i=0;i<a->totalRow;i++)
	{
		// sub row by row ,
		// increase min(col1,col2) until row a or row b is null
		// then insert the remaining of that row
		entryA = a->rowIndex[i]->rowLink;
		entryB = b->rowIndex[i]->rowLink;
		while(entryA!=NULL && entryB!=NULL)
		{
			if(entryA->col < entryB->col)
			{
				const QuadElement *targetA = entryA->data;
				setPidSparseQuadMatrix(cResult,targetA,i,entryA->col,pid);
				entryA = entryA->rowLink;
			}
			else if(entryA->col > entryB->col)
			{
				const QuadElement *targetA = NULL;
				const QuadElement *targetB = entryB->data;
				subQuadElement(entryC,targetA,targetB);
				setPidSparseQuadMatrix(cResult,entryC,i,entryB->col,pid);
				entryB = entryB->rowLink;
			}
			else
			{
				const QuadElement *targetA = entryA->data;
				const QuadElement *targetB = entryB->data;
				subQuadElement(entryC,targetA,targetB);
				setPidSparseQuadMatrix(cResult,entryC,i,entryA->col,pid);
				entryA = entryA->rowLink;
				entryB = entryB->rowLink;
			}
		}
		// insert the remaining of a
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			setPidSparseQuadMatrix(cResult,targetA,i,entryA->col,pid);
			entryA = entryA->rowLink;
		}
		// insert the remaining of b
		while(entryB!=NULL)
		{
			const QuadElement *targetA = NULL;
			const QuadElement *targetB = entryB->data;
			subQuadElement(entryC,targetA,targetB);
			setPidSparseQuadMatrix(cResult,entryC,i,entryB->col,pid);
			entryB = entryB->rowLink;
		}
	}
	copyPidSparseQuadMatrix(c,cResult,pid);
	
	freePidSparseQuadMatrix(cResult,pid);
	freePidQuadElement(entryC,pid);
}




// c = a * b
void mulSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a, const SparseQuadMatrix *b)
{
	SparseQuadMatrix *cResult = createSparseQuadMatrix(c->totalRow,c->totalCol,c->gvNum);
	QuadElement *entryC = createQuadElement(a->gvNum);
	SparseQuadElement *entryA;
	SparseQuadElement *entryB;
	int i,j;
	for(i=0;i<a->totalRow;i++)
	{
		for(j=0;j<b->totalCol;j++)
		{
			int flag = 0;
			QuadElement *tempSum = createQuadElement(a->gvNum);
			entryA = a->rowIndex[i]->rowLink;
			entryB = b->colIndex[j]->colLink;
			// calculate the product of row a(i,:) and col b(:,j) 
			// and then set to c(i,j)
			while(entryA!=NULL && entryB!=NULL)
			{
				if(entryA->col < entryB->row)
				{
					entryA = entryA->rowLink;
				}
				else if(entryA->col > entryB->row)
				{
					entryB = entryB->colLink;
				}
				else
				{
					const QuadElement *targetA = entryA->data;
					const QuadElement *targetB = entryB->data;
					mulQuadElement(entryC,targetA,targetB);
					addQuadElement(tempSum,tempSum,entryC);
					entryA = entryA->rowLink;
					entryB = entryB->colLink;
					flag = 1;
				}
			}
			if(flag==1) setSparseQuadMatrix(cResult,tempSum,i,j);
			freeQuadElement(tempSum);
		}
	}
	copySparseQuadMatrix(c,cResult);
	
	freeSparseQuadMatrix(cResult);
	freeQuadElement(entryC);
}


// c = k * a , k is a scale
void scaleSparseQuadMatrix(SparseQuadMatrix *c,const double k,const SparseQuadMatrix *a)
{
	int i;
	SparseQuadElement *entryA;
	QuadElement *entryC = createQuadElement(a->gvNum);
	SparseQuadMatrix *cResult = createSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			scaleQuadElement(entryC,k,targetA);
			setSparseQuadMatrix(cResult,entryC,i,entryA->col);
			entryA = entryA->rowLink;
		}
	}
	copySparseQuadMatrix(c,cResult);
	freeQuadElement(entryC);
	freeSparseQuadMatrix(cResult);
}




void scalePidSparseQuadMatrix(SparseQuadMatrix *c,const double k,const SparseQuadMatrix *a,const int pid)
{
	int i;
	SparseQuadElement *entryA;
	QuadElement *entryC = createPidQuadElement(a->gvNum,pid);
	SparseQuadMatrix *cResult = createPidSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum,pid);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			scaleQuadElement(entryC,k,targetA);
			setPidSparseQuadMatrix(cResult,entryC,i,entryA->col,pid);
			entryA = entryA->rowLink;
		}
	}
	copyPidSparseQuadMatrix(c,cResult,pid);
	freePidQuadElement(entryC,pid);
	freePidSparseQuadMatrix(cResult,pid);
}




// c = trans(A)
void transSparseQuadMatrix(SparseQuadMatrix *c, const SparseQuadMatrix *a)
{
	int i;
	SparseQuadElement *entryA;
	QuadElement *entryC = createQuadElement(a->gvNum);
	SparseQuadMatrix *cResult = createSparseQuadMatrix(a->totalRow,a->totalCol,a->gvNum);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			setSparseQuadMatrix(cResult,targetA,entryA->col,i);
			entryA = entryA->rowLink;
		}
	}
	copySparseQuadMatrix(c,cResult);
	freeQuadElement(entryC);
	freeSparseQuadMatrix(cResult);
}



// c = a * b, (sparse multiply dense)
// a is a sparse matrix
// b is a n*1 QuadMatrix
// c is a n*1 QuadMatrix
void mulVecSparseQuadMatrix(QuadMatrix *c,const SparseQuadMatrix *a, const QuadMatrix *b)
{
	int i,j;
	SparseQuadElement *entryA;
	QuadElement *entryC = createQuadElement(a->gvNum);
	QuadMatrix *cResult = createQuadMatrix(a->totalRow,1,a->gvNum);
	QuadElement *sum = createQuadElement(a->gvNum);

	setZeroQuadElement(sum,a->gvNum);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		resetQuadElement(sum);
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			const QuadElement *targetB = getPtrEntryQuadMatrix(b,entryA->col,0);
			mulQuadElement(entryC,targetA,targetB);
			addQuadElement(sum,sum,entryC);
			entryA = entryA->rowLink;
		}
		setQuadMatrix(cResult,sum,i,0);
	}
	copyQuadMatrix(c,cResult);
	
	freeQuadElement(entryC);
	freeQuadMatrix(cResult);
	freeQuadElement(sum);
}





void mulVecPidSparseQuadMatrix(QuadMatrix *c,const SparseQuadMatrix *a, const QuadMatrix *b,const int pid)
{
	int i,j;
	SparseQuadElement *entryA;
	QuadElement *entryC = createPidQuadElement(a->gvNum,pid);
	QuadMatrix *cResult = createPidQuadMatrix(a->totalRow,1,a->gvNum,pid);
	QuadElement *sum = createPidQuadElement(a->gvNum,pid);

	setZeroPidQuadElement(sum,a->gvNum,pid);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		resetQuadElement(sum);
		while(entryA!=NULL)
		{
			const QuadElement *targetA = entryA->data;
			const QuadElement *targetB = getPtrEntryQuadMatrix(b,entryA->col,0);
			mulPidQuadElement(entryC,targetA,targetB,pid);
			addQuadElement(sum,sum,entryC);
			entryA = entryA->rowLink;
		}
		setPidQuadMatrix(cResult,sum,i,0,pid);
	}
	copyPidQuadMatrix(c,cResult,pid);
	
	freePidQuadElement(entryC,pid);
	freePidQuadMatrix(cResult,pid);
	freePidQuadElement(sum,pid);
}



void identitySparseQuadMatrix(SparseQuadMatrix *a)
{
	identityPidSparseQuadMatrix(a,0);
}



// set a as identity matrix
void identityPidSparseQuadMatrix(SparseQuadMatrix *a,const int pid)
{
	clearPidSparseQuadMatrix(a,pid);
	QuadElement *element = createPidQuadElement(a->gvNum,pid);
	element->m = 1.0;
	int i;
	for(i=0;i<a->totalRow;i++) setPidSparseQuadMatrix(a,element,i,i,pid);

	freePidQuadElement(element,pid);
}



struct SparseLUParParallel
{
	int currentRow;
	// for(i = rowLow ; i<rowHigh+1 ; i++)
	int rowLow;
	int rowHigh; // op will include this row
	SparseQuadMatrix *aCopy;
	SparseQuadMatrix *partialL;
};


static void *eliminateRows(void *ptr)
{
	struct SparseLUParParallel *par = (struct SparseLUParParallel *) ptr;
	const SparseQuadElement *aii = par->aCopy->rowIndex[par->currentRow]->rowLink;
	QuadElement *scale = createQuadElement(par->aCopy->gvNum);
	QuadElement *element = createQuadElement(par->aCopy->gvNum);
	
	SparseQuadElement *eachRow = aii->colLink;
	while(eachRow != NULL) // update E
	{
		if(eachRow->row < par->rowLow )
		{
			eachRow = eachRow->colLink;
			continue;
		}
		else if( eachRow->row > par->rowHigh) break;
		else
		{
			divQuadElement(scale,eachRow->data,aii->data);
			
//			fprintf(stderr,"(%d,%d)\n",eachRow->row,par->currentRow);
			setSparseQuadMatrix(par->partialL,scale,eachRow->row,par->currentRow);
	
			SparseQuadElement *inRow = aii->rowLink;
			while(inRow != NULL)
			{
				mulQuadElement(element,scale,getSparseQuadMatrix(par->aCopy,par->currentRow,inRow->col));
				subQuadElement(element,getSparseQuadMatrix(par->aCopy,eachRow->row,inRow->col),element);
				setSparseQuadMatrix(par->aCopy,element,eachRow->row,inRow->col);
				inRow = inRow->rowLink;
			}
			delSparseQuadMatrix(par->aCopy,eachRow->row,par->currentRow);
			eachRow = eachRow->colLink;
		}
	}

	freeQuadElement(scale);
	freeQuadElement(element);

	fprintf(stderr,"end of eliminate\n");

	pthread_exit(NULL);
}





// a(i,j) = a(i,j) + element
// will allocate the memory if necessary
void incSparseQuadMatrix(SparseQuadMatrix *a,const QuadElement *element,const int row, const int col)
{
	const QuadElement *ptr = getSparseQuadMatrix(a,row,col);
	if(ptr == NULL)
	{
		setSparseQuadMatrix(a,element,row,col);
	}
	else
	{
		QuadElement *result = createQuadElement(a->gvNum);
		addQuadElement(result,ptr,element);
		setSparseQuadMatrix(a,result,row,col);
		freeQuadElement(result);
	}
}



// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
void decSparseQuadMatrix(SparseQuadMatrix *a,const QuadElement *element,const int row, const int col)
{
	const QuadElement *ptr = getSparseQuadMatrix(a,row,col);
	QuadElement *result = createQuadElement(a->gvNum);
	if(ptr == NULL)
	{
		scaleQuadElement(result,-1,element);	
	}
	else
	{
		subQuadElement(result,ptr,element);
	}
	setSparseQuadMatrix(a,result,row,col);
	freeQuadElement(result);
}



// dest and src should be allocated before calling this function
void quad2SparseQuadMatrix(SparseQuadMatrix *dest,const QuadMatrix *src)
{
	int i,j;
	clearSparseQuadMatrix(dest);
	for(i=0;i<src->row;i++)
	{
		for(j=0;j<src->col;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(src,i,j);
			if( isEmptyQuadElement(element) ) continue;
			else setSparseQuadMatrix(dest,element,i,j);
		}
	}
}




// dest and src should be allocated before calling this function
void quad2PidSparseQuadMatrix(SparseQuadMatrix *dest,const QuadMatrix *src,const int pid)
{
	int i,j;
	clearPidSparseQuadMatrix(dest,pid);
	for(i=0;i<src->row;i++)
	{
		for(j=0;j<src->col;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(src,i,j);
			if( isEmptyQuadElement(element) ) continue;
			else setPidSparseQuadMatrix(dest,element,i,j,pid);
		}
	}
}




void delQuad2PidSparseQuadMatrix(SparseQuadMatrix *dest, QuadMatrix *src,const int pid)
{
	int i,j;
	clearPidSparseQuadMatrix(dest,pid);
	for(i=0;i<src->row;i++)
	{
		for(j=0;j<src->col;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(src,i,j);
			if( isEmptyQuadElement(element) )
			{
				continue;
			}
			else
			{
				setPidSparseQuadMatrix(dest,element,i,j,pid);
				delPidQuadMatrix(src,i,j,pid);
			}
		}
	}
}




// dest and src should be allocated before calling this function
void toDenseSparseQuadMatrix(QuadMatrix *dest, const SparseQuadMatrix *src)
{
	int i,j;
	resetQuadMatrix(dest);
	for(i=0;i<src->totalRow;i++)
	{
		SparseQuadElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = ptr->row;
			const int col = ptr->col;
			const QuadElement *data = ptr->data;
			setQuadMatrix(dest,data,row,col);
			ptr = ptr->rowLink;
		}
	}
}




// dest and src should be allocated before calling this function
void toDensePidSparseQuadMatrix(QuadMatrix *dest, const SparseQuadMatrix *src,const int pid)
{
	int i,j;
	resetQuadMatrix(dest);
	for(i=0;i<src->totalRow;i++)
	{
		SparseQuadElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = ptr->row;
			const int col = ptr->col;
			const QuadElement *data = ptr->data;
			setPidQuadMatrix(dest,data,row,col,pid);
			ptr = ptr->rowLink;
		}
	}
}






void setDenseCol2SparseQuadMatrix(SparseQuadMatrix *dest, const QuadMatrix *src,const int col)
{
	int i,j;
	for(i=0;i<src->row;i++)
	{
		const QuadElement *element = getPtrEntryQuadMatrix(src,i,0);
		if( isEmptyQuadElement(element) ) delSparseQuadMatrix(dest,i,col);
		else setSparseQuadMatrix(dest,element,i,col);
	}
}



void setDenseColQuick2SparseQuadMatrix(SparseQuadMatrix *dest, const QuadMatrix *src,const int col)
{
	setDenseColQuick2SparseQuadMatrixPid(dest,src,col,0);
}





void setDenseColQuick2SparseQuadMatrixPid(SparseQuadMatrix *dest, const QuadMatrix *src,const int col,const int pid)
{
	int i,j;
	for(i=0;i<src->row;i++)
	{
		const QuadElement *element = getPtrEntryQuadMatrix(src,i,0);
		if( isEmptyQuadElement(element) ) continue;
		else setPidSparseQuadMatrix(dest,element,i,col,pid);
	}
}





// the boundaries row1,col1,row2,col2 are included
void clearBlockSparseQuadMatrix(SparseQuadMatrix *dest, const int row1, const int col1, const int row2, const int col2)
{
	int i;
	// del the old data in the block
	for(i=row1;i<=row2;i++)
	{
		SparseQuadElement *eachRow = dest->rowIndex[i]->rowLink;
		while(eachRow != NULL)
		{
			if(eachRow->col < col1) eachRow = eachRow->rowLink;
			else if(eachRow->col > col2) break;
			else
			{
				const int delRow = eachRow->row;
				const int delCol = eachRow->col;
				eachRow = eachRow->rowLink;
				delSparseQuadMatrix(dest,delRow,delCol);
			}
		}
	}
}




// the boundaries row1,col1,row2,col2 are included
void appendBlockSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int row1, const int col1, const int row2, const int col2)
{
	if(dest == src) return;

	int i;


	// insert new data to the block
	for(i=row1;i<=row2;i++)
	{
		SparseQuadElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow != NULL)
		{
			if(eachRow->col < col1) eachRow = eachRow->rowLink;
			else if(eachRow->col > col2) break;
			else
			{
				const QuadElement *data = eachRow->data;
				const int insRow = eachRow->row;
				const int insCol = eachRow->col;
				setSparseQuadMatrix(dest,data,insRow,insCol);
				eachRow = eachRow->rowLink;
			}
		}
	}


/*
	for(i=row1;i<=row2;i++)
	{
		for(j=col1;j<=col2;j++)
		{
			const QuadElement *data = getSparseQuadMatrix(src,i,j);
			if(data!=NULL)	setSparseQuadMatrix(dest,data,i,j);
			else
			{
//				fprintf(stderr,"start to del (%d,%d)\n",i,j);
				delSparseQuadMatrix(dest,i,j);
//				fprintf(stderr,"del complete\n");
			}
		}
	}
*/
}



// get the m part of sparseQuadMatrix
void mSparseQuadMatrix(SparseDoubleMatrix *dest, const SparseQuadMatrix *src)
{
	clearSparseDoubleMatrix(dest);
	int i;
	for(i=0;i<src->totalRow;i++)
	{
		const SparseQuadElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = i;
			const int col = ptr->col;
			const double val = ptr->data->m;
			setSparseDoubleMatrix(dest,val,row,col);
			ptr = ptr->rowLink;
		}
	}
}


// set the m part of sparseQuadMatrix
void setMSparseQuadMatrix(SparseQuadMatrix *dest, const SparseDoubleMatrix *src)
{
	clearSparseQuadMatrix(dest);
	int i;
	for(i=0;i<src->totalRow;i++)
	{
		const SparseDoubleElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = i;
			const int col = ptr->col;
			const double valDouble = ptr->data;
			QuadElement *val = createQuadElement(dest->gvNum);
			setQuadElement(val,valDouble,0.0,NULL,NULL);
			setSparseQuadMatrix(dest,val,row,col);
			ptr = ptr->rowLink;
			freeQuadElement(val);
		}
	}

}




// permutate 
void permutateSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *pRow, const SparseQuadMatrix *pCol, const SparseQuadMatrix *src)
{
	SparseQuadMatrix *destTemp = createSparseQuadMatrix(dest->totalRow,dest->totalCol,dest->gvNum);
	SparseQuadMatrix *destTemp2 = createSparseQuadMatrix(dest->totalRow,dest->totalCol,dest->gvNum);
	int i;

	// permutate row
	for(i=0;i<pRow->totalRow;i++)
	{
		const int colForRowI = pRow->rowIndex[i]->rowLink->col;
		const SparseQuadElement *ptr = src->rowIndex[colForRowI]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
			setSparseQuadMatrix(destTemp,ptr->data,i,col);			
			ptr = ptr->rowLink;
		}
	}


	int *iperm = getMempoolSet(sizeof(int)*pCol->totalRow);
	for(i=0;i<pCol->totalRow;i++)
	{
		const int col = pCol->rowIndex[i]->rowLink->col;
		iperm[i] = col;
	}

	for(i=0;i<destTemp->totalRow;i++)
	{
		SparseQuadElement *ptr = destTemp->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = i;
			const int col = ptr->col;
			setSparseQuadMatrix(destTemp2,ptr->data,row,iperm[col]);
			ptr = ptr->rowLink;
		}
	}
	retMempoolSet(iperm,sizeof(int)*pCol->totalRow);

	copySparseQuadMatrix(dest,destTemp2);
	freeSparseQuadMatrix(destTemp);
	freeSparseQuadMatrix(destTemp2);
}




void getSubPidSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc,const int pid)
{
	int i;

	clearPidSparseQuadMatrix(dest,pid);

	const int baseRow = ltRowSrc;
	const int baseCol = ltColSrc;

	for(i=ltRowSrc;i<=rbRowSrc;i++)
	{
		const SparseQuadElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			if(eachRow->col < ltColSrc) eachRow = eachRow->rowLink;
			else if(eachRow->col > rbColSrc) break;
			else
			{
				const QuadElement *data = eachRow->data;
				const int insRow = eachRow->row - baseRow;
				const int insCol = eachRow->col - baseCol;
				setPidSparseQuadMatrix(dest,data,insRow,insCol,pid);
				eachRow = eachRow->rowLink;
			}
		}
	}
}




void getSubSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc)
{
	getSubPidSparseQuadMatrix(dest,src,ltRowSrc,ltColSrc,rbRowSrc,rbColSrc,0);
}




void mergeSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest)
{
	mergePidSparseQuadMatrix(dest,src,destRow,destCol,ltRowDest,ltColDest,0);
}





void mergePidSparseQuadMatrix(SparseQuadMatrix *dest, const SparseQuadMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest,const int pid)
{
	int i;
	
	const int baseRow = ltRowDest;
	const int baseCol = ltColDest;

	SparseQuadElement **rowCache = getPidMempoolSet(sizeof(SparseQuadElement *)*dest->totalRow,pid);
	for(i=0;i<dest->totalRow;i++) rowCache[i] = NULL;
	SparseQuadElement **colCache = getPidMempoolSet(sizeof(SparseQuadElement *)*dest->totalCol,pid);
	for(i=0;i<dest->totalCol;i++) colCache[i] = NULL;

	for(i=0;i<src->totalRow;i++)
	{
		const SparseQuadElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			const QuadElement* data = eachRow->data;	
			const int insRow = baseRow + eachRow->row;
			const int insCol = baseCol + eachRow->col;
			setFastPidSparseQuadMatrix(dest,data,insRow,insCol,&rowCache[insRow],&colCache[insCol],pid);
			eachRow = eachRow->rowLink;
		}
	}

	retPidMempoolSet(rowCache,sizeof(SparseQuadElement *)*dest->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseQuadElement *)*dest->totalCol,pid);
}



