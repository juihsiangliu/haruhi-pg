#include "sparsedoublematrix.h"




SparseDoubleElement * createSparseDoubleElement(const double data)
{
	return createPidSparseDoubleElement(data,0);
}


SparseDoubleElement * createPidSparseDoubleElement(const double data,const int pid)
{
	SparseDoubleElement *ptr = getPidMempoolSet(sizeof(SparseDoubleElement),pid);

	ptr->row = ptr->col = -1;
	ptr->rowLink = ptr->colLink = NULL;
	ptr->data = data;
	
	return ptr;
}




void freeSparseDoubleElement(SparseDoubleElement *ptr)
{
	freePidSparseDoubleElement(ptr,0);
}




void freePidSparseDoubleElement(SparseDoubleElement *ptr,const int pid)
{
	retPidMempoolSet(ptr,sizeof(SparseDoubleElement),pid);
}


//========================================================

SparseDoubleMatrix *createSparseDoubleMatrix(const int row,const int col)
{
	return createPidSparseDoubleMatrix(row,col,0);
}




SparseDoubleMatrix *createPidSparseDoubleMatrix(const int row,const int col,const int pid)
{
	int i;
	SparseDoubleMatrix *ptr = (SparseDoubleMatrix *)getPidMempoolSet(sizeof(SparseDoubleMatrix),pid);
	ptr->totalRow = row;
	ptr->totalCol = col;
	ptr->nnz = 0;

	ptr->rowIndex = (SparseDoubleElement **)getPidMempoolSet(row*sizeof(SparseDoubleElement *),pid);
	ptr->colIndex = (SparseDoubleElement **)getPidMempoolSet(col*sizeof(SparseDoubleElement *),pid);
	for(i=0;i<row;i++)	ptr->rowIndex[i] = createPidSparseDoubleElement(0,pid);
	for(i=0;i<col;i++)	ptr->colIndex[i] = createPidSparseDoubleElement(0,pid);
	
	ptr->rowCache = (SparseDoubleElement **)getPidMempoolSet(row*sizeof(SparseDoubleElement *),pid);
	ptr->colCache = (SparseDoubleElement **)getPidMempoolSet(col*sizeof(SparseDoubleElement *),pid);
	for(i=0;i<row;i++)	ptr->rowCache[i] = NULL;
	for(i=0;i<col;i++)	ptr->colCache[i] = NULL;

	return ptr;
}





void freeSparseDoubleMatrix(SparseDoubleMatrix *ptr)
{
	freePidSparseDoubleMatrix(ptr,0);
}




void freePidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const int pid)
{
	// free the memory row by row
	int i;
	SparseDoubleElement *currentEntry;
	SparseDoubleElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry != NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freePidSparseDoubleElement(removeEntry,pid);
		}
	}

	for(i=0;i<ptr->totalRow;i++) freePidSparseDoubleElement(ptr->rowIndex[i],pid);
	for(i=0;i<ptr->totalCol;i++) freePidSparseDoubleElement(ptr->colIndex[i],pid);
	retPidMempoolSet(ptr->rowIndex,ptr->totalRow*sizeof(SparseDoubleElement *),pid);
	retPidMempoolSet(ptr->colIndex,ptr->totalCol*sizeof(SparseDoubleElement *),pid);

	retPidMempoolSet(ptr->rowCache,ptr->totalRow*sizeof(SparseDoubleElement *),pid);
	retPidMempoolSet(ptr->colCache,ptr->totalCol*sizeof(SparseDoubleElement *),pid);
	retPidMempoolSet(ptr,sizeof(SparseDoubleMatrix),pid);
}



// ===========================================================================

/*
static void setPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,const int pid)
{
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];
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
		SparseDoubleElement *insert = createPidSparseDoubleElement(element,pid);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data = element;
}
*/



static void setFastPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element ,const int rowIndex, const int colIndex,SparseDoubleElement **baseRow,SparseDoubleElement **baseCol,const int pid)
{
//	setPidSparseDoubleMatrix(ptr,element,rowIndex,colIndex,pid);
//	return;

	SparseDoubleElement *rowTarget = NULL;
	SparseDoubleElement *colTarget = NULL;
	

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
		SparseDoubleElement *insert = createPidSparseDoubleElement(element,pid);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	
		
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data = element;

	*baseRow = rowTarget;
	*baseCol = colTarget;
}



static void setFastSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element ,const int rowIndex, const int colIndex,SparseDoubleElement **baseRow,SparseDoubleElement **baseCol)
{
	setFastPidSparseDoubleMatrix(ptr,element,rowIndex,colIndex,baseRow,baseCol,0);
}



// important: element can NOT be null
// set ptr->data[rowIndex][colIndex] = element
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null && element!=null
void setSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex)
{
	setPidSparseDoubleMatrix(ptr,element,rowIndex,colIndex,0);
}




void setPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,const int pid)
{
	setFastPidSparseDoubleMatrix(ptr,element,rowIndex,colIndex,&ptr->rowCache[rowIndex],&ptr->colCache[colIndex],pid);
}



// ===========================================================================

/*
static double getSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex)
{
	SparseDoubleElement *target = ptr->rowIndex[rowIndex]->rowLink;
	// set target
	while(target != NULL)
	{
		if(target->col >= colIndex ) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=colIndex) return 0;
	else return target->data;
}
*/




static double getFastRowSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex, SparseDoubleElement **baseRow)
{
	SparseDoubleElement *target = NULL;
	if(*baseRow!=NULL)
	{
		if( (*baseRow)->row == rowIndex && (*baseRow)->col < colIndex )
		{
			target = *baseRow;
		}
	}

	if(target==NULL) target = ptr->rowIndex[rowIndex];

	// set target
	while(target->rowLink != NULL)
	{
		if(target->rowLink->col >= colIndex ) break;
		else target = target->rowLink;
	}

	if(target->rowLink == NULL || target->rowLink->col!=colIndex)
	{
		return 0;
	}
	else
	{
		*baseRow = target;
		return target->rowLink->data;
	}
}






static double getFastColSparseDoubleMatrix(const SparseDoubleMatrix *ptr,const int rowIndex, const int colIndex,SparseDoubleElement **baseCol)
{
//	return getSparseDoubleMatrix(ptr,rowIndex,colIndex);

	SparseDoubleElement *target = NULL;
	if(*baseCol!=NULL)
	{
		if((*baseCol)->col == colIndex  && (*baseCol)->row < rowIndex)
		{
			target = *baseCol;
		}
	}
	if(target == NULL) target = ptr->colIndex[colIndex];

	// set target
	while(target->colLink != NULL)
	{
		if(target->colLink->row >= rowIndex) break;
		else target = target->colLink;
	}

	if(target->colLink == NULL || target->colLink->row != rowIndex)
	{
		return 0;
	}
	else
	{
		*baseCol = target;
		return target->colLink->data;
	}
}





double getSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex,const char *type)
{
	if(strcmp(type,"col") == 0)
	{
		return getFastColSparseDoubleMatrix(ptr,rowIndex,colIndex,&ptr->colCache[colIndex]);
	}
	else
	{
		return getFastRowSparseDoubleMatrix(ptr,rowIndex,colIndex,&ptr->rowCache[rowIndex]);
	}
}




//==============================================================

/*
static void delPidSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex,const int pid)
{
	SparseDoubleElement *rowPrev = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colPrev = ptr->colIndex[colIndex];
	SparseDoubleElement *del;
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
		freePidSparseDoubleElement(del,pid);
		ptr->nnz--;
	}
}




static void delSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex)
{
	delPidSparseDoubleMatrix(ptr,rowIndex,colIndex,0);
}
*/


static void delFastPidSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex, SparseDoubleElement **baseRow, SparseDoubleElement **baseCol,const int pid)
{
//	delPidSparseDoubleMatrix(ptr,rowIndex,colIndex,pid);
//	return;

	SparseDoubleElement *rowTarget = NULL;
	SparseDoubleElement *colTarget = NULL;
	SparseDoubleElement *del;
	
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
		freePidSparseDoubleElement(del,pid);
		ptr->nnz--;
	}
	*baseRow = rowTarget;
	*baseCol = colTarget;
}




void delPidSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex,const int pid)
{
	delFastPidSparseDoubleMatrix(ptr,rowIndex,colIndex,&ptr->rowCache[rowIndex],&ptr->colCache[colIndex],pid);
}


void delSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex)
{
	delPidSparseDoubleMatrix(ptr,rowIndex,colIndex,0);
}

//==============================================================



void dumpSparseDoubleMatrix(FILE *fp,const SparseDoubleMatrix *ptr)
{
	int i;
//	fprintf(fp,"totalRow:%d ,totalCol:%d\n",ptr->totalRow,ptr->totalCol);
//	fprintf(fp,"nnz= %lld\n",ptr->nnz);
	fprintf(fp,"%d %d %lld\n",ptr->totalRow,ptr->totalCol,ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseDoubleElement *target = ptr->rowIndex[i]->rowLink;
		while(target!=NULL)
		{
//			fprintf(fp,"row: %d, col:%d, ",target->row,target->col);
//			fprintf(fp,"node: %p, rowLink: %p, colLink: %p\n",target,target->rowLink,target->colLink);
//			fprintf(fp,"data = %g\n",target->data);

//			fprintf(fp,"%d %d %g\n",target->row,target->col,target->data);
			target = target->rowLink;
		}
	}
}



void plotSparseDoubleMatrix(FILE *fp, const SparseDoubleMatrix *ptr)
{
	int i,j;
	fprintf(fp,"totalRow:%d ,totalCol:%d\n",ptr->totalRow,ptr->totalCol);
	fprintf(fp,"nnz= %lld\n",ptr->nnz);
	char buf[ptr->totalCol];
	memset(buf,'-',ptr->totalCol);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseDoubleElement *target = ptr->rowIndex[i]->rowLink;
		while(target!=NULL)
		{
			buf[target->col] = 'x';
			target = target->rowLink;
		}
		for(j=0;j<ptr->totalCol;j++) fprintf(fp,"%c",buf[j]);
		fprintf(fp,"\n");
		memset(buf,'-',ptr->totalCol);
	}
}





void swapRowSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int row1,const int row2)
{
	SparseDoubleMatrix *tmpMatrix = createSparseDoubleMatrix(ptr->totalRow,ptr->totalCol);

	while(ptr->rowIndex[row1]->rowLink != NULL)
	{
		const int rowNew = row2;
		const int col = ptr->rowIndex[row1]->rowLink->col;
		const double tmp = ptr->rowIndex[row1]->rowLink->data;
		setSparseDoubleMatrix(tmpMatrix,tmp,rowNew,col);
		delSparseDoubleMatrix(ptr,row1,col);
	}
	while(ptr->rowIndex[row2]->rowLink != NULL)
	{
		const int rowNew = row1;
		const int col = ptr->rowIndex[row2]->rowLink->col;
		const double tmp = ptr->rowIndex[row2]->rowLink->data;
		setSparseDoubleMatrix(tmpMatrix,tmp,rowNew,col);
		delSparseDoubleMatrix(ptr,row2,col);
	}
	// dump tempMatrix->row1 to ptr->row1
	while(tmpMatrix->rowIndex[row1]->rowLink != NULL)
	{
		const int col = tmpMatrix->rowIndex[row1]->rowLink->col;
		const double tmp = tmpMatrix->rowIndex[row1]->rowLink->data;
		setSparseDoubleMatrix(ptr,tmp,row1,col);
		delSparseDoubleMatrix(tmpMatrix,row1,col);
	}
	// dump tempMatrix->row2 to ptr->row2
	while(tmpMatrix->rowIndex[row2]->rowLink != NULL)
	{
		const int col = tmpMatrix->rowIndex[row2]->rowLink->col;
		const double tmp = tmpMatrix->rowIndex[row2]->rowLink->data;
		setSparseDoubleMatrix(ptr,tmp,row2,col);
		delSparseDoubleMatrix(tmpMatrix,row2,col);
	}
	freeSparseDoubleMatrix(tmpMatrix);
}




void clearSparseDoubleMatrix(SparseDoubleMatrix *ptr)
{
	clearPidSparseDoubleMatrix(ptr,0);
}




void clearPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const int pid)
{
	// free the memory row by row
	int i;
	SparseDoubleElement *currentEntry;
	SparseDoubleElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freePidSparseDoubleElement(removeEntry,pid);
		}
	}

	for(i=0;i<ptr->totalRow;i++) ptr->rowIndex[i]->rowLink = NULL;
	for(i=0;i<ptr->totalCol;i++) ptr->colIndex[i]->colLink = NULL;
	
	for(i=0;i<ptr->totalRow;i++) ptr->rowCache[i] = NULL;
	for(i=0;i<ptr->totalCol;i++) ptr->colCache[i] = NULL;
	ptr->nnz = 0;
}




void copySparseDoubleMatrix(SparseDoubleMatrix *dest,const SparseDoubleMatrix *src)
{
	copyPidSparseDoubleMatrix(dest,src,0);
}





void copyPidSparseDoubleMatrix(SparseDoubleMatrix *dest,const SparseDoubleMatrix *src,const int pid)
{
	if(dest == src) return;

	int i;
	SparseDoubleElement **rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*src->totalRow,pid);
	for(i=0;i<src->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*src->totalCol,pid);
	for(i=0;i<src->totalCol;i++) colCache[i] = NULL;

	clearPidSparseDoubleMatrix(dest,pid);
	
	for(i=0;i<src->totalRow;i++)
	{
		const SparseDoubleElement *currentEntry = src->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			const int insRow = currentEntry->row;
			const int insCol = currentEntry->col;
			const double data = currentEntry->data;
//			setPidSparseDoubleMatrix(dest,data,insRow,insCol,pid);
			setFastPidSparseDoubleMatrix(dest,data,insRow,insCol,&rowCache[i],&colCache[insCol],pid);
			currentEntry = currentEntry->rowLink;
		}
	}
	
	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*src->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*src->totalCol,pid);
}





void addSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a,const SparseDoubleMatrix *b)
{
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleElement *entryA;
	SparseDoubleElement *entryB;

	int i;
	SparseDoubleElement **rowCache = getMempoolSet(sizeof(SparseDoubleElement *)*c->totalRow);
	for(i=0;i<c->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getMempoolSet(sizeof(SparseDoubleElement *)*c->totalCol);
	for(i=0;i<c->totalCol;i++) colCache[i] = NULL;

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
				const double targetA = entryA->data;
//				setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
				setFastSparseDoubleMatrix(cResult,targetA,i,entryA->col,&rowCache[i],&colCache[entryA->col]);
				entryA = entryA->rowLink;
			}
			else if(entryA->col > entryB->col)
			{
				const double targetB = entryB->data;
//				setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
				setFastSparseDoubleMatrix(cResult,targetB,i,entryB->col,&rowCache[i],&colCache[entryB->col]);
				entryB = entryB->rowLink;
			}
			else
			{
				const double targetA = entryA->data;
				const double targetB = entryB->data;
//				setSparseDoubleMatrix(cResult,targetA+targetB,i,entryA->col);
				setFastSparseDoubleMatrix(cResult,targetA+targetB,i,entryA->col,&rowCache[i],&colCache[entryA->col]);
				entryA = entryA->rowLink;
				entryB = entryB->rowLink;
			}
		}
		// insert the remaining of a
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
//			setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
			setFastSparseDoubleMatrix(cResult,targetA,i,entryA->col,&rowCache[i],&colCache[entryA->col]);
			entryA = entryA->rowLink;
		}
		// insert the remaining of b
		while(entryB!=NULL)
		{
			const double targetB = entryB->data;
//			setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
			setFastSparseDoubleMatrix(cResult,targetB,i,entryB->col,&rowCache[i],&colCache[entryB->col]);
			entryB = entryB->rowLink;
		}
	}

	retMempoolSet(rowCache,sizeof(SparseDoubleElement *)*c->totalRow);
	retMempoolSet(colCache,sizeof(SparseDoubleElement *)*c->totalCol);

	copySparseDoubleMatrix(c,cResult);
	
	freeSparseDoubleMatrix(cResult);
}



// c = a - b
void subSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a,const SparseDoubleMatrix *b)
{
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleElement *entryA;
	SparseDoubleElement *entryB;
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
				const double targetA = entryA->data;
				setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
				entryA = entryA->rowLink;
			}
			else if(entryA->col > entryB->col)
			{
				const double targetB = -1*entryB->data;
				setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
				entryB = entryB->rowLink;
			}
			else
			{
				const double targetA = entryA->data;
				const double targetB = entryB->data;
				setSparseDoubleMatrix(cResult,targetA-targetB,i,entryA->col);
				entryA = entryA->rowLink;
				entryB = entryB->rowLink;
			}
		}
		// insert the remaining of a
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
			entryA = entryA->rowLink;
		}
		// insert the remaining of b
		while(entryB!=NULL)
		{
			const double targetB = -1*entryB->data;
			setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
			entryB = entryB->rowLink;
		}
	}
	copySparseDoubleMatrix(c,cResult);
	
	freeSparseDoubleMatrix(cResult);
}



// c = a * b
void mulSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a, const SparseDoubleMatrix *b)
{
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(c->totalRow,c->totalCol);
	SparseDoubleElement *entryA;
	SparseDoubleElement *entryB;
	int i,j;
	for(i=0;i<a->totalRow;i++)
	{
		for(j=0;j<b->totalCol;j++)
		{
			int flag = 0;
			double tempSum = 0;
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
					const double targetA = entryA->data;
					const double targetB = entryB->data;
					tempSum = tempSum + targetA * targetB;
					entryA = entryA->rowLink;
					entryB = entryB->colLink;
					flag = 1;
				}
			}
			if(flag==1) setSparseDoubleMatrix(cResult,tempSum,i,j);
		}
	}
	copySparseDoubleMatrix(c,cResult);
	
	freeSparseDoubleMatrix(cResult);
}


// c = k * a , k is a scale
void scaleSparseDoubleMatrix(SparseDoubleMatrix *c,const double k,const SparseDoubleMatrix *a)
{
	int i;
	SparseDoubleElement *entryA;
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			setSparseDoubleMatrix(cResult,k*targetA,i,entryA->col);
			entryA = entryA->rowLink;
		}
	}
	copySparseDoubleMatrix(c,cResult);
	freeSparseDoubleMatrix(cResult);
}




// c = trans(A)
void transSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a)
{
	int i;
	SparseDoubleElement *entryA;

	clearSparseDoubleMatrix(c);

//	SparseDoubleElement **rowCache = getMempoolSet(sizeof(SparseDoubleElement *)*c->totalRow);
//	for(i=0;i<c->totalRow;i++) rowCache[i] = NULL;
//	SparseDoubleElement **colCache = getMempoolSet(sizeof(SparseDoubleElement *)*c->totalCol);
//	for(i=0;i<c->totalCol;i++) colCache[i] = NULL;

	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
//			setFastSparseDoubleMatrix(c,targetA,entryA->col,i,&rowCache[entryA->col],&colCache[i]);
			setSparseDoubleMatrix(c,targetA,entryA->col,i);
			entryA = entryA->rowLink;
		}
	}
	
//	retMempoolSet(rowCache,sizeof(SparseDoubleElement *)*c->totalRow);
//	retMempoolSet(colCache,sizeof(SparseDoubleElement *)*c->totalCol);
}



static void *blockMulVec(void *ptr)
{
	int i;
	struct ParOfMul2 *par = (struct ParOfMul2 *) ptr;
	const SparseDoubleMatrix *a = par->a;
	const double *b = par->b;
	double *c = par->c;

	SparseDoubleElement *entryA;
	double sum = 0;
	for(i = par->rowBegin ; i<par->rowEnd ; i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		sum = 0;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			const double targetB = b[entryA->col];
			sum = sum + targetA * targetB;
			entryA = entryA->rowLink;
		}
		c[i] = sum;
	}
	pthread_exit(0);
}



void parallelMulVecSparseDoubleMatrix(double *c,const SparseDoubleMatrix *a, const double *b,const int threadNum)
{
	int i,j;
	memset(c,0,sizeof(double)*a->totalCol);

	const int offset = ceil((double)a->totalRow / (double) threadNum);
	int *start = getMempoolSet(sizeof(int)*(threadNum+1));
	for(i=0;i<threadNum;i++) start[i] = i*offset;
	start[i] = a->totalRow;

	struct ParOfMul2 *parList = getMempoolSet(sizeof(struct ParOfMul2)*threadNum);
	for(i=0;i<threadNum;i++)
	{
		parList[i].rowBegin = start[i];
		parList[i].rowEnd = start[i+1];
		parList[i].a = a;
		parList[i].b = b;
		parList[i].c = c;
	}

	pthread_t pid[threadNum];
	for(i=0;i<threadNum;i++) pthread_create(&pid[i],NULL,blockMulVec,&parList[i]);
	for(i=0;i<threadNum;i++) pthread_join(pid[i],NULL);

	retMempoolSet(parList,sizeof(struct ParOfMul2)*threadNum);
	retMempoolSet(start,sizeof(int)*(threadNum+1));
}





// c = a * b, (sparse multiply dense)
// a is a sparse matrix
// b is a n*1 QuadMatrix
// c is a n*1 QuadMatrix
void mulVecSparseDoubleMatrix(double *c,const SparseDoubleMatrix *a, const double *b)
{
	int i,j;
	SparseDoubleElement *entryA;
	double *cResult = getMempoolSet(sizeof(double)*a->totalCol);
	double sum = 0;

	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		sum = 0;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			const double targetB = b[entryA->col];
			sum = sum + targetA * targetB;
			entryA = entryA->rowLink;
		}
		cResult[i] = sum;
	}
	memcpy(c,cResult,sizeof(double)*a->totalRow);
	
	retMempoolSet(cResult,sizeof(double)*a->totalCol);
}





// set a as identity matrix
void identitySparseDoubleMatrix(SparseDoubleMatrix *a)
{
	identityPidSparseDoubleMatrix(a,0);
}



void identityPidSparseDoubleMatrix(SparseDoubleMatrix *a,const int pid)
{
	clearPidSparseDoubleMatrix(a,pid);
	int i;
	for(i=0;i<a->totalRow;i++) setPidSparseDoubleMatrix(a,1.0,i,i,pid);
}



// a(i,j) = a(i,j) + element
// will allocate the memory if necessary
void incSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex)
{
//	const double ptr = getSparseDoubleMatrix(a,row,col);
//	setSparseDoubleMatrix(a,ptr+element,row,col);
	
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];
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
		SparseDoubleElement *insert = createSparseDoubleElement(element);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data += element;
}


// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
void decSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex)
{
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];
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
		SparseDoubleElement *insert = createSparseDoubleElement(-1*element);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data -= element;
	
//	const double ptr = getSparseDoubleMatrix(a,row,col);
//	setSparseDoubleMatrix(a,ptr-element,row,col);
}



// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
SparseDoubleElement* decFastSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,SparseDoubleElement *baseRow,SparseDoubleElement *baseCol)
{
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];

	if(baseRow->col < colIndex) rowTarget = baseRow;
//	if(baseCol!=NULL) colTarget = baseCol;

	if(baseCol!=NULL)
	{
		colTarget = baseCol;
		if(baseCol->colLink!= NULL)
		{
			if(baseCol->colLink->row == rowIndex && baseCol->colLink->col ==  rowTarget->col)
			{
				baseCol->colLink->data -= element;
				return baseCol->colLink;
			}
		}
	}


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
		SparseDoubleElement *insert = createSparseDoubleElement(-1*element);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
		return insert;
	}
	else
	{
		rowTarget->rowLink->data -= element;
		return rowTarget->rowLink;
	}
	
//	const double ptr = getSparseDoubleMatrix(a,row,col);
//	setSparseDoubleMatrix(a,ptr-element,row,col);
}





void dense2SparseDoubleMatrix(SparseDoubleMatrix *dest,const double *src)
{
	dense2PidSparseDoubleMatrix(dest,src,0);
}






void dense2PidSparseDoubleMatrix(SparseDoubleMatrix *dest,const double *src,const int pid)
{
	clearPidSparseDoubleMatrix(dest,pid);
	const int row = dest->totalRow;
	const int col = dest->totalCol;
	int i,j,k;
	k = 0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			const double tmp = src[k];
			k++;
			if(tmp!=0)
			{
				setPidSparseDoubleMatrix(dest,tmp,i,j,pid);
			}
		}
	}
}





void sparse2DenseDoubleMatrix(double *dest,const SparseDoubleMatrix *src)
{
	memset(dest,0,sizeof(double)*src->totalRow*src->totalCol);
	int i,j;
	for(i=0;i<src->totalRow;i++)
	{
		SparseDoubleElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
			const double data = ptr->data;
			setMyMatrix(dest,data,src->totalRow,src->totalCol,i,col);
			ptr = ptr->rowLink;
		}
	}
}




void setDenseCol2SparseDoubleMatrix(SparseDoubleMatrix *dest, const double *src,const int totalRow, const int targetCol)
{
	int i,j;
	for(i=0;i<totalRow;i++)
	{
		const double element = src[i];
		if( element == 0 ) delSparseDoubleMatrix(dest,i,targetCol);
		else setSparseDoubleMatrix(dest,element,i,targetCol);
	}
}



void setDenseColQuick2SparseDoubleMatrix(SparseDoubleMatrix *dest, const double *src,const int totalRow, const int targetCol)
{
	int i,j;
	for(i=0;i<totalRow;i++)
	{
		const double element = src[i];
		if(element==0) continue;
		else setSparseDoubleMatrix(dest,element,i,targetCol);
	}
}





// the boundaries row1,col1,row2,col2 are included
void clearBlockSparseDoubleMatrix(SparseDoubleMatrix *dest, const int row1, const int col1, const int row2, const int col2)
{
	int i;
	// del the old data in the block
	for(i=row1;i<=row2;i++)
	{
		SparseDoubleElement *eachRow = dest->rowIndex[i]->rowLink;
		while(eachRow != NULL)
		{
			if(eachRow->col < col1) eachRow = eachRow->rowLink;
			else if(eachRow->col > col2) break;
			else
			{
				const int delRow = eachRow->row;
				const int delCol = eachRow->col;
				eachRow = eachRow->rowLink;
				delSparseDoubleMatrix(dest,delRow,delCol);
			}
		}
	}
}




// the boundaries row1,col1,row2,col2 are included
void appendBlockSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int row1, const int col1, const int row2, const int col2)
{
	if(dest == src) return;

	int i;


	// insert new data to the block
	for(i=row1;i<=row2;i++)
	{
		SparseDoubleElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow != NULL)
		{
			if(eachRow->col < col1) eachRow = eachRow->rowLink;
			else if(eachRow->col > col2) break;
			else
			{
				const double data = eachRow->data;
				const int insRow = eachRow->row;
				const int insCol = eachRow->col;
				setSparseDoubleMatrix(dest,data,insRow,insCol);
				eachRow = eachRow->rowLink;
			}
		}
	}
}




void mergeSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest)
{
	mergePidSparseDoubleMatrix(dest,src,destRow,destCol,ltRowDest,ltColDest,0);
}





void mergePidSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest,const int pid)
{
	int i;
	const int baseRow = ltRowDest;
	const int baseCol = ltColDest;
/*	
	SparseDoubleElement **rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	for(i=0;i<dest->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalCol,pid);
	for(i=0;i<dest->totalCol;i++) colCache[i] = NULL;
*/
	for(i=0;i<src->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = src->rowIndex[i]->rowLink;
		int insRow = 0;
		if(eachRow!=NULL) insRow = baseRow + eachRow->row; // must put here ... sometimes the whole row may be empty..
		while(eachRow!=NULL)
		{
			const double data = eachRow->data;	
			const int insCol = baseCol + eachRow->col;
//			setFastPidSparseDoubleMatrix(dest,data,insRow,insCol,&rowCache[insRow],&colCache[insCol],pid);
			setPidSparseDoubleMatrix(dest,data,insRow,insCol,pid);
			eachRow = eachRow->rowLink;
		}
	}
//	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*dest->totalRow,pid);
//	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*dest->totalCol,pid);
}




void getSubSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc)
{
	getPidSubSparseDoubleMatrix(dest,src,ltRowSrc,ltColSrc,rbRowSrc,rbColSrc,0);
}





void getPidSubSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc,const int pid)
{
	clearPidSparseDoubleMatrix(dest,pid);
	int i;
	const int baseRow = ltRowSrc;
	const int baseCol = ltColSrc;

	SparseDoubleElement **rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	for(i=0;i<dest->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalCol,pid);
	for(i=0;i<dest->totalCol;i++) colCache[i] = NULL;

	for(i=ltRowSrc;i<=rbRowSrc;i++)
	{
		const SparseDoubleElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			if(eachRow->col < ltColSrc) eachRow = eachRow->rowLink;
			else if(eachRow->col > rbColSrc) break;
			else
			{
				const double data = eachRow->data;
				const int insRow = eachRow->row - baseRow;
				const int insCol = eachRow->col - baseCol;
//				setPidSparseDoubleMatrix(dest,data,insRow,insCol,pid);
				setFastPidSparseDoubleMatrix(dest,data,insRow,insCol,&rowCache[insRow],&colCache[insCol],pid);
				eachRow = eachRow->rowLink;
			}
		}
	}
	
	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*dest->totalCol,pid);
}


/*
// permutate 
void permutateSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *pRow, const SparseDoubleMatrix *pCol, const SparseDoubleMatrix *src)
{
	clearSparseDoubleMatrix(dest);
	SparseDoubleMatrix *destTemp = createSparseDoubleMatrix(dest->totalRow,dest->totalCol);
	int i;

	// permutate row
	for(i=0;i<pRow->totalRow;i++)
	{
		const int colForRowI = pRow->rowIndex[i]->rowLink->col;
		const SparseDoubleElement *ptr = src->rowIndex[colForRowI]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
			setSparseDoubleMatrix(destTemp,ptr->data,i,col);			
			ptr = ptr->rowLink;
		}
	}

	// permutate col
	for(i=0;i<pCol->totalCol;i++)
	{
		const int rowForColI = pCol->colIndex[i]->colLink->row;
		const SparseDoubleElement *ptr = destTemp->colIndex[rowForColI]->colLink;
		while(ptr!=NULL)
		{
			const int row = ptr->row;
			setSparseDoubleMatrix(dest,ptr->data,row,i);
			ptr = ptr->colLink;
		}
	}
	
	freeSparseDoubleMatrix(destTemp);
	clearPidMempoolSet(0);
}
*/




// permutate 
void permutateSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *pRow, const SparseDoubleMatrix *pCol,const SparseDoubleMatrix *src)
{
	clearSparseDoubleMatrix(dest);
	
	int i;
	for(i=0;i<src->totalRow;i++)
	{
		SparseDoubleElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = pRow->colIndex[i]->colLink->row;
			const int col = pCol->rowIndex[ptr->col]->rowLink->col;
			setSparseDoubleMatrix(dest,ptr->data,row,col);
			ptr = ptr->rowLink;
		}
	}
}





void dpermutateSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *pRow, const SparseDoubleMatrix *pCol,SparseDoubleMatrix *src)
{
	clearSparseDoubleMatrix(dest);
	
	int i;
	for(i=0;i<src->totalRow;i++)
	{
		SparseDoubleElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int row = pRow->colIndex[i]->colLink->row;
			const int col = pCol->rowIndex[ptr->col]->rowLink->col;
			setSparseDoubleMatrix(dest,ptr->data,row,col);

			SparseDoubleElement *del = ptr;
			ptr = ptr->rowLink;
			freeSparseDoubleElement(del);
		}
	}

	for(i=0;i<src->totalRow;i++) freeSparseDoubleElement(src->rowIndex[i]);
	retMempoolSet(src->rowIndex,sizeof(SparseDoubleElement *)*(src->totalRow));
	for(i=0;i<src->totalCol;i++) freeSparseDoubleElement(src->colIndex[i]);
	retMempoolSet(src->colIndex,sizeof(SparseDoubleElement *)*(src->totalCol));
	retMempoolSet(src->rowCache,sizeof(SparseDoubleElement *)*(src->totalRow));
	retMempoolSet(src->colCache,sizeof(SparseDoubleElement *)*(src->totalCol));
	retMempoolSet(src,sizeof(SparseDoubleMatrix));
}





static void updateInvL(SparseDoubleMatrix *invL,const double scale,const int baseRow,const int targetRow)
{
	SparseDoubleElement *baseRowPtr = invL->rowIndex[baseRow]->rowLink;
	while(baseRowPtr!=NULL)
	{
		const double elementSrc = baseRowPtr->data;
		const double elementTarget = getSparseDoubleMatrix(invL,targetRow,baseRowPtr->col,"row");
		const double element = scale*elementSrc + elementTarget;
		setSparseDoubleMatrix(invL,element,targetRow,baseRowPtr->col);
		baseRowPtr = baseRowPtr->rowLink;
	}
}




// inverse lowerTriangular
void invLTSparseDoubleMatrix(SparseDoubleMatrix *invL, const SparseDoubleMatrix *lSrc)
{
	SparseDoubleMatrix *l = createSparseDoubleMatrix(lSrc->totalRow,lSrc->totalCol);
	copySparseDoubleMatrix(l,lSrc);
	identitySparseDoubleMatrix(invL);
	int i;
	for(i = 0; i<l->totalRow;i++)
	{
		SparseDoubleElement *eachRowL = l->rowIndex[i]->rowLink;
		const double base = eachRowL->data;
		eachRowL = eachRowL->colLink;
		while(eachRowL != NULL)
		{
			double element = eachRowL->data;
			double scale = -1.0 * (element/base);
//			fprintf(stderr,"scale:%g,baseRow:%d,targetRow:%d\n",scale,i,eachRowL->row);
			// update invL
			updateInvL(invL,scale,i,eachRowL->row);
//			dumpSparseDoubleMatrix(stderr,invL);
			// updata L
			if(scale!=0)
			{
				SparseDoubleElement *inRow = l->rowIndex[i]->rowLink->rowLink;
				while(inRow != NULL)
				{
					double elementIn = scale * getSparseDoubleMatrix(l,i,inRow->col,"row");
					elementIn = getSparseDoubleMatrix(l,eachRowL->row,inRow->col,"row") + elementIn;
					setSparseDoubleMatrix(l,elementIn,eachRowL->row,inRow->col);
					inRow = inRow->rowLink;
				}
			}
			SparseDoubleElement *next = eachRowL->colLink;
			delSparseDoubleMatrix(l,eachRowL->row,i);
			eachRowL = next;
		}
	}
	freeSparseDoubleMatrix(l);
}




/*
SparseDoubleMatrix *read_pid_SparseDoubleMatrix(const char *filename,const int pid)
{
	// init
	int i;
	FILE *fp = fopen(filename,"rb");
	int ret;
	int *row_ptr = NULL;
	int totalRow;
	int totalCol;
	long long nnz;
	// head
	ret = fread(&totalRow,sizeof(int),1,fp);	
	ret = fread(&totalCol,sizeof(int),1,fp);	
	ret = fread(&nnz,sizeof(long long),1,fp);	
	SparseDoubleMatrix *ptr = createPidSparseDoubleMatrix(totalRow,totalCol,pid);
	// row_ptr
	row_ptr = getPidMempoolSet(sizeof(int)*(ptr->totalRow+1),pid);
	ret = fread(row_ptr,sizeof(int),ptr->totalRow+1,fp);
	// col
	for(i=0;i<ptr->totalRow;i++)
	{
		int j = 0;
		const int length = row_ptr[i+1] - row_ptr[i];
		int *col_list = getPidMempoolSet(sizeof(int)*length,pid);
		ret = fread(col_list,sizeof(int),length,fp);
		for(j=0;j<length;j++)
		{
			const int row = i;
			const int col = col_list[j];
			setPidSparseDoubleMatrix(ptr,0.0,row,col,pid);
		}
		retPidMempoolSet(col_list,sizeof(int)*length,pid);
	}
	// values
	for(i=0;i<ptr->totalRow;i++)
	{
		int j = 0;
		const int length = row_ptr[i+1] - row_ptr[i];
		double *val_list = getPidMempoolSet(sizeof(double)*length,pid);
		ret = fread(val_list,sizeof(double),length,fp);
		SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			eachRow->data = val_list[j];
			j++;
			eachRow = eachRow->rowLink;
		}
		retPidMempoolSet(val_list,sizeof(double)*length,pid);
	}
	// final
	fclose(fp);
	retPidMempoolSet(row_ptr,sizeof(int)*(ptr->totalRow+1),pid);
	return ptr;
}
*/




SparseDoubleMatrix *read_pid_SparseDoubleMatrix(const char *prefix,const int pid)
{
	// init
	char sizeName[32] = {0};
	char rowName[32] = {0};
	char colName[32] = {0};
	char valName[32] = {0};
	sprintf(sizeName,"%s.size",prefix);
	sprintf(rowName,"%s.row",prefix);
	sprintf(colName,"%s.col",prefix);
	sprintf(valName,"%s.val",prefix);

	int i;
	int ret;
	int totalRow;
	int totalCol;
	long long nnz;
	// size
	FILE *fp_size = fopen(sizeName,"rb");
	ret = fread(&totalRow,sizeof(int),1,fp_size);	
	ret = fread(&totalCol,sizeof(int),1,fp_size);	
	ret = fread(&nnz,sizeof(long long),1,fp_size);	
	SparseDoubleMatrix *ptr = createPidSparseDoubleMatrix(totalRow,totalCol,pid);
	fclose(fp_size);
	// row_ptr
	FILE *fp_row = fopen(rowName,"rb");
	int *row_ptr = getPidMempoolSet(sizeof(int)*(ptr->totalRow+1),pid);
	ret = fread(row_ptr,sizeof(int),ptr->totalRow+1,fp_row);
	fclose(fp_row);
	// col
	FILE *fp_col = fopen(colName,"rb");
	for(i=0;i<ptr->totalRow;i++)
	{
		int j = 0;
		const int length = row_ptr[i+1] - row_ptr[i];
		int *col_list = getPidMempoolSet(sizeof(int)*length,pid);
		ret = fread(col_list,sizeof(int),length,fp_col);
		for(j=0;j<length;j++)
		{
			const int row = i;
			const int col = col_list[j];
			setPidSparseDoubleMatrix(ptr,0.0,row,col,pid);
		}
		retPidMempoolSet(col_list,sizeof(int)*length,pid);
	}
	fclose(fp_col);
	// values
	FILE *fp_val = fopen(valName,"rb");
	for(i=0;i<ptr->totalRow;i++)
	{
		int j = 0;
		const int length = row_ptr[i+1] - row_ptr[i];
		double *val_list = getPidMempoolSet(sizeof(double)*length,pid);
		ret = fread(val_list,sizeof(double),length,fp_val);
		SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			eachRow->data = val_list[j];
			j++;
			eachRow = eachRow->rowLink;
		}
		retPidMempoolSet(val_list,sizeof(double)*length,pid);
	}
	fclose(fp_val);
	// final
	retPidMempoolSet(row_ptr,sizeof(int)*(ptr->totalRow+1),pid);
	return ptr;
}





/*
void write_SparseDoubleMatrix(const char *filename,const SparseDoubleMatrix *ptr)
{
	// init
	int ret;
	int i;
	FILE *fp = fopen(filename,"wb");
	int *row_ptr = getMempoolSet(sizeof(int)*(ptr->totalRow+1));
	// head
	ret = fwrite(&ptr->totalRow,sizeof(int),1,fp);
	ret = fwrite(&ptr->totalCol,sizeof(int),1,fp);
	ret = fwrite(&ptr->nnz,sizeof(long long),1,fp);
	// row_ptr
	row_ptr[0] = 0;
	for(i=0;i<ptr->totalRow;i++)
	{
		row_ptr[i+1] = row_ptr[i];
		const SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			row_ptr[i+1]++;
			eachRow = eachRow->rowLink;
		}
	}
	ret = fwrite(row_ptr,sizeof(int),ptr->totalRow+1,fp);
	retMempoolSet(row_ptr,sizeof(int)*(ptr->totalRow+1));
	// col
	// values	
	long long col_index = 0;
	int *col_list = malloc(sizeof(int)*ptr->nnz);
	double *val_list = malloc(sizeof(double)*ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			col_list[col_index] = eachRow->col;
			val_list[col_index] = eachRow->data;
			col_index++;
			eachRow = eachRow->rowLink;
		}
	}
	ret = fwrite(col_list,sizeof(int),ptr->nnz,fp);
	ret = fwrite(val_list,sizeof(double),ptr->nnz,fp);
	free(col_list);
	free(val_list);
	// final
	fclose(fp);
}
*/




void write_SparseDoubleMatrix(const char *prefix,const SparseDoubleMatrix *ptr)
{
	// init
	char sizeName[32] = {0};
	char rowName[32] = {0};
	char colName[32] = {0};
	char valName[32] = {0};
	sprintf(sizeName,"%s.size",prefix);
	sprintf(rowName,"%s.row",prefix);
	sprintf(colName,"%s.col",prefix);
	sprintf(valName,"%s.val",prefix);

	int ret;
	int i;
	// size
	FILE *fp_size = fopen(sizeName,"wb");
	ret = fwrite(&ptr->totalRow,sizeof(int),1,fp_size);
	ret = fwrite(&ptr->totalCol,sizeof(int),1,fp_size);
	ret = fwrite(&ptr->nnz,sizeof(long long),1,fp_size);
	fclose(fp_size);
	// row_ptr
	FILE *fp_row = fopen(rowName,"wb");
	int *row_ptr = getMempoolSet(sizeof(int)*(ptr->totalRow+1));
	row_ptr[0] = 0;
	for(i=0;i<ptr->totalRow;i++)
	{
		row_ptr[i+1] = row_ptr[i];
		const SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			row_ptr[i+1]++;
			eachRow = eachRow->rowLink;
		}
	}
	ret = fwrite(row_ptr,sizeof(int),ptr->totalRow+1,fp_row);
	retMempoolSet(row_ptr,sizeof(int)*(ptr->totalRow+1));
	fclose(fp_row);
	// col
	// values
	FILE *fp_col = fopen(colName,"wb");
	FILE *fp_val = fopen(valName,"wb");
	long long col_index = 0;
	int *col_list = getMempoolSet(sizeof(int)*ptr->nnz);
	double *val_list = getMempoolSet(sizeof(double)*ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			col_list[col_index] = eachRow->col;
			val_list[col_index] = eachRow->data;
			col_index++;
			eachRow = eachRow->rowLink;
		}
	}
	ret = fwrite(col_list,sizeof(int),ptr->nnz,fp_col);
	ret = fwrite(val_list,sizeof(double),ptr->nnz,fp_val);
	retMempoolSet(col_list,sizeof(int)*ptr->nnz);
	retMempoolSet(val_list,sizeof(double)*ptr->nnz);
	fclose(fp_col);
	fclose(fp_val);
}





SparseDoubleMatrix *read_ind_SparseDoubleMatrix(const char *filename)
{
	// init
	long long i;
	int ret;
	FILE *fp = fopen(filename,"r");
	int totalRow;
	int totalCol;
	long long nnz;
	SparseDoubleElement **rowCache = NULL;
	SparseDoubleElement **colCache = NULL;
	// head
	ret = fscanf(fp,"%d %d %lld\n",&totalRow,&totalCol,&nnz);
	SparseDoubleMatrix *ptr = createSparseDoubleMatrix(totalRow,totalCol);
	
	int insRow;
	int insCol;
	double val;
	for(i=0;i<nnz;i++)
	{
		ret = fscanf(fp,"%d %d %lf\n",&insRow,&insCol,&val);
		setSparseDoubleMatrix(ptr,val,insRow,insCol);
	}

	fclose(fp);
	return ptr;
}



void write_ind_SparseDoubleMatrix(const char *filename,const SparseDoubleMatrix *ptr)
{
	int i;
	FILE *fp = fopen(filename,"w");
	fprintf(fp,"%d %d %lld\n",ptr->totalRow,ptr->totalCol,ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			fprintf(fp,"%d %d %lf\n",i,eachRow->col,eachRow->data);
			eachRow = eachRow->rowLink;
		}
	}
	fclose(fp);
}





double colNormSparseDoubleMatrix(const SparseDoubleMatrix *ptr,const int col)
{
	double sum = 0.0;
	const SparseDoubleElement *colPtr = ptr->colIndex[col]->colLink;
	while(colPtr!=NULL)
	{
//		sum += pow(colPtr->data,2);
		sum += colPtr->data * colPtr->data;
		colPtr = colPtr->colLink;
	}
	return sqrt(sum);
}



// ===================================================

CSR_SparseDoubleMatrix *create_CSR_SparseDoubleMatrix(int row, int col, long long nnz,int pid)
{
	CSR_SparseDoubleMatrix *mtx = getPidMempoolSet(sizeof(CSR_SparseDoubleMatrix),pid);

	mtx->pid = pid;
	mtx->totalRow = row;
	mtx->totalCol = col;
	mtx->nnz = nnz;

	mtx->rowPtr = getPidMempoolSet(sizeof(int)*(row+1),pid);
//	mtx->rowPtr = malloc(sizeof(int)*(row+1));
	memset(mtx->rowPtr,0,sizeof(int)*(row+1));
	
	mtx->col = getPidMempoolSet(sizeof(int)*nnz,pid);
//	mtx->col = malloc(sizeof(int)*nnz);
	memset(mtx->col,0,sizeof(int)*nnz);
	
	mtx->val = getPidMempoolSet(sizeof(double)*nnz,pid);
//	mtx->val = malloc(sizeof(double)*nnz);
	memset(mtx->val,0,sizeof(double)*nnz);

	return mtx;
}




void free_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx)
{
	retPidMempoolSet(mtx->rowPtr,sizeof(int)*(mtx->totalRow+1),mtx->pid);
	retPidMempoolSet(mtx->col,sizeof(int)*(mtx->nnz),mtx->pid);
	retPidMempoolSet(mtx->val,sizeof(double)*(mtx->nnz),mtx->pid);
//	free(mtx->rowPtr);
//	free(mtx->col);
//	free(mtx->val);

	retPidMempoolSet(mtx,sizeof(CSR_SparseDoubleMatrix),mtx->pid);
}



/*
CSR_SparseDoubleMatrix *read_to_CSR_SparseDoubleMatrix(const char *filename,const int pid)
{
	int i;
	int totalRow;
	int totalCol;
	int ret;
	long long nnz;
	FILE *fp = fopen(filename,"rb");
	// head
	ret = fread(&totalRow,sizeof(int),1,fp);	
	ret = fread(&totalCol,sizeof(int),1,fp);	
	ret = fread(&nnz,sizeof(long long),1,fp);	
	
	CSR_SparseDoubleMatrix *mtx = create_CSR_SparseDoubleMatrix(totalRow,totalCol,nnz,pid);
	ret = fread(mtx->rowPtr,sizeof(int),totalRow+1,fp);
	ret = fread(mtx->col,sizeof(int),nnz,fp);
	ret = fread(mtx->val,sizeof(double),nnz,fp);

	fclose(fp);
	return mtx;
}
*/




CSR_SparseDoubleMatrix *read_to_CSR_SparseDoubleMatrix(const char *prefix,const int pid)
{
	char sizeName[32] = {0};
	char rowName[32] = {0};
	char colName[32] = {0};
	char valName[32] = {0};
	sprintf(sizeName,"%s.size",prefix);
	sprintf(rowName,"%s.row",prefix);
	sprintf(colName,"%s.col",prefix);
	sprintf(valName,"%s.val",prefix);

	int i;
	int totalRow;
	int totalCol;
	int ret;
	long long nnz;
	
	FILE *fp_size = fopen(sizeName,"rb");
	ret = fread(&totalRow,sizeof(int),1,fp_size);	
	ret = fread(&totalCol,sizeof(int),1,fp_size);	
	ret = fread(&nnz,sizeof(long long),1,fp_size);	
	fclose(fp_size);
	
	CSR_SparseDoubleMatrix *mtx = create_CSR_SparseDoubleMatrix(totalRow,totalCol,nnz,pid);

	FILE *fp_row = fopen(rowName,"rb");
	ret = fread(mtx->rowPtr,sizeof(int),totalRow+1,fp_row);
	fclose(fp_row);

	FILE *fp_col = fopen(colName,"rb");
	ret = fread(mtx->col,sizeof(int),nnz,fp_col);
	fclose(fp_col);

	FILE *fp_val = fopen(valName,"rb");
	ret = fread(mtx->val,sizeof(double),nnz,fp_val);
	fclose(fp_val);

	return mtx;
}






/*
CSR_SparseDoubleMatrix *write_CSR_SparseDoubleMatrix(const char *filename,CSR_SparseDoubleMatrix *mtx)
{
	int ret;
	FILE *fp = fopen(filename,"wb");
	// head
	ret = fwrite(&mtx->totalRow,sizeof(int),1,fp);
	ret = fwrite(&mtx->totalCol,sizeof(int),1,fp);
	ret = fwrite(&mtx->nnz,sizeof(long long),1,fp);

	ret = fwrite(mtx->rowPtr,sizeof(int),mtx->totalRow+1,fp);
	ret = fwrite(mtx->col,sizeof(int),mtx->nnz,fp);
	ret = fwrite(mtx->val,sizeof(double),mtx->nnz,fp);

	fclose(fp);
}
*/





CSR_SparseDoubleMatrix *write_CSR_SparseDoubleMatrix(const char *prefix,CSR_SparseDoubleMatrix *mtx)
{
	char sizeName[32] = {0};
	char rowName[32] = {0};
	char colName[32] = {0};
	char valName[32] = {0};
	sprintf(sizeName,"%s.size",prefix);
	sprintf(rowName,"%s.row",prefix);
	sprintf(colName,"%s.col",prefix);
	sprintf(valName,"%s.val",prefix);

	int ret;

	FILE *fp_size = fopen(sizeName,"wb");
	ret = fwrite(&mtx->totalRow,sizeof(int),1,fp_size);
	ret = fwrite(&mtx->totalCol,sizeof(int),1,fp_size);
	ret = fwrite(&mtx->nnz,sizeof(long long),1,fp_size);
	fclose(fp_size);

	FILE *fp_row = fopen(rowName,"wb");
	ret = fwrite(mtx->rowPtr,sizeof(int),mtx->totalRow+1,fp_row);
	fclose(fp_row);

	FILE *fp_col = fopen(colName,"wb");
	ret = fwrite(mtx->col,sizeof(int),mtx->nnz,fp_col);
	fclose(fp_col);

	FILE *fp_val = fopen(valName,"wb");
	ret = fwrite(mtx->val,sizeof(double),mtx->nnz,fp_val);
	fclose(fp_val);
}






void expand_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx,const long long extra)
{
	const long long newNNZ = extra + mtx->nnz;
	//===============================================
	int *new_col = getPidMempoolSet(sizeof(int)*newNNZ,mtx->pid);
	memset(new_col,0,sizeof(int)*newNNZ);
	
	if(mtx->nnz < newNNZ) memcpy(new_col,mtx->col,sizeof(int)*mtx->nnz);
	else memcpy(new_col,mtx->col,sizeof(int)*newNNZ);
	
	retPidMempoolSet(mtx->col,sizeof(int)*mtx->nnz,mtx->pid);
	mtx->col = new_col;
	//===============================================
	double *new_val = getPidMempoolSet(sizeof(double)*newNNZ,mtx->pid);
	memset(new_val,0,sizeof(double)*newNNZ);
	
	if(mtx->nnz < newNNZ) memcpy(new_val,mtx->val,sizeof(double)*mtx->nnz);
	else memcpy(new_val,mtx->val,sizeof(double)*newNNZ);
	
	retPidMempoolSet(mtx->val,sizeof(double)*mtx->nnz,mtx->pid);
	mtx->val = new_val;
	//===============================================
	mtx->nnz = newNNZ;
}




void copy_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *dest, const CSR_SparseDoubleMatrix *src)
{
	if(dest == src) return;

	assert(dest->totalRow == src->totalRow);
	assert(dest->totalCol == src->totalCol);

	if(dest->nnz != src->nnz)
	{
		memcpy(dest->rowPtr,src->rowPtr,sizeof(int)*(src->totalRow+1));
		retPidMempoolSet(dest->col,sizeof(int)*dest->nnz,dest->pid);
		retPidMempoolSet(dest->val,sizeof(double)*dest->nnz,dest->pid);
		
		int *new_col = getPidMempoolSet(sizeof(int)*src->nnz,dest->pid);
		memcpy(new_col,src->col,sizeof(int)*src->nnz);

		double *new_val = getPidMempoolSet(sizeof(double)*src->nnz,dest->pid);
		memcpy(new_val,src->val,sizeof(double)*src->nnz);
		
		dest->col = new_col;
		dest->val = new_val;
		dest->nnz = src->nnz;
	}
	else
	{
		memcpy(dest->rowPtr,src->rowPtr,sizeof(int)*(src->totalRow+1));
		memcpy(dest->col,src->col,sizeof(int)*(src->nnz));
		memcpy(dest->val,src->val,sizeof(double)*(src->nnz));
	}
}





void dump_CSR_SparseDoubleMatrix(FILE *fp,const  CSR_SparseDoubleMatrix *mtx)
{
	fprintf(fp,"row = %d col = %d nnz = %lld true nnz = %d\n",mtx->totalRow,mtx->totalCol,mtx->nnz,mtx->rowPtr[mtx->totalRow]);
	int i;
	fprintf(fp,"rowPtr = ");
	for(i=0;i<mtx->totalRow+1;i++) fprintf(fp,"%d ",mtx->rowPtr[i]);
	fprintf(fp,"\n");

	fprintf(fp,"col = ");
	for(i=0;i<mtx->rowPtr[mtx->totalRow];i++) fprintf(fp,"%d ",mtx->col[i]);
	fprintf(fp,"\n");
	
	fprintf(fp,"vall = ");
	for(i=0;i<mtx->rowPtr[mtx->totalRow];i++) fprintf(fp,"%lf ",mtx->val[i]);
	fprintf(fp,"\n");
}




static void expand_if_necessary(CSR_SparseDoubleMatrix *mtx,const int new_dummy)
{
	clock_t t1,t2;
	t1 = clock();
	if( (mtx->nnz-mtx->rowPtr[mtx->totalRow]) < new_dummy) expand_CSR_SparseDoubleMatrix(mtx,new_dummy);
	t2 = clock();
//	fprintf(stderr,"\texpand: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
}




static void setLCol(CSR_SparseDoubleMatrix *l,const int i,const int *eliminate_row,const int nnzInBuf)
{
	int check = 0;
	int ind = l->rowPtr[i];
	int k,m;
	memmove(&l->col[l->rowPtr[i]+nnzInBuf],&l->col[l->rowPtr[i]],sizeof(int)*(l->rowPtr[l->totalRow]-l->rowPtr[i]));
	for(m=i;m<l->totalRow;m++)
	{
		memmove(&l->col[ind],&l->col[l->rowPtr[m]+nnzInBuf],sizeof(int)*(l->rowPtr[m+1]-l->rowPtr[m]));
		ind += (l->rowPtr[m+1] - l->rowPtr[m]);
		if(eliminate_row[check] == m && check < nnzInBuf)
		{
			l->col[ind] = i;
			ind++;
			check++;
		}
	}
}










static void setLVal(CSR_SparseDoubleMatrix *l,const int i,const int *eliminate_row,const double *eliminate_val,const int nnzInBuf)
{
	int check = 0;
	int ind = l->rowPtr[i];
	int k,m;
	memmove(&l->val[l->rowPtr[i]+nnzInBuf],&l->val[l->rowPtr[i]],sizeof(double)*(l->rowPtr[l->totalRow]-l->rowPtr[i]));
	for(m=i;m<l->totalRow;m++)
	{
		memcpy(&l->val[ind],&l->val[l->rowPtr[m]+nnzInBuf],sizeof(double)*(l->rowPtr[m+1]-l->rowPtr[m]));
		ind += (l->rowPtr[m+1]-l->rowPtr[m]);
		if(eliminate_row[check] == m && check < nnzInBuf)
		{
			l->val[ind] = eliminate_val[check];
			ind++;
			check++;
		}
	}
}







static void updateL(CSR_SparseDoubleMatrix *l,const int i,const int *eliminate_row, const double *eliminate_val,const int nnzInBuf)
{
	clock_t t1,t2;

	t1 = clock();
	expand_if_necessary(l,nnzInBuf);

	setLCol(l,i,eliminate_row,nnzInBuf);
	setLVal(l,i,eliminate_row,eliminate_val,nnzInBuf);

	int k;
	int add = 0;
	int check_ind = 0;
	for(k=0;k<l->totalRow;k++)
	{
		if(eliminate_row[check_ind] == k && check_ind < nnzInBuf)
		{
			add++;
			check_ind++;
		}
		l->rowPtr[k+1] += add;
	}
	t2 = clock();
//	fprintf(stderr,"\tupdateL: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
}








static void dumpBufDouble(FILE *fp, const char *name, const double *buf, const int bufSize)
{
	int i;
	fprintf(fp,"%s:",name);
	for(i=0;i<bufSize;i++) fprintf(fp,"%lf ",buf[i]);
	fprintf(fp,"\n");
}




static void dumpBufInt(FILE *fp, const char *name, const int *buf, const int bufSize)
{
	int i;
	fprintf(fp,"%s:",name);
	for(i=0;i<bufSize;i++) fprintf(fp,"%d ",buf[i]);
	fprintf(fp,"\n");
}





// 0  5 84  86
// 0  3  5  85
// 0  3  5  84  85  86
static int getFlattenSize(const CSR_SparseDoubleMatrix *mtx, const int m, const int n)
{
	clock_t t1,t2;
	t1 = clock();
	// get cross size
	int crossSize = 0;
	int i = mtx->rowPtr[m];
	int j = mtx->rowPtr[n];
	while(i!=mtx->rowPtr[m+1] && j!=mtx->rowPtr[n+1])
	{
		if(mtx->col[i] > mtx->col[j])
		{
			j++;
		}
		else if(mtx->col[i] < mtx->col[j])
		{
			i++;
		}
		else
		{
			crossSize++;
			i++;
			j++;
		}
	}
	t2 = clock();
//	fprintf(stderr,"\tget flatten size: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
	
	return (mtx->rowPtr[m+1]-mtx->rowPtr[m]) + (mtx->rowPtr[n+1]-mtx->rowPtr[n]) - crossSize;
}





static void getFlattenCols(CSR_SparseDoubleMatrix *mtx,const int m,const int n, int *ret_col_buf)
{
	clock_t t1,t2;
	t1 = clock();
	int *col_buf = ret_col_buf;
	
	int i = mtx->rowPtr[m]+1;
	int j = mtx->rowPtr[n]+1;
	int k = 0;
	while(i!=mtx->rowPtr[m+1] && j!=mtx->rowPtr[n+1])
	{
		if(mtx->col[i] > mtx->col[j])
		{
			col_buf[k] = mtx->col[j];
			j++;
		}
		else if(mtx->col[i] < mtx->col[j])
		{
			col_buf[k] = mtx->col[i];
			i++;
		}
		else
		{
			col_buf[k] = mtx->col[i];
			i++;
			j++;
		}
		k++;
	}
	// append the remaining...
	if(i==mtx->rowPtr[m+1]) // m is done
	{
		if(j!=mtx->rowPtr[n+1]) memcpy(&col_buf[k],&mtx->col[j],sizeof(int)*(mtx->rowPtr[n+1]-j));
	}
	else if(j==mtx->rowPtr[n+1]) // n is done
	{
		if(i!=mtx->rowPtr[m+1]) memcpy(&col_buf[k],&mtx->col[i],sizeof(int)*(mtx->rowPtr[m+1]-i));
	}
	t2 = clock();
//	fprintf(stderr,"\tget flatten col: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
}





static void eliminateRow(double *flatten_val,const int i, const int j,const CSR_SparseDoubleMatrix *u,const double scale)
{
	clock_t t1,t2;
	t1 = clock();
	int ind_i = u->rowPtr[i] + 1;
	int ind_j = u->rowPtr[j] + 1;
	int k = 0;
	while(ind_i!=u->rowPtr[i+1] && ind_j!=u->rowPtr[j+1])
	{
		if(u->col[ind_i] > u->col[ind_j])
		{
			flatten_val[k] = u->val[ind_j];
			ind_j++;
		}
		else if(u->col[ind_i] < u->col[ind_j])
		{
			flatten_val[k] = -1.0 * scale * u->val[ind_i];
			ind_i++;
		}
		else
		{
			flatten_val[k] = u->val[ind_j] - scale*u->val[ind_i];
			ind_i++;
			ind_j++;
		}
		k++;
	}
	// append the remaining
	if(ind_i==u->rowPtr[i+1]) // row i is done
	{
		if(ind_j!=u->rowPtr[j+1]) memcpy(&flatten_val[k],&u->val[ind_j],sizeof(double)*(u->rowPtr[j+1]-ind_j));
	}
	else if(ind_j==u->rowPtr[j+1]) // row j is done
	{
		if(ind_i!=u->rowPtr[i+1]) memcpy(&flatten_val[k],&u->val[ind_i],sizeof(double)*(u->rowPtr[i+1]-ind_i));
		while(ind_i!=u->rowPtr[i+1])
		{
			flatten_val[k] *= -scale;
			k++;
			ind_i++;
		}
	}
	t2 = clock();
//	fprintf(stderr,"\teliminate: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
}	




static int getNumToEliminate(const CSR_SparseDoubleMatrix *u,const int i)
{
	int num = 0;
	int j;
	for(j=i+1;j<u->totalRow;j++) 
		if( u->col[u->rowPtr[j]] == i ) num++; // eliminate if necessary
	
	return num;
}




static void getEliminateRows(const CSR_SparseDoubleMatrix *u,const int i,const int numToEliminate, int **ret_row_buf, double **ret_val_buf)
{
	clock_t t1,t2;
	t1 = clock();
	const int bufSize = numToEliminate + 1;
	int *row_buf = getPidMempoolSet(sizeof(int)*bufSize,u->pid);
	double *val_buf = getPidMempoolSet(sizeof(double)*bufSize,u->pid);
	memset(val_buf,0,sizeof(double)*bufSize);
	row_buf[0] = i;
	val_buf[0] = 1.0;

	int j;
	int k = 1;
	for(j=i+1;j<u->totalRow;j++)
	{
		if( u->col[u->rowPtr[j]] == i )  // eliminate if necessary
		{
			row_buf[k] = j;
			k++;
		}
	}
	assert(k==bufSize);
	
	*ret_row_buf = row_buf;
	*ret_val_buf = val_buf;
	t2 = clock();
//	fprintf(stderr,"\tget eliminateRows: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
}



static int getRowILength(const CSR_SparseDoubleMatrix *mtx, const int i)
{
	return (mtx->rowPtr[i+1] - mtx->rowPtr[i]);
}




static int setDummy(int *flattenSizeList,int list_size,const int i,const CSR_SparseDoubleMatrix *u)
{
	clock_t t1,t2;
	t1 = clock();
	int j;
	memset(flattenSizeList,0,sizeof(int)*list_size);
	int total_u_dummy = 0;
	for(j=i+1;j<u->totalRow;j++)
	{
		if( u->col[u->rowPtr[j]] == i )  // eliminate if necessary
		{
			flattenSizeList[j] = getFlattenSize(u,i,j);
			total_u_dummy += (-(flattenSizeList[j] - getRowILength(u,i) - getRowILength(u,j)) + 1);
		}
	}
	assert(total_u_dummy >= 0);

	t2 = clock();
//	fprintf(stderr,"\tset dummy: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
	return total_u_dummy;
}







// only col and val
static int copy_first_i(CSR_SparseDoubleMatrix *u_after,const CSR_SparseDoubleMatrix *u_current,const int i)
{
/*
	clock_t t1,t2;
	t1 = clock();
	const int len = u_current->rowPtr[i+1];
	memcpy(u_after->col,u_current->col,sizeof(int)*len);
	memcpy(u_after->val,u_current->val,sizeof(double)*len);
	t2 = clock();
//	fprintf(stderr,"\tcopy: %2.10f\n",(double)(t2-t1)/CLOCKS_PER_SEC);
	return len;
*/
	const int len = u_current->rowPtr[i+1] - u_current->rowPtr[i];
	memcpy(&u_after->col[u_current->rowPtr[i]],&u_current->col[u_current->rowPtr[i]],sizeof(int)*len);
	memcpy(&u_after->val[u_current->rowPtr[i]],&u_current->val[u_current->rowPtr[i]],sizeof(double)*len);

	return u_current->rowPtr[i+1];
}




void lu_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *l, CSR_SparseDoubleMatrix *u, const CSR_SparseDoubleMatrix *a)
{
	assert(a->totalRow == u->totalRow);
	assert(a->totalCol == u->totalCol);
	assert(a->totalRow == l->totalRow);
	assert(a->totalRow == l->totalCol);
	int i,j;
	CSR_SparseDoubleMatrix *u_temp[2];
	u_temp[0] = create_CSR_SparseDoubleMatrix(u->totalRow,u->totalCol,u->nnz,u->pid);
	u_temp[1] = create_CSR_SparseDoubleMatrix(u->totalRow,u->totalCol,u->nnz,u->pid);
	copy_CSR_SparseDoubleMatrix(u_temp[0],a);
	copy_CSR_SparseDoubleMatrix(u_temp[1],a);

	CSR_SparseDoubleMatrix *u_current = u_temp[1];
	CSR_SparseDoubleMatrix *u_after = u_temp[0];
	expand_if_necessary(u_current,8*a->nnz);
	expand_if_necessary(u_after,8*a->nnz);
	int *flattenSizeList = getPidMempoolSet(sizeof(int)*a->totalRow,a->pid);
	for(i=0;i<a->totalRow;i++) // for each row i in a
	{
//		fprintf(stderr,"i=%d\n",i);
		u_current = u_temp[(i+1)%2];
		u_after = u_temp[i%2];
		const int numToEliminate = getNumToEliminate(u_current,i);
		int *eliminate_row;
		double *eliminate_val; // size = numToEliminate + 1
		getEliminateRows(u_current,i,numToEliminate,&eliminate_row,&eliminate_val);
		// set dummy
		const int total_u_dummy = setDummy(flattenSizeList,a->totalRow,i,u_current);
		expand_if_necessary(u_after,total_u_dummy);
		int ins_u = copy_first_i(u_after,u_current,i);
		int ins_l = 1;
		memcpy(u_after->rowPtr,u_current->rowPtr,sizeof(int)*(i+2));
		for(j=i+1;j<a->totalRow;j++)
		{
//			printf("%d %d\n",i,j);
			if( u_current->col[u_current->rowPtr[j]] == i ) // eliminate if necessary
			{
				assert(u_current->val[u_current->rowPtr[i]]!=0);
				const double scale = u_current->val[u_current->rowPtr[j]] / u_current->val[u_current->rowPtr[i]];
				const int flattenSize = flattenSizeList[j];
				getFlattenCols(u_current,i,j,&u_after->col[ins_u]);
				eliminateRow(&u_after->val[ins_u],i,j,u_current,scale);
				ins_u += (flattenSize-1);
				eliminate_val[ins_l] = scale;
				ins_l++;
				u_after->rowPtr[j+1] = u_after->rowPtr[j] + (flattenSize-1);
			}
			else
			{
				// not eliminate ... copy row j to u ...
				memcpy(&u_after->col[ins_u],&u_current->col[u_current->rowPtr[j]],sizeof(int)*(u_current->rowPtr[j+1]-u_current->rowPtr[j]));
				memcpy(&u_after->val[ins_u],&u_current->val[u_current->rowPtr[j]],sizeof(double)*(u_current->rowPtr[j+1]-u_current->rowPtr[j]));
				ins_u += (u_current->rowPtr[j+1]-u_current->rowPtr[j]);
				u_after->rowPtr[j+1] = u_after->rowPtr[j] + (u_current->rowPtr[j+1]-u_current->rowPtr[j]);
			}
		}
		assert(ins_l==numToEliminate+1);
		updateL(l,i,eliminate_row,eliminate_val,numToEliminate+1);
		
		retPidMempoolSet(eliminate_row,sizeof(int)*(numToEliminate+1),u->pid);
		retPidMempoolSet(eliminate_val,sizeof(double)*(numToEliminate+1),u->pid);
	}
	retPidMempoolSet(flattenSizeList,sizeof(int)*a->totalRow,a->pid);

	// compress l,u
	free_CSR_SparseDoubleMatrix(u_current);
	copy_CSR_SparseDoubleMatrix(u,u_after);
	free_CSR_SparseDoubleMatrix(u_after);
	
	expand_CSR_SparseDoubleMatrix(u,u->rowPtr[u->totalRow] - u->nnz);
	expand_CSR_SparseDoubleMatrix(l,l->rowPtr[l->totalRow] - l->nnz);

//	dump_CSR_SparseDoubleMatrix(stdout,l);
//	printf("\n");
//	dump_CSR_SparseDoubleMatrix(stdout,u);
//	printf("\n");
//	exit(0);
}










CSR_SparseDoubleMatrix *sparse2CSR(const SparseDoubleMatrix *mtx, const int pid)
{
	CSR_SparseDoubleMatrix *ret = create_CSR_SparseDoubleMatrix(mtx->totalRow,mtx->totalCol,mtx->nnz,pid);
	long long ind = 0;
	int i;
	ret->rowPtr[0] = 0;
	for(i=0;i<mtx->totalRow;i++)
	{
		const SparseDoubleElement *ptr = mtx->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			ret->col[ind] = ptr->col;
			ret->val[ind] = ptr->data;
			ind++;
			ptr = ptr->rowLink;
		}
		ret->rowPtr[i+1] = ind;
	}
	return ret;
}



// just set to zero
void clear_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx)
{
	memset(mtx->rowPtr,0,sizeof(int)*(mtx->totalRow+1));
	memset(mtx->col,0,sizeof(int)*(mtx->nnz));
	memset(mtx->val,0,sizeof(double)*(mtx->nnz));
}




/*
CSR_SparseDoubleMatrix *linus_read_to_CSR_SparseDoubleMatrix(const char *filename,const int pid)
{
	int i;
	int totalRow;
	int totalCol;
	int ret;
	long long nnz;
	struct stat s;
	int fd = open(filename,O_RDONLY);
	i = stat(filename,&s);
	int len = s.st_size;
	// head
	ret = posix_fadvise(fd,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_WILLNEED);
	void *buf = mmap(0,len,PROT_READ,MAP_SHARED,fd,0);

	memcpy(&totalRow,buf,sizeof(int));
	memcpy(&totalCol,buf+sizeof(int),sizeof(int));
	memcpy(&nnz,buf+2*sizeof(int),sizeof(long long));

	CSR_SparseDoubleMatrix *mtx = getPidMempoolSet(sizeof(CSR_SparseDoubleMatrix),pid);
	mtx->pid = pid;
	mtx->totalRow = totalRow;
	mtx->totalCol = totalCol;
	mtx->nnz = nnz;
	mtx->buf = buf;

	const int headSize = sizeof(int) + sizeof(int) + sizeof(long long); 
	mtx->rowPtr = buf + headSize;
	mtx->col = buf + headSize + sizeof(int)*(totalRow+1);
	mtx->val = buf + headSize + sizeof(int)*(totalRow+1) + sizeof(int)*(nnz);

	ret = madvise(buf,len,MADV_WILLNEED);

	close(fd);
	return mtx;
}
*/



static void *linus_read_file(void *par)
{
	int ret;
	Linus_read *ptr = (Linus_read *)par;
	int fd = open(ptr->name,O_RDONLY);
	ret = posix_fadvise(fd,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_WILLNEED);
	ptr->buf = mmap(0,ptr->size,PROT_READ,MAP_SHARED,fd,0);
	ret = madvise(ptr->buf,ptr->size,MADV_WILLNEED|MADV_SEQUENTIAL);
//	ret = madvise(mtx->rowPtr,sizeof(int)*(totalRow+1),MADV_WILLNEED|MADV_SEQUENTIAL);
	close(fd);
	pthread_exit(0);
}



static Linus_read* initLinusRead(char *name,int size)
{
	Linus_read *ptr = malloc(sizeof(Linus_read));
	ptr->name = name;
	ptr->size = size;
	return ptr;
}




CSR_SparseDoubleMatrix *linus_read_to_CSR_SparseDoubleMatrix(const char *prefix,const int pid)
{
	char sizeName[32] = {0};
	char rowName[32] = {0};
	char colName[32] = {0};
	char valName[32] = {0};
	sprintf(sizeName,"%s.size",prefix);
	sprintf(rowName,"%s.row",prefix);
	sprintf(colName,"%s.col",prefix);
	sprintf(valName,"%s.val",prefix);

	int i,ret;
	int totalRow;
	int totalCol;
	long long nnz;

	FILE *fp_size = fopen(sizeName,"rb");
	ret = fread(&totalRow,sizeof(int),1,fp_size);	
	ret = fread(&totalCol,sizeof(int),1,fp_size);	
	ret = fread(&nnz,sizeof(long long),1,fp_size);
	fclose(fp_size);

	CSR_SparseDoubleMatrix *mtx = getPidMempoolSet(sizeof(CSR_SparseDoubleMatrix),pid);
	mtx->pid = pid;
	mtx->totalRow = totalRow;
	mtx->totalCol = totalCol;
	mtx->nnz = nnz;
/*
	int fd_row = open(rowName,O_RDONLY);
	ret = posix_fadvise(fd_row,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_WILLNEED);
	mtx->rowPtr = mmap(0,sizeof(int)*(totalRow+1),PROT_READ,MAP_SHARED,fd_row,0);
//	ret = madvise(mtx->rowPtr,sizeof(int)*(totalRow+1),MADV_WILLNEED|MADV_SEQUENTIAL);
	ret = madvise(mtx->rowPtr,sizeof(int)*(totalRow+1),MADV_WILLNEED);
	close(fd_row);

	int fd_col = open(colName,O_RDONLY);
	ret = posix_fadvise(fd_col,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_WILLNEED);
	mtx->col = mmap(0,sizeof(int)*(nnz),PROT_READ,MAP_SHARED,fd_col,0);
//	ret = madvise(mtx->col,sizeof(int)*(nnz),MADV_WILLNEED|MADV_SEQUENTIAL);
	ret = madvise(mtx->col,sizeof(int)*(nnz),MADV_WILLNEED);
	close(fd_col);

	int fd_val = open(valName,O_RDONLY);
	ret = posix_fadvise(fd_val,0,0,POSIX_FADV_SEQUENTIAL|POSIX_FADV_WILLNEED);
	mtx->val = mmap(0,sizeof(double)*(nnz),PROT_READ,MAP_SHARED,fd_val,0);
//	ret = madvise(mtx->val,sizeof(double)*(nnz),MADV_WILLNEED|MADV_SEQUENTIAL);
	ret = madvise(mtx->val,sizeof(double)*(nnz),MADV_WILLNEED);
	close(fd_val);
*/
	pthread_t thread[3];
	Linus_read *ptr[3];
	ptr[0] = initLinusRead(rowName,sizeof(int)*(totalRow+1));
	ptr[1] = initLinusRead(colName,sizeof(int)*(nnz));
	ptr[2] = initLinusRead(valName,sizeof(double)*(nnz));
	for(i=0;i<3;i++) pthread_create(&thread[i],NULL,linus_read_file,ptr[i]);
	for(i=0;i<3;i++) pthread_join(thread[i],NULL);
	mtx->rowPtr = ptr[0]->buf;
	mtx->col = ptr[1]->buf;
	mtx->val = ptr[2]->buf;
	for(i=0;i<3;i++) free(ptr[i]);

	return mtx;
}





void linus_free_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx)
{
	int ret;
	ret = munmap(mtx->rowPtr,sizeof(int)*(mtx->totalRow+1));
	ret = munmap(mtx->col,sizeof(int)*(mtx->nnz));
	ret = munmap(mtx->val,sizeof(double)*(mtx->nnz));

	retPidMempoolSet(mtx,sizeof(CSR_SparseDoubleMatrix),mtx->pid);
}





void mulVec_CSR_SparseDoubleMatrix(double *res,const CSR_SparseDoubleMatrix *mtx,const double *vec) // res and vec can not be the same
{
	assert(res!=vec);
	memset(res,0,sizeof(double)*mtx->totalCol);
	int i,j;
	for(i=0;i<mtx->totalRow;i++)
	{
		double sum = 0;
		for(j=mtx->rowPtr[i];j<mtx->rowPtr[i+1];j++)
		{
			sum += mtx->val[j] * vec[mtx->col[j]];
		}
		res[i] = sum;
	}
}




static void *mulVecBlock_CSR_SparseDoubleMatrix(void *par) // res and vec can not be the same
{
	Parallel_mulVec *ptr = (Parallel_mulVec *)par;
	const int low = block_low(ptr->p,ptr->threadNum,ptr->mtx->totalRow);
	const int high = block_high(ptr->p,ptr->threadNum,ptr->mtx->totalRow);
	int i,j;

	for(i=low;i<=high;i++)
	{
		double sum = 0;
		for(j=ptr->mtx->rowPtr[i];j<ptr->mtx->rowPtr[i+1];j++)
		{
			sum += ptr->mtx->val[j] * ptr->vec[ptr->mtx->col[j]];
		}
		ptr->res[i] = sum;
	}


	pthread_exit(0);
}




// res and vec can not be the same
void parallel_mulVec_CSR_SparseDoubleMatrix(double *res,const CSR_SparseDoubleMatrix *mtx,const double *vec,const int threadNum)
{
	assert(res!=vec);
	if(threadNum == 1) mulVec_CSR_SparseDoubleMatrix(res,mtx,vec);
	else
	{
		int i,j;
		int *p = getMempoolSet(sizeof(int)*threadNum);
		for(i=0;i<threadNum;i++) p[i] = i;
		memset(res,0,sizeof(double)*mtx->totalCol);
		Parallel_mulVec *parList = getMempoolSet(sizeof(Parallel_mulVec)*threadNum);
		
		pthread_t *pid = getMempoolSet(threadNum*sizeof(pthread_t));
		for(i=0;i<threadNum;i++)
		{
			parList[i].p = i;
			parList[i].threadNum = threadNum;
			parList[i].res = res;
			parList[i].mtx = mtx;
			parList[i].vec = vec;
			pthread_create(&pid[i],NULL,mulVecBlock_CSR_SparseDoubleMatrix,&parList[i]);
		}
		for(i=0;i<threadNum;i++) pthread_join(pid[i],NULL);

		retMempoolSet(pid,threadNum*sizeof(pthread_t));
		retMempoolSet(parList,sizeof(Parallel_mulVec)*threadNum);
		retMempoolSet(p,sizeof(int)*threadNum);
	}
}




//==============================================



CSC_SparseDoubleMatrix *create_CSC_SparseDoubleMatrix(int row, int col, long long nnz,int pid)
{
	CSC_SparseDoubleMatrix *ptr = getPidMempoolSet(sizeof(CSC_SparseDoubleMatrix),pid);
	ptr->pid = pid;
	ptr->totalRow = row;
	ptr->totalCol = col;
	ptr->nnz = nnz;

	ptr->row = getPidMempoolSet(sizeof(int)*nnz,pid);
	ptr->colPtr = getPidMempoolSet(sizeof(int)*(col+1),pid);
	ptr->val = getPidMempoolSet(sizeof(double)*nnz,pid);

	memset(ptr->row,0,sizeof(int)*nnz);
	memset(ptr->colPtr,0,sizeof(int)*(col+1));
	memset(ptr->val,0,sizeof(double)*nnz);

	return ptr;
}





void free_CSC_SparseDoubleMatrix(CSC_SparseDoubleMatrix *ptr)
{
	retPidMempoolSet(ptr->row,sizeof(int)*ptr->nnz,ptr->pid);
	retPidMempoolSet(ptr->colPtr,sizeof(int)*(ptr->totalCol+1),ptr->pid);
	retPidMempoolSet(ptr->val,sizeof(double)*ptr->nnz,ptr->pid);

	retPidMempoolSet(ptr,sizeof(CSC_SparseDoubleMatrix),ptr->pid);
}






CSC_SparseDoubleMatrix *sparse2CSC(const SparseDoubleMatrix *mtx, const int pid)
{
	CSC_SparseDoubleMatrix *ptr = create_CSC_SparseDoubleMatrix(mtx->totalRow,mtx->totalCol,mtx->nnz,pid);

	int i;
	int j = 0;
	for(i=0;i<mtx->totalCol;i++)
	{
		const SparseDoubleElement *eachCol = mtx->colIndex[i]->colLink;
		int len = 0;
		while(eachCol!=NULL)
		{
			ptr->row[j] = eachCol->row;
			ptr->val[j] = eachCol->data;
			j++;
			len++;
			eachCol = eachCol->colLink;
		}
		ptr->colPtr[i+1] = ptr->colPtr[i] + len;
	}

	return ptr;
}






