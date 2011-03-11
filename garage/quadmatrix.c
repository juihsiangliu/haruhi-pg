#include "quadmatrix.h"
#include <stdlib.h>


static int indexConvert(const QuadMatrix *ptr,const int rowIndex,const int colIndex)
{
	return (rowIndex*ptr->col + colIndex);
}


// swap a->data[row1][:] and a->data[row2][:] 
// direct swap the pointer
// (use swapQuadMatrix)
void static swapRowQuadMatrix(QuadMatrix *a,const int row1,const int row2)
{
	const int total = a->col;
	int i;
	for(i=0;i<total;i++) swapQuadMatrix(a,row1,i,row2,i);
}


// find out the pivot from currentRow and update the pivot "p" and the goal matrix "a"
static void updateQuadMatrix(QuadMatrix *a, QuadMatrix *p, const int currentRow)
{
	int i;
	int index = currentRow;
	// find out the pivot row
	for(i=currentRow;a->row;i++)
	{
		if(!isEmptyQuadElement(getPtrEntryQuadMatrix(a,i,currentRow)));
		{
			index = i;
			break;
		}
	}
	// update the permutation matrix
	swapQuadMatrix(p,currentRow,currentRow,index,currentRow);
	swapQuadMatrix(p,currentRow,index,index,index);
	// swap currentRow and the pivot row
	swapRowQuadMatrix(a,currentRow,index);	
}



// setup the quadmatrix to an identity matrix
static void identityQuadMatrix(QuadMatrix *a)
{
	const int total = a->row;
	int i,j;
	// reset
	for(i=0;i<total;i++)
	{
		for(j=0;j<total;j++)
		{
			const int index = indexConvert(a,i,j);
			QuadElement *current = getPtrEntryQuadMatrix(a,i,j);
			freeQuadElement(current);
			a->data[index] = NULL;
		}
	}
	// set diagonal
	QuadElement *unit = createQuadElement(a->gvNum);
	unit->m = 1;
	for(i=0;i<total;i++)
	{
		setQuadMatrix(a,unit,i,i);	
	}
	freeQuadElement(unit);
}





// setup the quadmatrix to an identity matrix
static void identityPidQuadMatrix(QuadMatrix *a,const int pid)
{
	const int total = a->row;
	int i,j;
	// reset
	for(i=0;i<total;i++)
	{
		for(j=0;j<total;j++)
		{
			const int index = indexConvert(a,i,j);
			QuadElement *current = getPtrEntryQuadMatrix(a,i,j);
			freePidQuadElement(current,pid);
			a->data[index] = NULL;
		}
	}
	// set diagonal
	QuadElement *unit = createPidQuadElement(a->gvNum,pid);
	unit->m = 1;
	for(i=0;i<total;i++)
	{
		setPidQuadMatrix(a,unit,i,i,pid);	
	}
	freePidQuadElement(unit,pid);
}




//======================================================================


// c = a op b , op the nth row of c
struct QuadOpParallel
{
	const QuadMatrix *a;
	const QuadMatrix *b;
	QuadMatrix *c;
	int n;
};

typedef struct QuadOpParallel QuadOpParallel;



//======================================================================

QuadMatrix *createQuadMatrix(const int row,const int col,const int gvNum)
{
	int i;
	QuadMatrix *ptr = (QuadMatrix *)getMempoolSet(sizeof(QuadMatrix));
	ptr->row = row;
	ptr->col = col;
	ptr->gvNum = gvNum;
	ptr->data = (QuadElement **)getMempoolSet( row*col*sizeof(QuadElement *));
	for(i=0;i<row*col;i++) ptr->data[i] = NULL;

	return ptr;
}




QuadMatrix *createPidQuadMatrix(const int row,const int col,const int gvNum,const int pid)
{
	int i;
	QuadMatrix *ptr = (QuadMatrix *)getPidMempoolSet(sizeof(QuadMatrix),pid);
	ptr->row = row;
	ptr->col = col;
	ptr->gvNum = gvNum;
	ptr->data = (QuadElement **)getPidMempoolSet( row*col*sizeof(QuadElement *),pid);
	for(i=0;i<row*col;i++) ptr->data[i] = NULL;

	return ptr;
}





void freeQuadMatrix(QuadMatrix *ptr)
{
	int i;
	const int total = ptr->row * ptr->col;
	for(i=0;i<total;i++)
	{
		if( ptr->data[i] != NULL)
			freeQuadElement(ptr->data[i]);
	}
	retMempoolSet(ptr->data,total*sizeof(QuadElement *));
	retMempoolSet(ptr,sizeof(QuadMatrix));
}





void freePidQuadMatrix(QuadMatrix *ptr,const int pid)
{
	int i;
	const int total = ptr->row * ptr->col;
	for(i=0;i<total;i++)
	{
		if( ptr->data[i] != NULL)
			freePidQuadElement(ptr->data[i],pid);
	}
	retPidMempoolSet(ptr->data,total*sizeof(QuadElement *),pid);
	retPidMempoolSet(ptr,sizeof(QuadMatrix),pid);
}




void setQuadMatrix(QuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex)
{
	const int gvNum = ptr->gvNum;
	const int index = indexConvert(ptr,rowIndex,colIndex);
	if(isEmptyQuadElement(element))
	{
		freeQuadElement(ptr->data[index]);
		ptr->data[index] = NULL;
	}
	else
	{
		QuadElement *dest = getPtrEntryQuadMatrix(ptr,rowIndex,colIndex);
		if( dest == NULL )
		{
			ptr->data[index] = createQuadElement(ptr->gvNum);
			copyQuadElement(ptr->data[index],element);
		}
		else
		{
			copyQuadElement(dest,element);
		}
	}
}





void setPidQuadMatrix(QuadMatrix *ptr,const QuadElement *element,const int rowIndex, const int colIndex,const int pid)
{
	const int gvNum = ptr->gvNum;
	const int index = indexConvert(ptr,rowIndex,colIndex);
	if(isEmptyQuadElement(element))
	{
		freePidQuadElement(ptr->data[index],pid);
		ptr->data[index] = NULL;
	}
	else
	{
		QuadElement *dest = getPtrEntryQuadMatrix(ptr,rowIndex,colIndex);
		if( dest == NULL )
		{
			ptr->data[index] = createPidQuadElement(ptr->gvNum,pid);
			copyQuadElement(ptr->data[index],element);
		}
		else
		{
			copyQuadElement(dest,element);
		}
	}
}






QuadElement *getPtrEntryQuadMatrix(const QuadMatrix *ptr,const int rowIndex,const int colIndex)
{
	const int index = indexConvert(ptr,rowIndex,colIndex);
	return ptr->data[index];
}




QuadElement *getCopyEntryQuadMatrix(QuadElement *copy,const QuadMatrix *ptr,const int rowIndex,const int colIndex)
{
	const int index = indexConvert(ptr,rowIndex,colIndex);
	copyQuadElement(copy,ptr->data[index]);
}



static void *addQuadMatrixRow(void *ptr)
{
	QuadOpParallel *par = (QuadOpParallel *) ptr;
	int i;
	const int col = par->a->col;
	QuadElement *entryC = createQuadElement(par->a->gvNum);
	for(i=0;i<col;i++)
	{
		const QuadElement *a = getPtrEntryQuadMatrix(par->a,par->n,i);
		const QuadElement *b = getPtrEntryQuadMatrix(par->b,par->n,i);
		addQuadElement(entryC,a,b);
		setQuadMatrix(par->c,entryC,par->n,i);
	}
	freeQuadElement(entryC);	
	pthread_exit(NULL);
}



// c = a + b
void addQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b)
{
/*
	int i,j;
	QuadElement *entryC = createQuadElement(a->gvNum);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *entryA = getPtrEntryQuadMatrix(a,i,j);
			const QuadElement *entryB = getPtrEntryQuadMatrix(b,i,j);
			addQuadElement(entryC,entryA,entryB);
			setQuadMatrix(c,entryC,i,j);
		}
	}	
	freeQuadElement(entryC);
*/
	addPidQuadMatrix(c,a,b,0);
}




void addPidQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b,const int pid)
{
	int i,j;
	QuadElement *entryC = createPidQuadElement(a->gvNum,pid);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *entryA = getPtrEntryQuadMatrix(a,i,j);
			const QuadElement *entryB = getPtrEntryQuadMatrix(b,i,j);
			addQuadElement(entryC,entryA,entryB);
			setPidQuadMatrix(c,entryC,i,j,pid);
		}
	}	
	freePidQuadElement(entryC,pid);
}






// c = a - b
void subQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b)
{
/*
	int i,j;
	QuadElement *entryC = createQuadElement(a->gvNum);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *entryA = getPtrEntryQuadMatrix(a,i,j);
			const QuadElement *entryB = getPtrEntryQuadMatrix(b,i,j);
			subQuadElement(entryC,entryA,entryB);
			setQuadMatrix(c,entryC,i,j);
		}
	}
	freeQuadElement(entryC);
*/
	subPidQuadMatrix(c,a,b,0);
}





// c = a - b
void subPidQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b,const int pid)
{
	int i,j;
	QuadElement *entryC = createPidQuadElement(a->gvNum,pid);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *entryA = getPtrEntryQuadMatrix(a,i,j);
			const QuadElement *entryB = getPtrEntryQuadMatrix(b,i,j);
			subQuadElement(entryC,entryA,entryB);
			setPidQuadMatrix(c,entryC,i,j,pid);
		}
	}
	freePidQuadElement(entryC,pid);
}





// c = a * b
void mulQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b)
{
	mulPidQuadMatrix(c,a,b,0);
}




void mulPidQuadMatrix(QuadMatrix *c, const QuadMatrix *a, const QuadMatrix *b,const int pid)
{
	int i,j,k;
	QuadElement *temp = createPidQuadElement(a->gvNum,pid);
	QuadElement *sum = createPidQuadElement(a->gvNum,pid);
	QuadMatrix *result;
	
	// to avoid the dummy copy and allocation
	if(c==a || c==b)
		result = createPidQuadMatrix(a->row,b->col,a->gvNum,pid);
	else
		result = c;

	for(i=0;i<a->row;i++)
	{
		for(j=0;j<b->col;j++)
		{
			for(k=0;k<a->col;k++)
			{
				const QuadElement *entryA = getPtrEntryQuadMatrix(a,i,k);
				const QuadElement *entryB = getPtrEntryQuadMatrix(b,k,j);
				mulPidQuadElement(temp,entryA,entryB,pid);
				addQuadElement(sum,sum,temp);
			}
			setPidQuadMatrix(result,sum,i,j,pid);
			resetQuadElement(sum);
		}
	}

	if(c==a || c==b)
	{
		copyPidQuadMatrix(c,result,pid);
		freePidQuadMatrix(result,pid);
	}

	freePidQuadElement(temp,pid);
	freePidQuadElement(sum,pid);
}






// c = k * a, // k is a scale
void scaleQuadMatrix(QuadMatrix *c, const double k, const QuadMatrix *a)
{
	int i,j;
	QuadElement *dst = createQuadElement(a->gvNum);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *src = getPtrEntryQuadMatrix(a,i,j);
			scaleQuadElement(dst,k,src);
			setQuadMatrix(c,dst,i,j);
		}
	}
	freeQuadElement(dst);
}






// c = k * a, // k is a scale
void scalePidQuadMatrix(QuadMatrix *c, const double k, const QuadMatrix *a,const int pid)
{
	int i,j;
	QuadElement *dst = createPidQuadElement(a->gvNum,pid);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *src = getPtrEntryQuadMatrix(a,i,j);
			scaleQuadElement(dst,k,src);
			setPidQuadMatrix(c,dst,i,j,pid);
		}
	}
	freePidQuadElement(dst,pid);
}





// c = q * a // q is a quadelement
void scaleqQuadMatrix(QuadMatrix *c, const QuadElement *q, const QuadMatrix *a)
{
	int i,j;
	QuadElement *dst = createQuadElement(a->gvNum);
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *src = getPtrEntryQuadMatrix(a,i,j);
			mulQuadElement(dst,q,src);
			setQuadMatrix(c,dst,i,j);
		}
	}
	freeQuadElement(dst);
}



// b = a'
void transposeQuadMatrix(QuadMatrix *b,const QuadMatrix *a)
{
	transposePidQuadMatrix(b,a,0);
}





void transposePidQuadMatrix(QuadMatrix *b,const QuadMatrix *a,const int pid)
{
	int i,j;
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *ptrInA = getPtrEntryQuadMatrix(a,i,j);
			setPidQuadMatrix(b,ptrInA,j,i,pid);
		}
	}
}






// swap ptr->data[row1][col1] and ptr->data[row2][col2]
void swapQuadMatrix(QuadMatrix *ptr,const int row1,const int col1,const int row2, const int col2)
{
	int index;
	// temp = &2
	QuadElement *temp = getPtrEntryQuadMatrix(ptr,row2,col2);
	// &2 = &1
	index = row2 * ptr->col + col2;
	ptr->data[index] = getPtrEntryQuadMatrix(ptr,row1,col1);
	// &1 = temp
	index = row1 * ptr->col + col1;
	ptr->data[index] = temp;
}




void copyQuadMatrix(QuadMatrix *dest, const QuadMatrix *src)
{
	copyPidQuadMatrix(dest,src,0);
}




void copyPidQuadMatrix(QuadMatrix *dest, const QuadMatrix *src,const int pid)
{
	const int row = src->row;
	const int col = src->col;
	int i,j;

	dest->row = src->row;
	dest->col = src->col;
	dest->gvNum = src->gvNum;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			const QuadElement *srcElement = getPtrEntryQuadMatrix(src,i,j);
			setPidQuadMatrix(dest,srcElement,i,j,pid);
		}
	}
}






// dump only 1 entry
void dumpEntryQuadMatrix(const QuadMatrix *ptr,const int rowIndex,const int colIndex)
{
	printf("row: %d , col: %d\n",rowIndex,colIndex);
	const QuadElement *entry = getPtrEntryQuadMatrix(ptr,rowIndex,colIndex);
	dumpQuadElement(entry);
}




// dump the whole matrix
void dumpQuadMatrix(const QuadMatrix *ptr)
{
	int i,j;
	for(i=0;i<ptr->row;i++)
	{
		for(j=0;j<ptr->col;j++)
		{
			dumpEntryQuadMatrix(ptr,i,j);
		}
	}
}



// a = plu
void luQuadMatrix(QuadMatrix *p,QuadMatrix *l,QuadMatrix *u,const QuadMatrix *a)
{
	QuadMatrix *aTemp = createQuadMatrix(a->row,a->col,a->gvNum);
	copyQuadMatrix(aTemp,a);
	QuadMatrix *E = createQuadMatrix(a->row,a->col,a->gvNum);
	identityQuadMatrix(l);
	identityQuadMatrix(p);
	int i,j,k;
	for(i=0;i<a->row;i++)
	{
		updateQuadMatrix(aTemp,p,i);
		for(j=i+1;j<a->col;j++)
		{
			if( isEmptyQuadElement( getPtrEntryQuadMatrix(aTemp,j,i)))
			{
				continue;
			}
			QuadElement *scale = createQuadElement(a->gvNum);
			QuadElement *target = getPtrEntryQuadMatrix(aTemp,j,i);
			divQuadElement(scale,target,getPtrEntryQuadMatrix(aTemp,i,i));
			identityQuadMatrix(E);
			setQuadMatrix(E,scale,j,i);
			mulQuadMatrix(l,l,E);
			freeQuadElement(target);
			const int index = indexConvert(aTemp,j,i);
			aTemp->data[index] = NULL;

			for(k=i+1;k<a->col;k++)
			{
				const QuadElement *rowIElement = getPtrEntryQuadMatrix(aTemp,i,k);
				QuadElement *rowJElement = getPtrEntryQuadMatrix(aTemp,j,k);
				QuadElement *rowTemp = createQuadElement(a->gvNum);
				QuadElement *resultTemp = createQuadElement(a->gvNum);
				mulQuadElement(rowTemp,scale,rowIElement);
				// use set function to perform fill-in operation
				subQuadElement(resultTemp,rowJElement,rowTemp);
				setQuadMatrix(aTemp,resultTemp,j,k);
				freeQuadElement(rowTemp);
				freeQuadElement(resultTemp);
			}
			freeQuadElement(scale);
		}
	}
	
	copyQuadMatrix(u,aTemp);
	
	freeQuadMatrix(aTemp);
	freeQuadMatrix(E);
}




// a = lu
void luPidQuadMatrix(QuadMatrix *l,QuadMatrix *u,const QuadMatrix *a,const int pid)
{
	QuadMatrix *aTemp = createPidQuadMatrix(a->row,a->col,a->gvNum,pid);
	copyPidQuadMatrix(aTemp,a,pid);
	QuadMatrix *E = createPidQuadMatrix(a->row,a->col,a->gvNum,pid);
	identityPidQuadMatrix(l,pid);
	int i,j,k;
	for(i=0;i<a->row;i++)
	{
		for(j=i+1;j<a->col;j++)
		{
			if( isEmptyQuadElement( getPtrEntryQuadMatrix(aTemp,j,i)))
			{
				continue;
			}
			QuadElement *scale = createPidQuadElement(a->gvNum,pid);
			QuadElement *target = getPtrEntryQuadMatrix(aTemp,j,i);
			divPidQuadElement(scale,target,getPtrEntryQuadMatrix(aTemp,i,i),pid);
			identityPidQuadMatrix(E,pid);
			setPidQuadMatrix(E,scale,j,i,pid);
			mulPidQuadMatrix(l,l,E,pid);
			freePidQuadElement(target,pid);
			const int index = indexConvert(aTemp,j,i);
			aTemp->data[index] = NULL;

			for(k=i+1;k<a->col;k++)
			{
				const QuadElement *rowIElement = getPtrEntryQuadMatrix(aTemp,i,k);
				QuadElement *rowJElement = getPtrEntryQuadMatrix(aTemp,j,k);
				QuadElement *rowTemp = createPidQuadElement(a->gvNum,pid);
				QuadElement *resultTemp = createPidQuadElement(a->gvNum,pid);
				mulPidQuadElement(rowTemp,scale,rowIElement,pid);
				// use set function to perform fill-in operation
				subQuadElement(resultTemp,rowJElement,rowTemp);
				setPidQuadMatrix(aTemp,resultTemp,j,k,pid);
				freePidQuadElement(rowTemp,pid);
				freePidQuadElement(resultTemp,pid);
			}
			freePidQuadElement(scale,pid);
		}
	}
	
	copyPidQuadMatrix(u,aTemp,pid);
	
	freePidQuadMatrix(aTemp,pid);
	freePidQuadMatrix(E,pid);
}






//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveQuadMatrix(QuadMatrix *x,const QuadMatrix *p,const QuadMatrix *l,const QuadMatrix *u,const QuadMatrix *b)
{
	int i,k;
	QuadMatrix *y = createQuadMatrix(x->row,1,x->gvNum);
	QuadMatrix *pb = createQuadMatrix(x->row,1,x->gvNum);
/*
	// to fill out x,y,pb
	QuadElement *e = createQuadElement(x->gvNum);
	setQuadElement(e,1,0,NULL,NULL);
	for(i=0;i<x->row;i++)
	{
		setQuadMatrix(x,e,i,0);
		setQuadMatrix(y,e,i,0);
		setQuadMatrix(pb,e,i,0);
	}
*/
	mulQuadMatrix(pb,p,b);
	// solve ly = pb
	for(i=0;i<x->row;i++)
	{
		if(!isEmptyQuadElement( getPtrEntryQuadMatrix(l,i,i) ))
		{
			QuadElement *temp1 = createQuadElement(b->gvNum);
			for(k=0;k<i;k++)
			{
				QuadElement *temp2 = createQuadElement(b->gvNum);
				mulQuadElement(temp2,getPtrEntryQuadMatrix(l,i,k),getPtrEntryQuadMatrix(y,k,0));
				addQuadElement(temp1,temp1,temp2);
				freeQuadElement(temp2);
			}
			// check the null dest
			const int indexPB = indexConvert(pb,i,0);
			if(pb->data[indexPB] == NULL)
				pb->data[indexPB] = createQuadElement(pb->gvNum);
			subQuadElement(getPtrEntryQuadMatrix(pb,i,0),getPtrEntryQuadMatrix(pb,i,0),temp1);
			// check the null dest
			const int indexY = indexConvert(y,i,0);
			if(y->data[indexY] == NULL)
				y->data[indexY] = createQuadElement(y->gvNum);
			divQuadElement(getPtrEntryQuadMatrix(y,i,0),getPtrEntryQuadMatrix(pb,i,0),getPtrEntryQuadMatrix(l,i,i));
			freeQuadElement(temp1);
		}
	}


	// solve ux = y
	for(i=x->row-1;i>-1;i--)
	{
		if(!isEmptyQuadElement( getPtrEntryQuadMatrix(u,i,i)  ))
		{
			QuadElement *temp1 = createQuadElement(b->gvNum);
			for(k=i+1;k<x->row;k++)
			{
				QuadElement *temp2 = createQuadElement(b->gvNum);	
				mulQuadElement(temp2,getPtrEntryQuadMatrix(u,i,k),getPtrEntryQuadMatrix(x,k,0));
				addQuadElement(temp1,temp1,temp2);
				freeQuadElement(temp2);
			}
			subQuadElement(getPtrEntryQuadMatrix(y,i,0),getPtrEntryQuadMatrix(y,i,0),temp1);
			// check the null dest
			const int indexX = indexConvert(x,i,0);
			if(x->data[indexX] == NULL)
				x->data[indexX] = createQuadElement(x->gvNum);
			divQuadElement(getPtrEntryQuadMatrix(x,i,0),getPtrEntryQuadMatrix(y,i,0),getPtrEntryQuadMatrix(u,i,i));
			freeQuadElement(temp1);
		 }
	}

	freeQuadMatrix(y);
	freeQuadMatrix(pb);
}





//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolvePidQuadMatrix(QuadMatrix *x,const QuadMatrix *l,const QuadMatrix *u,const QuadMatrix *b,const int pid)
{
	int i,k;
	QuadMatrix *y = createPidQuadMatrix(x->row,1,x->gvNum,pid);
	const QuadMatrix *pb = b;
	
	// solve ly = pb
	for(i=0;i<x->row;i++)
	{
		if(!isEmptyQuadElement( getPtrEntryQuadMatrix(l,i,i) ))
		{
			QuadElement *temp1 = createPidQuadElement(b->gvNum,pid);
			for(k=0;k<i;k++)
			{
				QuadElement *temp2 = createPidQuadElement(b->gvNum,pid);
				mulPidQuadElement(temp2,getPtrEntryQuadMatrix(l,i,k),getPtrEntryQuadMatrix(y,k,0),pid);
				addQuadElement(temp1,temp1,temp2);
				freePidQuadElement(temp2,pid);
			}
			// check the null dest
			const int indexPB = indexConvert(pb,i,0);
			if(pb->data[indexPB] == NULL)
				pb->data[indexPB] = createPidQuadElement(pb->gvNum,pid);
			subQuadElement(getPtrEntryQuadMatrix(pb,i,0),getPtrEntryQuadMatrix(pb,i,0),temp1);
			// check the null dest
			const int indexY = indexConvert(y,i,0);
			if(y->data[indexY] == NULL)
				y->data[indexY] = createPidQuadElement(y->gvNum,pid);
			divPidQuadElement(getPtrEntryQuadMatrix(y,i,0),getPtrEntryQuadMatrix(pb,i,0),getPtrEntryQuadMatrix(l,i,i),pid);
			freePidQuadElement(temp1,pid);
		}
	}


	// solve ux = y
	for(i=x->row-1;i>-1;i--)
	{
		if(!isEmptyQuadElement( getPtrEntryQuadMatrix(u,i,i)  ))
		{
			QuadElement *temp1 = createPidQuadElement(b->gvNum,pid);
			for(k=i+1;k<x->row;k++)
			{
				QuadElement *temp2 = createPidQuadElement(b->gvNum,pid);	
				mulPidQuadElement(temp2,getPtrEntryQuadMatrix(u,i,k),getPtrEntryQuadMatrix(x,k,0),pid);
				addQuadElement(temp1,temp1,temp2);
				freePidQuadElement(temp2,pid);
			}
			subQuadElement(getPtrEntryQuadMatrix(y,i,0),getPtrEntryQuadMatrix(y,i,0),temp1);
			// check the null dest
			const int indexX = indexConvert(x,i,0);
			if(x->data[indexX] == NULL)
				x->data[indexX] = createPidQuadElement(x->gvNum,pid);
			divPidQuadElement(getPtrEntryQuadMatrix(x,i,0),getPtrEntryQuadMatrix(y,i,0),getPtrEntryQuadMatrix(u,i,i),pid);
			freePidQuadElement(temp1,pid);
		 }
	}

	freePidQuadMatrix(y,pid);
}






// directly use luQuadMatrix() and triSolveQuadMatrix() to solve Ax = b
void solveQuadMatrix(QuadMatrix *x,const QuadMatrix *A,const QuadMatrix *b)
{
	const int nodeNum = A->row;
	const int gvNum = A->gvNum;
	QuadMatrix *P = createQuadMatrix(nodeNum,nodeNum,gvNum);
	QuadMatrix *L = createQuadMatrix(nodeNum,nodeNum,gvNum);
	QuadMatrix *U = createQuadMatrix(nodeNum,nodeNum,gvNum);
	luQuadMatrix(P,L,U,A);
	triSolveQuadMatrix(x,P,L,U,b);
	freeQuadMatrix(P);
	freeQuadMatrix(L);
	freeQuadMatrix(U);
}




// no permutation inside ...
void solvePidQuadMatrix(QuadMatrix *x,const QuadMatrix *A, const QuadMatrix *b,const int pid)
{
	const int nodeNum = A->row;
	const int gvNum = A->gvNum;
	QuadMatrix *L = createPidQuadMatrix(nodeNum,nodeNum,gvNum,pid);
	QuadMatrix *U = createPidQuadMatrix(nodeNum,nodeNum,gvNum,pid);
	luPidQuadMatrix(L,U,A,pid);
	triSolvePidQuadMatrix(x,L,U,b,pid);


	freePidQuadMatrix(L,pid);
	freePidQuadMatrix(U,pid);

}



// return the copy of the nth column in a
void getColCopyQuadMatrix(QuadMatrix *col,const int n,const QuadMatrix *a)
{
	int i = 0;
	for(i=0;i<a->row;i++)
	{
		const QuadElement *ptr = getPtrEntryQuadMatrix(a,i,n);
		setQuadMatrix(col,ptr,i,0);
	}
}


// return the copy of the nth row in a
void getRowCopyQuadMatrix(QuadMatrix *row,const int n,const QuadMatrix *a)
{
	int i = 0;
	for(i=0;i<a->col;i++)
	{
		const QuadElement *ptr = getPtrEntryQuadMatrix(a,n,i);
		setQuadMatrix(row,ptr,0,i);
	}
}



//set the nth col in a as col
void setColQuadMatrix(QuadMatrix *a,const QuadMatrix *col,const int n)
{
	int i = 0;
	for(i=0;i<a->row;i++)
	{
		const QuadElement *ptr = getPtrEntryQuadMatrix(col,i,0);
		setQuadMatrix(a,ptr,i,n);
	}
}


//set the nth row in a as row
void setRowQuadMatrix(QuadMatrix *a,const QuadMatrix *row,const int n)
{
	int i = 0;
	for(i=0;i<a->col;i++)
	{
		const QuadElement *ptr = getPtrEntryQuadMatrix(row,0,i);
		setQuadMatrix(a,ptr,n,i);
	}
}



void meanQuadMatrix(gsl_matrix *meanMatrix,const QuadMatrix *a,const double *r)
{
	int i,j;
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(a,i,j);
			if(!isEmptyQuadElement(element))
			{
				const double val = meanQuadElement(element,r);
				gsl_matrix_set(meanMatrix,i,j,val);
			}
			else gsl_matrix_set(meanMatrix,i,j,0.0);
		}
	}
}


void varQuadMatrix(gsl_matrix *varMatrix,const QuadMatrix *a,const double *r)
{
	int i,j;
	for(i=0;i<a->row;i++)
	{
		for(j=0;j<a->col;j++)
		{
			const QuadElement *element = getPtrEntryQuadMatrix(a,i,j); 
			if(!isEmptyQuadElement(element))
			{
				const double val = varQuadElement(element,r);
				gsl_matrix_set(varMatrix,i,j,val);
			}
			else gsl_matrix_set(varMatrix,i,j,0.0);
		}
	}
}




void resetQuadMatrix(QuadMatrix *ptr)
{
	const int row = ptr->row;
	const int col = ptr->col;
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			resetQuadElement(getPtrEntryQuadMatrix(ptr,i,j));
		}
	}
}




// ptr need be allocated before
// all entry will be set to zero
// the null entry will be also filled
void setZeroQuadMatrix(QuadMatrix *ptr)
{
	const int row = ptr->row;
	const int col = ptr->col;
	const int gvNum = ptr->gvNum;
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			QuadElement * result = setZeroQuadElement(getPtrEntryQuadMatrix(ptr,i,j),gvNum);
			QuadElement *dest = getPtrEntryQuadMatrix(ptr,i,j);
			if( dest == NULL )
			{
				const int index = indexConvert(ptr,i,j);
				ptr->data[index] = createQuadElement(ptr->gvNum);
				copyQuadElement(dest,result);
			}
			else
			{
				copyQuadElement(dest,result);
			}
			freeQuadElement(result);
		}
	}
}



// c = a'*b,  a and b are n x 1 vector
void innerQuadMatrix(QuadElement *c,const QuadMatrix *a, const QuadMatrix *b)
{
	int i;
	const int row = a->row;
	QuadElement *temp = createQuadElement(a->gvNum);
	resetQuadElement(c);
	for(i=0;i<row;i++)
	{
		const QuadElement *op1 = getPtrEntryQuadMatrix(a,i,0);
		const QuadElement *op2 = getPtrEntryQuadMatrix(b,i,0);
		mulQuadElement(temp,op1,op2);
		addQuadElement(c,c,temp);
	}
	freeQuadElement(temp);
}



void delPidQuadMatrix(QuadMatrix *ptr,const int rowIndex, const int colIndex,const int pid)
{
	const int index = indexConvert(ptr,rowIndex,colIndex);
	if(ptr->data[index]!=NULL)
	{
		freePidQuadElement(ptr->data[index],pid);
		ptr->data[index] = NULL;
	}
}




void delQuadMatrix(QuadMatrix *ptr,const int rowIndex, const int colIndex)
{
	delPidQuadMatrix(ptr,rowIndex,colIndex,0);
}




