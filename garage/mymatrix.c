#include "mymatrix.h"
#include <stdio.h>



static int indexConvert(const int totalRow,const int totalCol, const int rowIndex,const int colIndex)
{
	return (rowIndex*totalCol + colIndex);
}



void setMyMatrix(double *a,double val,const int rowA,const int colA,const int i, const int j)
{
//	assert(i>=0);
//	assert(i<rowA);
//	assert(j>=0);
//	assert(j<colA);
	const int index = indexConvert(rowA,colA,i,j);
//	assert(index < rowA*colA);

	a[index] = val;
	return;
}


double getMyMatrix(const double *a,const int rowA,const int colA,const int i, const int j)
{
//	assert(i>=0);
//	assert(i<rowA);
//	assert(j>=0);
//	assert(j<colA);
	const int index = indexConvert(rowA,colA,i,j);
//	assert(index < rowA*colA);
	const double val = a[index];
	return val;
}



void addMyMatrix(double *c,const double *a,const double *b, const int rowA,const int colA)
{
	int i;
	const int totalNum = rowA*colA;
	for(i=0;i<totalNum;i++) c[i] = a[i] + b[i];
}



void subMyMatrix(double *c,const double *a,const double *b, const int rowA,const int colA)
{
	int i;
	const int totalNum = rowA*colA;
	for(i=0;i<totalNum;i++) c[i] = a[i] - b[i];
}


// the address of a and c  or  b and c can not be the same
void mulMyMatrix(double *c,const double *a,const double *b, const int rowA,const int colA, const int rowB, const int colB)
{
	int i,j,k,index;
	double sum,elementA,elementB;
	const int rowC = rowA;
	const int colC = colB;
	for(i=0;i<rowC;i++)
	{
		for(j=0;j<colC;j++)
		{
			sum = 0.0;
			for(k=0;k<colA;k++)
			{
				elementA = getMyMatrix(a,rowA,colA,i,k); 
				elementB = getMyMatrix(b,rowB,colB,k,j);
				sum += elementA * elementB;
			}
			index = indexConvert(rowC,colC,i,j);
			c[index] = sum;
		}
	}
}


int isNullMyMatrix(const double *a,const int rowA,const int colA)
{
	int i;
	double val;
	const totalNum = rowA*colA;
	for(i=0;i<totalNum;i++)
	{
		if(a[i]!=0) return 0;
	}
	return 1;
}


void scaleMyMatrix(double *c, const double k, const double *a,const int rowA,const int colA)
{
	int i;
	const int totalNum = rowA*colA;
	for(i=0;i<totalNum;i++) c[i] = k*a[i];
}


void transMyMatrix(double *c, const double *a,const int rowA,const int colA)
{
	int i,j;
	int index;
	for(i=0;i<rowA;i++)
	{
		for(j=0;j<colA;j++)
		{
			index = indexConvert(colA,rowA,j,i);
			c[index] = a[indexConvert(rowA,colA,i,j)];
		}
	}
}


double traceMyMatrix(double *a,const int rowA,const int colA)
{
	int i;
	int index;
	double sum = 0.0;
	for(i=0;i<rowA;i++)
	{
		sum += a[indexConvert(rowA,colA,i,i)];
	}
	return sum;
}



void copyMyMatrix(double *dest, const double *src, const int row, const int col)
{
	memcpy(dest,src,sizeof(double)*row*col);
	/*
	int i,j;
	int index;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			index = indexConvert(row,col,i,j);
			dest[index] = src[index];
		}
	}
	*/
}



double absMaxMyMatrix(const double *src, const int row, const int col)
{
	double val;
	double ret;
	int i,j;

	ret = fabs(getMyMatrix(src,row,col,0,0));
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			val = fabs(getMyMatrix(src,row,col,i,j));
			if(val > ret) ret = val;	
		}
	}

	return ret;
}



/*
int main(void)
{
	const int row = 30;
	const int col = 900;
	double *x = malloc(sizeof(double)*row*col);

	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			const double tmp = i+j + 0.5;
			setMyMatrix(x,tmp,row,col,i,j);
		}
	}
	
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			printf("%lf ",getMyMatrix(x,row,col,i,j));
		}
		printf("\n\n");
	}

}
*/
