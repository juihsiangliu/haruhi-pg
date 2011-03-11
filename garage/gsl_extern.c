#include "gsl_extern.h"


static int indexConvert(const gsl_matrix *src,const int rowIndex,const int colIndex)
{
	const int col = gsl_matrix_col_num(src);
	return (rowIndex*col + colIndex);
}

//========================================



int gsl_matrix_row_num(const gsl_matrix *a)
{
	return (int)a->size1;
}

int gsl_matrix_col_num(const gsl_matrix *a)
{
	return (int)a->size2;
}


double gsl_matrix_trace(const gsl_matrix *a)
{

	double ret = 0;
	const int num = (int)a->size1;
	int i;
	for(i=0;i<num;i++)
	{
		ret += gsl_matrix_get(a,i,i);
	}
	return ret;
}

double gsl_matrix_abs_max(const gsl_matrix *a)
{
	const int row = gsl_matrix_row_num(a);
	const int col = gsl_matrix_col_num(a);
	double ret = 0;
	int i,j;
	double tmp;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			tmp = gsl_matrix_get(a,i,j);
			if(abs(tmp) > ret) ret = tmp;
		}
	}

	return ret;
}


double gsl_vector_abs_max(const gsl_vector *a)
{
	const int num = (int)a->size;
	double val;
	int i;
	gsl_vector *aCopy = gsl_vector_alloc(num);
	gsl_vector_memcpy(aCopy,a);
	for(i=0;i<num;i++)
	{
		val = gsl_vector_get(aCopy,i);
		if(val < 0)
		{
			val = -1*val;
			gsl_vector_set(aCopy,i,val);
		}
	}
	return gsl_vector_max(aCopy);
}


void gsl_matrix_solve(gsl_vector *x,const gsl_matrix *a,const gsl_vector *b)
{
	int s;
	const int size = gsl_matrix_row_num(a);
	gsl_matrix *aCopy = gsl_matrix_alloc(size,size);
	gsl_permutation * p = gsl_permutation_alloc(size);
	gsl_matrix_memcpy(aCopy,a);
	gsl_linalg_LU_decomp(aCopy, p, &s);
	gsl_linalg_LU_solve(aCopy,p,b,x);
	gsl_permutation_free(p);
	gsl_matrix_free(aCopy);
}

     
			     

// convert a matrix from 2d gsl_matrix to 1d double *dest,
// the data save in dest is row-major placement
void gsl_matrix_to_double(double *dest,const gsl_matrix *src)
{
	const int row = gsl_matrix_row_num(src);
	const int col = gsl_matrix_col_num(src);
	int i,j,index;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			index = indexConvert(src,i,j);
			dest[index] = gsl_matrix_get(src,i,j);
		}
	}
}



// convert a matrix from 1d dobule *dest to 2d gsl_matrix
// the data save in src is row-major placement
// note: dest should be allocated before
void double_to_gsl_matrix(gsl_matrix *dest, const double *src)
{
	const int row = gsl_matrix_row_num(dest);
	const int col = gsl_matrix_col_num(dest);
	int i,j,index;
	double val;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			index = col*i + j;
			val = src[index];
			gsl_matrix_set(dest,i,j,val);
		}
	}
}



