#ifndef GSL_EXTERN_H
#define GSL_EXTERN_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>


int gsl_matrix_row_num(const gsl_matrix *a);
int gsl_matrix_col_num(const gsl_matrix *a);
double gsl_matrix_trace(const gsl_matrix *a);
double gsl_matrix_abs_max(const gsl_matrix *a);
double gsl_vector_abs_max(const gsl_vector *a);
void gsl_matrix_solve(gsl_vector *x,const gsl_matrix *a,const gsl_vector *b);

// convert a matrix from 2d gsl_matrix to 1d double *dest,
// the data save in dest is row-major placement
void gsl_matrix_to_double(double *dest,const gsl_matrix *src);



// convert a matrix from 1d dobule *dest to 2d gsl_matrix
// the data save in src is row-major placement
// note: dest should be allocated before
void double_to_gsl_matrix(gsl_matrix *dest, const double *src);


#endif
