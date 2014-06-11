/*

	Lapack_wrapper.h
	
	Header file for Lapack wrapper functions. See lapack_wrapper.c for more details
	
	Rob Heylen, aug 2006
	http://itf.fys.kuleuven.be/~rob/computer/lapack_wrapper/index.htm

*/

#ifndef __lapack_wrapper_h__
#define __lapack_wrapper_h__


#ifdef __cplusplus
extern "C"
{
#endif
   
void mat2cvec(int m, int n, double** mat, double* cvec);

void mat2fvec(int m, int n, double** mat, double* fvec);

void cvec2mat(int m, int n, double* cvec, double** mat);

void fvec2mat(int m, int n, double* fvec, double** mat);

void cvec2fvec(int m, int n, double *cvec, double* fvec);

void fvec2cvec(int m, int n, double *fvec, double* cvec);

void matrix_matrix_mult(int m, int n, int k, double alpha, double beta, double* A, double* B, double* C);

void matrix_add(int m, int n, double a, double **X, double **Y);

void vector_add(int n, double a, double *X, double *Y);

double dotprod(int n, double *X, double *Y);

void vector_copy(int n, double* X, double* Y);

void vector_scale(int n, double a, double* X);

void matrix_vector_mult(int m, int n, double a, double b, double *A, double *x, double *y);

void matrix_transpose(int m, int n, double *X, double *Y);

int eigen_decompositionf(int n, float* X, float *eigvec, float *eigval);
int eigen_decomposition(int n, double* X, double *eigvec, double *eigval);

int matrix_square_root_n(int n, double *X, double *I, double *Y);

int matrix_square_root(int n, double *X, double *Y);

int matrix_invert(int n, double *X, double *Y);

int linsolve(int n, double *A, double *B, double *x);

#ifdef __cplusplus
} //extern "C"
#endif

#endif // __lapack_wrapper_h__
