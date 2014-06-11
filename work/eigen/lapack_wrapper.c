/*

	Lapack_wrapper.c	
	
	This file contains an interface towards some commonly used matrix algebra 
	functions.
	
	See individual functions for more info on usage and functioning.
	
	A few things to keep in mind:
	
	- Input parameters always come first, output parameters last. e.g.
		
		some_function(m, n, a, X, Y)
		m, n, a, X are input, result will be stored in Y.
	
	- All floating point types are "double". All integer types are "int".
	
	- All matrices passed as argument need to be in Fortran vector format.
	  This means the argument is a pointer to an array of doubles. This array
	  stores the matrix columnwise.
	  
	- All input and output parameters need to be properly allocated. The 
	  functions do not dynamically allocate memory for output variables. 
	  
	- The functions use the Lapack and the Atlas-blas libraries. You need to
	  have these installed in order to use this file.
	
	
	More info on
	http://itf.fys.kuleuven.be/~rob/computer/lapack_wrapper/index.htm
	
	
	Rob Heylen, aug 2006
	
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include </usr/local/atlas/include/cblas.h>
#include </usr/local/atlas/include/clapack.h>
#include "lapack_wrapper.h"


// Wrapper for Lapack eigenvalue function
int ssyevr (char JOBZ, char RANGE, char UPLO, int N,
        float *A, int LDA, float VL, float VU,
        int IL, int IU, float ABSTOL, int *M,
        float *W, float *Z, int LDZ, int *ISUPPZ,
        float *WORK, int LWORK, int *IWORK, int LIWORK) 
{
        extern  void  ssyevr_ (char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                 float *A, int *LDAp, float *VLp, float *VUp,
                 int *ILp, int *IUp, float *ABSTOLp, int *Mp,
                 float *W, float *Z, int *LDZp, int *ISUPPZ,
                 float *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                 int *INFOp);
        int  INFO;
        ssyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
           &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
           WORK, &LWORK, IWORK, &LIWORK, &INFO);

  return  INFO;
}

// Wrapper for Lapack eigenvalue function
int dsyevr (char JOBZ, char RANGE, char UPLO, int N,
	double *A, int LDA, double VL, double VU,
	int IL, int IU, double ABSTOL, int *M,
	double *W, double *Z, int LDZ, int *ISUPPZ,
	double *WORK, int LWORK, int *IWORK, int LIWORK) 
{
	extern  void  dsyevr_ (char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
		 double *A, int *LDAp, double *VLp, double *VUp,
		 int *ILp, int *IUp, double *ABSTOLp, int *Mp,
		 double *W, double *Z, int *LDZp, int *ISUPPZ,
		 double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
		 int *INFOp);
	int  INFO;
	dsyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
	   &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
	   WORK, &LWORK, IWORK, &LIWORK, &INFO);

  return  INFO;
}

// Wrapper for Lapack machine precision function
static double dlamch (char CMACH)
{
  extern  double  dlamch_ (char *CMACHp);
  return  dlamch_ (&CMACH);
}

int check_input(double* x, double*y, char* c) {
	int out=0;
	if (x==y) {
		printf("Error in %s: input equals output \n", c);
		out=1;
	}
	return out;
}

int check_inputf(float* x, float*y, char* c) {
        int out=0;
        if (x==y) {
                printf("Error in %s: input equals output \n", c);
                out=1;
        }
        return out;
}

void mat2cvec(int m, int n, double** mat, double* cvec) {
	/*
		Turns a matrix (given by a double pointer) into its C vector format 
		(single vector, rowwise). The matrix "mat" needs to be an n*m matrix
		The vector "vec" is a vector of lenght nm
	*/
	int i,j;
	for (i=0; i<m; i++) for (j=0; j<n; j++) cvec[i*n+j] = mat[i][j];
}

void mat2fvec(int m, int n, double** mat, double* fvec) {
	/*
		Turns a matrix (given by a double pointer) into its Fortran vector format 
		(single vector, colunmwise). The matrix "mat" needs to be an m*n matrix
		The vector "vec" is a vector of lenght mn
	*/
	int i,j;
	for (i=0; i<m; i++) for (j=0; j<n; j++) fvec[i+m*j] = mat[i][j];
}

void cvec2mat(int m, int n, double* cvec, double** mat) {
	/*
		Turns a matrix in C vector format into a matrix given by a double pointer
		The matrices have dimension m*n
	*/
	int i,j;
	for (i=0; i<m; i++) for (j=0; j<n; j++)  mat[i][j] = cvec[i*n+j];
}

void fvec2mat(int m, int n, double* fvec, double** mat) {
	/*
		Turns a matrix in C vector format into a matrix given by a double pointer
		The matrices have dimension m*n
	*/
	int i,j;
	for (i=0; i<m; i++) for (j=0; j<n; j++)  mat[i][j] = fvec[i+m*j];
}

void cvec2fvec(int m, int n, double *cvec, double* fvec) {
	/*
		Turn matrix in C vector format into matrix in Fortran vector format
		Matrix has dimension m*n
	*/
	int i,j;
	check_input(cvec, fvec, "cvec2fvec");
	for (i=0; i<m; i++) for (j=0; j<n; j++) fvec[i+m*j] = cvec[n*i+j];
}

void fvec2cvec(int m, int n, double *fvec, double* cvec) {
	/*
		Turn matrix in Fortran vector format into matrix in C vector format
		Matrix has dimension m*n
	*/
	int i,j;
	check_input(cvec, fvec, "fvec2cvec");
	for (i=0; i<m; i++) for (j=0; j<n; j++) cvec[n*i+j] = fvec[i+m*j];
}

void matrix_matrix_mult(int m, int n, int k, double alpha, double beta, double* A, double* B, double* C){
	/*
		Calculates C = a*A*B + b*C
		A, B and C need to be matrices in Fortran vector format (single vector, columnwise)
		dimensions: A is m*k, B is k*n, so C has to be m*n
		This is an O(n^3) operation
	*/
	check_input(C, A, "matrix_matrix_mult");
	check_input(C, B, "matrix_matrix_mult");
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, m, B, k, beta, C, m);
}

void matrix_add(int m, int n, double a, double **X, double **Y) {
	/*
		Calculates Y = a*X + Y
		X and Y are m*n matrices in double pointer format	
	*/
	int i,j;
	for (i=0; i<m; i++) for (j=0; j<n; j++) Y[i][j]+=a*X[i][j];
}

void vector_add(int n, double a, double *X, double *Y) {
	/*
		Calculates Y = a*X + Y
		X and Y are vectors of length n
	*/
	check_input(X, Y, "vector_add");
	cblas_daxpy(n, a, X, 1, Y, 1);
}

double dotprod(int n, double *X, double *Y) {
	/*
		Returns the dot-product of X and Y
	*/
	return cblas_ddot(n, X, 1, Y, 1);
}

void vector_copyf(int n, float* X, float* Y) {
        /*
                Copies a vector X of lenght n to vector Y
        */
        int i;
        for (i=0; i<n; i++) Y[i]=X[i];
}

void vector_copy(int n, double* X, double* Y) {
	/*
		Copies a vector X of lenght n to vector Y
	*/
	int i;
	for (i=0; i<n; i++) Y[i]=X[i];
}

void vector_scale(int n, double a, double* X) {
	/*
		Scale vector of length n: X = a*X
	*/
	int i;
	for (i=0; i<n; i++) X[i]=a*X[i];
}

void matrix_vector_mult(int m, int n, double a, double b, double *A, double *x, double *y) {
	/*
		Calculates y = a*A*x + b*y
		A is m*n matrix in Fortran vector format, so x and y need to have length m
	*/
	check_input(x, y, "matrix_vector_mult");
	cblas_dgemv (CblasColMajor, CblasNoTrans, m, n, a, A, m, x, 1, b, y, 1);
}

void matrix_transpose(int m, int n, double *X, double *Y) {
	/*
		Transposes an m*n matrix: Y = X^T
		Matrix can be in either C or Fortran vector format
	*/
	check_input(X,Y, "matrix_transpose"); 
	int i,j;
	for (i=0; i<m; i++) for (j=0; j<n; j++) Y[n*i+j] = X[i+m*j];
}

int eigen_decompositionf(int n, float* X, float *eigvec, float *eigval) {
        /*
                This function calculates the eigenvalues and eigenvectors of 
                the n*n symmetric matrix X. 
                The matrices have to be in Fortran vector format.
                The eigenvectors will be put columnwise in the n*n matrix eigvec,
                where the corresponding eigenvalues will be put in the vector 
                eigval (length n of course). Only the lower triangle of the matrix
                X is used. The content of X is not changed.
                
                This function first queries the Lapack routines for optimal workspace 
                sizes. These memoryblocks are then allocated and the decomposition is 
                calculated using the Lapack function "dsyevr". The allocated memory 
                is then freed. 
        */
        
        float *WORK, *Xc;
        int *ISUPPZ, *IWORK;
        int  numeig, info, sizeWORK, sizeIWORK;
        
        if (check_inputf(X, eigvec, "eigen_decomposition")) return 1;
        
        /*  Use a copy of X so we don't need to change its value or use its memoryblock */
        Xc=malloc(n*n*sizeof(float));
                
        /*  The support of the eigenvectors. We will not use this but the routine needs it  */
        ISUPPZ = malloc (2*n*sizeof(float));
        
        /*  Allocate temporarily minimally allowed size for workspace arrays */
        WORK = malloc (26*n*sizeof(float));
        IWORK = malloc (10*n*sizeof(int));
                
        /*  Check for NULL-pointers.  */
        if ((Xc==NULL)||(ISUPPZ==NULL)||(WORK==NULL)||(IWORK==NULL)) {
                printf("malloc failed in eigen_decomposition\n"); 
                return 2;
        }
        
        vector_copyf(n*n, X, Xc);
        
        /*  Query the Lapack routine for optimal sizes for workspace arrays  */
        info=ssyevr ('V', 'A', 'L', n, Xc, n, 0, 0, 0, 0, dlamch('S'), &numeig, eigval, eigvec, n, ISUPPZ, WORK, -1, IWORK, -1);
        sizeWORK = (int)WORK[0]; 
        sizeIWORK = IWORK[0]; 
        
        /*  Free previous allocation and reallocate preferable workspaces, Check result  */
        free(WORK);free(IWORK);
        WORK = malloc (sizeWORK*sizeof(float));
        IWORK = malloc (sizeIWORK*sizeof(int));
        if ((WORK==NULL)||(IWORK==NULL)) {
                printf("malloc failed in eigen_decomposition\n"); 
                return 2;
        }
        printf("starting\n");
        fflush(stdout);
        
        /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
        info=ssyevr ('V', 'A', 'L', n, Xc, n, 0, 0, 0, 0, dlamch('S'), &numeig, eigval, eigvec, n, ISUPPZ, WORK, sizeWORK, IWORK, sizeIWORK);
        
        /*  Cleanup and exit  */
        free(WORK); free(IWORK); free(ISUPPZ); free(Xc);
        return info;
}       

int eigen_decomposition(int n, double* X, double *eigvec, double *eigval) {
	/*
		This function calculates the eigenvalues and eigenvectors of 
		the n*n symmetric matrix X. 
		The matrices have to be in Fortran vector format.
		The eigenvectors will be put columnwise in the n*n matrix eigvec,
		where the corresponding eigenvalues will be put in the vector 
		eigval (length n of course). Only the lower triangle of the matrix
		X is used. The content of X is not changed.
		
		This function first queries the Lapack routines for optimal workspace 
		sizes. These memoryblocks are then allocated and the decomposition is 
		calculated using the Lapack function "dsyevr". The allocated memory 
		is then freed. 
	*/
	
	double *WORK, *Xc;
	int *ISUPPZ, *IWORK;
	int  numeig, info, sizeWORK, sizeIWORK;
	
	if (check_input(X, eigvec, "eigen_decomposition")) return 1;
	
	/*  Use a copy of X so we don't need to change its value or use its memoryblock */
	Xc=malloc(n*n*sizeof(double));
		
	/*  The support of the eigenvectors. We will not use this but the routine needs it  */
	ISUPPZ = malloc (2*n*sizeof(int));
	
	/*  Allocate temporarily minimally allowed size for workspace arrays */
	WORK = malloc (26*n*sizeof(double));
	IWORK = malloc (10*n*sizeof(int));
		
	/*  Check for NULL-pointers.  */
	if ((Xc==NULL)||(ISUPPZ==NULL)||(WORK==NULL)||(IWORK==NULL)) {
		printf("malloc failed in eigen_decomposition\n"); 
		return 2;
	}
	
	vector_copy(n*n, X, Xc);
	
	/*  Query the Lapack routine for optimal sizes for workspace arrays  */
	info=dsyevr ('V', 'A', 'L', n, Xc, n, 0, 0, 0, 0, dlamch('S'), &numeig, eigval, eigvec, n, ISUPPZ, WORK, -1, IWORK, -1);
	sizeWORK = (int)WORK[0]; 
	sizeIWORK = IWORK[0]; 
	
	/*  Free previous allocation and reallocate preferable workspaces, Check result  */
	free(WORK);free(IWORK);
	WORK = malloc (sizeWORK*sizeof(double));
	IWORK = malloc (sizeIWORK*sizeof(int));
	if ((WORK==NULL)||(IWORK==NULL)) {
		printf("malloc failed in eigen_decomposition\n"); 
		return 2;
	}
	
	/*  Now calculate the eigenvalues and vectors using optimal workspaces  */
	info=dsyevr ('V', 'A', 'L', n, Xc, n, 0, 0, 0, 0, dlamch('S'), &numeig, eigval, eigvec, n, ISUPPZ, WORK, sizeWORK, IWORK, sizeIWORK);
	
	/*  Cleanup and exit  */
	free(WORK); free(IWORK); free(ISUPPZ); free(Xc);
	return info;
}	

int matrix_square_root_n(int n, double *X, double *I, double *Y) {

	/*  
		This function calculates one of the square roots of the matrix X and stores it in Y:
		X = sqrt(Y);
		X needs to be a symmetric positive definite matrix of dimension n*n in Fortran vector
		format. Y is of course a vector of similar size, and will contain the result on exit.
		Y will be used as a workspace variable in the meantime.
		The variable I is a vector of length n containing +1 and -1 elements. It can be used 
		to select one of the 2^n different square roots of X
		
		This function first calculates the eigenvalue decomposition of X: X = U*D*U^T
		A new matrix F is then calculated with on the diagonal the square roots of D, with
		signs taken from I.
		The square root is then obtained by calculating U*F*U^T
	*/
	
	if (check_input(X, Y, "matrix_square_root_n")) return 1;
	
	double *eigval, *eigvec, *temp;
	int info = 0, i, j;
	
	/*  Memory allocation  */
	eigval=malloc(n*sizeof(double));
	eigvec=malloc(n*n*sizeof(double));
	temp=malloc(n*n*sizeof(double));
	if ((eigval==NULL)||(eigvec==NULL)||(temp==NULL)) {
		printf("malloc failed in matrix_square_root_n\n"); 
		return 2;
	}
	
	/*  Eigen decomposition  */
	info=eigen_decomposition(n, X, eigvec, eigval);
	if (info != 0) return info;
	
	/*  Check for positive definitiveness*/
	for (i=0; i<n; i++) if (eigval[i]<0) {
		fprintf(stderr, "In matrix_square_root_n: Matrix is not positive definite.\n");
		return 1;
	}
	
	/*  Build square rooted diagonal matrix, with sign signature I  */
	for (i=0; i<n; i++) for (j=0; j<n; j++) Y[i*n+j] = 0.0;
	for (i=0; i<n; i++) Y[i*n+i] = I[i]*sqrt(eigval[i]);
	
	/*  Multiply the eigenvectors with this diagonal matrix Y. Store back in Y  */
	matrix_matrix_mult(n, n, n, 1.0, 0, eigvec, Y, temp);
	vector_copy(n*n, temp, Y);
	
	/*  Transpose eigenvectors. Store in temp  */
	matrix_transpose(n, n, eigvec, temp);
	
	/*  Multiply Y with this temp. Store in eigvec which is no longer required. Copy to Y  */
	matrix_matrix_mult(n, n, n, 1.0, 0, Y, temp, eigvec);
	vector_copy(n*n, eigvec, Y);
			
	/*  Cleanup and exit  */
	free(eigvec); free(eigval); free(temp);
	return info;
}

int matrix_square_root(int n, double *X, double *Y) {
	/*
		Calculates the unique positive semi-definite square root of the n*n positive definite
		symmetric matrix X. This function calls matrix_square_root_n() with I=[1 1 ... 1]. 
		See matrix_square_root_n() for more details
	*/
	
	if (check_input(X, Y, "matrix_square_root")) return 1;
	
	int i, info=0;
	double *I;
	I = malloc(n*sizeof(double));
	if (I==NULL) {
		printf("malloc failed in matrix_square_root\n"); 
		return 2;
	}
	for (i=0; i<n; i++) I[i]=1.0;
	info=matrix_square_root_n(n, X, I, Y);
	free(I);
	return info;
}

int matrix_invert(int n, double *X, double *Y) {
	/*  
		Calculates the inverse of the n*n matrix X: Y = X^-1
		Does not change the value of X, unless Y=X
	*/
	
	int info=0;
	
	/*  When X!=Y we want to keep X unchanged. Copy to Y and use this as working variable  */
	if (X!=Y) vector_copy(n*n, X, Y);
	
	/*  We need to store the pivot matrix obtained by the LU factorisation  */
	int *ipiv;
	ipiv=malloc(n*sizeof(int));
	if (ipiv==NULL) {
		printf("malloc failed in matrix_invert\n"); 
		return 2;
	}
	
	/*  Turn Y into its LU form, store pivot matrix  */
	info = clapack_dgetrf (CblasColMajor, n, n, Y, n, ipiv);
	
	/*  Don't bother continuing when illegal argument (info<0) or singularity (info>0) occurs  */
	if (info!=0) return info;
		
	/*  Feed this to the lapack inversion routine.  */
	info = clapack_dgetri (CblasColMajor, n, Y, n, ipiv);
	
	/*  Cleanup and exit  */
 	free(ipiv);
	return info;
}

int linsolve(int n, double *A, double *B, double *x) {
	/*
		Solves the matrix equation A*x = B for x
		A is an n*n matrix in Fortran vector notation, B and x are vectors of lenght n
		A and B are not changed. x will contain the solution.
		x can equal B if B may be overwritten.
	*/
	
	int info=0;
	
	/*  When x!=B we want to keep B unchanged. Copy to x and use this as working variable  */
	if (x!=B) vector_copy(n, B, x);
	
	/*  We need space for the pivot matrix  */
	int *ipiv;
	ipiv=malloc(n*sizeof(int));
	if (ipiv==NULL) {
		printf("malloc failed in linsolve\n"); 
		return 2;
	}
	
	/*  Now we can call the Lapack function  */
	info = clapack_dgesv (CblasColMajor, n, 1, A, n, ipiv, x, n);
	
	/*  Cleanup and exit  */
	free(ipiv);
	return info;
}

