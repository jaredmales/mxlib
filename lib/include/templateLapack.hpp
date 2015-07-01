/** \file templateLapack.hpp
  * \brief Declares and defines templatized wrappers for the Lapack library
  * \ingroup genutils
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __templateLapack_hpp__
#define __templateLapack_hpp__

/// Determine machine parameters. 
/** Wrapper for Lapack xLAMCH
  *
  * See more details at http://www.netlib.org/lapack/lapack-3.1.1/html/slamch.f.html.
  * 
  * * \tparam dataT is the data type of the arrays, and determines which underlying Lapack routine is called.
  *
  * \param[in] CMACH Specifies the value to be returned by SLAMCH: <pre>
  *     = 'E' or 'e',   returns eps, the relative machine precision  
  *     = 'S' or 's ,   returns sfmin, the safe minimum, such that 1/sfmin does not overflow  
  *     = 'B' or 'b',   returns base, the base of the machine 
  *     = 'P' or 'p',   returns eps*base  
  *     = 'N' or 'n',   returns t, the number of (base) digits in the mantissa  
  *     = 'R' or 'r',   returns rnd, 1.0 when rounding occurs in addition, 0.0 otherwise  
  *     = 'M' or 'm',   returns emin, the minimum exponent before (gradual) underflow  
  *     = 'U' or 'u',   returns rmin, the underflow threshold - base^(emin-1)  
  *     = 'L' or 'l',   returns emax, the largest exponent before overflow  
  *     = 'O' or 'o',   returns rmax, the overflow threshold  - (base^emax)*(1-eps)  
  * </pre>
  * 
  * \returns the value of the specified machine parameters for the specified precision
  */

template<typename dataT>
dataT lamch (char CMACH)
{
   return -1;
}


extern "C"
{
   extern  float  slamch_ (char *CMACHp);
   extern  float  dlamch_ (char *CMACHp);
}
      
// Float specialization of lamch, a wrapper for Lapack SLAMCH
template<>
inline
float lamch<float>(char CMACH)
{ 
   return  slamch_ (&CMACH);
}

// Double specialization of lamch, a wrapper for Lapack DLAMCH
template<>
inline
double lamch<double>(char CMACH)
{ 
   return  dlamch_ (&CMACH);
}

/// Reduce a real symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation
/** xSYTRD reduces a real symmetric matrix A to real symmetric
  *  tridiagonal form T by an orthogonal similarity transformation:
  * 
  * \f$ Q^T * A * Q = T. \f$
  * 
  * For more see: http://www.netlib.org/lapack/lapack-3.1.1/html/ssytrd.f.html
  * 
  * \tparam dataT is the data type
  *  
  * \param UPLO 'U':  Upper triangle of A is stored, 'L':  Lower triangle of A is stored.
  * \param N The order of the matrix A.  N >= 0.
  * \param A  array, dimension (LDA,N)
  * \param LDA The leading dimension of the array A.  LDA >= max(1,N).
  * \param D (output) array, dimension (N), the diagonal elements of the tridiagonal matrix T:
  * \param E (output) array, dimension (N-1), the off-diagonal elements of the tridiagonal matrix T:
  * \param TAU (output) array, dimension (N-1), the scalar factors of the elementary reflectors (see Further Details).
  * \param WORK (workspace/output) array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  * \param LWORK (input)  The dimension of the array WORK.  LWORK >= 1. For optimum performance LWORK >= N*NB, where NB is the optimal blocksize.
  * \param INFO (output)  0:  successful exit < 0:  if INFO = -i, the i-th argument had an illegal value
  * 
  * \returns the value of INFO from the LAPACK routine
  */
template<typename dataT>
int sytrd( char UPLO, int N, dataT * A, int LDA, dataT *D, dataT *E, dataT *TAU, dataT *WORK, int LWORK, int INFO)
{
   return -1;
}

//Declarations of the actual LAPACK functions
extern "C"
{
   void ssytrd_( char * UPLO, int * N, float * A, int * LDA, float *D, float *E, float* TAU, float *WORK, int *LWORK, int *INFO );
   void dsytrd_( char * UPLO, int * N, double * A, int * LDA, double *D, double *E, double* TAU, double *WORK, int *LWORK, int *INFO );
}

template<>
inline
int sytrd<float>( char UPLO, int N, float * A, int LDA, float *D, float *E, float *TAU, float *WORK, int LWORK, int INFO)
{
  
   ssytrd_(&UPLO, &N, A, &LDA, D, E, TAU, WORK, &LWORK, &INFO);
   
   return INFO;
}

template<>
inline
int sytrd<double>( char UPLO, int N, double * A, int LDA, double *D, double *E, double *TAU, double *WORK, int LWORK, int INFO)
{
  
   dsytrd_(&UPLO, &N, A, &LDA, D, E, TAU, WORK, &LWORK, &INFO);
   
   return INFO;
}




/// Compute selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix
/** xSYEVR computes selected eigenvalues and, optionally, eigenvectors
  * of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
  * selected by specifying either a range of values or a range of
  * indices for the desired eigenvalues.
  * 
  * See more details: http://www.netlib.org/lapack/lapack-3.1.1/html/ssyevr.f.html.
  */
template<typename dataT>
int syevr (  char JOBZ, char RANGE, char UPLO, int N, dataT *A, int LDA, dataT VL, dataT VU, int IL, 
             int IU, dataT ABSTOL, int *M, dataT *W, dataT *Z, int LDZ, int *ISUPPZ, dataT *WORK, int
             LWORK, int *IWORK, int LIWORK ) 
{
   return -1;
}


extern "C"
{
   void  ssyevr_ (char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                             float *A, int *LDAp, float *VLp, float *VUp,
                              int *ILp, int *IUp, float *ABSTOLp, int *Mp,
                               float *W, float *Z, int *LDZp, int *ISUPPZ,
                                float *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                                 int *INFOp);
   
   void  dsyevr_ (char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                             double *A, int *LDAp, double *VLp, double *VUp,
                              int *ILp, int *IUp, double *ABSTOLp, int *Mp,
                               double *W, double *Z, int *LDZp, int *ISUPPZ,
                                double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                                 int *INFOp);
}

      
// Float specialization of syevr, a wrapper for Lapack SSYEVR
template<>
inline
int syevr<float> ( char JOBZ, char RANGE, char UPLO, int N, float *A, int LDA, float VL, float VU,
                    int IL, int IU,  float ABSTOL, int *M, float *W, float *Z, int LDZ, int *ISUPPZ,
                     float *WORK, int LWORK, int *IWORK, int LIWORK ) 
{

   int  INFO;
   
   ssyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
           &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
           WORK, &LWORK, IWORK, &LIWORK, &INFO);

   return  INFO;
}

// Double specialization of syevr, a wrapper for Lapack DSYEVR
template<>
inline
int syevr<double> ( char JOBZ, char RANGE, char UPLO, int N, double *A, int LDA, double VL, double VU,
                    int IL, int IU,  double ABSTOL, int *M, double *W, double *Z, int LDZ, int *ISUPPZ,
                     double *WORK, int LWORK, int *IWORK, int LIWORK ) 
{

   int  INFO;
   
   dsyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
           &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
           WORK, &LWORK, IWORK, &LIWORK, &INFO);

   return  INFO;
}


/// Compute the singular value decomposition (SVD) of a real matrix
/** xGESVD computes the singular value decomposition (SVD) of a real
  * M-by-N matrix A, optionally computing the left and/or right singular
  * vectors. The SVD is written
  *
  *  \f$  A = U * \Sigma * V^T \f$
  *
  * where \f$ \Sigma \f$ is an M-by-N matrix which is zero except for its
  * min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  * V is an N-by-N orthogonal matrix.  The diagonal elements of \f$ \Sigma \f$
  * are the singular values of A; they are real and non-negative, and
  * are returned in descending order.  The first min(m,n) columns of
  * U and V are the left and right singular vectors of A.
  *
  * Note that the routine returns \f$ V^T \f$, not \f$ V \f$.
  * 
  * See more details: http://www.netlib.org/lapack/explore-html/d8/d49/sgesvd_8f.html.
  * This documentation is taken from there.
  * 
  * \tparam dataT is the data type of the arrays, and determines which underlying Lapack routine is called.
  * 
  * \param[in] JOBU
  * \parblock
    (char) Specifies options for computing all or part of the matrix U: <br />
     = 'A':  all M columns of U are returned in array U  <br />
     = 'S':  the first min(m,n) columns of U (the left singular vectors) 
               are returned in the array U   <br /> 
     = 'O':  the first min(m,n) columns of U (the left singular  <br />
           vectors) are overwritten on the array A  <br />
     = 'N':  no columns of U (no left singular vectors) are
           computed.
    
  * \param[in] JOBVT
  * \parblock
    (char) Specifies options for computing all or part of the matrix  V**T:
    = 'A':  all N rows of V**T are returned in the array VT <br />
    = 'S':  the first min(m,n) rows of V**T (the right singular
            vectors) are returned in the array VT <br />
    = 'O':  the first min(m,n) rows of V**T (the right singular
            vectors) are overwritten on the array A<br />
    = 'N':  no rows of V**T (no right singular vectors) are
            computed.
  
    JOBVT and JOBU cannot both be 'O'.
  *
  * \param[in] M
  * \parblock
    (int)
    The number of rows of the input matrix A.  M >= 0.
  *
  * \param[in] N
  * \parblock
    (int)
    The number of columns of the input matrix A.  N >= 0.
  *
  * \param[in,out] A
  * \parblock
     (dataT *, dimension (LDA,N))
     On entry: the M-by-N matrix A.<br />
     On exit:<br />
     if JOBU = 'O',  A is overwritten with the first min(m,n)
                     columns of U (the left singular vectors,
                     stored columnwise)<br />
     if JOBVT = 'O', A is overwritten with the first min(m,n)
                     rows of V**T (the right singular vectors,
                     stored rowwise)<br />
     if JOBU != 'O' and JOBVT .ne. 'O', the contents of A
                     are destroyed.
  *
  * \param[in] LDA
  * \parblock
    (int)
    The leading dimension of the array A.  LDA >= max(1,M).
  *
  * \param[out] S
  * \parblock
     (dataT *, dimension (min(M,N))
     The singular values of A, sorted so that S(i) >= S(i+1).
  * 
  * \param[out] U
  * \parblock
     (dataT *, dimension (LDU,UCOL))
     (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.<br />
     If JOBU = 'A', U contains the M-by-M orthogonal matrix U<br />
     if JOBU = 'S', U contains the first min(m,n) columns of U
     (the left singular vectors, stored columnwise)<br />
     if JOBU = 'N' or 'O', U is not referenced.
  *
  * \param[in] LDU
  * \parblock
     (int)
     The leading dimension of the array U.  <br />
     LDU >= 1<br />
     if JOBU == 'S' or 'A', LDU >= M.
  *
  * \param[out] VT
  * \parblock
     (dataT *, dimension (LDVT,N))
     If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
     V**T <br />
     if JOBVT = 'S', VT contains the first min(m,n) rows of
     V**T (the right singular vectors, stored rowwise) <br />
     if JOBVT = 'N' or 'O', VT is not referenced.
  *
  * \param[in] LDVT
  * \parblock
     (int)
     The leading dimension of the array VT. <br /> 
     LDVT >= 1<br /> 
     if JOBVT = 'A', LDVT >= N<br /> 
     if JOBVT = 'S', LDVT >= min(M,N).
  *
  * \param[out] WORK
   \parblock
     (dataT *, dimension (MAX(1,LWORK)) )<br />
     On exit, if INFO = 0, WORK[1] returns the optimal LWORK<br />
     if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
     superdiagonal elements of an upper bidiagonal matrix B
     whose diagonal is in S (not necessarily sorted). B
     satisfies A = U * B * VT, so it has the same singular values
     as A, and singular vectors related by U and VT.
  *
  * \param[in] LWORK
  \parblock
     (int)
     The dimension of the array WORK.<br />
     LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
        - PATH 1  (M much larger than N, JOBU='N') 
        - PATH 1t (N much larger than M, JOBVT='N')
     LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths 
     For good performance, LWORK should generally be larger.
     
     If LWORK = -1, then a workspace query is assumed; the routine
     only calculates the optimal size of the WORK array, returns
     this value as the first entry of the WORK array, and no error
     message related to LWORK is issued by XERBLA.
  * \endparblock
  * \returns
  * \parblock
     =0:  successful exit. <br />
     <0:  if INFO = -i, the i-th argument had an illegal value. <br />
     >0:  if SBDSQR did not converge, INFO specifies how many
           superdiagonals of an intermediate bidiagonal form B
           did not converge to zero. See the description of WORK
          above for details.
    \endparblock
  */
template<typename dataT>
int gesvd( char JOBU, char JOBVT, int M, int N, dataT * A, int LDA, dataT * S, dataT *U, int LDU, 
                dataT * VT, int LDVT, dataT * WORK, int LWORK)       
{
   return -1;
}

//Declarations of the lapack calls
extern "C"
{
   void sgesvd_( char *JOBUp, char *JOBVTp, int *Mp, int *Np, float * A, int *LDAp, float * S, float *U, int *LDUp, 
                float * VT, int *LDVTp, float * WORK, int *LWORKp, int *INFOp);

   void dgesvd_( char *JOBUp, char *JOBVTp, int *Mp, int *Np, double * A, int *LDAp, double * S, double *U, int *LDUp, 
                double * VT, int *LDVTp, double * WORK, int *LWORKp, int *INFOp);
}

//float specialization of gesvd
template<>
inline
int gesvd<float>( char JOBU, char JOBVT, int M, int N, float * A, int LDA, float * S, float *U, int LDU, 
                float * VT, int LDVT, float * WORK, int LWORK)       
{
   int INFO;
   
   sgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU,VT, &LDVT, WORK, &LWORK, &INFO);
   
   return INFO;
}

//double specialization of gesvd
template<>
inline
int gesvd<double>( char JOBU, char JOBVT, int M, int N, double * A, int LDA, double * S, double *U, int LDU, 
                double * VT, int LDVT, double * WORK, int LWORK)       
{
   int INFO;
   
   dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU,VT, &LDVT, WORK, &LWORK, &INFO);
   
   return INFO;
}

// 

#endif //__templateLapack_hpp__
