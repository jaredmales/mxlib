/** \file templateLapack.hpp
  * \brief Declares and defines templatized wrappers for the Lapack library
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef math_templateLapack_hpp
#define math_templateLapack_hpp

#include <complex>

extern "C"
{
#ifdef MXLIB_MKL

   #define MKL_Complex8 std::complex<float>
   #define MKL_Complex16 std::complex<double>
        

   #include <mkl.h>

#else

   #include <lapack.h>

#endif
}


//MKL can use 64-bit integer types, standard BLAS and LAPACK use int
//MKL declares some pointer args const.  Compiling with ATLAS doesn't seem to care, but just to be safe...
#ifdef MXLIB_MKL
   typedef MKL_INT MXLAPACK_INT;
   #define MKL_CONST_PTR const
#else
   typedef int MXLAPACK_INT;
   #define MKL_CONST_PTR
#endif

//typedef int MXLAPACK_INT;

namespace mx
{
namespace math
{


   
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
  * 
  * \ingroup template_lapack
  */
template<typename dataT>
dataT lamch (char CMACH)
{
   return -1;
}
 
/*
extern "C"
{
   extern  float   slamch_ (MKL_CONST_PTR char *CMACHp);
   extern  double  dlamch_ (MKL_CONST_PTR char *CMACHp);
}*/

// Float specialization of lamch, a wrapper for Lapack SLAMCH
template<>
float lamch<float>(char CMACH);

// Double specialization of lamch, a wrapper for Lapack DLAMCH
template<>
double lamch<double>(char CMACH);

/// Compute the Cholesky factorization of a real symmetric positive definite matrix A.
/**
  * The factorization has the form
  *  A = U**T * U,  if UPLO = 'U', or
  *  A = L  * L**T,  if UPLO = 'L',
  * where U is an upper triangular matrix and L is lower triangular.
  */ 
template<typename dataT>
MXLAPACK_INT potrf( char UPLO,          ///< [in] 'U' if upper triangle of A is stored, 'L' if lower triangle of A is stored.
                    MXLAPACK_INT N,     ///< [in] The order of the matrix A,  \>= 0.
                    dataT * A,          ///< [in/out] Symmetric matrix of dimension (LDA,N), stored as specified in UPLO.  Note that the opposite half is not referenced.
                    MXLAPACK_INT LDA,   ///< [in] The leading dimension of A.
                    MXLAPACK_INT & INFO ///< [out] 0 on success, \< 0 -INFO means the i-th argument had an illegal value, \>0 the leading minor of order INFO is not positive definite, and the factorization could not be completed.
                  );

template<>
MXLAPACK_INT potrf<float> ( char UPLO, MXLAPACK_INT N, float * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO );

template<>
MXLAPACK_INT potrf<double> ( char UPLO, MXLAPACK_INT N, double * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO );

template<>
MXLAPACK_INT potrf<std::complex<float>> ( char UPLO, MXLAPACK_INT N, std::complex<float> * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO );

template<>
MXLAPACK_INT potrf<std::complex<double>> ( char UPLO, MXLAPACK_INT N, std::complex<double> * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO );

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
  * 
  * \ingroup template_lapack
  */
template<typename dataT>
MXLAPACK_INT sytrd( char UPLO, 
                    MXLAPACK_INT N, 
                    dataT * A, 
                    MXLAPACK_INT LDA, 
                    dataT *D, 
                    dataT *E, 
                    dataT *TAU, 
                    dataT *WORK, 
                    MXLAPACK_INT LWORK, 
                    MXLAPACK_INT INFO
                  )
{
   return -1;
}

/*
//Declarations of the actual LAPACK functions
extern "C"
{
   void ssytrd_( MKL_CONST_PTR char * UPLO, MKL_CONST_PTR MXLAPACK_INT * N, float * A, MKL_CONST_PTR MXLAPACK_INT * LDA, float *D, float *E, float* TAU, float *WORK, MKL_CONST_PTR MXLAPACK_INT *LWORK, MXLAPACK_INT *INFO );
   void dsytrd_( MKL_CONST_PTR char * UPLO, MKL_CONST_PTR MXLAPACK_INT * N, double * A, MKL_CONST_PTR MXLAPACK_INT * LDA, double *D, double *E, double* TAU, double *WORK, MKL_CONST_PTR MXLAPACK_INT *LWORK, MXLAPACK_INT *INFO );
}*/

template<>
MXLAPACK_INT sytrd<float>( char UPLO, MXLAPACK_INT N, float * A, MXLAPACK_INT LDA, float *D, float *E, float *TAU, float *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT INFO);

template<>
MXLAPACK_INT sytrd<double>( char UPLO, MXLAPACK_INT N, double * A, MXLAPACK_INT LDA, double *D, double *E, double *TAU, double *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT INFO);

/// Compute selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix
/** xSYEVR computes selected eigenvalues and, optionally, eigenvectors
  * of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
  * selected by specifying either a range of values or a range of
  * indices for the desired eigenvalues.
  * 
  * See more details: http://www.netlib.org/lapack/lapack-3.1.1/html/ssyevr.f.html.
  * 
  * \ingroup template_lapack
  */
template<typename dataT>
MXLAPACK_INT syevr (  char JOBZ, char RANGE, char UPLO, MXLAPACK_INT N, dataT *A, MXLAPACK_INT LDA, dataT VL, dataT VU, MXLAPACK_INT IL, 
             MXLAPACK_INT IU, dataT ABSTOL, MXLAPACK_INT *M, dataT *W, dataT *Z, MXLAPACK_INT LDZ, MXLAPACK_INT *ISUPPZ, dataT *WORK, MXLAPACK_INT
             LWORK, MXLAPACK_INT *IWORK, MXLAPACK_INT LIWORK ) 
{
   return -1;
}
   
/*
extern "C"
{
   void  ssyevr_ (MKL_CONST_PTR char *JOBZp, MKL_CONST_PTR char *RANGEp, MKL_CONST_PTR char *UPLOp, MKL_CONST_PTR MXLAPACK_INT *Np,
                     float *A, MKL_CONST_PTR MXLAPACK_INT *LDAp, MKL_CONST_PTR float *VLp, MKL_CONST_PTR float *VUp,
                       MKL_CONST_PTR MXLAPACK_INT *ILp, MKL_CONST_PTR MXLAPACK_INT *IUp, MKL_CONST_PTR float *ABSTOLp, MXLAPACK_INT *Mp,
                         float *W, float *Z, MKL_CONST_PTR MXLAPACK_INT *LDZp, MXLAPACK_INT *ISUPPZ,
                            float *WORK, MKL_CONST_PTR MXLAPACK_INT *LWORKp, MXLAPACK_INT *IWORK, MKL_CONST_PTR MXLAPACK_INT *LIWORKp,
                              MXLAPACK_INT *INFOp);
   
   
   void  dsyevr_ (MKL_CONST_PTR char *JOBZp, MKL_CONST_PTR char *RANGEp, MKL_CONST_PTR char *UPLOp, MKL_CONST_PTR MXLAPACK_INT *Np,
                     double *A, MKL_CONST_PTR MXLAPACK_INT *LDAp, MKL_CONST_PTR double *VLp, MKL_CONST_PTR double *VUp,
                        MKL_CONST_PTR MXLAPACK_INT *ILp, MKL_CONST_PTR MXLAPACK_INT *IUp, MKL_CONST_PTR double *ABSTOLp, MXLAPACK_INT *Mp,
                               double *W, double *Z, MKL_CONST_PTR MXLAPACK_INT *LDZp, MXLAPACK_INT *ISUPPZ,
                                double *WORK, MKL_CONST_PTR MXLAPACK_INT *LWORKp, MXLAPACK_INT *IWORK, MKL_CONST_PTR MXLAPACK_INT *LIWORKp,
                                 MXLAPACK_INT *INFOp);
}
*/

// Float specialization of syevr, a wrapper for Lapack SSYEVR
template<>
MXLAPACK_INT syevr<float> ( char JOBZ, char RANGE, char UPLO, MXLAPACK_INT N, float *A, MXLAPACK_INT LDA, float VL, float VU,
                    MXLAPACK_INT IL, MXLAPACK_INT IU,  float ABSTOL, MXLAPACK_INT *M, float *W, float *Z, MXLAPACK_INT LDZ, MXLAPACK_INT *ISUPPZ,
                     float *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT *IWORK, MXLAPACK_INT LIWORK );

// Double specialization of syevr, a wrapper for Lapack DSYEVR
template<>
MXLAPACK_INT syevr<double> ( char JOBZ, char RANGE, char UPLO, MXLAPACK_INT N, double *A, MXLAPACK_INT LDA, double VL, double VU,
                    MXLAPACK_INT IL, MXLAPACK_INT IU,  double ABSTOL, MXLAPACK_INT *M, double *W, double *Z, MXLAPACK_INT LDZ, MXLAPACK_INT *ISUPPZ,
                     double *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT *IWORK, MXLAPACK_INT LIWORK );

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
    (MXLAPACK_INT)
    The number of rows of the input matrix A.  M >= 0.
  *
  * \param[in] N
  * \parblock
    (MXLAPACK_INT)
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
    (MXLAPACK_INT)
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
     (MXLAPACK_INT)
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
     (MXLAPACK_INT)
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
     (MXLAPACK_INT)
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
           superdiagonals of an MXLAPACK_INTermediate bidiagonal form B
           did not converge to zero. See the description of WORK
          above for details.
    \endparblock
  *
  * \ingroup template_lapack
  */
template<typename dataT>
MXLAPACK_INT gesvd( char JOBU, char JOBVT, MXLAPACK_INT M, MXLAPACK_INT N, dataT * A, MXLAPACK_INT LDA, dataT * S, dataT *U, MXLAPACK_INT LDU, 
                dataT * VT, MXLAPACK_INT LDVT, dataT * WORK, MXLAPACK_INT LWORK)       
{
   return -1;
}

/*
//Declarations of the lapack calls
extern "C"
{
   void sgesvd_( MKL_CONST_PTR char *JOBUp, MKL_CONST_PTR char *JOBVTp, MKL_CONST_PTR MXLAPACK_INT *Mp, MKL_CONST_PTR MXLAPACK_INT *Np, 
                  float * A, MKL_CONST_PTR MXLAPACK_INT *LDAp, float * S, float *U, MKL_CONST_PTR MXLAPACK_INT *LDUp, 
                    float * VT, MKL_CONST_PTR MXLAPACK_INT *LDVTp, float * WORK, MKL_CONST_PTR MXLAPACK_INT *LWORKp, MXLAPACK_INT *INFOp);

   void dgesvd_( MKL_CONST_PTR char *JOBUp, MKL_CONST_PTR char *JOBVTp, MKL_CONST_PTR MXLAPACK_INT *Mp, MKL_CONST_PTR MXLAPACK_INT *Np, 
                  double * A, MKL_CONST_PTR MXLAPACK_INT *LDAp, double * S, double *U, MKL_CONST_PTR MXLAPACK_INT *LDUp, 
                    double * VT, MKL_CONST_PTR MXLAPACK_INT *LDVTp, double * WORK, MKL_CONST_PTR MXLAPACK_INT *LWORKp, MXLAPACK_INT *INFOp);
}*/

//float specialization of gesvd
template<>
MXLAPACK_INT gesvd<float>( char JOBU, char JOBVT, MXLAPACK_INT M, MXLAPACK_INT N, float * A, MXLAPACK_INT LDA, float * S, float *U, MXLAPACK_INT LDU, 
                float * VT, MXLAPACK_INT LDVT, float * WORK, MXLAPACK_INT LWORK);

//double specialization of gesvd
template<>
MXLAPACK_INT gesvd<double>( char JOBU, char JOBVT, MXLAPACK_INT M, MXLAPACK_INT N, double * A, MXLAPACK_INT LDA, double * S, double *U, MXLAPACK_INT LDU, 
                double * VT, MXLAPACK_INT LDVT, double * WORK, MXLAPACK_INT LWORK);

/// Compute the singular value decomposition (SVD) of a real matrix with GESDD
/** 
 This documentation copied shamelessly from the LAPACK source at <a href="http://www.netlib.org/lapack/explore-html/d4/dca/group__real_g_esing.html#gac2cd4f1079370ac908186d77efcd5ea8">netlib</a>.
  
 SGESDD computes the singular value decomposition (SVD) of a real
 M-by-N matrix A, optionally computing the left and right singular
 vectors.  If singular vectors are desired, it uses a
 divide-and-conquer algorithm.

 The SVD is written

\f[ A = U \Sigma V^T \f]

 where \f$ \Sigma \f$ is an M-by-N matrix which is zero except for its
 min(m,n) diagonal elements, \f$ U \f$ is an M-by-M orthogonal matrix, and
 \f$ V  is an N-by-N orthogonal matrix.  The diagonal elements of \f$ \Sigma \f$
 are the singular values of A; they are real and non-negative, and
 are returned in descending order.  The first min(m,n) columns of
 U and V are the left and right singular vectors of A.

 Note that the routine returns \f$ V^T \f$, not \f$ V \f$.

 The divide and conquer algorithm makes very mild assumptions about
 floating point arithmetic. It will work on machines with a guard
 digit in add/subtract, or on those binary machines without guard
 digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 Cray-2. It could conceivably fail on hexadecimal or decimal machines
 without guard digits, but the authors know of none.
 \endparblock

 \param[in] JOBZ
 \parblock
          JOBZ is CHARACTER*1
          Specifies options for computing all or part of the matrix U:
          = 'A':  all M columns of U and all N rows of V**T are
                  returned in the arrays U and VT;
          = 'S':  the first min(M,N) columns of U and the first
                  min(M,N) rows of V**T are returned in the arrays U
                  and VT;
          = 'O':  If M >= N, the first N columns of U are overwritten
                  on the array A and all rows of V**T are returned in
                  the array VT;
                  otherwise, all columns of U are returned in the
                  array U and the first M rows of V**T are overwritten
                  in the array A;
          = 'N':  no columns of U or rows of V**T are computed.
 \param[in] M
 \parblock
          M is INTEGER
          The number of rows of the input matrix A.  M >= 0.
 \param[in] N
 \parblock
          N is INTEGER
          The number of columns of the input matrix A.  N >= 0.
 \param[in,out] A
 \parblock
          A is REAL array, dimension (LDA,N)
          On entry, the M-by-N matrix A.
          On exit,
          if JOBZ = 'O',  A is overwritten with the first N columns
                          of U (the left singular vectors, stored
                          columnwise) if M >= N;
                          A is overwritten with the first M rows
                          of V**T (the right singular vectors, stored
                          rowwise) otherwise.
          if JOBZ .ne. 'O', the contents of A are destroyed.
 \param[in] LDA
 \parblock
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 

 \param[out] S
 \parblock
          S is REAL array, dimension (min(M,N))
          The singular values of A, sorted so that S(i) >= S(i+1).

 \param[out] U
 \parblock
          U is REAL array, dimension (LDU,UCOL)
          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
          UCOL = min(M,N) if JOBZ = 'S'.
          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
          orthogonal matrix U;
          if JOBZ = 'S', U contains the first min(M,N) columns of U
          (the left singular vectors, stored columnwise);
          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.

 \param[in] LDU
 \parblock
          LDU is INTEGER
          The leading dimension of the array U.  LDU >= 1; if
          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.

 \param[out] VT
 \parblock
          VT is REAL array, dimension (LDVT,N)
          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
          N-by-N orthogonal matrix V**T;
          if JOBZ = 'S', VT contains the first min(M,N) rows of
          V**T (the right singular vectors, stored rowwise);
          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.

 \param[in] LDVT
 \parblock
          LDVT is INTEGER
          The leading dimension of the array VT.  LDVT >= 1; if
          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
          if JOBZ = 'S', LDVT >= min(M,N).

 \param[out] WORK
 \parblock
          WORK is REAL array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;

 \param[in] LWORK
 \parblock
          LWORK is INTEGER
          The dimension of the array WORK. LWORK >= 1.
          If JOBZ = 'N',
            LWORK >= 3*min(M,N) + max(max(M,N),6*min(M,N)).
          If JOBZ = 'O',
            LWORK >= 3*min(M,N) + 
                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
          If JOBZ = 'S' or 'A'
            LWORK >= min(M,N)*(7+4*min(M,N))
          For good performance, LWORK should generally be larger.
          If LWORK = -1 but other input arguments are legal, WORK(1)
          returns the optimal LWORK.

 \param[out] IWORK
 \parblock
          IWORK is INTEGER array, dimension (8*min(M,N))
 \endparblock

 \returns
 \parblock
          = 0:  successful exit.<br />
          < 0:  if INFO = -i, the i-th argument had an illegal value.<br />
          > 0:  SBDSDC did not converge, updating process failed.
 \endparblock
 
\ingroup template_lapack

*/
template<typename dataT>
MXLAPACK_INT gesdd(char JOBZ, MXLAPACK_INT M, MXLAPACK_INT N, dataT *A, MXLAPACK_INT LDA, dataT *S, dataT * U, MXLAPACK_INT LDU, dataT * VT, MXLAPACK_INT LDVT, dataT *WORK, MXLAPACK_INT  LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT INFO)
{
   return -1;
}

/*
//Declarations of the lapack calls
extern "C"
{
   void sgesdd_(MKL_CONST_PTR char * JOBZ, MKL_CONST_PTR MXLAPACK_INT *M, MKL_CONST_PTR MXLAPACK_INT *N, float *A, MKL_CONST_PTR MXLAPACK_INT *LDA, float *S, float * U, MKL_CONST_PTR MXLAPACK_INT *LDU, float * VT, MKL_CONST_PTR MXLAPACK_INT *LDVT, float *WORK, MKL_CONST_PTR MXLAPACK_INT * LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT *INFO);
   
   void dgesdd_(MKL_CONST_PTR char * JOBZ, MKL_CONST_PTR MXLAPACK_INT *M, MKL_CONST_PTR MXLAPACK_INT *N, double *A, MKL_CONST_PTR MXLAPACK_INT *LDA, double *S, double * U, MKL_CONST_PTR MXLAPACK_INT *LDU, double * VT, MKL_CONST_PTR MXLAPACK_INT *LDVT, double *WORK, MKL_CONST_PTR MXLAPACK_INT * LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT *INFO);
   
}*/

//float specialization of gesdd
template<>
MXLAPACK_INT gesdd<float>(char JOBZ, MXLAPACK_INT M, MXLAPACK_INT N, float *A, MXLAPACK_INT LDA, float *S, float * U, MXLAPACK_INT LDU, float * VT, MXLAPACK_INT LDVT, float *WORK, MXLAPACK_INT  LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT INFO);

//double specialization of gesdd
template<>
MXLAPACK_INT gesdd<double>(char JOBZ, MXLAPACK_INT M, MXLAPACK_INT N, double *A, MXLAPACK_INT LDA, double *S, double * U, MXLAPACK_INT LDU, double * VT, MXLAPACK_INT LDVT, double *WORK, MXLAPACK_INT  LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT INFO);

} //namespace math
} //namespace mx

// 

#endif //math_templateLapack_hpp
