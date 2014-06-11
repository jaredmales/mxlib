
#ifndef __templateLapack_hpp__
#define __templateLapack_hpp__

/// Wrapper for Lapack xLAMCH
/** xLAMCH determines x-precision machine parameters.
  *
  * See more details at http://www.netlib.org/lapack/lapack-3.1.1/html/slamch.f.html.
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
      
/// Float specialization of lamch, a wrapper for Lapack SLAMCH
template<>
float lamch<float>(char CMACH)
{ 
   return  slamch_ (&CMACH);
}

/// Double specialization of lamch, a wrapper for Lapack DLAMCH
template<>
double lamch<double>(char CMACH)
{ 
   return  dlamch_ (&CMACH);
}


/// Wrapper for Lapack xSYEVR
/** xSYEVR computes selected eigenvalues and, optionally, eigenvectors
  * of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
  * selected by specifying either a range of values or a range of
  * indices for the desired eigenvalues.
  * 
  * See more details: http://www.netlib.org/lapack/lapack-3.1.1/html/ssyevr.f.html.
  */
template<typename dataT>
int syevr (  char JOBZ, 
             char RANGE, 
             char UPLO, 
             int N,
             dataT *A, 
             int LDA, 
             dataT VL, 
             dataT VU,
             int IL, 
             int IU, 
             dataT ABSTOL, 
             int *M,
             dataT *W, 
             dataT *Z, 
             int LDZ, 
             int *ISUPPZ,
             dataT *WORK, 
             int LWORK, 
             int *IWORK, 
             int LIWORK ) 
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

      
/// Float specialization of syevr, a wrapper for Lapack SSYEVR
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

/// Double specialization of syevr, a wrapper for Lapack DSYEVR
template<>
inline
int syevr<double> ( char JOBZ, char RANGE, char UPLO, int N, double *A, int LDA, double VL, double VU,
                    int IL, int IU,  double ABSTOL, int *M, double *W, double *Z, int LDZ, int *ISUPPZ,
                     double *WORK, int LWORK, int *IWORK, int LIWORK ) 
{

   int  INFO;
   
   ssyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
           &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
           WORK, &LWORK, IWORK, &LIWORK, &INFO);

   return  INFO;
}

// 

#endif //__templateLapack_hpp__
