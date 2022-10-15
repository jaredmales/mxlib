/** \file templateLapack.cpp
  * \brief Implementation of templatized wrappers for the Lapack library
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "math/templateLapack.hpp"


namespace mx
{
namespace math
{

      
template<>
float lamch<float>(char CMACH)
{ 
   return  slamch_ (&CMACH
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
}

// Double specialization of lamch, a wrapper for Lapack DLAMCH
template<>
double lamch<double>(char CMACH)
{ 
   return  dlamch_ (&CMACH
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
}

template<>
MXLAPACK_INT potrf<float> ( char UPLO, MXLAPACK_INT N, float * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO )
{
   spotrf_(&UPLO, &N, A, &LDA, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

template<>
MXLAPACK_INT potrf<double> ( char UPLO, MXLAPACK_INT N, double * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO )
{
   dpotrf_(&UPLO, &N, A, &LDA, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

template<>
MXLAPACK_INT potrf<std::complex<float>> ( char UPLO, MXLAPACK_INT N, std::complex<float> * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO )
{
   cpotrf_(&UPLO, &N,
   (float _Complex*)A,
   &LDA, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

template<>
MXLAPACK_INT potrf<std::complex<double>> ( char UPLO, MXLAPACK_INT N, std::complex<double> * A, MXLAPACK_INT LDA, MXLAPACK_INT &INFO )
{
   zpotrf_(&UPLO, &N,
   (double _Complex*)A,
   &LDA, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

template<>
MXLAPACK_INT sytrd<float>( char UPLO, MXLAPACK_INT N, float * A, MXLAPACK_INT LDA, float *D, float *E, float *TAU, float *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT INFO)
{
  
   ssytrd_(&UPLO, &N, A, &LDA, D, E, TAU, WORK, &LWORK, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

template<>
MXLAPACK_INT sytrd<double>( char UPLO, MXLAPACK_INT N, double * A, MXLAPACK_INT LDA, double *D, double *E, double *TAU, double *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT INFO)
{
  
   dsytrd_(&UPLO, &N, A, &LDA, D, E, TAU, WORK, &LWORK, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}
      
// Float specialization of syevr, a wrapper for Lapack SSYEVR
template<>
MXLAPACK_INT syevr<float> ( char JOBZ, char RANGE, char UPLO, MXLAPACK_INT N, float *A, MXLAPACK_INT LDA, float VL, float VU,
                    MXLAPACK_INT IL, MXLAPACK_INT IU,  float ABSTOL, MXLAPACK_INT *M, float *W, float *Z, MXLAPACK_INT LDZ, MXLAPACK_INT *ISUPPZ,
                     float *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT *IWORK, MXLAPACK_INT LIWORK ) 
{

   MXLAPACK_INT  INFO;
   
   ssyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
           &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
           WORK, &LWORK, IWORK, &LIWORK, &INFO
           #ifdef LAPACK_FORTRAN_STRLEN_END
           , 1, 1, 1
           #endif
           );

   return  INFO;
}

// Double specialization of syevr, a wrapper for Lapack DSYEVR
template<>
MXLAPACK_INT syevr<double> ( char JOBZ, char RANGE, char UPLO, MXLAPACK_INT N, double *A, MXLAPACK_INT LDA, double VL, double VU,
                    MXLAPACK_INT IL, MXLAPACK_INT IU,  double ABSTOL, MXLAPACK_INT *M, double *W, double *Z, MXLAPACK_INT LDZ, MXLAPACK_INT *ISUPPZ,
                     double *WORK, MXLAPACK_INT LWORK, MXLAPACK_INT *IWORK, MXLAPACK_INT LIWORK ) 
{

   MXLAPACK_INT  INFO;
   
   dsyevr_ (&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
           &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
           WORK, &LWORK, IWORK, &LIWORK, &INFO
           #ifdef LAPACK_FORTRAN_STRLEN_END
           , 1, 1, 1
           #endif
   );

   return  INFO;
}

//float specialization of gesvd
template<>
MXLAPACK_INT gesvd<float>( char JOBU, char JOBVT, MXLAPACK_INT M, MXLAPACK_INT N, float * A, MXLAPACK_INT LDA, float * S, float *U, MXLAPACK_INT LDU, 
                float * VT, MXLAPACK_INT LDVT, float * WORK, MXLAPACK_INT LWORK)       
{
   MXLAPACK_INT INFO;
   
   sgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU,VT, &LDVT, WORK, &LWORK, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1, 1
   #endif
   );
   
   return INFO;
}

//double specialization of gesvd
template<>
MXLAPACK_INT gesvd<double>( char JOBU, char JOBVT, MXLAPACK_INT M, MXLAPACK_INT N, double * A, MXLAPACK_INT LDA, double * S, double *U, MXLAPACK_INT LDU, 
                double * VT, MXLAPACK_INT LDVT, double * WORK, MXLAPACK_INT LWORK)       
{
   MXLAPACK_INT INFO;
   
   dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU,VT, &LDVT, WORK, &LWORK, &INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1, 1
   #endif
   );
   
   return INFO;
}

//float specialization of gesdd
template<>
MXLAPACK_INT gesdd<float>(char JOBZ, MXLAPACK_INT M, MXLAPACK_INT N, float *A, MXLAPACK_INT LDA, float *S, float * U, MXLAPACK_INT LDU, float * VT, MXLAPACK_INT LDVT, float *WORK, MXLAPACK_INT  LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT INFO)
{
   sgesdd_(&JOBZ,&M,&N,A,&LDA,S,U,&LDU,VT,&LDVT,WORK,&LWORK,IWORK,&INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

//double specialization of gesdd
template<>
MXLAPACK_INT gesdd<double>(char JOBZ, MXLAPACK_INT M, MXLAPACK_INT N, double *A, MXLAPACK_INT LDA, double *S, double * U, MXLAPACK_INT LDU, double * VT, MXLAPACK_INT LDVT, double *WORK, MXLAPACK_INT  LWORK, MXLAPACK_INT * IWORK, MXLAPACK_INT INFO)
{
   dgesdd_(&JOBZ,&M,&N,A,&LDA,S,U,&LDU,VT,&LDVT,WORK,&LWORK,IWORK,&INFO
   #ifdef LAPACK_FORTRAN_STRLEN_END
   , 1
   #endif
   );
   
   return INFO;
}

} //namespace math
} //namespace mx


