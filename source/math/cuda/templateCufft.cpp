/** \file templateCufft.cpp
  * \author Jared R. Males
  * \brief Implementation of a template interface to cufft
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2019,2020 Jared R. Males (jaredmales@gmail.com)
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

#include "math/cuda/templateCufft.hpp"

namespace mx
{
namespace cuda
{

template<>
cufftResult cufftPlan2d<std::complex<float>, std::complex<float>>( cufftHandle *plan, 
                                                                  int nx, 
                                                                  int ny
                                                                )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_C2C);
}

template<>
cufftResult cufftPlan2d<cuComplex, cuComplex>( cufftHandle *plan, 
                                               int nx, 
                                               int ny
                                             )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_C2C);
}

template<>
cufftResult cufftPlan2d<std::complex<double>, std::complex<double>>( cufftHandle *plan, 
                                                                    int nx, 
                                                                    int ny
                                                                  )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_Z2Z);
}

template<>
cufftResult cufftPlan2d<cuDoubleComplex, cuDoubleComplex>( cufftHandle *plan, 
                                                           int nx, 
                                                           int ny
                                                         )
{
   return ::cufftPlan2d(plan, nx, ny, CUFFT_Z2Z);
}

template<>
cufftResult cufftExec<std::complex<float>, std::complex<float>>( cufftHandle plan, 
                                                                 std::complex<float> *idata, 
                                                                 std::complex<float> *odata, 
                                                                 int direction
                                                               )
{
   return ::cufftExecC2C(plan, (cuComplex *) idata, (cuComplex *) odata, direction);
}

template<>
cufftResult cufftExec<cuComplex, cuComplex>( cufftHandle plan, 
                                             cuComplex *idata, 
                                             cuComplex *odata, 
                                             int direction
                                           )
{
   return ::cufftExecC2C(plan, idata, odata, direction);
}

template<>
cufftResult cufftExec<std::complex<double>, std::complex<double>>( cufftHandle plan, 
                                                                   std::complex<double> *idata, 
                                                                   std::complex<double> *odata, 
                                                                   int direction
                                                                 )
{
   return ::cufftExecZ2Z(plan, (cuDoubleComplex *) idata, (cuDoubleComplex *) odata, direction);
}

template<>
cufftResult cufftExec<cuDoubleComplex, cuDoubleComplex>( cufftHandle plan, 
                                                         cuDoubleComplex *idata, 
                                                         cuDoubleComplex *odata, 
                                                         int direction
                                                       )
{
   return ::cufftExecZ2Z(plan, idata, odata, direction);
}


}//namespace cuda 
}//namespace mx

