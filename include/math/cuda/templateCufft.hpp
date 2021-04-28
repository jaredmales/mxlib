/** \file templateCufft.hpp
  * \author Jared R. Males
  * \brief A template interface to cufft
  * \ingroup cuda_files
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

#ifndef math_templateCufft_hpp
#define math_templateCufft_hpp

#include <complex>
#include <cuda_runtime.h>
#include <cufft.h>


namespace mx
{
namespace cuda
{

template<typename inputT, typename outputT>
cufftResult cufftPlan2d( cufftHandle *plan, 
                         int nx, 
                         int ny
                       );

template<>
cufftResult cufftPlan2d<std::complex<float>, std::complex<float>>( cufftHandle *plan, 
                                                                  int nx, 
                                                                  int ny
                                                                );

template<>
cufftResult cufftPlan2d<cuComplex, cuComplex>( cufftHandle *plan, 
                                               int nx, 
                                               int ny
                                             );

template<>
cufftResult cufftPlan2d<std::complex<double>, std::complex<double>>( cufftHandle *plan, 
                                                                    int nx, 
                                                                    int ny
                                                                  );

template<>
cufftResult cufftPlan2d<cuDoubleComplex, cuDoubleComplex>( cufftHandle *plan, 
                                                           int nx, 
                                                           int ny
                                                         );

template<typename inputT, typename outputT>
cufftResult cufftExec( cufftHandle plan, 
                       inputT *idata, 
                       inputT *odata, 
                       int direction
                     );

template<>
cufftResult cufftExec<std::complex<float>, std::complex<float>>( cufftHandle plan, 
                                                                 std::complex<float> *idata, 
                                                                 std::complex<float> *odata, 
                                                                 int direction
                                                               );

template<>
cufftResult cufftExec<cuComplex, cuComplex>( cufftHandle plan, 
                                             cuComplex *idata, 
                                             cuComplex *odata, 
                                             int direction
                                           );

template<>
cufftResult cufftExec<std::complex<double>, std::complex<double>>( cufftHandle plan, 
                                                                   std::complex<double> *idata, 
                                                                   std::complex<double> *odata, 
                                                                   int direction
                                                                 );

template<>
cufftResult cufftExec<cuDoubleComplex, cuDoubleComplex>( cufftHandle plan, 
                                                         cuDoubleComplex *idata, 
                                                         cuDoubleComplex *odata, 
                                                         int direction
                                                       );

}//namespace cuda 
}//namespace mx

#endif // math_templateCufft_hpp
