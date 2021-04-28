/** \file templateCurand.hpp
  * \author Jared R. Males
  * \brief A template interface to curand
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

#ifndef math_templateCurand_hpp
#define math_templateCurand_hpp


#include <cuda_runtime.h>
#include <curand.h>


namespace mx
{
namespace cuda
{


template<typename realT>
curandStatus_t curandGenerateNormal( curandGenerator_t generator, 
                                     realT *outputPtr, 
                                     size_t n, 
                                     realT mean, 
                                     realT stddev
                                   );

template<>
curandStatus_t curandGenerateNormal<float>( curandGenerator_t generator, 
                                            float *outputPtr, 
                                            size_t n, 
                                            float mean, 
                                            float stddev
                                          );

template<>
curandStatus_t curandGenerateNormal<double>( curandGenerator_t generator, 
                                             double *outputPtr, 
                                             size_t n, 
                                             double mean, 
                                             double stddev
                                           );


}//namespace cuda 
}//namespace mx

#endif //math_templateCurand_hpp
