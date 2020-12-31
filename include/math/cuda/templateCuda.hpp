/** \file templateCuda.hpp
  * \author Jared R. Males
  * \brief Utilities for a template interface to cuda
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

#ifndef math_templateCuda_hpp
#define math_templateCuda_hpp

#include <complex>

#include <cuda_runtime.h>
#include <cuComplex.h>

namespace mx
{
namespace cuda
{

template<typename realT>
struct complex;

template<>
struct complex<float>
{
   typedef cuComplex cudaType;
};

template<>
struct complex<double>
{
   typedef cuDoubleComplex cudaType;
};


template<typename cppType>
struct cpp2cudaType
{
   typedef cppType cudaType;
};

template<>
struct cpp2cudaType<std::complex<float>>
{
   typedef complex<float>::cudaType cudaType;
};

template<>
struct cpp2cudaType<std::complex<double>>
{
   typedef complex<double>::cudaType cudaType;
};


}//namespace cuda 
}//namespace mx

#endif // templateCuda_hpp
