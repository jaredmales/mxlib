/** \file templateCurand.cpp
 * \author Jared R. Males
 * \brief Implementation of a template interface to curand
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

#include "math/cuda/templateCurand.hpp"

namespace mx
{
namespace cuda
{

template <>
curandStatus_t
curandGenerateNormal<float>( curandGenerator_t generator, float *outputPtr, size_t n, float mean, float stddev )
{
    return ::curandGenerateNormal( generator, outputPtr, n, mean, stddev );
}

template <>
curandStatus_t
curandGenerateNormal<double>( curandGenerator_t generator, double *outputPtr, size_t n, double mean, double stddev )
{
    return ::curandGenerateNormalDouble( generator, outputPtr, n, mean, stddev );
}

} // namespace cuda
} // namespace mx
