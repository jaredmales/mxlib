/** \file fftwEnvironment.cpp
 * \brief Definitions for the fftwEnvironment manager
 * \ingroup fft_files
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

#include "math/fft/fftwEnvironment.hpp"

namespace mx
{
namespace math
{
namespace fft
{

template <>
std::string fftw_typename<float>()
{
    return "float";
}

template <>
std::string fftw_typename<double>()
{
    return "double";
}

template <>
std::string fftw_typename<long double>()
{
    return "long_double";
}

#ifdef HASQUAD
template <>
std::string fftw_typename<__float128>()
{
    return "quad";
}
#endif

} // namespace fft
} // namespace math
} // namespace mx
