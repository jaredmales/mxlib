/** \file gslInterpolator.cpp
 * \brief Implementation of Class for managing 1-D interpolation using the GNU Scientific Library
 * \ingroup gen_math_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
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

#include "math/gslInterpolator.hpp"

namespace mx
{
namespace math
{

template class gslInterpolator<gsl_interp_linear<double>>;
template class gslInterpolator<gsl_interp_steffen<double>>;

} // namespace math
} // namespace mx
