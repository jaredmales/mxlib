/** \file jinc.cpp
  * \brief Instantiations of the Jinc and Jinc2 functions
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#include "math/func/jinc.hpp"

namespace mx
{
namespace math
{
namespace func 
{


template
float jinc<float>(const float & x);

template
double jinc<double>(const double & x);

template
long double jinc<long double>(const long double & x);

#ifdef HASQUAD
template
__float128 jinc<__float128>(const __float128 & x);
#endif


template
float jinc2<float>(const float & x);

template
double jinc2<double>(const double & x);

template
long double jinc2<long double>(const long double & x);

#ifdef HASQUAD
template
__float128 jinc2<__float128>(const __float128 & x);
#endif

} //namespace func 
} //namespace math
} //namespace mx

