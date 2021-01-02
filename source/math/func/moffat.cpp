/** \file moffat.cpp
  * \author Jared R. Males
  * \brief Implementation utilities related to the Moffat function.
  * \ingroup gen_math_files
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

#include "math/func/moffat.hpp"

namespace mx
{
namespace math 
{
namespace func 
{
   

template 
float moffat<float>(const float x, const float I0, const float Ipk, const float x0, const float alpha, const float beta);

template 
double moffat<double>(const double x, const double I0, const double Ipk, const double x0, const double alpha, const double beta);

template 
long double moffat<long double>(const long double x, const long double I0, const long double Ipk, const long double x0, const long double alpha, const long double beta);

#ifdef HASQUAD
template 
__float128 moffat<__float128>(const __float128 x, const __float128 I0, const __float128 Ipk, const __float128 x0, const __float128 alpha, const __float128 beta);
#endif


template 
float moffat2D<float>(const float x, const float y, const float I0, const float Ipk, const float x0, const float y0, const float alpha, const float beta);

template 
double moffat2D<double>(const double x, const double y, const double I0, const double Ipk, const double x0, const double y0, const double alpha, const double beta);

template 
long double moffat2D<long double>(const long double x, const long double y, const long double I0, const long double Ipk, const long double x0, const long double y0, const long double alpha, const long double beta);

#ifdef HASQUAD
template 
__float128 moffat2D<__float128>(const __float128 x, const __float128 y, const __float128 I0, const __float128 Ipk, const __float128 x0, const __float128 y0, const __float128 alpha, const __float128 beta);
#endif


template
float moffatFWHM(float alpha, float beta);

template
double moffatFWHM(double alpha, double beta);

template
long double moffatFWHM(long double alpha, long double beta);

#ifdef HASQUAD
template
__float128 moffatFWHM(__float128 alpha, __float128 beta);
#endif
} //namespace func 
} //namespace math
} //namespace mx




