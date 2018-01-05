/** \file roots.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines functions for finding roots
  * \ingroup gen_math_files
  * 
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __roots_hpp__
#define __roots_hpp__

#include <vector>
#include <complex>
#include <cmath>

namespace mx
{


///Find the roots of the general quartic equation
/** Finds the roots of
  * \f[
    f(x) = a x^4 + b x^3 + c x^2 + d x + e
    \f]
  * using the general formula for quartic roots. See https://en.wikipedia.org/wiki/Quartic_function.
  *
  * \param [out] x is resized to length 4, and on exit contains the 4 roots.
  * \param [in] a is the coefficient of the \f$x^4\f$ term. 
  * \param [in] b is the coefficient of the \f$x^3\f$ term.
  * \param [in] c is the coefficient of the \f$x^2\f$ term.
  * \param [in] d is the coefficient of the \f$x^1\f$ term.
  * \param [in] e is the coefficient of the \f$x^0\f$ term.
  *
  * \tparam realT is the floating point type used for calculations.
  * 
  * \ingroup gen_math
  */
template<typename realT> 
void quarticRoots( std::vector<std::complex<realT> > & x, 
                   realT a, 
                   realT b, 
                   realT c, 
                   realT d, 
                   realT e )
{
   std::complex<realT> p, q, S, Q, Delta0, Delta1;
   
   p = (8.0*a*c - 3.0*b*b)/(8.0*a*a);
   
   q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d)/(8.0*a*a*a);
   
   Delta0 = c*c - 3.0*b*d + 12.0*a*e;
   Delta1 = 2.0*c*c*c - 9.0*b*c*d + 27.0*b*b*e + 27.0*a*d*d - 72.0*a*c*e;
   
   Q = pow( static_cast<realT>(0.5)* (Delta1 + sqrt( Delta1*Delta1 - static_cast<realT>(4)*Delta0*Delta0*Delta0)), 1./3.);
   
   S = static_cast<realT>(0.5)*sqrt( - static_cast<realT>(2.0/3.0)*p + static_cast<realT>(1.0/(3.0*a))*(Q + Delta0/Q));
   
   x.resize(4);
   
   x[0] = -b/(4*a) - S + static_cast<realT>(0.5)*sqrt( static_cast<realT>(-4)*S*S - static_cast<realT>(2)*p + q/S);
   x[1] = -b/(4*a) - S - static_cast<realT>(0.5)*sqrt( static_cast<realT>(-4)*S*S - static_cast<realT>(2)*p + q/S);
   x[2] = -b/(4*a) + S + static_cast<realT>(0.5)*sqrt( static_cast<realT>(-4)*S*S - static_cast<realT>(2)*p - q/S);
   x[3] = -b/(4*a) + S - static_cast<realT>(0.5)*sqrt( static_cast<realT>(-4)*S*S - static_cast<realT>(2)*p - q/S);
   
} //quarticRoots

} //namespace mx

#endif //__roots_hpp__
