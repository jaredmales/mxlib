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

#ifndef math_roots_hpp
#define math_roots_hpp

#include <vector>
#include <complex>
#include <cmath>

namespace mx
{
namespace math 
{

   /// Calculate the descriminant for the general cubic
/**
  * \returns the descriminant
  */
template<typename realT>
realT cubicDescriminant( const realT & a, ///< [in] the 3rd order coefficient of the general cubic
                         const realT & b, ///< [in] the 2nd order coefficient of the general cubic
                         const realT & c, ///< [in] the 1st order coefficient of the general cubic
                         const realT & d  ///< [in] the 0th order coefficient of the general cubic
                       )
{
   return 18*a*b*c*d
             - 4*b*b*b*d
                  + b*b*c*c
                       - 4*a*c*c*c
                           - 27*a*a*d*d;
}

/// Calculate the descriminant for the depressed cubic
/**
  * \returns the descriminant
  */
template<typename realT>
realT cubicDescriminant( const realT & p, ///< [in] the 0th order coefficient of the depressed cubic
                         const realT & q  ///< [in] the 0th order coefficient of the depressed cubic
                       )
{
   return -1*(4*p*p*p + 27*q*q);
}

/// Convert a general cubic equation to depressed form
/** The general cubic is
  * \f[
  * ax^3 + bx^2 + cx + d = 0
  * \f]
  * which can be converted to compressed form
  * \f[
  * t^3 + pt + q = 0
  * \f]
  *
  */
template<typename realT>
void cubicDepressed( realT & p,       ///< [out] the 1st order coefficient of the depressed cubic
                     realT & q,       ///< [out] the 0th order coefficient of the depressed cubic
                     const realT & a, ///< [in] the 3rd order coefficient of the general cubic
                     const realT & b, ///< [in] the 2nd order coefficient of the general cubic
                     const realT & c, ///< [in] the 1st order coefficient of the general cubic
                     const realT & d  ///< [in] the 0th order coefficient of the general cubic
                   )
{
   p = (3*a*c-b*b)/(3*a*a);
   q = (2*b*b*b-9*a*b*c+27*a*a*d)/(27*a*a*a);
}

/// Calculate the real root for a depressed cubic with negative descriminant
/**
  * \returns the real root for a depressed cubic defined by \p p and \p q
  *
  * \throws std::runtime_error if there are 3 real roots
  *
  */
template<typename realT>
realT cubicRealRoot( const realT & p, ///< [in] the 1st order coefficient of the depressed cubic
                     const realT & q  ///< [in] the 0th order coefficient of the depressed cubic
                   )
{
   realT D = q*q/4 + p*p*p/27;

   if(D < 0) throw std::runtime_error("cannot apply Cardano's formula, find roots using....");

   D = sqrt(D);

   realT u1 = -q/2 + D;
   realT u2 = -q/2 - D;

   return std::cbrt(u1) + std::cbrt(u2);
}

/// Find the roots of the general quartic equation
/** Finds the roots of
  * \f[
    f(x) = a x^4 + b x^3 + c x^2 + d x + e
    \f]
  * using the general formula for quartic roots. See https://en.wikipedia.org/wiki/Quartic_function.
  *
  * \tparam realT is the floating point type used for calculations.
  * 
  * \ingroup gen_math
  */
template<typename realT> 
void quarticRoots( std::vector<std::complex<realT> > & x, ///< [out] On exit contains the 4 roots, is resized to length 4.
                   realT a, ///< [in] the coefficient of the \f$x^4\f$ term.
                   realT b, ///< [in] the coefficient of the \f$x^3\f$ term.
                   realT c, ///< [in] the coefficient of the \f$x^2\f$ term.
                   realT d, ///< [in] the coefficient of the \f$x^1\f$ term. 
                   realT e  ///< [in] the coefficient of the \f$x^0\f$ term. 
                 )
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

} //namespace math
} //namespace mx

#endif //math_roots_hpp
