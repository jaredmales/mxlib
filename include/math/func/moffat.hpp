/** \file moffat.hpp
  * \author Jared R. Males
  * \brief Declarations for utilities related to the Moffat function.
  * \ingroup gen_math_files
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

#ifndef math_func_moffat_hpp
#define math_func_moffat_hpp

#include <cmath>


namespace mx
{
namespace math 
{
namespace func 
{
   
/** \addtogroup gen_math_moffats
  * The Moffat Function\cite moffat_1969, a.k.a. the Moffat Profile, a.k.a. the Moffat Distribution, has the form
  * \f[
    I(x) = I_{pk}\left[ 1 + \frac{x^2}{\alpha^2}\right]^{-\beta}
  * \f] 
  * With \f$\beta=1\f$ it is the 
  * Lorentzian or Cauchy distribution.  See also https://en.wikipedia.org/wiki/Moffat_distribution and
  * https://en.wikipedia.org/wiki/Cauchy_distribution.
  *
  * 1-D and 2-D symmetric forms are provided.  Utilities are provided for normalizing and calculating the full-width at half-maximum.
  */
  
///Find value at position (x) of the 1D arbitrarily-centered symmetric unnormalized Moffat function
/** The Moffat distribution is due to \cite moffat_1969.  With \f$\beta=1\f$ it is the 
  * Lorentzian or Cauchy distribution.  See also https://en.wikipedia.org/wiki/Moffat_distribution and
  * https://en.wikipedia.org/wiki/Cauchy_distribution.
  * 
  * Here we use the unnormalized general form, most useful for peak fitting.
  * 
  * This function computes:
  * \f[
    I(x) = I_0 + I_{pk}\left[ 1 + \frac{(x-x_0)^2}{\alpha^2}\right]^{-\beta} 
  * \f] 
  * 
  * \returns the value of the 1D arbitrarily-centered unnormalized Moffat at (x)
  * 
  * \tparam realT is type to use for arithmetic
  * 
  * \test Scenario: compiling 1D Moffat function \ref tests_math_func_moffat1D "[test doc]"
  *
  * \ingroup gen_math_moffats
  */ 
template<typename realT>
realT moffat( const realT x,     ///< [in] is the x-position at which to evaluate the Moffat function
              const realT I0,    ///< [in] is the constant to add to the Moffat function
              const realT Ipk,   ///< [in] is the scaling factor (peak = A)
              const realT x0,    ///< [in] is the x-coordinate of the center
              const realT alpha, ///< [in] is the width parameter of the Moffat function.
              const realT beta   ///< [in] is the shape parameter of the Moffat function
            )
{ 
   return I0 + Ipk * pow( static_cast<realT>(1) + pow(x-x0,2)/pow(alpha,2), -beta);
}

extern template 
float moffat<float>(const float x, const float I0, const float Ipk, const float x0, const float alpha, const float beta);

extern template 
double moffat<double>(const double x, const double I0, const double Ipk, const double x0, const double alpha, const double beta);

extern template 
long double moffat<long double>(const long double x, const long double I0, const long double Ipk, const long double x0, const long double alpha, const long double beta);

#ifdef HASQUAD
extern template 
__float128 moffat<__float128>(const __float128 x, const __float128 I0, const __float128 Ipk, const __float128 x0, const __float128 alpha, const __float128 beta);
#endif

/// Find value at position (x,y) of the 2D arbitrarily-centered unnormalized symmetric Moffat function
/** The Moffat distribution is due to \cite moffat_1969.  With \f$\beta=1\f$ it is the 
  * Lorentzian or Cauchy distribution.  See also https://en.wikipedia.org/wiki/Moffat_distribution and
  * https://en.wikipedia.org/wiki/Cauchy_distribution.
  * 
  * Here we use the unnormalized general form, most useful for peak fitting.
  * 
  * This function omputes:
  * \f[ 
    I(x) = I_0 + I_{pk}\left[ 1 + \frac{(x-x_0)^2 + (y-y_0)^2}{\alpha^2}\right]^{-\beta} 
  * \f]
  * 
  * \returns the value of the 2D arbitrarily-centered symmetric Moffat function at (x,y)
  * 
  * \tparam realT is type to use for arithmetic
  * 
  * \test Scenario: compiling 2D Moffat function \ref tests_math_func_moffat2D "[test doc]"
  *
  * \ingroup gen_math_moffats
  */ 
template<typename realT>
realT moffat2D( const realT x,     ///< [in] the x-position at which to evaluate the Moffat function
                const realT y,     ///< [in] the y-positoin at which to evaluate the Moffat function
                const realT I0,    ///< [in] the constant to add to the Moffat function
                const realT Ipk,   ///< [in] the scaling factor (peak height is  A-G0)
                const realT x0,    ///< [in] the x-coordinate of the center
                const realT y0,    ///< [in] the y-coordinate of the center
                const realT alpha, ///< [in] the width parameter of the Moffat function.
                const realT beta   ///< [in] the shape parameter of the Moffat function.
              )
{ 
   return I0 + Ipk * pow( static_cast<realT>(1) + (pow(x-x0,2)+pow(y-y0,2))/pow(alpha,2), -beta);
}

extern template 
float moffat2D<float>(const float x, const float y, const float I0, const float Ipk, const float x0, const float y0, const float alpha, const float beta);

extern template 
double moffat2D<double>(const double x, const double y, const double I0, const double Ipk, const double x0, const double y0, const double alpha, const double beta);

extern template 
long double moffat2D<long double>(const long double x, const long double y, const long double I0, const long double Ipk, const long double x0, const long double y0, const long double alpha, const long double beta);

#ifdef HASQUAD
extern template 
__float128 moffat2D<__float128>(const __float128 x, const __float128 y, const __float128 I0, const __float128 Ipk, const __float128 x0, const __float128 y0, const __float128 alpha, const __float128 beta);
#endif

/// Compute the full-width at half-maximum of a Moffat profile
/** This returns the value of
  * \f[
  * FWHM = 2 \alpha \sqrt{2^{1/\beta} - 1}
  * \f]
  *
  * \returns the FWHM of the Moffat profile
  * 
  * \tparam realT is the type to use for arithmetic
  * 
  * \test Scenario: compiling Moffat FWHM \ref tests_math_func_moffatFWHM "[test doc]"
  * 
  * \ingroup gen_math_moffats
  */
template<typename realT>
realT moffatFWHM( realT alpha, ///< [in] the width parameter of the Moffat function.
                  realT beta   ///< [in] the shape parameter of the Moffat function.
                )
{
   return 2*alpha*sqrt( pow(static_cast<realT>(2), static_cast<realT>(1)/beta) - 1);
}

extern template
float moffatFWHM(float alpha, float beta);

extern template
double moffatFWHM(double alpha, double beta);

extern template
long double moffatFWHM(long double alpha, long double beta);

#ifdef HASQUAD
extern template
__float128 moffatFWHM(__float128 alpha, __float128 beta);
#endif

} //namespace func 
} //namespace math
} //namespace mx

#endif //math_func_moffat_hpp



