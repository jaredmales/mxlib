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

#ifndef moffat_hpp
#define moffat_hpp

#include <cmath>


namespace mx
{
namespace math 
{
namespace func 
{
   
   
///Find value at position (x) of the 1D arbitrarily-centered symmetric unnormalized Moffat function
/** The Moffat distribution is due to \cite moffat_1969.  With \f$\beta=1\f$ it is the 
  * Lorentzian or Cauchy distribution.  See also https://en.wikipedia.org/wiki/Moffat_distribution and
  * https://en.wikipedia.org/wiki/Cauchy_distribution.
  * 
  * Here we use the unnormalized general form, most useful for peak fitting.
  * 
  * This function computes:
  * \f$ I(x) = I_0 + I\left[ 1 + \frac{(x-x_0)^2}{\alpha^2}\right]^{-\beta} \f$
  * 
  * 
  * \returns the value of the 1D arbitrarily-centered unnormalized Moffat at (x)
  * 
  * \tparam realT is type to use for arithmetic
  * 
  * \ingroup gen_math_moffats
  */ 
template<typename realT>
realT moffat( const realT x,   ///< [in] is the x-position at which to evaluate the Moffat function
              const realT I0,    ///< [in] is the constant to add to the Moffat function
              const realT I,     ///< [in] is the scaling factor (peak = A)
              const realT x0,    ///< [in] is the x-coordinate of the center
              const realT alpha, ///< [in] is the width parameter of the Moffat function.
              const realT beta   ///< [in] is the shape parameter of the Moffat function
            )
{ 
   return I0 + I * pow( static_cast<realT>(1) + pow(x-x0,2)/pow(alpha,2), -beta);
}

///Find value at position (x,y) of the 2D arbitrarily-centered unnormalized symmetric Moffat function
/** The Moffat distribution is due to \cite moffat_1969.  With \f$\beta=1\f$ it is the 
  * Lorentzian or Cauchy distribution.  See also https://en.wikipedia.org/wiki/Moffat_distribution and
  * https://en.wikipedia.org/wiki/Cauchy_distribution.
  * 
  * Here we use the unnormalized general form, most useful for peak fitting.
  * 
  * This function omputes:
  * \f[ 
    I(x) = I_0 + I\left[ 1 + \frac{(x-x_0)^2 + (y-y_0)^2}{\alpha^2}\right]^{-\beta} 
  * \f]
  * 
  * \returns the value of the 2D arbitrarily-centered symmetric Moffat function at (x,y)
  * 
  * \tparam realT is type to use for arithmetic
  * 
  * \ingroup gen_math_moffats
  */ 
template<typename realT>
realT moffat2D( const realT x,     ///< [in] the x-position at which to evaluate the Moffat function
                const realT y,     ///< [in] the y-positoin at which to evaluate the Moffat function
                const realT I0,    ///< [in] the constant to add to the Moffat function
                const realT I,     ///< [in] the scaling factor (peak height is  A-G0)
                const realT x0,    ///< [in] the x-coordinate of the center
                const realT y0,    ///< [in] the y-coordinate of the center
                const realT alpha, ///< [in] the width parameter of the Moffat function.
                const realT beta   ///< [in] the shape parameter of the Moffat function.
              )
{ 
   return I0 + I * pow( static_cast<realT>(1) + (pow(x-x0,2)+pow(y-y0,2))/pow(alpha,2), -beta);
}

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
  * \ingroup gen_math_moffats
  */
template<typename realT>
realT moffatFWHM( realT alpha, ///< [in] the width parameter of the Moffat function.
                  realT beta   ///< [in] the shape parameter of the Moffat function.
                )
{
   return 2*alpha*sqrt( pow(static_cast<realT>(2), static_cast<realT>(1)/beta) - 1);
}

} //namespace func 
} //namespace math
} //namespace mx

#endif //moffat_hpp



