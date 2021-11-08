/** \file airyPattern.hpp
  * \author Jared R. Males
  * \brief Utilities related to the Airy pattern point spread function.
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

#ifndef math_func_airyPattern_hpp
#define math_func_airyPattern_hpp


#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../constants.hpp"
#include "bessel.hpp"
#include "jinc.hpp"


namespace mx
{
   
namespace math
{

namespace func
{
   
///The classical Airy pattern
/** Returns the intensity distribution of the Airy pattern at a given \f$ \lambda/D \f$
  * 
  * References: \cite born_and_wolf, \cite mahajan_1986, https://en.wikipedia.org/wiki/Airy_disk.
  * 
  * \tparam realT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT airyPattern( realT x /**< [in] is the separation in units of \f$ \lambda/D \f$. */)
{
   return pow(2*jinc(pi<realT>()*x),2); 
}

///The centrally obscured Airy pattern
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$
  * 
  * References: \cite born_and_wolf, \cite mahajan_1986, https://en.wikipedia.org/wiki/Airy_disk.
  * 
  * \tparam realT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT airyPattern( realT x,  ///< [in] is the separation in units of  \f$ \lambda/D \f$.
                   realT eps ///< [in] is the ratio of the circular central obscuration diameter to the diameter.
                 )
{
   return (1./pow(1.-eps*eps,2))*pow(2.*jinc( pi<realT>()*x)-eps*eps*2.*jinc(eps*pi<realT>()*x),2);
}


///The general centrally obscured Airy pattern, with arbitrary center and platescale.
/** Returns the intensity distribution of the centrally obscured Airy pattern at a given \f$ \lambda/D \f$.  This version allows
  * for an arbitrary center, scaling, and platescale.
  * 
  * References: \cite born_and_wolf, \cite mahajan_1986, https://en.wikipedia.org/wiki/Airy_disk.
  * 
  * \tparam realT is the floating point type used for arithmetic.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT airyPattern( realT x,  ///< [in] is the x-coordinate in units of pixels
                   realT y,  ///< [in] is the y-coordinate in units of pixels
                   realT A0, ///< [in] constant value added to the Airy pattern
                   realT A, ///< [in] peak scale of the Airy pattern.
                   realT x0, ///< [in] is the x-center in units of pixels
                   realT y0, ///< [in] is the y-center in units of pixels
                   realT ps, ///< [in] the platescale in \f$ (\lambda/D)/pixel  \f$
                   realT eps ///< [in] is the ratio of the circular central obscuration diameter to the diameter.
                 )
{
   realT r = sqrt( pow(x-x0,2) + pow(y-y0,2)) * ps;
   
   return A0 + A*airyPattern(r, eps);
}

/// Fill in an array with the 2D arbitrarily-centered classical Airy pattern
/**
  * At each pixel (x,y) of the array this computes:
  * 
  * \f$ f(x,y) = A_0 + A\times Airy(\sqrt(x^2+y^2) \f$
  * 
  * \tparam realT is the type to use for arithmetic
  * 
  * \ingroup gen_math_gaussians
  */ 
template<typename realT>
void airyPattern2D( realT * arr,    ///< [out] is the allocated array to fill in
                    size_t nx,      ///< [in] is the size of the x dimension of the array (rows)
                    size_t ny,      ///< [in] is the size of the y dimension of the array (columns)
                    const realT A0, ///< [in] is the constant to add to the Gaussian
                    const realT A,  ///< [in] is the scaling factor (peak height = A-A0)
                    const realT x0, ///< [in] is the x-coordinate of the center
                    const realT y0, ///< [in] is the y-coordinate of the center
                    realT ps        ///< [in] the platescale in \f$ (\lambda/D)/pixel  \f$
                  )
{
   size_t idx;
   
   for(size_t j=0;j<ny; ++j)
   {
      for(size_t i=0;i<nx; ++i)
      {
         idx = i + j*nx;
         
         realT rad = sqrt( pow(i-x0,2) + pow(j-y0,2) ) * ps;

         arr[idx] = A0 + A* airyPattern(rad);
      }
   }
}

///Seeing Halo Profile
/** A Moffat profile due to Roddier (1981)\cite roddier_1981, see also Racine et al (1999)\cite racine_1999, which
  * can be used to describe the seeing limited point-spread function of a telescope imaging through turbulence, or the halo in partially corrected imaging. 
  *  
  * \returns the value of the profile at x, in units of fractional flux per unit area.
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */
template<typename realT>
realT seeingHalo( realT x, ///< [in] the separation in the same units as fwhm.
                   realT fwhm ///< [in] the fwhm in arbitrary units.  Note that this defines the area units of the density.
                 )
{
   return (0.488/(fwhm*fwhm))*pow(1. + (11./6.)*pow(x/fwhm,2), -11./6.);
}

/// Calculate the fraction of enclosed power at a given radius for the unobscured Airy Pattern 
/** See Mahajan (1986) \cite mahajan_1986 and https://en.wikipedia.org/wiki/Airy_disk.
  *
  * \returns the fraction of enclosed power at radius x
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */ 
template<typename realT>
realT airyPatternEnclosed( realT x /**< [in] the radius */)
{
   realT x1 = x*pi<realT>();
   
   realT b0 = bessel_j<realT>(0,x1);
   b0=b0*b0;
   
   realT b1 = bessel_j<realT>(1,x1);
   b1=b1*b1;
   
   realT encp = static_cast<realT>(1) - b0 - b1;
   
   return encp;
}

template<typename realT>
realT apeInt( realT x,
              void * params
            )
{
   realT eps = *static_cast<realT*>(params);
   
   return bessel_j<realT>(1,x)*bessel_j<realT>(1,eps*x)/x;
}

/// Calculate the fraction of enclosed power at a given radius for the centrally obscured Airy Pattern 
/** See Mahajan (1986) \cite mahajan_1986 and https://en.wikipedia.org/wiki/Airy_disk.
  * If eps = 0, this calls the unobscured version.
  * 
  * \returns the fraction of enclosed power at radius x
  * 
  * \ingroup psfs
  * \ingroup gen_math_airy_pattern
  */ 
template<typename realT>
realT airyPatternEnclosed( realT x,  ///< [in] the radius
                           realT eps ///< [in] the central obscuration fraction
                         )
{
   if(eps == 0) return airyPatternEnclosed(x);
   
   gsl_function func;
   func.function = apeInt<realT>;
   func.params = &eps;
   
   realT jint;
   realT abserr;
   size_t neval;
   
   realT x1 = x*pi<realT>();
   
   gsl_set_error_handler_off();
   gsl_integration_qng( &func, 0, x1, 1e-7, 1e-7, &jint,&abserr, &neval);
   
   realT b0 = bessel_j<realT>(0,x1);
   b0=b0*b0;
   realT b1 = bessel_j<realT>(1,x1);
   b1=b1*b1;
   
   realT b0e = bessel_j<realT>(0,eps*x1);
   b0e=b0e*b0e;
   realT b1e = bessel_j<realT>(1,eps*x1);
   b1e=b1e*b1e;
   
   realT eps2 = pow(eps,2);
   
   realT encp = static_cast<realT>(1) - b0 - b1 + eps2*(static_cast<realT>(1) - b0e - b1e);
   
   encp = encp - 4*eps*jint;
   encp = encp/(static_cast<realT>(1)-eps2);
   
   return encp;
}

} //namespace func
}//namespace math
}//namespace mx

#endif //math_func_airyPattern_hpp

