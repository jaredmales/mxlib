/** \file arpsd.hpp
  * \brief Tools for working with autogressive PSDs
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2023 Jared R. Males (jaredmales@gmail.com)
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

#ifndef arpsd_hpp
#define arpsd_hpp

#include <vector>
#include <complex>

#include <mx/math/constants.hpp>

namespace mx
{
namespace sigproc 
{

/// Get the phase of a given frequency
/**
  * \returns \f$ 2\pi \frac{f}{f_s} \f$
  *
  * \ingroup arpsds
  */
template<typename realT>
realT ar1FreqPhase( realT f,    ///< [in] the desired frequency, usually Hz
                    realT fsamp ///< [in] the samplingfrequency, same units as \p freq
                  )
{
    return math::two_pi<realT>()*f/fsamp;
}

/// Calculate the norm of an AR1 PSD
/** For the PSD of the form
  * \f[
  * P(f) = \frac{\beta}{ \left| 1 - \alpha e ^{-i 2\pi f/f_s}\right|^2 }
  * \f]
  * where \f$ f_s \f$ is the sampling frequency and \f$ \alpha \f$ is complex with
  * the pole frequency set by
  * \f[
  * f_p = \arg(\alpha)\frac{f_s}{2\pi}
  * \f]
  * this returns the one-sided norm
  * \f[
    \int_0^{f_s/2} P(f) df
  * \f]
  * 
  * Uses 2.553 \#3 of Gradshteyn (2007) \cite gradshteyn_2007
  *
  * \ingroup arpsds
  */
template<typename realT>
realT ar1PSDNorm( realT beta,  ///< the normalization parameter
                  realT alpha, ///< the magnitude of the AR1 constant 
                  realT fpole, ///< the pole frequency, which sets the phase of complex \f$\alpha\f$
                  realT fsamp  ///< the sampling frequency
                )
{
    realT a = 1 + pow(alpha,2);
    realT b = -2*alpha;

    realT sqab = sqrt(a*a - b*b);

    realT wf = math::two_pi<realT>()*(fsamp/2.0-fpole)/fsamp;
    realT w0 = math::two_pi<realT>()*(0.0-fpole)/fsamp;

    return (beta*fsamp/math::two_pi<realT>()) *  2/sqab*(atan( ((a-b)*tan(wf/2) / sqab)) - atan( ((a-b)*tan(w0/2)) / sqab));
    


}

/// Generate an AR1 PSD
/** Fills in the the PSD of the form
  * \f[
  * P(f) = \frac{\beta}{ \left| 1 - \alpha e ^{-i 2\pi f/f_s}\right|^2 }
  * \f]
  * where \f$ f_s \f$ is the sampling frequency and the AR1 constant \f$ \alpha \f$ is complex 
  * with phase set by the pole frequency 
  * \f[
  * \alpha = \left|\alpha\right|^2 e^{-i 2\pi f_p/ f_s}
  * \f]
  *
  * \ingroup arpsds
  */
template<typename realT>
void ar1PSD( std::vector<realT> & psd,     ///< [out] the PSD, will be resized and filled in with the values of the AR1 PSD
             const std::vector<realT> & f, ///< [in] the frequency scale.  2*f.back() is the sampling frequency.
             realT beta,                   ///< [in] the normalization parameter
             realT alpha,                  ///< [in] magnitude of the AR1 constant
             realT fpole                   ///< [in] the pole frequency, which sets the phase of complex \f$ \alpha \f$.
           )
{
    psd.resize(f.size());
    
    realT fmax = 2*f.back();

    for(size_t n = 0; n < f.size(); ++n)
    {
        realT denom = 1 + alpha*alpha - 2*alpha*cos(math::two_pi<realT>()*(f[n] - fpole)/fmax);

        psd[n] = beta/denom;
    }

}

/// Generate an AR1 PSD
/** Fills in the the PSD of the form
  * \f[
  * P(f) = \frac{\beta}{ \left| 1 - \alpha e ^{-i 2\pi f/f_s}\right|^2 }
  * \f]
  * where \f$ f_s \f$ is the sampling frequency and the AR1 constant \f$ \alpha \f$.
  *
  * \overload
  * \ingroup arpsds
  */
template<typename realT>
void ar1PSD( std::vector<realT> & psd,     ///< [out] the PSD, will be resized and filled in with the values of the AR1 PSD
             const std::vector<realT> & f, ///< [in] the frequency scale.  2*f.back() is the sampling frequency.
             realT beta,                   ///< [in] the normalization parameter
             std::complex<realT> alpha     ///< [in] the complex AR1 constant
           )
{
    realT fpole = atan2(alpha.imag(), alpha.real()) * (2*f.back())/math::two_pi<realT>();

    return ar1PSD(psd, f, beta, abs(alpha), fpole);
}

template<typename realT>
void arNPSD( std::vector<realT> & psd,
             const std::vector<realT> & f,
             realT fmax,
             realT var,
             std::vector<realT> alphaMag,
             std::vector<realT> alphaPhase
           )
{
    psd.resize(f.size());
    
    std::complex j = std::complex<realT>(0,1);
    std::vector<std::complex<realT>> alpha(alphaMag.size());
    
    for(size_t n=0; n < alphaMag.size(); ++n)
    {
        alpha[n] = alphaMag[n]*exp(j*alphaPhase[n]);
    }

    for(size_t n = 0; n < f.size(); ++n)
    {

        std::complex<realT> denom = 1.0;
        
        for(size_t m=0; m < alpha.size(); ++m)
        {
            denom -= alpha[m] * exp(-j*math::two_pi<realT>()*f[n]*std::complex<realT>((m+1),0)/fmax);
        }

        psd[n] = var / std::norm(denom);
    }

}

} //namespace sigproc
} //namespace mx

#endif //arpsd_hpp
