/** \file psdVarMean.hpp
  * \brief Tools for calculating the variance of the mean of a PSD
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2018,2019,2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef psdVarMean_hpp
#define psdVarMean_hpp

#include <vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../math/func/bessel.hpp"
#include "../math/constants.hpp"
#include "../math/gslInterpolator.hpp"

#include "psdFilter.hpp"

namespace mx
{
namespace sigproc 
{
 
/// Parameters for calculating the variance of the mean from a numerical PSD.
/** Interpolates on the PSD for the integration.
  * 
  *  \tparam realT the real floating point type of the frequency and PSD vectors. [Must be double due to GSL] 
  */
template<typename _realT>
struct psdVarMeanParams
{
   typedef _realT realT;
   
   std::vector<realT> * freq {nullptr}; ///< Pointer to the frequency vector
   std::vector<realT> * psd {nullptr}; ///< Pointer to the psd vector
   realT T {0}; ///< The sample length (units appropriate for freq)
   
   realT minf {0}; ///< The minimum frequency
   realT maxf {0}; ///< The maximum frequency
   
   realT extrapAlpha{0};
   
   bool initialized {false}; ///< Flag controlling whether the interpolator is initialized
   
   #ifdef MX_OLD_GSL
   math::gslInterpolator<math::gsl_interp_linear<realT>> terp; ///< Interpolator for the PSD
   #else 
   math::gslInterpolator<math::gsl_interp_steffen<realT>> terp; ///< Interpolator for the PSD
   #endif
   
   /// Get the value of the PSD at k using interpolation
   /** Extrapolates outside the range of frequencies in freq.
     * Right now this is just constant. 
     *
     * \todo implement better extrapolation past max-freq.
     */ 
   realT psdVal( realT k /**< [in] The transformed frequency coordinate */)
   {
      if( !initialized )
      {
         if(freq  == nullptr || psd == nullptr)
         {
            std::cerr << "Vectors not set.\n";
            exit(-1);
         }
         if(T == 0)
         {
            std::cerr << "T not set.\n";
            exit(-1);
         }
         
         terp.setup(*freq, *psd);
         //#ifdef MX_OLD_GSL
         //terp.setup(gsl_interp_linear, *freq, *psd);
         //#else 
         //terp.setup(gsl_interp_steffen, *freq, *psd);
         //#endif

         minf = (*freq)[0];
         maxf = (*freq)[freq->size()-1];
         
         initialized = true;
      }
      
      realT f = (2.0/T)*k;
      
      if(f < minf) return (*psd)[0];
      
      if(f > maxf)
      {
         //Enter extrapolation here
         //return 0;//(*psd)[psd->size()-1];
         if(extrapAlpha == 0) return 0;
            
         return (*psd)[psd->size()-1] * pow( f/maxf, -1*extrapAlpha);
      }
      
      return terp(f);
   }
};

/// Integration worker function for psdVarMean::operator()
/**
  * \tparam paramsT the type of the parameter structure. See psdVarMean.
  * 
  * \returns the value of the integrand at k
  * 
  *
  */ 
template<typename paramsT>
typename paramsT::realT psdVarMeanFunc( typename paramsT::realT k, ///< [in] the scaled frequency coordinate
                                        void * params              ///< [in] pointer to the parameter structure
                                      )
{
   typedef typename paramsT::realT realT;
   
   paramsT * pvar = (paramsT *) params;
   
   realT J = math::func::bessel_j<realT, realT>(0.5, math::two_pi<realT>()*k);
   
   return J*J/k*pvar->psdVal(k);
   
}

/// Calculate the variance of the mean for a process given its PSD 
/** Functor which manages the GSL integration workspace memory.
  *
  * \todo document the paramsT interface
  * 
  * \tparam paramsT is a parameter structure with the correct interface.
  * 
  *  \ingroup psds
  */
template<typename paramsT>
struct psdVarMean
{
   typedef typename paramsT::realT realT;
   
   realT extrapAlpha {0}; ///< The extrapolation exponent.  No extrapolation if 0.
   
protected:
   gsl_integration_workspace * w {nullptr}; ///< The GSL integration workspace memory
      
   size_t wSize {10000}; ///< The size of the GSL workspace
   
   /// Initializes the gsl library.
   void init()
   {
      ///\todo We'll get the occasional failure to reach tolerance error, just ignore them all for now.
      gsl_set_error_handler_off();
   }
   
public:
   /// Default c'tor
   psdVarMean() 
   {
      init();
   } 

   /// Constructor which initializes workspace size
   psdVarMean( size_t wm /**< [in] Size of the GSL integration workspace*/)
   {
      init();
      wSize = wm;
   }

   /// d'tor cleans up the workspace memory
   ~psdVarMean()
   {
      if(w != nullptr)
      {
         gsl_integration_workspace_free(w);
      }
   }

   /// Calculate the variance of the mean for a psd over a sample length 
   /** 
     * \returns the variance of the mean for the process governed by the PSD.
     */ 
   realT operator()( realT & error,             ///< [out] The error estimate returned by gsl_integration_qagiu 
                     std::vector<realT> & freq, ///< [in] The frequency vector
                     std::vector<realT> & psd,  ///< [in] The psd vector
                     realT T                    ///< [in] The sample length (units approprate for freq)
                   )
   {
      if(w == nullptr)
      {
         w = gsl_integration_workspace_alloc(wSize);
      }
      
      //Setup the integration parameters structure
      paramsT pvar;
   
      pvar.freq = &freq;
      pvar.psd = &psd;
      pvar.T = T;
   
      pvar.extrapAlpha = extrapAlpha;
      
      //Setup the function
      gsl_function func;

      func.function = &psdVarMeanFunc<paramsT>;
      func.params = &pvar;
   
      realT result;

      int ec = gsl_integration_qagiu (&func, 0.0, 1e-14, 1e-8, wSize, w, &result, &error);
      static_cast<void>(ec);
      
      return result/(2*T);
   }
};



} //namespace sigproc 
} //namespace mx

#endif //psdVarMean_hpp
