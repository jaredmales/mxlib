/** \file fitExpModGaussian.hpp
 * \author Jared R. Males
 * \brief Tools for fitting the GammaDistribution distribution to data.
 * \ingroup fitting_files
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

#ifndef fitExpModGaussian_hpp
#define fitExpModGaussian_hpp

#include "levmarInterface.hpp"
#include "../func/expModGaussian.hpp"

namespace mx
{
namespace math
{
namespace fit
{

template<typename realT>
struct array2FitExpModGaussian;


/** \defgroup expModGauss_peak_fit Exponetially Modified Gaussian Distribution
  * \brief Fitting the Exponetially Modified Gaussian Distribution to data.
  * 
  * The Exponetially Modified Gaussian Distribution is fit to data.
  * 
  * \ingroup peak_fit 
  */

///Class to manage fitting the Exponetially Modified Gaussian Distribution  to data via the \ref levmarInterface
/** In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires this definition in fitterT
  * \code
  * static const int nparams = 3;
  * \endcode
  * where the number 3 is replaced by the number of parameters that fitterT expects to fit.
  *
  * \tparam fitterT a type meeting the above requirements.
  *
  * \ingroup expModGauss_peak_fit
  *
  */
template<typename fitterT>
class fitExpModGaussian : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

   static const int nparams = fitterT::nparams;

   array2FitExpModGaussian<realT> arr;
   
   void initialize()
   {
      this->allocate_params(nparams);
      this->adata = &arr;      
   }
   
   fitExpModGaussian()
   {
      initialize();
   }
      
   ~fitExpModGaussian()
   {
   }
   
   ///Set the initial guess.
   void setGuess( realT G0,    ///< [in] the location parameter
                  realT mu,     ///< [in] the shape parameter
                  realT sigma, ///< [in] the scale parameter
                  realT lambda  ///< [in] the denominator or 1/peak-scale
                )
   {
      static_assert( nparams==4 , "fitExpModGaussian: Wrong setGuess called for no location parameter.");

      this->p[3] = G0; //this is p[2] to make it easy to leave out
      this->p[0] = mu;
      this->p[1] = sigma;
      this->p[2] = lambda;
   }

   ///Set the initial guess.
   void setGuess( realT mu,   ///< [in] the location parameter
                  realT sigma,    ///< [in] the shape parameter
                  realT lambda ///< [in] the scale parameter
                )
   {
      static_assert( nparams==3 , "fitExpModGaussian: Wrong setGuess called for no location parameter.");
      
      this->p[0] = mu; //this is p[2] to make it easy to leave out
      this->p[1] = sigma;
      this->p[2] = lambda;
   }

   void setArray(realT *data, int n)
   {
      arr.data = data;
      arr.n = n;

      this->n = n;
      
   }
   
   void G0( realT nG0 )
   {
      arr.G0 = nG0;
      if(nparams == 4)
      {
         this->p[3] = nG0;  //this is p[2] to make it easy to leave out
      }
   }
   
   void mu( realT nmu )
   {
      arr.mu = nmu;
      this->p[0] = nmu;
   }
   
   void sigma( realT ns )
   {
      arr.sigma = ns;
      this->p[1] = ns;
   }

   void lambda( realT nl)
   {
      arr.lambda = nl;
      this->p[2] = nl;
   }

   int fit()
   {
      return levmarInterface<fitterT>::fit();      
   }
      
   realT G0()
   {
      if(nparams == 4)
      {
         return this->p[3];
      }
      else
      {
         return 0;
      }
   }
   
   realT mu()
   {
      return this->p[0];
   }
   
   realT sigma()
   {
      return this->p[1];
   }
   
   realT lambda()
   {
      return this->p[2];
   }

   
};


///Wrapper for a native array to pass to \ref levmarInterface, with Exponentially Modified Gaussian details.
/** \ingroup gammaDist_peak_fit
  */
template<typename realT>
struct array2FitExpModGaussian
{
   realT * data {nullptr}; ///< Pointer to the array
   size_t n {0}; ///< dimension of the array
   
   realT G0 {0}; ///< the location parameter.
   realT mu {0}; ///< the shape parameter
   realT sigma {0}; ///< the scale parameter
   realT lambda {0}; ///< the denominator or 1/peak-scale
};

///\ref levmarInterface fitter structure for the Exponentially Modified Gaussian with fixed constant level.
/**
  *
  * \ingroup gammaDist_peak_fit
  *
  */
template<typename _realT>
struct fitExpModGaussian_3param_fitter
{
   typedef _realT realT;

   static const int nparams = 3;

   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitExpModGaussian<realT> * arr = (array2FitExpModGaussian<realT> *) adata;

      for(int i=0;i<arr->n; i++)
      {
         hx[i] = func::expModGaussian<realT>(i, p[0], p[1], p[2]) - arr->data[i];
      }
   }
};

///\ref levmarInterface fitter structure for the Exponentially Modified Gaussian with arbitrary constant level.
/**
  *
  * \ingroup gammaDist_peak_fit
  *
  */
template<typename _realT>
struct fitExpModGaussian_4param_fitter
{
   typedef _realT realT;

   static const int nparams = 4;

   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitExpModGaussian<realT> * arr = (array2FitExpModGaussian<realT> *) adata;

      for(int i=0;i<arr->n; i++)
      {
         hx[i] = p[3] + func::expModGaussian<realT>(i, p[0], p[1], p[2]) - arr->data[i];
      }
   }
};



} //namespace fit
} //namespace math
} //namespace mx

#endif //fitExpModGaussian_hpp

