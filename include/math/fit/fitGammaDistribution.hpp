/** \file fitGammaDistribution.hpp
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

#ifndef fitGammaDistribution_hpp
#define fitGammaDistribution_hpp

#include "levmarInterface.hpp"
#include "../func/gammaDistribution.hpp"

namespace mx
{
namespace math
{
namespace fit
{

template<typename realT>
struct array2FitGammaDistribution;


/** \defgroup gammaDist_peak_fit Gamma Distribution
  * \brief Fitting the Gamma Distribution to data.
  * 
  * The Gamma Distribution is fit to data.
  * 
  * \ingroup peak_fit 
  */

///Class to manage fitting the GammaDistribution Distribution to data via the \ref levmarInterface
/** In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires this definition in fitterT
  * \code
  * static const int nparams = 3;
  * \endcode
  * where the number 3 is replaced by the number of parameters that fitterT expects to fit.
  *
  * \tparam fitterT a type meeting the above requirements.
  *
  * \ingroup gammaDist_peak_fit
  *
  */
template<typename fitterT>
class fitGammaDistribution : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

   static const int nparams = fitterT::nparams;

   array2FitGammaDistribution<realT> arr;
   
   void initialize()
   {
      this->allocate_params(nparams);
      this->adata = &arr;      
   }
   
   fitGammaDistribution()
   {
      initialize();
   }
      
   ~fitGammaDistribution()
   {
   }
   
   ///Set the initial guess.
   void setGuess( realT x0,    ///< [in] the location parameter
                  realT k,     ///< [in] the shape parameter
                  realT theta, ///< [in] the scale parameter
                  realT denom  ///< [in] the denominator or 1/peak-scale
                )
   {
      static_assert( nparams==4 , "fitGammaDistribution: Wrong setGuess called for no location parameter.");

      this->p[2] = x0; //this is p[2] to make it easy to leave out
      this->p[0] = k;
      this->p[1] = theta;
      this->p[3] = denom;
   }

   ///Set the initial guess.
   void setGuess( realT x0,   ///< [in] the location parameter
                  realT k,    ///< [in] the shape parameter
                  realT theta ///< [in] the scale parameter
                )
   {
      static_assert( nparams==3 , "fitGammaDistribution: Wrong setGuess called for no location parameter.");
      
      this->p[2] = x0; //this is p[2] to make it easy to leave out
      this->p[0] = k;
      this->p[1] = theta;
   }
   
   ///Set the initial guess when no location parameter is used.
   void setGuess( realT k,    ///< [in] the shape parameter
                  realT theta ///< [in] the scale parameter
                )
   {
      static_assert( nparams==2 , "fitGammaDistribution: Wrong setGuess called for location parameter.");

      this->p[0] = k;
      this->p[1] = theta;
   }

   void setArray(realT *data, int n)
   {
      arr.data = data;
      arr.n = n;

      this->n = n;
      
   }
   
   void x0( realT nx0 )
   {
      arr.x0 = nx0;
      if(nparams == 3)
      {
         this->p[2] = nx0;  //this is p[2] to make it easy to leave out
      }
   }
   
   void k( realT nk )
   {
      arr.k = nk;
      this->p[0] = nk;
   }
   
   void theta( realT na )
   {
      arr.lambda = na;
      this->p[1] = na;
   }

   void denom( realT nd)
   {
      arr.denom = nd;
      if(nparams == 4)
      {
         this->p[3] = nd;
      }
   }

   int fit()
   {
      return levmarInterface<fitterT>::fit();      
   }
      
   realT x0()
   {
      if(nparams == 3)
      {
         return this->p[2];
      }
      else
      {
         return 0;
      }
   }
   
   realT k()
   {
      return this->p[0];
   }
   
   realT theta()
   {
      return this->p[1];
   }
   
   realT denom()
   {
      if(nparams == 4)
      {
         return this->p[3];
      }
      else
      {
         return func::gammaDistributionDenom<realT>(this->p[0], this->p[1]);
      }
   }

   
};


///Wrapper for a native array to pass to \ref levmarInterface, with GammaDistribution details.
/** \ingroup gammaDist_peak_fit
  */
template<typename realT>
struct array2FitGammaDistribution
{
   realT * data {nullptr}; ///< Pointer to the array
   size_t n {0}; ///< dimension of the array
   
   realT x0 {0}; ///< the location parameter.
   realT k {0}; ///< the shape parameter
   realT theta {0}; ///< the scale parameter
   realT denom {0}; ///< the denominator or 1/peak-scale
};

///\ref levmarInterface fitter structure for the shifted Gamma Distribution with arbitrary peak scaling
/**
  *
  * \ingroup gammaDist_peak_fit
  *
  */
template<typename _realT>
struct gammaDistribution_4param_fitter
{
   typedef _realT realT;

   static const int nparams = 4;

   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitGammaDistribution<realT> * arr = (array2FitGammaDistribution<realT> *) adata;

      for(int i=0;i<arr->n; i++)
      {
         hx[i] = func::gammaDistribution<realT>(i, p[2], p[0], p[1], p[3]) - arr->data[i];
      }
   }
};

///\ref levmarInterface fitter structure for the shifted Gamma Distribution
/**
  * 
  * \ingroup gammaDist_peak_fit
  *
  */
template<typename _realT>
struct gammaDistribution_3param_fitter
{
   typedef _realT realT;
   
   static const int nparams = 3;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitGammaDistribution<realT> * arr = (array2FitGammaDistribution<realT> *) adata;

      realT denom = func::gammaDistributionDenom<realT>(p[0], p[1]);

      for(int i=0;i<arr->n; i++)
      {
         hx[i] = func::gammaDistribution<realT>(i, p[2], p[0], p[1], denom) - arr->data[i];
      }
   }
};

///\ref levmarInterface fitter structure for the Gamma Distribution
/**
  *
  * \ingroup gammaDist_peak_fit
  *
  */
template<typename _realT>
struct gammaDistribution_2param_fitter
{
   typedef _realT realT;

   static const int nparams = 2;

   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitGammaDistribution<realT> * arr = (array2FitGammaDistribution<realT> *) adata;

      realT denom = func::gammaDistributionDenom<realT>(p[0], p[1]);

      for(int i=0;i<arr->n; i++)
      {
         hx[i] = func::gammaDistribution<realT>(i, 0.0, p[0], p[1], denom) - arr->data[i];
      }
   }
};


} //namespace fit
} //namespace math
} //namespace mx

#endif //fitGammaDistribution_hpp

