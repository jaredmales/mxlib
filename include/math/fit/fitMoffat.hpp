/** \file fitMoffat.hpp
 * \author Jared R. Males
 * \brief Tools for fitting Moffat functions to data.
 * \ingroup fitting_files
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

#ifndef fitMoffat_hpp
#define fitMoffat_hpp

#include "levmarInterface.hpp"
#include "../func/moffat.hpp"

#include "../../improc/eigenImage.hpp"

namespace mx
{
namespace math
{
namespace fit
{

/** \defgroup moffat_peak_fit Moffat Functions
  * \brief Fitting Moffat functions to data.
  * 
  * The Moffat Function\cite moffat_1969, a.k.a. the Moffat Profile, a.k.a. the Moffat Distribution, has the form
  * \f[
    I(x) = I_{pk}\left[ 1 + \frac{x^2}{\alpha^2}\right]^{-\beta}
  * \f] 
  * With \f$\beta=1\f$ it is the 
  * Lorentzian or Cauchy distribution.  See also https://en.wikipedia.org/wiki/Moffat_distribution and
  * https://en.wikipedia.org/wiki/Cauchy_distribution.
  * 
  * For utilities for working with Moffat functions see \ref gen_math_moffats.
  *
  * \ingroup peak_fit 
  */

//forward
template<typename realT>
struct array2FitMoffat;

//forward
template<typename _realT>
struct moffat2D_sym_fitter;

///Class to manage fitting a 2D Moffat to data via the \ref levmarInterface
/** Fits \ref gen_math_moffats to a 2-dimensional array of data.
  *
  * This class allows for treating any of the parameters as fixed.
  * 
  * \tparam fitterT a type meeting the requirements specified in \ref levmarInterface.
  *
  * \ingroup moffat_peak_fit
  *
  */
template<typename fitterT>
class fitMoffat2D : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

protected:
   
   array2FitMoffat<realT> arr;
   
   void initialize()
   {
      this->allocate_params(arr.nparams());
      this->adata = &arr;      
   }

public:
   
   fitMoffat2D()
   {
      initialize();
   }
      
   ~fitMoffat2D()
   {
   }
   
   /// Set whether each parameter is fixed.
   /** Sets the parameter indices appropriately.
     */
   void setFixed( bool I0,    ///< [in] if true, then I0 will be not be part of the fit
                  bool I,     ///< [in] if true, then I will be not be part of the fit
                  bool x0,    ///< [in] if true, then x0 will be not be part of the fit
                  bool y0,    ///< [in] if true, then y0 will be not be part of the fit
                  bool alpha, ///< [in] if true, then alpha will be not be part of the fit
                  bool beta   ///< [in] if true, then beta will be not be part of the fit
                )
   {
      arr.setFixed(I0, I, x0, y0, alpha, beta);
      this->allocate_params(arr.nparams());
   }
   
   ///Set the initial guess for a symmetric Moffat.
   void setGuess( realT I0,    ///< [in] the constant background level
                  realT I,     ///< [in] the peak scaling
                  realT x0,    ///< [in] the center x-coordinate
                  realT y0,    ///< [in] the center y-coordinate
                  realT alpha, ///< [in] the width parameter
                  realT beta   ///< [in] the shape parameter
                )
   {
      arr.I0(this->p,I0);
      arr.I(this->p,I);
      arr.x0(this->p,x0);
      arr.y0(this->p,y0);
      arr.alpha(this->p,alpha);
      arr.beta(this->p,beta);
   }
   
   ///Set the data aray.
   void setArray( realT *data, ///< [in] pointer to an nx X ny array of data to be fit
                  int nx,      ///< [in] the number of pixels in the x direction of the data array
                  int ny       ///< [in] the number of pixels in the y direction of the data array
                )
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = ny;
      
      this->n = nx*ny;
      
   }
   
   ///Get the current value of I0, the constant.
   /**
     * \returns the current value of I0
     */ 
   realT I0()
   {
      return arr.I0( this->p );
   }

   ///Get the current value of I, the peak scaling.
   /**
     * \returns the current value of I
     */
   realT I()
   {
      return arr.I( this->p );
   }

   ///Get the center x-
   /**
     * \returns the current value of x0
     */ 
   realT x0()
   {
      return arr.x0( this->p );
   }
   
   ///Get the center y-coordinate
   /**
     * \returns the current value of y0
     */ 
   realT y0()
   {
      return arr.y0( this->p );
   }
   
   ///Return the width parameter
   /**
     * \returns the current value of alpha
     */ 
   realT alpha()
   {
      return arr.alpha( this->p );
   }
   
   /// Return the shape parameter
   /**
     * \returns the current value of beta
     */ 
   realT beta()
   {
      return arr.beta( this->p );
   }
   
   ///Return the full-width at half maximum
   /**
     * \returns the FWHM calculated from alpha and beta
     */ 
   realT fwhm()
   {
      return func::moffatFWHM(alpha(), beta());
   }
   
   
   
};

///Wrapper for a native array to pass to \ref levmarInterface, with Moffat details.
/** \ingroup moffat_peak_fit
  */
template<typename realT>
struct array2FitMoffat
{
   realT * data {nullptr}; ///< Pointer to the array
   size_t nx {0}; ///< X dimension of the array
   size_t ny {0}; ///< Y dimension of the array
   
   realT m_I0 {0};
   realT m_I {0};
   realT m_x0 {0};
   realT m_y0 {0};
   realT m_alpha {0};
   realT m_beta {0};
   
   int m_I0_idx {0};
   int m_I_idx {1};
   int m_x0_idx {2};
   int m_y0_idx {3};
   int m_alpha_idx {4};
   int m_beta_idx {5};
   
   int m_nparams = 6;
   
   /// Set whether each parameter is fixed.
   /** Sets the parameter indices appropriately.
     */
   void setFixed( bool I0, 
                  bool I, 
                  bool x0, 
                  bool y0, 
                  bool alpha, 
                  bool beta 
                )
   {
      int idx = 0;
      
      if(I0) m_I0_idx = -1;
      else m_I0_idx = idx++;
      
      if(I) m_I_idx = -1;
      else m_I_idx = idx++;
   
      if(x0) m_x0_idx = -1;
      else m_x0_idx = idx++;
      
      if(y0) m_y0_idx = -1;
      else m_y0_idx = idx++;
      
      if(alpha) m_alpha_idx = -1;
      else m_alpha_idx = idx++;
      
      if(beta) m_beta_idx = -1;
      else m_beta_idx = idx++;
      
      m_nparams = idx;
   }
   
   realT I0( realT * p )
   {
      if( m_I0_idx < 0 )
      {
         return m_I0;
      }
      else
      {
         return p[m_I0_idx];
      }
   }
   
   void I0( realT * p,
            realT nI0
           )
   {
      if( m_I0_idx < 0 )
      {
         m_I0 = nI0;
      }
      else
      {
         p[m_I0_idx]=nI0;
      }
   }
   
   realT I( realT * p )
   {
      if( m_I_idx < 0 )
      {
         return m_I;
      }
      else
      {
         return p[m_I_idx];
      }
   }
   
   void I( realT * p,
           realT nI
         )
   {
      if( m_I_idx < 0 )
      {
         m_I = nI;
      }
      else
      {
         p[m_I_idx] = nI;
      }
   }
   
   realT x0( realT * p )
   {
      if( m_x0_idx < 0 )
      {
         return m_x0;
      }
      else
      {
         return p[m_x0_idx];
      }
   }
   
   void x0( realT * p,
            realT nx0
          )
   {
      if( m_x0_idx < 0 )
      {
         m_x0 = nx0;
      }
      else
      {
         p[m_x0_idx] = nx0;
      }
   }
   
   realT y0( realT * p )
   {
      if( m_y0_idx < 0 )
      {
         return m_y0;
      }
      else
      {
         return p[m_y0_idx];
      }
   }
   
   void y0( realT * p,
            realT ny0
          )
   {
      if( m_y0_idx < 0 )
      {
         m_y0 = ny0;
      }
      else
      {
         p[m_y0_idx] = ny0;
      }
   }
   
   
   realT alpha( realT * p )
   {
      if( m_alpha_idx < 0 )
      {
         return m_alpha;
      }
      else
      {
         return p[m_alpha_idx];
      }
   }
   
   void alpha( realT * p,
               realT nalpha
             )
   {
      if( m_alpha_idx < 0 )
      {
         m_alpha = nalpha;
      }
      else
      {
         p[m_alpha_idx] = nalpha;
      }
   }
   
   
   realT beta( realT * p )
   {
      if( m_beta_idx < 0 )
      {
         return m_beta;
      }
      else
      {
         return p[m_beta_idx];
      }
   }
   
   void beta( realT * p,
              realT nbeta
            )
   {
      if( m_beta_idx < 0 )
      {
         m_beta = nbeta;
      }
      else
      {
         p[m_beta_idx] = nbeta;
      }
   }
   
   int nparams()
   {
      return m_nparams;
   }
   
   
};

///\ref levmarInterface fitter structure for the symmetric Moffat.
/** \ingroup moffat_peak_fit
  *
  */
template<typename _realT>
struct moffat2D_sym_fitter
{
   typedef _realT realT;
   
   static const int nparams = 6;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitMoffat<realT> * arr = (array2FitMoffat<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      realT I0 = arr->I0(p);
      realT I = arr->I(p);
      realT x0 = arr->x0(p);
      realT y0 = arr->y0(p);
      realT alpha = arr->alpha(p);
      realT beta = arr->beta(p);
      
      for(int i=0; i<arr->nx; ++i)
      {
         for(int j=0; j<arr->ny; ++j)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = func::moffat2D<realT>(i,j,I0,I, x0, y0, alpha, beta) - arr->data[idx_mat];
            
            ++idx_dat;
         }
      }
   }   
};


///Alias for the fitMoffat2D type fitting the symmetric Moffat profile.
/** \ingroup moffat_peak_fit
  */
template<typename realT>
using fitMoffat2Dsym = mx::math::fit::fitMoffat2D<mx::math::fit::moffat2D_sym_fitter<realT>>;


} //namespace fit
} //namespace math

} //namespace mx

#endif //fitMoffat_hpp

