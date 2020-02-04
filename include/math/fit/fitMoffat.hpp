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

/** \defgroup moffat_peak_fit Moffat functions
  * \brief Fitting Moffat functions to data.
  * 
  * The Moffat function is fit to data.  
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
/** In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires the following in fitterT
  * \code
    static const int nparams = 6; //where the number 7 is replaced by the number of parameters that fitterT expects to fit. 
   
    static realT I0( array2FitMoffat<realT> & arr,
                     realT * p
                   );
                  
    static void I0( array2FitMoffat<realT> & arr,
                    realT * p,
                    realT nI0
                  );
  * \endcode
  * 
  *
  * \tparam fitterT a type meeting the above requirements.
  *
  * \ingroup moffat_peak_fit
  *
  */
template<typename fitterT>
class fitMoffat2D : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

   static const int nparams = fitterT::nparams;

   array2FitMoffat<realT> arr;
   
   void initialize()
   {
      this->allocate_params(nparams);
      this->adata = &arr;      
   }
   
   fitMoffat2D()
   {
      initialize();
   }
      
   ~fitMoffat2D()
   {
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
      fitterT::I0(arr,this->p,I0);
      fitterT::I(arr,this->p,I);
      fitterT::x0(arr,this->p,x0);
      fitterT::y0(arr,this->p,y0);
      fitterT::alpha(arr,this->p,alpha);
      fitterT::beta(arr,this->p,beta);
   }
   
   
   
   ///Set the data aray.
   void setArray( realT *data, 
                  int nx, 
                  int ny
                )
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = ny;
      
      this->n = nx*ny;
      
   }
   
   ///Do the fit.
   int fit()
   {
      fitterT fitter;
      levmarInterface<fitterT>::fit();      
   }
     
   ///Get the current value of I0, the constant.
   /**
     * \returns the current value of I0, which is p[0].
     */ 
   realT I0()
   {
      return fitterT::I0( this->arr, this->p );
   }

   ///Get the peak scaling.
   realT I()
   {
      return fitterT::I( this->arr, this->p );
   }

   ///Get the center x-coordinate
   realT x0()
   {
      return fitterT::x0( this->arr, this->p );
   }
   
   ///Get the center y-coordinate
   realT y0()
   {
      return fitterT::y0( this->arr, this->p );
   }
   
   ///Return the width parameter
   realT alpha()
   {
      return fitterT::alpha( this->arr, this->p );
   }
   
   realT beta()
   {
      return fitterT::beta( this->arr, this->p );
   }
   
   ///Return the full-width at half maximum
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
   
      for(int i=0;i<arr->nx; i++)
      {
         for(int j=0;j<arr->ny;j++)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = func::moffat2D<realT>(i,j,p[0],p[1], p[2], p[3], p[4], p[5]) - arr->data[idx_mat];
            
            idx_dat++;
         }
      }
   }
   
   static realT I0( array2FitMoffat<realT> & arr,
                    realT * p
                  )
   {
      static_cast<void>(arr);
      return p[0];
   }
   
   static void I0( array2FitMoffat<realT> & arr,
                   realT * p,
                   realT nI0
                 )
   {
      static_cast<void>(arr);
      p[0] = nI0;
   }
   
   static realT I( array2FitMoffat<realT> & arr,
                   realT * p
                 )
   {
      static_cast<void>(arr);
      return p[1];
   }
   
   static void I( array2FitMoffat<realT> & arr,
                  realT * p,
                  realT nI
                )
   {
      static_cast<void>(arr);
      p[1] = nI;
   }
   
   static realT x0( array2FitMoffat<realT> & arr,
                    realT * p
                  )
   {
      static_cast<void>(arr);
      return p[2];
   }
   
   static void x0( array2FitMoffat<realT> & arr,
                   realT * p,
                   realT nx0
                 )
   {
      static_cast<void>(arr);
      p[2] = nx0;
   }
   
   static realT y0( array2FitMoffat<realT> & arr,
                    realT * p
                  )
   {
      static_cast<void>(arr);
      return p[3];
   }

   static void y0( array2FitMoffat<realT> & arr,
                   realT * p,
                   realT ny0
                 )
   {
      static_cast<void>(arr);
      p[3] = ny0;
   }
   
   static realT alpha( array2FitMoffat<realT> & arr,
                       realT * p
                     )
   {
      static_cast<void>(arr);
      return p[4];
   }
   
   static void alpha( array2FitMoffat<realT> & arr,
                      realT * p,
                      realT na
                    )
   {
      static_cast<void>(arr);
      p[4] = na;
   }
   
   static realT beta( array2FitMoffat<realT> & arr,
                      realT * p
                    )
   {
      static_cast<void>(arr);
      return p[5];
   }
   
   static void beta( array2FitMoffat<realT> & arr,
                     realT * p,
                     realT nb
                   )
   {
      static_cast<void>(arr);
      p[5] = nb;
   }
};

///\ref levmarInterface fitter structure for the symmetric Moffat.
/** \ingroup moffat_peak_fit
  *
  */
template<typename _realT>
struct moffat2D_sym_fitter_alpha_beta_only
{
   typedef _realT realT;
   
   static const int nparams = 6;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitMoffat<realT> * arr = (array2FitMoffat<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      for(int i=0;i<arr->nx; i++)
      {
         for(int j=0;j<arr->ny;j++)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = func::moffat2D<realT>(i,j,arr->m_I0,arr->m_I, arr->m_x0, arr->m_y0, p[0], p[1]) - arr->data[idx_mat];
            
            idx_dat++;
         }
      }
   }
   
   static realT I0( array2FitMoffat<realT> & arr,
                    realT * p
                  )
   {
      static_cast<void>(p);
      return arr->m_I0;
   }
   
   static void I0( array2FitMoffat<realT> & arr,
                   realT * p,
                   realT nI0
                 )
   {
      static_cast<void>(p);
      arr.m_I0 = nI0;
   }
   
   static realT I( array2FitMoffat<realT> & arr,
                   realT * p
                 )
   {
      static_cast<void>(p);
      return arr.m_I;
   }
   
   static void I( array2FitMoffat<realT> & arr,
                  realT * p,
                  realT nI
                )
   {
      static_cast<void>(p);
      arr.m_I = nI;
   }
   
   static realT x0( array2FitMoffat<realT> & arr,
                    realT * p
                  )
   {
      static_cast<void>(p);
      return arr.m_x0;
   }
   
   static void x0( array2FitMoffat<realT> & arr,
                   realT * p,
                   realT nx0
                 )
   {
      static_cast<void>(p);
      arr.m_x0 = nx0;
   }
   
   static realT y0( array2FitMoffat<realT> & arr,
                    realT * p
                  )
   {
      static_cast<void>(p);
      return arr.m_y0;
   }

   static void y0( array2FitMoffat<realT> & arr,
                   realT * p,
                   realT ny0
                 )
   {
      static_cast<void>(p);
      arr.m_y0 = ny0;
   }
   
   static realT alpha( array2FitMoffat<realT> & arr,
                       realT * p
                     )
   {
      static_cast<void>(arr);
      return p[0];
   }
   
   static void alpha( array2FitMoffat<realT> & arr,
                      realT * p,
                      realT a
                    )
   {
      static_cast<void>(arr);
      p[0] = a;
   }
   
   static realT beta( array2FitMoffat<realT> & arr,
                      realT * p
                    )
   {
      static_cast<void>(arr);
      return p[1];
   }
   
   static void beta( array2FitMoffat<realT> & arr,
                     realT * p,
                     realT b
                   )
   {
      static_cast<void>(arr);
      p[1] = b;
   }
   
};

///Alias for the fitMoffat2D type fitting the symmetric Moffat profile.
/** \ingroup moffat_peak_fit
  */
template<typename realT>
using fitMoffat2Dsym = mx::math::fit::fitMoffat2D<mx::math::fit::moffat2D_sym_fitter<realT>>;

///Alias for the fitMoffat2D type fitting the symmetric Moffat profile with only \f$\alpha\f$ and \f$\beta\f$ free.
/** \ingroup moffat_peak_fit
  */
template<typename realT>
using fitMoffat2Dsym_alpha_beta_only = mx::math::fit::fitMoffat2D<mx::math::fit::moffat2D_sym_fitter_alpha_beta_only<realT>>;

} //namespace fit
} //namespace math

} //namespace mx

#endif //fitMoffat_hpp

