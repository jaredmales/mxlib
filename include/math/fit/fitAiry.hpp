/** \file fitAiry.hpp
 * \author Jared R. Males
 * \brief Tools for fitting the Airy pattern to PSF images.
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

#ifndef fitAiry_hpp
#define fitAiry_hpp

#include "levmarInterface.hpp"
#include "../func/airyPattern.hpp"

namespace mx
{
namespace math
{
namespace fit
{

template<typename realT>
struct array2FitAiry;

template<typename _realT>
struct airy2D_obs_fitter;

template<typename _realT>
struct airy2D_obs_fitter_ps;

/** \defgroup airy_peak_fit Airy Patterns
  * \brief Fitting Airy patterns to data.
  * 
  * The Airy pattern is fit to data, typically diffraction-limited images of point sources.  
  * 
  * \ingroup peak_fit 
  */

///Class to manage fitting a 2D Airy pattern to data via the \ref levmarInterface
/** In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires this definition in fitterT
  * \code
  * static const int nparams = 5;
  * \endcode
  * where the number 5 is replaced by the number of parameters that fitterT expects to fit. 
  *
  * \tparam fitterT a type meeting the above requirements.
  *
  * \ingroup airy_peak_fit
  *
  */
template<typename fitterT>
class fitAiry2D : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

   static const int nparams = fitterT::nparams;

   array2FitAiry<realT> arr;
   
   void initialize()
   {
      this->allocate_params(nparams);
      this->adata = &arr;      
   }
   
   fitAiry2D()
   {
      initialize();
   }
      
   ~fitAiry2D()
   {
   }
   
   ///Set the initial guess when platescale and central obscuration are fixed.
   void setGuess( realT A0, ///< [in] the constant background
                  realT A,  ///< [in] the peak scaling
                  realT x0,  ///< [in] the center x-coordinate [pixels]
                  realT y0  ///< [in] the center y-coordinate [pixels]
                )
   {
      static_assert( nparams==4 , "fitAiry2D: Wrong setGuess called for platescale and/or cen-obs variable.");
      
      this->p[0] = A0;
      this->p[1] = A;
      this->p[2] = x0;
      this->p[3] = y0;      
   }
   
   ///Set the initial guess when platescale is variable, and central obscuration is fixed.
   void setGuess( realT A0, ///< [in] the constant background
                  realT A,  ///< [in] the peak scaling
                  realT x0,  ///< [in] the center x-coordinate [pixels]
                  realT y0,  ///< [in] the center y-coordinate [pixels]
                  realT ps  ///< [in] the platescale [ (lambda/D) / pixel]
                )
   {
      static_assert( nparams==5 , "fitAiry2D: Wrong setGuess called for only platescale variable.");
      
      this->p[0] = A0;
      this->p[1] = A;
      this->p[2] = x0;
      this->p[3] = y0;      
      this->p[4] = ps;
   }
   
   ///Set the initial guess when central-obscuration is variable.
   void setGuess( realT A0, ///< [in] the constant background
                  realT A,  ///< [in] the peak scaling
                  realT x0,  ///< [in] the center x-coordinate [pixels]
                  realT y0,  ///< [in] the center y-coordinate [pixels]
                  realT ps,  ///< [in] the platescale [ (lambda/D) / pixel]
                  realT co ///< [in] the central obscuration [ratio]
                )
   {
      static_assert( nparams==6 , "fitAiry2D: Wrong setGuess called for cen-obs not variable.");
      
      this->p[0] = A0;
      this->p[1] = A;
      this->p[2] = x0;
      this->p[3] = y0;      
      this->p[4] = ps;
      this->p[5] = co;
   }
   
   void setArray(realT *data, int nx, int ny)
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = ny;
      
      this->n = nx*ny;
      
   }
   
   void cenObs( realT co )
   {
      arr.cenObs = co;
   }
   
   void ps( realT ps )
   {
      arr.ps = ps;
   }
   
   int fit()
   {
      return levmarInterface<fitterT>::fit();      
   }
      
   realT A0()
   {
      return this->p[0];
   }
   
   realT A()
   {
      return this->p[1];
   }
   
   realT x0()
   {
      return this->p[2];
   }
   
   realT y0()
   {
      return this->p[3];
   }
   
   realT ps()
   {
      if(nparams < 5)
      {
         return arr.ps;
      }
      else
      {
         return this->p[4];
      }
   }
   
   realT cenObs()
   {
      if(nparams < 6)
      {
         return arr.cenObs;
      }
      else
      {
         return this->p[5];
      }
   }

   
};


///Wrapper for a native array to pass to \ref levmarInterface, with Airy details.
/** \ingroup airy_peak_fit
  */
template<typename realT>
struct array2FitAiry
{
   realT * data {nullptr}; ///< Pointer to the array
   size_t nx {0}; ///< X dimension of the array
   size_t ny {0}; ///< Y dimension of the array
   
   realT cenObs {0}; ///< the platescale in \f$ (\lambda/D)/pixel  \f$
   realT ps {0}; ///< is the ratio of the circular central obscuration diameter to the diameter.
};

///\ref levmarInterface fitter structure for the centrally obscured Airy pattern.
/**
  * Platescale and central obscuration are fixed.
  * 
  * \ingroup airy_peak_fit
  *
  */
template<typename _realT>
struct airy2D_obs_fitter
{
   typedef _realT realT;
   
   static const int nparams = 4;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitAiry<realT> * arr = (array2FitAiry<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      /*
        p[0] = floor
        p[1] = scale
        p[2] = x0
        p[3] = y0
      */
      
      realT r;
      
      for(int i=0;i<arr->nx; i++)
      {
         for(int j=0;j<arr->ny;j++)
         { 
            idx_mat = i+j*arr->nx;
   
            r = sqrt( pow( i-p[2],2) + pow(j-p[3],2));
            
            hx[idx_dat] = func::airyPattern( static_cast<realT>(i), static_cast<realT>(j), p[0], p[1], p[2], p[3], arr->ps, arr->cenObs)  - arr->data[idx_mat];
            
            //hx[idx_dat] *= fabs(arr->data[idx_mat]);
            
            idx_dat++;
         }
      }
   }
   

};

///\ref levmarInterface fitter structure for the obstructed Airy pattern, including platescale.
/** Central obscuration is fixed.
  * 
  * \ingroup airy_peak_fit
  *
  */
template<typename _realT>
struct airy2D_obs_fitter_ps
{
   typedef _realT realT;
   
   static const int nparams = 5;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitAiry<realT> * arr = (array2FitAiry<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      /*
        p[0] = floor
        p[1] = scale
        p[2] = x0
        p[3] = y0
        p[4] = ps [(lam/D)/pix]
      */
      
      realT r;
      
      for(int i=0;i<arr->nx; i++)
      {
         for(int j=0;j<arr->ny;j++)
         { 
            idx_mat = i+j*arr->nx;
   
            r = sqrt( pow( i-p[2],2) + pow(j-p[3],2));
            
            hx[idx_dat] = func::airyPattern( static_cast<realT>(i), static_cast<realT>(j), p[0], p[1], p[2], p[3], p[4], arr->cenObs)  - arr->data[idx_mat];
            
            //hx[idx_dat] *= fabs( pow(arr->data[idx_mat],2));
            
            idx_dat++;
         }
      }
   }
   

};


///\ref levmarInterface fitter structure for the obstructed Airy pattern, including fitting platescale and central obscuration.
/** \ingroup airy_peak_fit
  *
  */
template<typename _realT>
struct airy2D_obs_fitter_ps_eps
{
   typedef _realT realT;
   
   static const int nparams = 6;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitAiry<realT> * arr = (array2FitAiry<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      /*
        p[0] = floor
        p[1] = scale
        p[2] = x0
        p[3] = y0
        p[4] = ps [(lam/D)/pix]
        p[5] = cen-obs
      */
      
      //realT r;
      
      for(int i=0;i<arr->nx; ++i)
      {
         for(int j=0;j<arr->ny; ++j)
         { 
            idx_mat = i+j*arr->nx;
   
            //r = sqrt( pow( i-p[2],2) + pow(j-p[3],2));
            
            hx[idx_dat] = func::airyPattern( static_cast<realT>(i), static_cast<realT>(j), p[0], p[1], p[2], p[3], p[4], p[5])  - arr->data[idx_mat];
            
            //hx[idx_dat] *= fabs( pow(arr->data[idx_mat],2));
            
            idx_dat++;
         }
      }
   }
   

};

///Alias for the fitAiry2D type with both platescale and cen-obs fixed.
/** \ingroup airy_peak_fit
  */
template<typename realT>
using fitAiry2DbothFixed = mx::math::fit::fitAiry2D<mx::math::fit::airy2D_obs_fitter<realT>>;

///Alias for the fitAiry2D type with cen-obs fixed.
/** \ingroup airy_peak_fit
  */
template<typename realT>
using fitAiry2DcenObsFixed = mx::math::fit::fitAiry2D<mx::math::fit::airy2D_obs_fitter_ps<realT>>;

///Alias for the fitAiry2D type with none fixed.
/** \ingroup airy_peak_fit
  */
template<typename realT>
using fitAiry2DnoneFixed = mx::math::fit::fitAiry2D<mx::math::fit::airy2D_obs_fitter_ps_eps<realT>>;


} //namespace fit
} //namespace math
} //namespace mx

#endif //fitAiry_hpp

