/** \file fitGaussian.hpp
 * \author Jared R. Males
 * \brief Tools for fitting Gaussians to data.
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

#ifndef fitGaussian_hpp
#define fitGaussian_hpp

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include "levmarInterface.hpp"
#include "../func/gaussian.hpp"

#include "../../improc/eigenImage.hpp"

namespace mx
{
namespace math
{
namespace fit
{

/** \defgroup gaussian_peak_fit Gaussians
  * \brief Fitting Gaussians to data.
  * 
  * The Gaussian function is fit to data.  
  * 
  * \ingroup peak_fit 
  */


//forward
template<typename _realT>
struct gaussian1D_fitter;


///Class to manage fitting a 1D Gaussian to data via the \ref levmarInterface
/** In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires the following in fitterT
  * \code
  * static const int nparams = 4; 
  * 
  * \endcode
  * 
  *
  * \tparam fitterT a type meeting the above requirements.
  *
  * \ingroup gaussian_peak_fit
  *
  */
//template<typename fitterT>
template<typename _realT>
class fitGaussian1D : public levmarInterface<gaussian1D_fitter<_realT>> //fitterT>
{
   
public:
   
   typedef gaussian1D_fitter<_realT> fitterT;
   
   typedef typename fitterT::realT realT;

   static const int nparams = fitterT::nparams;

   array2Fit<realT> arr;
   
   void initialize()
   {
      this->allocate_params(nparams);
      this->adata = &arr;      
   }
   
   fitGaussian1D()
   {
      initialize();
   }
      
   ~fitGaussian1D()
   {
   }
   
   ///Set the initial guess for a symmetric Gaussian.
   /** Also works for the general case, setting the same width in both directions.
     */
   void setGuess( realT G0,    ///< [in] the constant background level
                  realT A,     ///< [in] the peak scaling
                  realT x0,     ///< [in] the center x-coordinate
                  realT sigma   ///< [in] the width parameter
                )
   {
      this->p[0] = G0;
      this->p[1] = A;
      this->p[2] = x0;
      this->p[3] = sigma;
   }
   
   ///Set the data aray.
   void setArray( realT *data, 
                  int nx
                )
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = 1;
      
      this->n = nx;
   }
   
   ///Do the fit.
   int fit()
   {
      fitterT fitter;
      
      levmarInterface<fitterT>::fit();
      
   }
     
   ///Get the current value of G0, the constant.
   /**
     * \returns the current value of G0, which is p[0].
     */ 
   realT G0()
   {
      return this->p[0];
   }

   ///Get the peak scaling.
   realT A()
   {
      return this->p[1];
   }

   ///Get the center x-coordinate
   realT x0()
   {
      return this->p[2];
   }
      
   ///Return the width parameter
   /** As described for the symmetric Gaussian.
     *
     * For the general Gaussian, this returns \f$ \sigma = \sqrt{ \sigma_x^2 + \sigma_y^2} \f$.
     */ 
   realT sigma()
   {
      return this->p[3];
   }
      
};


///\ref levmarInterface fitter structure for the symmetric Gaussian.
/** \ingroup gaussian_peak_fit
  *
  */
template<typename _realT>
struct gaussian1D_fitter
{
   typedef _realT realT;
   
   static const int nparams = 4;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2Fit<realT> * arr = (array2Fit<realT> *) adata;
   
      for(size_t i=0;i<arr->nx; i++)
      {
         hx[i] = func::gaussian<realT>(i,p[0],p[1], p[2], p[3]) - arr->data[i];
      }
      
   }
   
};

///Alias for the fitGaussian1D type fitting the gaussian.
/** \ingroup gaussian_peak_fit
  */
//template<typename realT>
//using fitGaussian1D = mx::math::fit::fitGaussian1D<mx::math::fit::gaussian1D_fitter<realT>>;













//forward
template<typename realT>
struct array2FitGaussian2D;

//forward
template<typename _realT>
struct gaussian2D_sym_fitter;

//forward
template<typename _realT>
struct gaussian2D_gen_fitter;

template<typename _realT>
struct gaussian2D_gen_fitter_bgfixed;

///Class to manage fitting a 2D Gaussian to data via the \ref levmarInterface
/** Can fit either the symmetric Gaussian, or the rotated asymmetric Gaussian.  Individual parameters can be fixed, 
  * except in the asymmetric case \f$\sigma_x, \sigma_y, \theta\f$ must currently be fixed together.
  * 
  * In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires the following in fitterT
  * \code
  * static const int nparams = 7; //where the number 7 is replaced by the number of parameters that fitterT expects to fit. 
  * 
  * void paramNormalizer( array2FitGaussian2D * arr,
  *                       realT *, 
  *                       int); //A function which converts from input parameters to fitting parameters.  May do nothing.
  * 
  * \endcode
  * 
  *
  * \tparam fitterT a type meeting the above requirements.
  *
  * \ingroup gaussian_peak_fit
  *
  */
template<typename fitterT>
class fitGaussian2D : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

   //static const int nparams = fitterT::nparams;

   array2FitGaussian2D<realT> arr;
   
   void initialize()
   {
      if(fitterT::maxNparams == 5)
      {
         arr.setSymmetric();
      }
      else
      {
         arr.setGeneral();
      }
      this->allocate_params(arr.nparams());
      this->adata = &arr;      
   }
   
   fitGaussian2D()
   {
      initialize();
   }
      
   ~fitGaussian2D()
   {
   }
   
   /// Set whether each parameter is fixed.
   /** Sets the parameter indices appropriately.
     */
   void setFixed( bool G0,      ///< [in] if true, then G0 will be not be part of the fit
                  bool G,       ///< [in] if true, then G will be not be part of the fit
                  bool x0,      ///< [in] if true, then x0 will be not be part of the fit
                  bool y0,      ///< [in] if true, then y0 will be not be part of the fit
                  bool sigma_x, ///< [in] if true, then sigma_x will be not be part of the fit
                  bool sigma_y, ///< [in] if true, then sigma_y will be not be part of the fit
                  bool theta    ///< [in] if true, then theta will be not be part of the fit
                )
   {
      arr.setFixed(G0, G, x0, y0, sigma_x, sigma_y, theta);
      this->allocate_params(arr.nparams());
   }
   
   /// Set whether each parameter is fixed.
   /** Sets the parameter indices appropriately.
     */
   void setFixed( bool G0,   ///< [in] if true, then G0 will be not be part of the fit
                  bool G,    ///< [in] if true, then G will be not be part of the fit
                  bool x0,   ///< [in] if true, then x0 will be not be part of the fit
                  bool y0,   ///< [in] if true, then y0 will be not be part of the fit
                  bool sigma ///< [in] if true, then sigma will be not be part of the fit
                )
   {
      arr.setFixed(G0, G, x0, y0, sigma, sigma, sigma);
      this->allocate_params(arr.nparams());
   }
   
   ///Set the initial guess for a symmetric Gaussian.
   /** Also works for the general case, setting the same width in both directions.
     */
   void setGuess( realT G0,    ///< [in] the constant background level
                  realT G,     ///< [in] the peak scaling
                  realT x0,    ///< [in] the center x-coordinate
                  realT y0,    ///< [in] the center y-coordinate
                  realT sigma  ///< [in] the width parameter
                )
   {
      arr.G0(this->p,G0);
      arr.G(this->p,G);
      arr.x0(this->p,x0);
      arr.y0(this->p,y0);
      arr.sigma(this->p,sigma);
      
      
   }
   
   ///Set the initial guess for the general Gaussian.
   void setGuess( realT G0,      ///< [in] the constant background level
                  realT G,       ///< [in] the peak scaling
                  realT x0,      ///< [in] the center x-coordinate
                  realT y0,      ///< [in] the center y-coordinate
                  realT sigma_x, ///< [in] the width parameter in the rotated x-direction (the long axis) 
                  realT sigma_y, ///< [in] the width parameter in the rotated y-direction
                  realT theta    ///< [in] the angle of the long axis (always sigma_x)
                )
   {
      
      arr.G0(this->p,G0);
      arr.G(this->p,G);
      arr.x0(this->p,x0);
      arr.y0(this->p,y0);
      arr.sigma_x(this->p,sigma_x);
      arr.sigma_y(this->p,sigma_y);
      arr.theta(this->p,theta);      
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
      fitter.paramNormalizer(&arr, this->p, 1);
      
      levmarInterface<fitterT>::fit();
      
      fitter.paramNormalizer(&arr, this->p, -1);
      fitter.paramNormalizer(&arr, this->init_p, -1); //The normalized version is stored, so fix it before possible output.
      
      return 0;
   }
     
   ///Get the current value of G0, the constant.
   /**
     * \returns the current value of G0.
     */ 
   realT G0()
   {
      return arr.G0(this->p);
   }

   ///Get the peak scaling.
   realT G()
   {
      return arr.G(this->p);
   }

   ///Get the center x-coordinate
   realT x0()
   {
      return arr.x0(this->p);
   }
   
   ///Get the center y-coordinate
   realT y0()
   {
      return arr.y0(this->p);
   }
   
   ///Return the width parameter
   /** As described for the symmetric Gaussian.
     *
     * For the general Gaussian, this returns \f$ \sigma = \sqrt{ \sigma_x^2 + \sigma_y^2} \f$.
     */ 
   realT sigma()
   {
      return arr.sigma(this->p);     
   }
   
   ///Return the full-width at half maximum
   /** This is a simple scaling of the sigma() result 
     */
   realT fwhm()
   {
      return func::sigma2fwhm(arr.sigma(this->p));
   }
   
   ///Return the width parameter on the long axis.
   realT sigma_x()
   {
      return arr.sigma_x(this->p);
   }
   
   ///Return the width parameter on the short axis.
   realT sigma_y()
   {
      return arr.sigma_y(this->p);
   }
   
   ///Return the orientation of the long axis.
   realT theta()
   {
      return arr.theta(this->p);
   }
   
};

///Wrapper for a native array to pass to \ref levmarInterface, with 2D Gaussian details.
/** Supports fixing G0, G, x0, and y0 independently.  The shape and orientation can be fixed, but
  * for the general form, sigma_x, sigma_y, and theta can only be fixed together. 
  * \ingroup gaussian_peak_fit
  */
template<typename realT>
struct array2FitGaussian2D
{
   realT * data {nullptr}; ///< Pointer to the array
   size_t nx {0}; ///< X dimension of the array
   size_t ny {0}; ///< Y dimension of the array
   
   realT m_G0 {0};
   realT m_G {0};
   realT m_x0 {0};
   realT m_y0 {0};
   realT m_sigma_x {0};
   realT m_sigma_y {0};
   realT m_theta {0};
   
   realT m_a {0};
   realT m_b {0};
   realT m_c {0};
   
   realT m_sigma {0};
   
   int m_G0_idx {0};
   int m_G_idx {1};
   int m_x0_idx {2};
   int m_y0_idx {3};
   int m_sigma_x_idx {4}; ///< Index of sigma_x in the parameters.  Re-used for a.
   int m_sigma_y_idx {5}; ///< Index of sigma_y in the parameters.  Re-used for b.
   int m_theta_idx {6};  ///< Index of theta in the parameters.  Re-used for c.
   
   realT m_sigma_idx {4};  ///< Index of sigma in the symmetric case
   
   int m_nparams {7};
   int m_maxNparams {7};
   
   void setSymmetric()
   {
      m_maxNparams = 5;
   }
   
   void setGeneral()
   {
      m_maxNparams = 7;
   }
   
   /// Set whether each parameter is fixed.
   /** Sets the parameter indices appropriately.
     */
   void setFixed( bool G0, 
                  bool G, 
                  bool x0, 
                  bool y0, 
                  bool sigma_x, 
                  bool sigma_y,
                  bool theta
                )
   {
      int idx = 0;
      
      if(G0) m_G0_idx = -1;
      else m_G0_idx = idx++;
      
      if(G) m_G_idx = -1;
      else m_G_idx = idx++;
   
      if(x0) m_x0_idx = -1;
      else m_x0_idx = idx++;
      
      if(y0) m_y0_idx = -1;
      else m_y0_idx = idx++;
      
      if(m_maxNparams == 5)
      {
         if(sigma_x) m_sigma_idx = -1;
         else m_sigma_idx = idx++;
      }
      else if(sigma_x && sigma_y && theta)
      {
         m_sigma_x_idx = -1;
         m_sigma_y_idx = -1;
         m_theta_idx = -1;
      }
      else
      {
         if(sigma_x || sigma_y || theta)
         {
            mxError("array2FitGaussian2D::setFixed", MXE_NOTIMPL, "cannot fix sigma_x, sigma_y, and theta separately");
         }
         
         m_sigma_x_idx = idx++;
         m_sigma_y_idx = idx++;
         m_theta_idx = idx++;
      }
      
      m_nparams = idx;
   }
   
   realT G0( realT * p )
   {
      if( m_G0_idx < 0 )
      {
         return m_G0;
      }
      else
      {
         return p[m_G0_idx];
      }
   }
   
   void G0( realT * p,
            realT nG0
          )
   {
      if( m_G0_idx < 0 )
      {
         m_G0 = nG0;
      }
      else
      {
         p[m_G0_idx]=nG0;
      }
   }
   
   realT G( realT * p )
   {
      if( m_G_idx < 0 )
      {
         return m_G;
      }
      else
      {
         return p[m_G_idx];
      }
   }
   
   void G( realT * p,
           realT nG
         )
   {
      if( m_G_idx < 0 )
      {
         m_G = nG;
      }
      else
      {
         p[m_G_idx] = nG;
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
   
   
   realT sigma_x( realT * p )
   {
      if( m_sigma_x_idx < 0 )
      {
         return m_sigma_x;
      }
      else
      {
         return p[m_sigma_x_idx];
      }
   }
   
   void sigma_x( realT * p,
                 realT nsigma_x
               )
   {
      if( m_sigma_x_idx < 0 )
      {
         m_sigma_x = nsigma_x;
      }
      else
      {
         p[m_sigma_x_idx] = nsigma_x;
      }
   }
   
   realT sigma_y( realT * p )
   {
      if( m_sigma_y_idx < 0 )
      {
         return m_sigma_y;
      }
      else
      {
         return p[m_sigma_y_idx];
      }
   }
   
   void sigma_y( realT * p,
                 realT nsigma_y
               )
   {
      if( m_sigma_y_idx < 0 )
      {
         m_sigma_y = nsigma_y;
      }
      else
      {
         p[m_sigma_y_idx] = nsigma_y;
      }
   }
   
   realT theta( realT * p )
   {
      if( m_theta_idx < 0 )
      {
         return m_theta;
      }
      else
      {
         return p[m_theta_idx];
      }
   }
   
   void theta( realT * p,
               realT ntheta
             )
   {
      if( m_theta_idx < 0 )
      {
         m_theta = ntheta;
      }
      else
      {
         p[m_theta_idx] = ntheta;
      }
   }
   
   realT sigma( realT * p )
   {
      if( m_sigma_idx < 0 )
      {
         return m_sigma;
      }
      else
      {
         return p[m_sigma_idx];
      }
   }
   
   void sigma( realT * p,
               realT nsigma
             )
   {
      if( m_sigma_idx < 0 )
      {
         m_sigma = nsigma;
      }
      else
      {
         p[m_sigma_idx] = nsigma;
      }
   }
   
   
   realT a( realT * p )
   {
      //This aliases sigma_x after param normalization
      if( m_sigma_x_idx < 0 )
      {
         return m_a;
      }
      else
      {
         return p[m_sigma_x_idx];
      }
   }
   
   void a( realT * p,
           realT na
         )
   {
      //This aliases sigma_x after param normalization
      if( m_sigma_x_idx < 0 )
      {
         m_a = na;
      }
      else
      {
         p[m_sigma_x_idx] = na;
      }
   }
   
   realT b( realT * p )
   {
      //This aliases sigma_y after param normalization
      if( m_sigma_y_idx < 0 )
      {
         return m_b;
      }
      else
      {
         return p[m_sigma_y_idx];
      }
   }
   
   void b( realT * p,
           realT nb
         )
   {
      //This aliases sigma_y after param normalization
      if( m_sigma_y_idx < 0 )
      {
         m_b = nb;
      }
      else
      {
         p[m_sigma_y_idx] = nb;
      }
   }
   
   realT c( realT * p )
   {
      //This aliases theta after param normalization
      if( m_theta_idx < 0 )
      {
         return m_c;
      }
      else
      {
         return p[m_theta_idx];
      }
   }
   
   void c( realT * p,
           realT nc
         )
   {
      //This aliases theta after param normalization
      if( m_theta_idx < 0 )
      {
         m_c = nc;
      }
      else
      {
         p[m_theta_idx] = nc;
      }
   }
   
   int nparams()
   {
      return m_nparams;
   }
   
   
};

///\ref levmarInterface fitter structure for the symmetric Gaussian.
/** \ingroup gaussian_peak_fit
  *
  */
template<typename _realT>
struct gaussian2D_sym_fitter
{
   typedef _realT realT;
   
   static const int maxNparams = 5;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2Fit<realT> * arr = (array2Fit<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      for(int j=0;j<arr->ny; j++)
      {
         for(int i=0;i<arr->nx;i++)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = func::gaussian2D<realT>(i,j,p[0],p[1], p[2], p[3], p[4]) - arr->data[idx_mat];
            
            idx_dat++;
         }
      }
      
   }
   
   ///Does nothing in this case.
   void paramNormalizer(realT * p, int dir)
   {
      return;
   }
};

///\ref levmarInterface fitter structure for the general elliptical Gaussian.
/** \ingroup gaussian_peak_fit
  *
  */
template<typename _realT>
struct gaussian2D_gen_fitter
{
   typedef _realT realT;
   
   static const int maxNparams = 7;
   
   static void func(realT *p, realT *hx, int m, int n, void *adata)
   {
      array2FitGaussian2D<realT> * arr = (array2FitGaussian2D<realT> *) adata;
   
      size_t idx_mat, idx_dat;

      realT G0 = arr->G0(p); //0
      realT G = arr->G(p); //1
      realT x0 = arr->x0(p); //2
      realT y0 = arr->y0(p); //3
      realT a = arr->a(p); //4
      realT b = arr->b(p); //5
      realT c = arr->c(p); //6
      
      //Check for positive-definiteness of {{a b}{b c}}
      if( a*c - b*b <= 0 || a <= 0 || c <= 0 || a+c <= 2*fabs(b))
      {
         idx_dat = 0;
         //If it's not positive-definite, then we just fill in with the value of the image itself.
         for(int j=0;j<arr->ny; ++j)
         {
            for(int i=0;i<arr->ny; ++i)
            { 
               idx_mat = i+j*arr->nx;
   
               hx[idx_dat] = arr->data[idx_mat];
                        
               ++idx_dat;
            }
         }         
         return;
      }
      
      //If positive-definite, now actually calculate
      idx_dat = 0;
   
      for(int j=0;j<arr->ny; ++j)
      {
         for(int i=0;i<arr->nx; ++i)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = func::gaussian2D<realT>(i,j, G0, G, x0, y0, a, b, c) - arr->data[idx_mat];
                        
            ++idx_dat;
         }
      }
   }
      
   void paramNormalizer( array2FitGaussian2D<realT> * arr,
                         realT * p, 
                         int dir
                       )
   {
      //Prepare for fit
      if(dir == 1)
      {
         realT na, nb, nc;
         realT sx = arr->sigma_x(p);
         realT sy = arr->sigma_y(p);
         realT th = arr->theta(p);
         
         func::gaussian2D_rot2gen(na,nb,nc, sx, sy, th);
         arr->a(p, na);
         arr->b(p, nb);
         arr->c(p, nc);
         return;
      }
      
      //Convert after fit
      if(dir == -1)
      {
         realT sx, sy, th;
         realT na = arr->a(p);
         realT nb = arr->b(p);
         realT nc = arr->c(p);
         func::gaussian2D_gen2rot(sx, sy, th, na, nb, nc);

         arr->sigma_x(p, sx);
         arr->sigma_y(p, sy);
         arr->theta(p, th);
         return;
      }
   }
};



///Alias for the fitGaussian2D type fitting the symmetric gaussian.
/** \ingroup gaussian_peak_fit
  */
template<typename realT>
using fitGaussian2Dsym = mx::math::fit::fitGaussian2D<mx::math::fit::gaussian2D_sym_fitter<realT>>;

///Alias for the fitGaussian2D type fitting the general elliptical gaussian.
/** \ingroup gaussian_peak_fit
  */
template<typename realT>
using fitGaussian2Dgen = mx::math::fit::fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<realT>>;


///Form an estimate of the parameters of an elliptical Gaussian from a 2D image.
/** Note that this assumes that there is no constant value (i.e. it is zero-ed).
  *
  *  \ingroup gaussian_peak_fit 
  */
template<typename realT>
int guessGauss2D_ang( realT & Ag, ///< [out] estimate of the peak
                      realT & xg, ///< [out] estimate of the x-coordinate of the peak
                      realT & yg, ///< [out] estimate of the y-coordinate of the peak
                      realT & xFWHM, ///< [out] estimate of the x-FWHM
                      realT & yFWHM, ///< [out] estimate of the y-FWHM
                      realT & angG, ///< [out] estimate of the angle of the ellipse
                      mx::improc::eigenImage<realT> & im,  ///< [in] the image with an elliptical gaussian
                      realT maxWidth,   ///< [in] the width of the box to search for the maximum
                      realT widthWidth,  ///< [in] the radius of the circle to search for the widths
                      realT nAngs,  ///< [in] the number of angles at which to search for the widths
                      realT xg0,  ///< [in] an initial guess at the x-coordinate of the peak
                      realT yg0  ///< [in] an initial guess at the y-coordinate of the peak
                    )
{
   xg = xg0;
   yg = yg0;
   Ag = im( (int) xg, (int) yg);
   for(int i =0; i< 2*maxWidth+1; ++i)
   {
      for(int j=0; j< 2*maxWidth+1; ++j)
      {
         if(  im( (int)(xg0 - maxWidth + i), (int) (yg0 - maxWidth + j)) > Ag )
         {
            Ag = im( (int) (xg0 - maxWidth + i), (int) (yg0 - maxWidth + j));
            xg = xg0 - maxWidth + i;
            yg = yg0 - maxWidth + j;
         }
      }
   }
   
   realT dAng = two_pi<realT>()/nAngs;
   //std::vector<realT> dist(nAngs);
   
   realT c, s;
   
   realT maxD = 0;
   int maxDidx = 0;
   realT minD = widthWidth;
   int minDidx = 0;
   
   for(int i=0; i < nAngs; ++i)
   {
      c = cos(i*dAng);
      s = sin(i*dAng);
      
      for(int j=0; j < widthWidth; ++j)
      {
         if( im( (int)(xg + j*c), (int)(yg + j*s)) <= 0.5*Ag )
         {
            //dist[i] = j;
            
            if(j > maxD) 
            {
               maxD = j;
               maxDidx = i;
            }
            
            if(j < minD)
            {
               minD = j;
               minDidx = i;
            }
            break;
         }
      }
   }
   
   //Take minang and move it by 90 degrees
   realT minang = fmod(minDidx * dAng - 0.5*pi<realT>(), pi<realT>());
   if(minang < 0) minang = fmod(minang + two_pi<realT>(), pi<realT>());
   
   realT maxang = fmod(maxDidx * dAng, pi<realT>());
   
   //Now average
   angG = 0.5*(minang + maxang);
   
   xFWHM = 2*maxD;
   yFWHM = 2*minD;
   
   return 0; ///\returns 0 if successful
}


} //namespace fit
} //namespace math

} //namespace mx

#endif //fitGaussian_hpp

