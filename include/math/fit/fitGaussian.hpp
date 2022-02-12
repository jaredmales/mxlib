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

#include "array2FitGaussian2D.hpp"

#include "../../mxError.hpp"
#include "levmarInterface.hpp"
#include "../func/gaussian.hpp"
#include "../constants.hpp"

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

extern template struct gaussian1D_fitter<float>;
extern template struct gaussian1D_fitter<double>;


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
      
      return levmarInterface<fitterT>::fit();
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

extern template class fitGaussian1D<float>;
extern template class fitGaussian1D<double>;


///Alias for the fitGaussian1D type fitting the gaussian.
/** \ingroup gaussian_peak_fit
  */
//template<typename realT>
//using fitGaussian1D = mx::math::fit::fitGaussian1D<mx::math::fit::gaussian1D_fitter<realT>>;












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
  * \test Scenario: Verify direction and accuracy of various image shifts \ref tests_improc_imageTransforms_imageShift "[test doc]"
  *
  */
template<typename fitterT>
class fitGaussian2D : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::realT realT;

   array2FitGaussian2D<realT> arr; ///< Data array to pass to the levmar library.  Contains the actual data plus the parameters if fixed.
   
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
   void setArray( realT *data, ///< [in] The 2D array of data to fit.
                  int nx,      ///< [in] the x size of the data
                  int ny       ///< [in] the y size of the data
                )
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = ny;
      arr.mask = nullptr;
      this->n = nx*ny;
      
   }
   
   ///Set the data aray, with a mask.
   void setArray( realT *data, ///< [in] The 2D array of data to fit.
                  int nx,      ///< [in] the x size of the data
                  int ny,      ///< [in] the y size of the data
                  realT *mask  ///< [in] Array of same size as data.  Any 0 pixels in this array will be excluded from the fit.
                )
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = ny;
      arr.mask = mask;
      
      this->n = 0;
      int idx_mat;
      for(int j=0;j<nx; ++j)
      {
         for(int i=0;i<ny; ++i)
         { 
            idx_mat = i+j*nx;
   
            this->n += arr.mask[idx_mat];
         }
      }
            
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



///\ref levmarInterface fitter structure for the symmetric Gaussian.
/** \ingroup gaussian_peak_fit
  *
  * \test Scenario: Verify direction and accuracy of various image shifts \ref tests_improc_imageTransforms_imageShift "[test doc]"
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
   void paramNormalizer( array2FitGaussian2D<realT> * arr, 
                         realT * p, 
                         int dir
                       )
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
         if(arr->mask == nullptr)
         {
            for(int j=0;j<arr->ny; ++j)
            {
               for(int i=0;i<arr->ny; ++i)
               { 
                  idx_mat = i+j*arr->nx;
                  hx[idx_dat] = arr->data[idx_mat];
                  ++idx_dat;
               }
            }
         }
         else
         {
            for(int j=0;j<arr->ny; ++j)
            {
               for(int i=0;i<arr->ny; ++i)
               { 
                  idx_mat = i+j*arr->nx;
                  if(arr->mask[idx_mat]==0) continue;
                  hx[idx_dat] = arr->data[idx_mat];
                  ++idx_dat;
               }
            }
         }   
         return;
      }
      
      
      //If positive-definite, now actually calculate
      idx_dat = 0;
   
      if(arr->mask == nullptr)
      {
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
      else
      {
         for(int j=0;j<arr->ny; ++j)
         {
            for(int i=0;i<arr->nx; ++i)
            { 
               idx_mat = i+j*arr->nx;
               
               if(arr->mask[idx_mat] == 0) 
               {
                  continue;
               }
               else
               {
                  hx[idx_dat] = func::gaussian2D<realT>(i,j, G0, G, x0, y0, a, b, c) - arr->data[idx_mat];
               
                  ++idx_dat;
               }
            }
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

extern template class fitGaussian2D<mx::math::fit::gaussian2D_sym_fitter<float>>;
extern template class fitGaussian2D<mx::math::fit::gaussian2D_sym_fitter<double>>;
extern template class fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<float>>;
extern template class fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<double>>;

///Alias for the fitGaussian2D type fitting the symmetric gaussian.
/** \ingroup gaussian_peak_fit
  */
template<typename realT>
using fitGaussian2Dsym = fitGaussian2D<mx::math::fit::gaussian2D_sym_fitter<realT>>;

///Alias for the fitGaussian2D type fitting the general elliptical gaussian.
/** \ingroup gaussian_peak_fit
  */
template<typename realT>
using fitGaussian2Dgen = fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<realT>>;


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
   
   realT dAng = math::two_pi<realT>()/nAngs;
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
   realT minang = fmod(minDidx * dAng - 0.5*math::pi<realT>(), math::pi<realT>());
   if(minang < 0) minang = fmod(minang + math::two_pi<realT>(), math::pi<realT>());
   
   realT maxang = fmod(maxDidx * dAng, math::pi<realT>());
   
   //Now average
   angG = 0.5*(minang + maxang);
   
   xFWHM = 2*maxD;
   yFWHM = 2*minD;
   
   return 0; ///\returns 0 if successful
}

extern template
int guessGauss2D_ang<float>( float & Ag, 
                             float & xg,        
                             float & yg, 
                             float & xFWHM,
                             float & yFWHM,
                             float & angG,  
                             mx::improc::eigenImage<float> & im,  
                             float maxWidth,  
                             float widthWidth,
                             float nAngs, 
                             float xg0,  
                             float yg0  
                           );

extern template
int guessGauss2D_ang<double>( double & Ag, 
                              double & xg,        
                              double & yg, 
                              double & xFWHM,
                              double & yFWHM,
                              double & angG,  
                              mx::improc::eigenImage<double> & im,  
                              double maxWidth,  
                              double widthWidth,
                              double nAngs, 
                              double xg0,  
                              double yg0  
                            );

} //namespace fit
} //namespace math

} //namespace mx

#endif //fitGaussian_hpp

