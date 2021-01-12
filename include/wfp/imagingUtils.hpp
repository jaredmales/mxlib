/** \file imagingUtils.hpp
  * \brief Utilities for modeling image formation
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
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

#ifndef imagingUtils_hpp
#define imagingUtils_hpp

#include <cmath>

#include "../math/constants.hpp"
#include "../mxError.hpp"

#include "imagingArray.hpp"


namespace mx
{

namespace wfp
{
   
///Calculate the angular plate scale (radians per pixel) of an image after propagation by FFT.
/** 
  *
  * \returns the platescale of the wavefront after propagation by FFT, in radians per pixel.
  * 
  * \tparam realT a real floating point type 
  *
  * \ingroup imaging
  */   
template<typename realT>
realT fftPlateScale( size_t pixels, ///< [in] the linear dimension of the FFT (including 0 pad, etc.)
                     realT metersPerPixel, ///< [in] the scale of the input wavefront [m/pix] 
                     realT lambda ///< [in] the wavelength of the wavefront [m]
                   ) 
{
   return (lambda/metersPerPixel) * (1./pixels);
}
   
/// Fill in an Eigen-like array with a circular pupil mask.
/** Sets any pixel which is at rad \<= r \< rad+(1.0/overscan) pixels to rho = 1,
  * 
  * 
  * \retval 0 on success
  * \retval -1 on error 
  * 
  * \ingroup imaging
  */  
template<class arrayT> 
int circularPupil( arrayT & m,  ///< [in/out] is the allocated Array.  Dimensions are used to create the pupil.
                   typename arrayT::Scalar eps=0, ///< [in] [optional] is the central obscuration.  0-1, default is 0. 
                   typename arrayT::Scalar rad=0, ///< [in] [optional] is the desired radius. If rad \<= 0, then the maximum radius based on dimensions of m is used. Default is 0.
                   typename arrayT::Scalar overscan = 0  ///< [in] [optional] overscan in fractional pixels, to include partial pixels on the edge. Default is 0.
                 )
{
   
   if( eps < 0)
   {
      mxError("circularPupil", MXE_INVALIDARG, "Central obscuration can not be < 0." );
      return -1;
   }
   
   if(eps > 1)
   {
      mxError("circularPupil", MXE_INVALIDARG, "Central obscuration can not be > 1." );
      return -1;
   }
   
   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   typename arrayT::Scalar r;
   typename arrayT::Scalar xc = 0.5*(l0-1);
   typename arrayT::Scalar yc = 0.5*(l1-1);
   
   if(rad <= 0) rad = 0.5*std::min(l0-1, l1-1);
   
   for(size_t i=0; i < l0; i++)
   {
      for(size_t j=0; j < l1; j++)
      {
         r = std::sqrt( std::pow(i-xc, 2) + std::pow(j-yc, 2) );
         
         if(r <= rad+( overscan ) && r >= eps*rad) m(i,j) = 1;
         else m(i,j) = 0;
      }
   }
   
   return 0;
}

///Draw a line in an image
/** \todo should handle width much more intelligently, this only works for ~45 degree lines.
  * 
  * \tparam arrayT is an Eigen-like array with public typedef Scalar
  */ 
template<class arrayT> 
void drawLine( arrayT & im,  ///< [in/out] The input image, modified.
               typename arrayT::Scalar x0, ///< [in] the x value, relative to image center, of the starting point
               typename arrayT::Scalar y0, ///< [in] the y value, relative to image center, of the starting point
               typename arrayT::Scalar x1,  ///< [in] the x value, relative to image center, of the end point
               typename arrayT::Scalar y1, ///< [in] the y value, relative to image center, of the end point
               typename arrayT::Scalar halfWidth  ///< [in] the half-width of the line.
             )
{

   int d1 = im.rows();
   int d2 = im.cols();
   
   typename arrayT::Scalar xc, yc;
   
   xc = 0.5*(im.rows()-1);
   yc = 0.5*(im.cols()-1);
   
   
   typename arrayT::Scalar m = (y1-y0)/(x1-x0);
   int y;
   
   if(x1 > x0)
   {
      for(int x = x0; x<=x1; ++x)
      {
         y = y0 + (x-x0)*m;
      
         for(int i=0; i<= halfWidth; ++i)
         {
            if( x+xc >= 0 && x+xc < d1 && y+yc+i >= 0 && y+yc+i < d2) im(x+xc,y+yc+i) = 0;
            if( x+xc >= 0 && x+xc < d1 && y+yc-i >= 0 && y+yc-i < d2) im(x+xc,y+yc-i) = 0;
         }
      }
   }
   else
   {
      for(int x = x1; x<=x0; ++x)
      {
         y = y0 + (x-x0)*m;
      
         for(int i=0; i<= halfWidth; ++i)
         {
            if( x+xc >= 0 && x+xc < d1 && y+yc+i >= 0 && y+yc+i < d2) im(x+xc,y+yc+i) = 0;
            if( x+xc >= 0 && x+xc < d1 && y+yc-i >= 0 && y+yc-i < d2) im(x+xc,y+yc-i) = 0;
         }
      }
      
   }
   
}
   

   
///Create a complex pupil plane wavefront from a real amplitude mask.
/** The input real amplitude mask is placed in the center of a 0-padded complex array.
  *
  * \param [out] complexPupil the complex pupil plane wavefront
  * \param [in] realPupil a real amplitude mask.
  * \param [in] wavefrontSizePixels the desired size of the ouput wavefront, should be at least as big as realPupil
  * 
  * \ingroup imaging
  */ 
template<typename arrayOutT, typename arrayInT>
void makeComplexPupil( arrayOutT & complexPupil, 
                       const arrayInT & realPupil, 
                       int wavefrontSizePixels)
{
   
   complexPupil.resize(wavefrontSizePixels, wavefrontSizePixels);
   complexPupil.set(typename arrayOutT::Scalar(0,0));
     
   //Lower-left corner of insertion region
   int bl = 0.5*(complexPupil.rows()-1) - 0.5*(realPupil.rows()-1.);
   
   for(int i=0; i< realPupil.cols(); ++i)
   {
      for(int j=0; j < realPupil.rows(); ++j)
      {
         complexPupil(bl+i, bl+j) = typename arrayOutT::Scalar(realPupil(i,j),0); //*exp( typename arrayOutT::Scalar(0,1)); 
      }
   }
   //complexPupil.block(bl, bl, realPupil.rows(), realPupil.rows()) = realPupil*std::complex<realT>(1,0);

}


///Create a complex wavefront from a real amplitude and a real phase.
/** The wavefront is placed in the center of a 0-padded complex array.
  *
  * \param [out] complexWavefront the complex pupil plane wavefront
  * \param [in] realAmplitude is the real-valued amplitude.
  * \param [in] realPhase is the real-valued phase in radians, same size as realAmplitude
  * \param [in] wavefrontSizePixels the desired size of the ouput wavefront, should be at least as big as the real arrays
  * 
  * \ingroup imaging
  */ 
template<typename arrayOutT, typename arrayInT>
void makeComplexPupil( arrayOutT & complexWavefront, 
                       const arrayInT & realAmplitude, 
                       const arrayInT & realPhase, 
                       int wavefrontSizePixels)
{
   
   complexWavefront.resize(wavefrontSizePixels, wavefrontSizePixels);
   complexWavefront.set(typename arrayOutT::Scalar(0,0));
     
   //Lower-left corner of insertion region
   int bl = 0.5*(complexWavefront.rows()-1) - 0.5*(realAmplitude.rows()-1.);
   
   for(int i=0; i< realAmplitude.cols(); ++i)
   {
      for(int j=0; j < realAmplitude.rows(); ++j)
      {
         complexWavefront(bl+i, bl+j) = realAmplitude(i,j)*exp(  (typename arrayOutT::Scalar(0,1)) * realPhase(i,j)); 
      }
   }
}


///Apply a tilt to a wavefront 
/**
  * \param complexWavefront [in/out] the complex wavefront to tilt, will be modified on output
  * \param xTilt [input] the amount of tilt in the x direction, in pixels 
  * \param yTilt [input] the amount of tilt in the y direction, in pixels 
  * 
  * \ingroup imaging
  */
template<typename wavefrontT>
void tiltWavefront( wavefrontT & complexWavefront, 
                    typename wavefrontT::Scalar::value_type xTilt, 
                    typename wavefrontT::Scalar::value_type yTilt)
{
   typedef typename wavefrontT::Scalar complexT;
   typedef typename wavefrontT::Scalar::value_type realT;
   
   realT pi = math::pi<realT>();
   
   int wfsSizeX = complexWavefront.cols();
   int wfsSizeY = complexWavefront.rows();
   
   realT xCen = 0.5*wfsSizeX;
   realT yCen = 0.5*wfsSizeY;

   realT argX = 2.0*pi/(wfsSizeX-1.0);
   realT argY = 2.0*pi/(wfsSizeY-1.0);
  
   for(int ii=0; ii < wfsSizeX; ++ii)
   {
      for(int jj=0; jj < wfsSizeY; ++jj)
      {     
         complexWavefront(ii,jj) = complexWavefront(ii,jj)*exp( complexT( (realT)0., argX*xTilt*(ii-xCen)+argY*yTilt*(jj-yCen)));
      }
   }
}

template< typename imageT1, typename imageT2>
void extractBlock(imageT1 & im,
                  int imX0,
                  int imXsz,
                  int imY0,
                  int imYsz,
                  imageT2 & wf,
                  int wfX0,
                  int wfY0)
{
   int im_rows = im.cols();
   
   int wf_rows = wf.cols();
   
   typedef typename imageT1::Scalar dataT;
   
   dataT * im_data;
   dataT * wf_data;
   
   

   for(int j =0; j< imYsz; ++j)
   {
      im_data = &im.data()[imX0 + (imY0+j)*im_rows];
      wf_data = &wf.data()[wfX0 + (wfY0+j)*wf_rows];
      
      memcpy( im_data, wf_data, sizeof(dataT)*imXsz);
   }
}

template< typename realImageT,
          typename complexImageT>
void extractIntensityImage(realImageT & im,
                           int imX0,
                           int imXsz,
                           int imY0,
                           int imYsz,
                           complexImageT & wf,
                           int wfX0,
                           int wfY0)
{
   int im_rows = im.cols();
   
   int wf_rows = wf.cols();
   
   typename realImageT::Scalar * im_data;
   typename complexImageT::Scalar * wf_data;
   
   for(int j =0; j< imXsz; ++j)
   {
      im_data = &im.data()[imX0 + (imY0+j)*im_rows];
      wf_data = &wf.data()[wfX0 + (wfY0+j)*wf_rows];
      for(int i=0; i<imYsz; ++i)
      {
         im_data[i] = norm(wf_data[i]);
      }
   }
}

template< typename realImageT,
          typename complexImageT>
void extractIntensityImageAccum(realImageT & im,
                           int imX0,
                           int imXsz,
                           int imY0,
                           int imYsz,
                           complexImageT & wf,
                           int wfX0,
                           int wfY0)
{
   int im_rows = im.cols();
   
   int wf_rows = wf.cols();
   
   typename realImageT::Scalar * im_data;
   typename complexImageT::Scalar * wf_data;
   
   for(int j =0; j< imXsz; ++j)
   {
      im_data = &im.data()[imX0 + (imY0+j)*im_rows];
      wf_data = &wf.data()[wfX0 + (wfY0+j)*wf_rows];
      for(int i=0; i<imYsz; ++i)
      {
         im_data[i] += norm(wf_data[i]);
      }
   }
}

} //namespace wfp
} //namespace mx

#endif //__imagingUtils_hpp__

