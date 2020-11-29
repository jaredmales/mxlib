/** \file imageXCorrDiscrete.hpp
  * \brief A class to register images using the discrete cross correlation.
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef imageXCorrDiscrete_hpp
#define imageXCorrDiscrete_hpp

#include "../mxError.hpp"
#include "../math/fit/fitGaussian.hpp"
#include "imageUtils.hpp"

namespace mx
{
namespace improc
{
   
enum class xcorrPeakMethod{ centroid, gaussfit, interp };

/// Find the optimum shift to align two images using the discrete cross correlation.
/** The reference image must be smaller than the target image.  Both the reference image and the section of
  * target image being analyzed are mean subtracted and variance normalized.  An optional mask can be supplied,
  * which limits the pixels used for the mean and variance calculation.
  * 
  * Typical usage will be to set the mask, then the reference, then repeatedly call operator() to 
  * determine the shifts for a sequence of imaages.  No new heap allocations take place on these calls
  * to operator(), and the reference image is not re-normalized on each call.
  * 
  * The shift is reported in pixels such that if the mxlib imageShift function is used
  * to shift the input image by the negative of the shifts, it will align with the 
  * reference at the center of the array.
  *
  * Three peak finding methods are provided.  xcorrPeakMethod::centroid uses center of light, 
  * xcorrPeakMethod::centroid uses Gaussian centroiding, and xcorrPeakMethod::interp uses interpolation
  * to find the peak to a given tolerance.
  * 
  * \tparam _ccImT is the Eigen-like array type used for image processing.  See typedefs.
  * 
  * \ingroup image_reg
  */ 
template<class _ccImT>
class imageXCorrDiscrete
{
public:
   typedef _ccImT ccImT; ///< the Eigen-like array type used for image processing
   typedef typename _ccImT::Scalar Scalar; ///< the scalar type of the image type
   
protected:
   
   /** \name Working Memory
     * @{
     */ 
   
   ccImT m_refIm;  ///< The normalized reference image.
   
   ccImT m_maskIm; ///< Mask image to use, may be needed for proper normalization even if refIm has 0 mask applied.
   
   bool m_haveMask {false};
   
   ccImT m_normIm; ///< The normalized image.
   
   ccImT m_ccIm;   ///< The cross-correlation image
   
   ccImT m_magIm; ///< The magnified image, used if m_peakMethod == xcorrPeakMethod::interp
   
   ///@}
   
   math::fit::fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<Scalar>> m_fitter;
   
   int m_maxLag {0}; ///< The maximum lag to consider in the initial cross-correlation.  
   
   Scalar m_tol {0.1}; ///< The tolerance of the interpolated-magnified image, in pixels.
   
   Scalar m_magSize {0}; ///< Magnified size of the ccIm when using interp.  Set as function of m_tol and m_maxLag.
   
public:
   xcorrPeakMethod m_peakMethod {xcorrPeakMethod::centroid};
   
public:
   
   /// Default c'tor
   imageXCorrDiscrete();
   
   /// Construct seeting maxLag.
   explicit imageXCorrDiscrete(int maxLag);
   
   /// Get the current maximum lag 
   /**
     * \returns the current value of m_maxLag 
     */
   int maxLag();
   
   /// Set the maximum lag
   void maxLag( int ml /**< [in] the new maximum lag */);
   
   /// Get the tolerance of the interpolated-magnified image, in pixels.
   /**
     * \returns the current value of m_tol.
     */ 
   Scalar tol();
   
   /// Set the tolerance of the interpolated-magnified image, in pixels.
   void tol( Scalar nt /**< [in] The new value of the interpolation tolerance. */ );
   
   /// Set the size of the cross-correlation images.
   /** This resizes all working memory.
     *
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int resize( int nrows, ///< [in] the number of rows in the images to register
               int ncols  ///< [in] the number of columns in the images to register
             );
   
   /// Set the mask image
   /** 
     * \returns 0 on success
     * \returns -1 on error
     */
   int maskIm( const ccImT & mask /**< [in] the new mask image */ );
   
   /// Get a reference to the mask image.
   /**
     * \returns a const referance to the mask image.
     */ 
   const ccImT & maskIm();
   
   /// Set the reference image
   /** Normalizes the reference image by mean subtraction and variance division.  Applies
     * the mask first if supplied.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */
   int refIm( const ccImT & im0 );
      
   /// Get a reference to the reference image.
   /**
     * \returns a const referent to m_refIm.
     */  
   const ccImT& refIm();
   
   /// Get a reference to the normalized image.
   /**
     * \returns a const referent to m_normIm.
     */
   const ccImT & normIm();
   
   /// Get a reference to the cross correlation image.
   /**
     * \returns a const referent to m_ccIm.
     */
   const ccImT & ccIm();

   /// Get a reference to the magnified image.
   /**
     * \returns a const referent to m_magIm.
     */
   const ccImT & magIm();
   
   /// Conduct the cross correlation to a specified tolerance
   /** 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   template<class imT>
   int operator()( Scalar & xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                   Scalar & yShift, ///< [out] the y shift of im w.r.t. im0, in pixels
                   const imT & im   ///< [in] the image to cross-correlate with the reference
                 );
};

template< class ccImT>
imageXCorrDiscrete<ccImT>::imageXCorrDiscrete()
{
}

template< class ccImT>
imageXCorrDiscrete<ccImT>::imageXCorrDiscrete(int mL)
{
   maxLag(mL);
}

template< class ccImT>
int imageXCorrDiscrete<ccImT>::maxLag()
{
   return m_maxLag;
}

template< class ccImT>
void imageXCorrDiscrete<ccImT>::maxLag( int ml )
{
   m_maxLag = ml;
   tol(m_tol);
}

template< class ccImT>
typename ccImT::Scalar imageXCorrDiscrete<ccImT>::tol()
{
   return m_tol;
}


template< class ccImT>
void imageXCorrDiscrete<ccImT>::tol( Scalar nt )
{
   
   m_magSize = ceil(((2.*m_maxLag + 1) - 1.0)/ nt)+1;
   
   Scalar mag = (m_magSize-1.0)/((2.*m_maxLag + 1) - 1.0);
   
   m_tol = 1.0/mag;
}

template< class ccImT>
int imageXCorrDiscrete<ccImT>::maskIm( const ccImT & mask )
{
   m_maskIm = mask;
   m_haveMask = true;
   
   return 0;
}

template< class ccImT>
const ccImT & imageXCorrDiscrete<ccImT>::maskIm()
{ 
   return m_maskIm; 
}

template< class ccImT>
int imageXCorrDiscrete<ccImT>::refIm( const ccImT & im )
{
   ccImT im0;
   if(m_haveMask)
   {
      if(im.rows()!=m_maskIm.rows() && im.cols() != m_maskIm.cols())
      {
         mxError("imageXCorFit::setReference", MXE_SIZEERR, "reference and mask are not the same size");
         return -1;
      }
      
      im0 = im*m_maskIm;
   }
   else
   {
      im0 = im;
   }
   
   Scalar m = imageMean(im0);
   Scalar v = imageVariance(im0, m);
   m_refIm = (im0 - m)/sqrt(v);
      
   return 0;
}

template< class ccImT>
const ccImT & imageXCorrDiscrete<ccImT>::refIm()
{ 
   return m_refIm; 
}

template< class ccImT>
const ccImT & imageXCorrDiscrete<ccImT>::normIm()
{ 
   return m_normIm; 
}

template< class ccImT>
const ccImT & imageXCorrDiscrete<ccImT>::ccIm()
{ 
   return m_ccIm; 
}

template< class ccImT>
const ccImT & imageXCorrDiscrete<ccImT>::magIm()
{ 
   return m_magIm; 
}

template< class ccImT>
template< class imT>
int imageXCorrDiscrete<ccImT>::operator()( Scalar & xShift,
                                           Scalar & yShift,
                                           const imT & im
                                         )
{
   if( im.rows() <= m_refIm.rows() )
   {
      mxError("imageXCorrDiscrete", MXE_SIZEERR, "reference must be smaller than target image (rows)");
      return -1;
   }
   
   if( im.cols() <= m_refIm.cols() )
   {
      mxError("imageXCorrDiscrete", MXE_SIZEERR, "reference must be smaller than target image (cols)");
      return -1;
   }
                                                  // 16/4    15/4   16/3   15/3
   int maxLag_r = 0.5*(im.rows()-m_refIm.rows()); // 6       5       6      6
   int maxLag_c = 0.5*(im.cols()-m_refIm.cols());
   
   
   
   if(maxLag_r > m_maxLag && m_maxLag != 0) maxLag_r = m_maxLag;
   if(maxLag_c > m_maxLag && m_maxLag != 0) maxLag_c = m_maxLag;
      
   m_ccIm.resize(2*maxLag_r + 1, 2*maxLag_c + 1);
   

   for( int rL = -maxLag_r; rL <= maxLag_r; ++rL)
   {
      for( int cL = -maxLag_c; cL <= maxLag_c; ++cL)
      {                          
                                                         //   16/4          15/4         16/3    15/3
         int r0 = 0.5*im.rows() + rL-0.5*m_refIm.rows(); //   0-15          0-13         0-14    0-14  
         int c0 = 0.5*im.cols() + cL-0.5*m_refIm.cols();
         
         
         if(m_haveMask)
         {
            m_normIm = im.block(r0,c0, m_refIm.rows(), m_refIm.cols())*m_maskIm;
         }
         else
         {  
            m_normIm = im.block(r0,c0, m_refIm.rows(), m_refIm.cols());
         }
         
         Scalar m = imageMean(m_normIm);
         Scalar sv = sqrt(imageVariance(m_normIm,m));
         if(sv <= 0) 
         {
            m_ccIm(maxLag_r +  rL, maxLag_c + cL) = 0;
         }
         else
         {
            m_normIm = (m_normIm-m)/sv;
            m_ccIm(maxLag_r +  rL, maxLag_c + cL) = (m_normIm*m_refIm).sum();    
            //ds9Interface ds9( m_normIm, 1);
         }
         
         
      }
   }
   
   if(m_peakMethod == xcorrPeakMethod::gaussfit)
   {
      int xLag0,yLag0;
      Scalar pk = m_ccIm.maxCoeff(&xLag0, &yLag0);
      Scalar mn = m_ccIm.minCoeff();
      m_fitter.setArray(m_ccIm.data(), m_ccIm.rows(), m_ccIm.cols());
      m_fitter.setGuess(mn, pk, xLag0, yLag0, 0.2*m_ccIm.rows(), 0.2*m_ccIm.cols(), 0); 
      m_fitter.fit();
   
      xShift = m_fitter.x0() - maxLag_r;
      yShift = m_fitter.y0() - maxLag_c;
   }
   else if(m_peakMethod == xcorrPeakMethod::interp)
   {
      m_magIm.resize(m_magSize, m_magSize);
      imageMagnify(m_magIm, m_ccIm, cubicConvolTransform<Scalar>());
      int x, y;
      m_magIm.maxCoeff(&x,&y);
      xShift = x*m_tol - maxLag_r;
      yShift = y*m_tol - maxLag_c;
   }
   else
   {
      Scalar x,y;
      imageCenterOfLight(x,y,m_ccIm);
      xShift = x - maxLag_r;
      yShift = y - maxLag_c;
   }
   
   return 0;
   
}
   
} //improc
} //mx 

#endif //imageXCorrDiscrete_hpp
