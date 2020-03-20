/** \file imageXCorrFFT.hpp
  * \brief A class to register images using the Fourier cross correlation with a peak fit.
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

#ifndef imageXCorrFFT_hpp
#define imageXCorrFFT_hpp

#include "../mxError.hpp"
#include "../fft/fft.hpp"
#include "../math/fit/fitGaussian.hpp"

#include "imageUtils.hpp"
#include "imageTransforms.hpp"

namespace mx
{
namespace improc
{
   
/// Find the optimum shift to align two images using the Fourier cross correlation and a peak fit.
/** 
  * 
  * The shift is reported in pixels such that if the mxlib imageShift function is used
  * to shift the input image by the negative of the shifts, it will align with the 
  * reference.
  *
  * \tparam _ccImT is the Eigen-like array type used for image processing.  See typedefs.
  * 
  * \ingroup image_reg
  */ 
template<class _ccImT>
class imageXCorrFFT
{
public:
   typedef _ccImT ccImT; ///< the Eigen-like array type used for image processing
   
   typedef typename _ccImT::Scalar Scalar; ///< the scalar type of the image type
   
   typedef typename _ccImT::Scalar realT; ///< the scalar type of the image type
   
   typedef std::complex<realT> complexT; ///< Complex floating point type.
   
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> realArrayT; ///< Real eigen array type
   
   typedef Eigen::Array<complexT, Eigen::Dynamic, Eigen::Dynamic> complexArrayT; ///< Complex eigen array type with Scalar==complexT
   
protected:
   
   int m_rows {0};
   
   int m_cols {0};
   
   /** \name Working Memory
     * @{
     */ 
   ccImT m_ccIm;     ///< The cross-correlation image
   
   complexArrayT m_ftIm0;  ///< Working memory for the FT of the reference image.
   
   complexArrayT m_ftWork; ///< Working memory for the FFT.
   
   fftT< realT, complexT,2,0> m_fft_fwd; ///< FFT object for the forward transform.
   
   fftT< complexT, realT,2,0> m_fft_back; ///< FFT object for the backward transfsorm.
   
   ///@}
   
   math::fit::fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<Scalar>> m_fitter;
   
   int m_maxLag {0}; ///< The maximum lag to consider in the initial cross-correlation.  
   
   
public:
   
   /// Default c'tor
   imageXCorrFFT();
   
   /// Construct seeting maxLag.
   explicit imageXCorrFFT(int maxLag);
   
   /// Get the current maximum lag 
   /**
     * \returns the current value of m_maxLag 
     */
   int maxLag();
   
   /// Set the maximum lag
   void maxLag( int ml /**< [in] the new maximum lag */);
   
   /// Set the size of the cross-correlation images.
   /** This resizes all working memory and conducts fftw planning.
     *
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int resize( int nrows, ///< [in] the number of rows in the images to register
               int ncols  ///< [in] the number of columns in the images to register
             );
   
   /// Set the reference image
   /** Normalizes, fourier transforms, and conjugates the reference image
     * 
     * \returns 0 on success
     * \returns -1 on error
     */
   int setReference( const ccImT & im0 );
   
protected:


      
public:
   
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
   
   /// Conduct the cross correlation to a specified tolerance
   /** 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   template<class im0T, class imT>
   int operator()( Scalar & xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                   Scalar & yShift, ///< [out] the y shift of im w.r.t. im0, in pixels
                   im0T & im0,      ///< [in] the reference image
                   imT & im         ///< [in] the image to cross-correlate with the reference
                 );
};

template< class ccImT>
imageXCorrFFT<ccImT>::imageXCorrFFT()
{
}

template< class ccImT>
imageXCorrFFT<ccImT>::imageXCorrFFT(int maxLag)
{
   m_maxLag = maxLag;
}

template< class ccImT>
int imageXCorrFFT<ccImT>::maxLag()
{
   return m_maxLag;
}

template< class ccImT>
void imageXCorrFFT<ccImT>::maxLag( int ml )
{
   m_maxLag = ml;
}

template< class ccImT>
int imageXCorrFFT<ccImT>::resize( int nrows,
                                  int ncols
                                )
{
   if( m_rows == nrows && m_cols == ncols)
   {
      return 0;
   }
   
   m_rows = nrows;
   
   m_cols = ncols;
   
   m_ccIm.resize(m_rows, m_cols);
   
   m_ftIm0.resize( (int) (0.5*m_rows) + 1, m_cols);
   
   m_ftWork.resize( (int) (0.5*m_rows) + 1, m_cols);

   //fftw is row-major, eigen defaults to column-major
   m_fft_fwd.plan(m_cols, m_rows, MXFFT_FORWARD, false);
      
   m_fft_back.plan(m_cols, m_rows, MXFFT_BACKWARD, false);   
   
   return 0;
}

template< class ccImT>
int imageXCorrFFT<ccImT>::setReference( const ccImT & im0 )
{
   resize(im0.rows(), im0.cols());
   
   imageShiftWP( m_ccIm, im0, 0.5*m_rows, 0.5*m_cols );
   realT m = imageMean(m_ccIm);
   realT v = imageVariance(m_ccIm, m);
   m_ccIm = (m_ccIm - m)/sqrt(v);
   
   m_fft_fwd(m_ftIm0.data(), m_ccIm.data());
   
   //Conjugate and normalize.
   for(int c =0; c < m_ftIm0.cols(); ++c)
   {
      for(int r=0; r < m_ftIm0.rows(); ++r)
      {
         complexT val = std::conj(m_ftIm0(r,c));
         val /= (m_rows*m_cols);
         m_ftIm0(r,c) = val;
      }
   }
   
   return 0;
}

template< class ccImT>
template< class imT>
int imageXCorrFFT<ccImT>::operator()( Scalar & xShift,
                                      Scalar & yShift,
                                      const imT & im
                                    )
{
   if( im.rows() != m_rows )
   {
      mxError("imageXCorrFFT", MXE_SIZEERR, "image must be same size as reference (rows)");
      return -1;
   }
   
   if( im.cols() != m_cols )
   {
      mxError("imageXCorrFFT", MXE_SIZEERR, "image must be same size as reference (cols)");
      return -1;
   }
   
   int maxLag = m_maxLag;
   if(maxLag == 0) 
   {
      maxLag = 0.25*m_rows-1;
   }
   
   realT m = imageMean(im);
   realT v = imageVariance(im, m);
   m_ccIm = (im - m)/sqrt(v);
   m_fft_fwd(m_ftWork.data(), m_ccIm.data());

   //m_ftIm0 is the conjugated, fftw-normalized, reference image in the Fourier domain
   //So this is the FT of the cross-correlation:
   m_ftWork *= m_ftIm0;
   
   m_fft_back(m_ccIm.data(), m_ftWork.data());
   
   int xLag0,yLag0;
   Scalar pk = m_ccIm.maxCoeff(&xLag0, &yLag0);
   Scalar mn = m_ccIm.minCoeff();
 
   ccImT subim = m_ccIm.block(xLag0-maxLag, yLag0-maxLag, 2*maxLag+1, 2*maxLag+1);
   
   m_fitter.setArray(subim.data(), subim.rows(), subim.cols());
   m_fitter.setGuess(mn, pk, maxLag, maxLag, 0.2*subim.rows(), 0.2*subim.cols(), 0); 
   m_fitter.fit();
   
   xShift = m_fitter.x0() + (xLag0-maxLag) - (int)(0.5*m_rows);
   yShift = m_fitter.y0() + (yLag0-maxLag) - (int)(0.5*m_cols);
   
   return 0;
   
}

template< class ccImT>
template< class im0T, class imT>
int imageXCorrFFT<ccImT>::operator()( Scalar & xShift,
                                      Scalar & yShift,
                                      im0T & im0,
                                      imT & im
                                    )
{
   setReference(im0);
   return operator()(xShift, yShift, im);
   
}
   
} //improc
} //mx 

#endif //imageXCorrFFT_hpp
