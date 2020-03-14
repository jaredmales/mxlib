/** \file imageXCorrFit.hpp
  * \brief A class to register images using the cross correlation with a peak fit.
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

#ifndef imageXCorrFit_hpp
#define imageXCorrFit_hpp

#include "../mxError.hpp"
#include "../math/fit/fitGaussian.hpp"

namespace mx
{
namespace improc
{
   
/// Find the optimum shift to align two images using the discrete cross correlation and a peak fit.
/** This is generally much faster than imageXCorrGrid, and much more precise for high S/N and 
  * closely related images.
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
class imageXCorrFit
{
public:
   typedef _ccImT ccImT; ///< the Eigen-like array type used for image processing
   typedef typename _ccImT::Scalar Scalar; ///< the scalar type of the image type
   
protected:
   
   /** \name Working Memory
     * @{
     */ 
   ccImT m_ccIm;     ///< The coarses cross-correlation image
   
   ///@}
   
   math::fit::fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<Scalar>> m_fitter;
   
   int m_maxLag {0}; ///< The maximum lag to consider in the initial cross-correlation.  
   
public:
   
   /// Default c'tor
   imageXCorrFit();
   
   /// Construct seeting maxLag.
   imageXCorrFit(int maxLag);
   
   /// Get the current maximum lag 
   /**
     * \returns the current value of m_maxLag 
     */
   int maxLag();
   
   /// Set the maximum lag
   void maxLag( int ml /**< [in] the new maximum lag */);
   
protected:


      
public:
   
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
imageXCorrFit<ccImT>::imageXCorrFit()
{
}

template< class ccImT>
imageXCorrFit<ccImT>::imageXCorrFit(int maxLag)
{
   m_maxLag = maxLag;
}

template< class ccImT>
int imageXCorrFit<ccImT>::maxLag()
{
   return m_maxLag;
}

template< class ccImT>
void imageXCorrFit<ccImT>::maxLag( int ml )
{
   m_maxLag = ml;
}

template< class ccImT>
template< class im0T, class imT>
int imageXCorrFit<ccImT>::operator()( Scalar & xShift,
                                   Scalar & yShift,
                                   im0T & im0,
                                   imT & im
                                 )
{
   if( im0.rows() < im.rows() )
   {
      mxError("imageXCorrFit", MXE_SIZEERR, "image must be same size or smaller than reference (rows)");
      return -1;
   }
   
   if( im0.cols() < im.cols() )
   {
      mxError("imageXCorrFit", MXE_SIZEERR, "image must be same size or smaller than reference (cols)");
      return -1;
   }
   
   int maxLag = m_maxLag;
   if(maxLag == 0) 
   {
      maxLag = 0.25*im0.rows()-1;
   }
   
   m_ccIm.resize(2*maxLag + 1, 2*maxLag + 1);
   
   int x0 = 0.5*im0.rows() - 0.5*im.rows() + maxLag;
   int y0 = 0.5*im0.cols() - 0.5*im.cols() + maxLag;

   int cMin = maxLag;
   int cMax = im.cols()-cMin;
   
   int rMin = maxLag;
   int rMax = im.rows()-rMin;
      
   for( int xL = -maxLag; xL <= maxLag; ++xL)
   {
      for( int yL = -maxLag; yL <= maxLag; ++yL)
      {
         m_ccIm(maxLag +  xL, maxLag + yL) = 0;
         
         for(int c=cMin; c < cMax; ++c)
         {
            for(int r=rMin; r < rMax; ++r)
            {
               m_ccIm(maxLag +  xL, maxLag + yL) += im(r,c) * im0(x0+xL + r-rMin, y0+yL+c - cMin);
            }
         }
      }
   }
         
   int xLag0,yLag0;
   Scalar pk = m_ccIm.maxCoeff(&xLag0, &yLag0);
   Scalar mn = m_ccIm.minCoeff();
 
   m_fitter.setArray(m_ccIm.data(), m_ccIm.rows(), m_ccIm.cols());
   m_fitter.setGuess(mn, pk, xLag0, yLag0, 0.2*m_ccIm.rows(), 0.2*m_ccIm.cols(), 0); 
   m_fitter.fit();
   
   xShift = -m_fitter.x0() + maxLag;
   yShift = -m_fitter.y0() + maxLag;
   return 0;
   
}
   
} //improc
} //mx 

#endif //imageXCorrFit_hpp
