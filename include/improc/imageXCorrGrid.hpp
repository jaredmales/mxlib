/** \file imageXCorrGrid.hpp
  * \brief A class to register images using the cross correlation and a collapsing grid search.
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

#ifndef imageXCorrGrid_hpp
#define imageXCorrGrid_hpp

#include "../mxError.hpp"
#include "imageTransforms.hpp"

namespace mx
{
namespace improc
{
   
/// Find the optimum shift to align two images using the discrete cross correlation and a refining grid search.
/** This is generally much slower for the same precision compared to imageXCorrFit.
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
class imageXCorrGrid
{
public:
   typedef _ccImT ccImT; ///< the Eigen-like array type used for image processing
   typedef typename _ccImT::Scalar Scalar; ///< the scalar type of the image type
   
protected:
   
   /** \name Working Memory
     * @{
     */ 
   ccImT m_ccIm;     ///< The coarses cross-correlation image
   ccImT m_subGrid0; ///< One of the fine cross-correlation images, always 5x5
   ccImT m_subGrid1; ///< One of the fine cross-correlation images, always 5x5
   ccImT m_shiftIm;  ///< Holds the shifted image during the fine cross correlation.
   
   ///@}
   
   int m_maxLag {0}; ///< The maximum lag to consider in the initial cross-correlation.  
   Scalar m_gridTol {0.0625}; ///< The grid tolerance to iterate to.
   
public:
   
   /// Default c'tor
   imageXCorrGrid();
   
   /// Construct seeting maxLag.
   explicit imageXCorrGrid(int maxLag);
   
   /// Get the current maximum lag 
   /**
     * \returns the current value of m_maxLag 
     */
   int maxLag();
   
   /// Set the maximum lag
   void maxLag( int ml /**< [in] the new maximum lag */);
   
   /// Get the current grid tolerance
   /**
     * \returns the current value of m_gridTol 
     */
   int gridTol();
   
   /// Set the grid tolerance
   void gridTol( Scalar gt /**< [in] the new grid tolerance */);
   
protected:

   /// Calculate the cross-correlation of a 5x5 subgrid of a sub-pixel spacing
   /** This is called iteratively with ever smaller spacing by the function operator.
     * New values of the c.c. are only calculated for pixels not sampled by the previous grid, which
     * is nominally 14/25.
     * 
     * \returns 0 on success
     * \returns -1 on an error 
     */
   template<class im0T, class imT>
   int subGrid( Scalar & xShift,     ///< [out] the total pixel shift corresponding to the peak c.c., in x
                Scalar & yShift,     ///< [out] the total pixel shift corresponding to the peak c.c., in y
                int & xLag,          ///< [out] the lag in pixels on the 5x5 grid corresponding to the peak c.c., in x
                int & yLag,          ///< [out] the lag in pixels on the 5x5 grid corresponding to the peak c.c., in y
                ccImT & nextSubGrid, ///< [out] the 5x5 grid of c.c. values on this new sampling
                Scalar xShift0,      ///< [in] the previous total pixel shift corresponding to the peak c.c., in x
                Scalar yShift0,      ///< [in] the previous total pixel shift corresponding to the peak c.c., in y
                int xLag0,           ///< [in] the previous lag in pixels on the 5x5 grid corresponding to the peak c.c., in x
                int yLag0,           ///< [in] the previous the lag in pixels on the 5x5 grid corresponding to the peak c.c., in y
                ccImT & lastSubGrid, ///< [in] 5x5 grid of c.c. values on the previous sampling
                Scalar dLag,         ///< [in] the lag per pixel of lastSubGrid
                im0T & im0,          ///< [in] the reference image
                imT & im,            ///< [in] the image to cross-correlate with the reference
                int maxLag           ///< [in] the maximum lag considered
              );
      
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
imageXCorrGrid<ccImT>::imageXCorrGrid()
{
   m_subGrid0.resize(5,5);
   m_subGrid1.resize(5,5);
}

template< class ccImT>
imageXCorrGrid<ccImT>::imageXCorrGrid(int maxLag)
{
   m_maxLag = maxLag;
   m_subGrid0.resize(5,5);
   m_subGrid1.resize(5,5);
}

template< class ccImT>
int imageXCorrGrid<ccImT>::maxLag()
{
   return m_maxLag;
}

template< class ccImT>
void imageXCorrGrid<ccImT>::maxLag( int ml )
{
   m_maxLag = ml;
}

template< class ccImT>
int imageXCorrGrid<ccImT>::gridTol()
{
   return m_gridTol;
}

template< class ccImT>
void imageXCorrGrid<ccImT>::gridTol( Scalar gt )
{
   m_gridTol = gt;
}
   
template< class ccImT>
template< class im0T, class imT>
int imageXCorrGrid<ccImT>::subGrid( Scalar & xShift, 
                                Scalar & yShift,
                                int & xLag,
                                int & yLag,
                                ccImT & nextSubGrid,
                                Scalar xShift0,
                                Scalar yShift0,
                                int xLag0,
                                int yLag0,
                                ccImT & lastSubGrid,
                                Scalar dLag,
                                im0T & im0,
                                imT & im,
                                int maxLag
                              )
{
   if( im0.rows() < im.rows() )
   {
      mxError("imageXCorrGrid", MXE_SIZEERR, "image must be same size or smaller than reference (rows)");
      return -1;
   }
   
   if( im0.cols() < im.cols() )
   {
      mxError("imageXCorrGrid", MXE_SIZEERR, "image must be same size or smaller than reference (cols)");
      return -1;
   }
   
   int cMin = maxLag;
   int cMax = im.cols()-cMin;

   int rMin = maxLag;
   int rMax = im.rows()-rMin;
      
   for( int xL = -2; xL <= 2; ++xL)
   {
      for( int yL = -2; yL <= 2; ++yL)
      {
         //If we don't need to recalculate we don't
         if( xL%2 == 0 && yL%2 == 0 && xLag0 +xL/2 > 0 && xLag0 +xL/2 < lastSubGrid.rows() && yLag0 +yL/2 > 0 && yLag0 +yL/2 < lastSubGrid.cols())
         {
            nextSubGrid(2 +  xL, 2 + yL) = lastSubGrid(xLag0 + xL/2, yLag0 + yL/2);
         }
         else
         {
            nextSubGrid( 2 +  xL, 2 + yL) = 0;
         
            imageShift(m_shiftIm, im0, xShift0 + xL*0.5*dLag, yShift0 + yL*0.5*dLag, cubicConvolTransform<Scalar>(-0.5));
         
            for(int c=cMin; c < cMax; ++c)
            {
               for(int r=rMin; r < rMax; ++r)
               {               
                  nextSubGrid( 2 +  xL, 2 + yL) += im(r,c) * m_shiftIm(r,c);
               }
            }
         }
      }
   }
   
   nextSubGrid.maxCoeff(&xLag, &yLag);
   xShift = xShift0 + (xLag-2)*0.5 * dLag;
   yShift = yShift0 + (yLag-2)*0.5 * dLag;
   
   return 0;
}


template< class ccImT>
template< class im0T, class imT>
int imageXCorrGrid<ccImT>::operator()( Scalar & xShift,
                                   Scalar & yShift,
                                   im0T & im0,
                                   imT & im
                                 )
{
   if( im0.rows() < im.rows() )
   {
      mxError("imageXCorrGrid", MXE_SIZEERR, "image must be same size or smaller than reference (rows)");
      return -1;
   }
   
   if( im0.cols() < im.cols() )
   {
      mxError("imageXCorrGrid", MXE_SIZEERR, "image must be same size or smaller than reference (cols)");
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
   m_ccIm.maxCoeff(&xLag0, &yLag0);
   Scalar xShift0 = -xLag0 + maxLag;
   Scalar yShift0 = -yLag0 + maxLag;
    
   if(m_gridTol >= 1) 
   {
      xShift=xShift0;
      yShift=yShift0;
      return 0;
   }
   
   int nSub = ceil(-log2(m_gridTol));
   
   int xLag, yLag;
   
   //First halving:
   subGrid( xShift, yShift, xLag, yLag, m_subGrid0, xShift0, yShift0, xLag0, yLag0, m_ccIm, 1, im0, im, maxLag);
   xShift0=xShift;
   yShift0=yShift;
   
   //Subsequent halvings, alternating sub-grids:
   for(int n=1; n< nSub; ++n)
   {
      if(n%2 == 1)
      {
         subGrid( xShift, yShift, xLag, yLag, m_subGrid1, xShift0, yShift0, xLag0, yLag0, m_subGrid0, pow(0.5, n), im0, im, maxLag);
      }
      else
      {
         subGrid( xShift, yShift, xLag, yLag, m_subGrid0, xShift0, yShift0, xLag0, yLag0, m_subGrid1, pow(0.5, n), im0, im, maxLag);
      }
      
      //std::cerr << n << " " << pow(0.5, n+1) << " " << xShift << " " << yShift << "\n";
   
      xShift0=xShift;
      yShift0=yShift;
      xLag0 = xLag;
      yLag0 = yLag;
   }
   
   return 0;
   
}
   
} //improc
} //mx 

#endif //imageXCorrGrid_hpp
