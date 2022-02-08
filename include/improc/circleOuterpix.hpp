/** \file circleOuterpix.hpp
  * \brief Declares and defines a class for finding the edge of a circle mask
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#ifndef improc_circleOuterpix_hpp
#define improc_circleOuterpix_hpp


#include "eigenImage.hpp"

#include "imageUtils.hpp"

namespace mx
{
namespace improc 
{
   
/// Find the center and the outermost pixels of a circular mask, giving an estimate of the radius
/** Takes a 1/0 mask of a filled circle or annulus, and calculates the center coordinate by 
  * center-of-light, then finds the outermost 1 pixel in each row and column.  Calculates the average radius
  * of these pixels w.r.t. to the c.o.l. coords, as well as w.r.t. the average coords of the outer pixels
  * themselves. 
  * 
  * Also produces an image containing the outer pixels.
  * 
  * \returns 0 on success
  * 
  * \tparam realT the real floating-point type for calculations
  * 
  * \ingroup image_utils
  * 
  */ 
template<typename realT>
int circleOuterpix( realT & x0,                      /// < [out] x coord of center of the mask by center-of-light
                    realT & y0,                      /// < [out] y coord of center of the mask by center-of-light
                    realT & avgr0,                   /// < [out] average radius of edge pixels w.r.t. center-of-light
                    realT & avgx,                    /// < [out] x coord of center of the mask by average edge pixel
                    realT & avgy,                    /// < [out] y coord of center of the mask by average edge pixel
                    realT & avgr,                    /// < [out] average radius of edge pixels w.r.t. to average edge pixel
                    eigenImage<realT> & circ,        /// < [out] image showing the edge pixels with value 1, all other pixels 0.  Resized.
                    const eigenImage<realT> & masked /// < [in] an image with a (roughly) circular mask of 1s 
                  )
{
   
   circ.resize(masked.rows(), masked.cols());
   circ.setZero();

   mx::improc::imageCenterOfLight(x0, y0, masked);

   for(size_t i=0; i< masked.rows(); ++i)
   {
      realT ndy =0;
      realT pdy = 0;

      int nj = -1;
      int pj = -1;
      
      for(size_t j=0; j<masked.cols(); ++j)
      {
         if(masked(i,j) == 1)
         {
            if(j-y0 < ndy)
            {
               ndy = j-y0;
               nj = j;
            }
            if(j-y0 > pdy)
            {
               pdy = j-y0;
               pj = j;
            }
         }
      }
      
      if(nj != -1) circ(i,nj) = 1;
      if(pj != -1) circ(i,pj) = 1;
   }
   
   for(size_t j=0; j< masked.cols(); ++j)
   {
      realT ndx =0;
      realT pdx = 0;

      int ni = -1;
      int pi = -1;
      
      for(size_t i=0; i<masked.rows(); ++i)
      {
         if(masked(i,j) == 1)
         {
            if(i-x0 < ndx)
            {
               ndx = i-x0;
               ni = i;
            }
            if(i-x0 > pdx)
            {
               pdx = i-x0;
               pi = i;
            }
         }
      }
      
      if(ni != -1) circ(ni,j) = 1;
      if(pi != -1) circ(pi,j) = 1;
   }
   
   int npix = 0;
   avgx = 0;
   avgy = 0;
   
   for(size_t i = 0; i < masked.rows(); ++i)
   {
      for(size_t j = 0; j < (size_t) masked.cols(); ++j)
      {
         if(circ(i,j) == 1)
         {
            ++npix;
            avgx += i;
            avgy += j;
         }
      }
   }
      
   avgx /= npix;
   avgy /= npix;

   avgr = 0;
   avgr0 = 0;
   for(size_t i = 0; i < (size_t) circ.rows(); ++i)
   {
      for(size_t j = 0; j < (size_t) circ.cols(); ++j)
      {
         if(circ(i,j) == 1)
         {
            avgr += sqrt( pow(i-avgx,2) + pow(j-avgy,2));
            avgr0 += sqrt( pow(i-x0,2) + pow(j-y0,2));
         }
      }
   }

   avgr /= npix;
   avgr0 /= npix;

   return 0;
   
}

} //namespace improc 
} //namespace mx

#endif // improc_circleOuterpix_hpp
