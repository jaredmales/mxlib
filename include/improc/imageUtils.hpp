/** \file imageUtils.hpp
  * \author Jared R. Males
  * \brief  Header for the image processing utilities
  * \ingroup image_processing_files
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

#ifndef improc_imageUtils_hpp
#define improc_imageUtils_hpp



namespace mx
{
namespace improc
{

/** \ingroup image_utils
  *@{
  */

/// Calculate the mean value of an image
/**
  *
  * \returns the mean of the input image
  */ 
template< class imageT>
typename imageT::Scalar imageMean(imageT & im /**< [in] the image of which to calculate the mean*/)
{
   typename imageT::Scalar m = 0;
   
   for(int c=0;c<im.cols();++c)
   {
      for(int r=0;r<im.rows();++r)
      {
         m += im(r,c);
      }
   }
   
   m/=(im.rows()*im.cols());
   
   return m;
}

/// Calculate the variance of an image w.r.t. a given value
/**
  *
  * \returns the variance of the input image
  */ 
template< class imageT>
typename imageT::Scalar imageVariance( imageT & im,                 ///< [in] the image of which to calculate the variance
                                       typename imageT::Scalar mean ///< [in] the value to calculate the variance with respect to
                                     )
{
   typename imageT::Scalar v = 0;
   
   for(int c=0;c<im.cols();++c)
   {
      for(int r=0;r<im.rows();++r)
      {
         v += pow(im(r,c)-mean,2);
      }
   }
   
   v /= (im.rows()*im.cols());
   
   return v;
}

template<typename imageT>
int centerOfLight( typename imageT::Scalar & x,
                   typename imageT::Scalar & y,
                   imageT & im
                 )
{
   x =0;
   y= 0;
   
   typename imageT::Scalar sum = im.sum();
   
   for(int i=0; i< im.rows(); ++i)
   {
      for(int j=0; j < im.cols(); ++j)
      {
         x += (i+1)*im(i,j);
         y += (j+1)*im(i,j);
      }
   }
   
   x = x/sum - 1;
   y = y/sum - 1;
   
   return 0;
}

///@}

} //namespace math
} //namespace mx



#endif // improc_imageUtils_hpp
