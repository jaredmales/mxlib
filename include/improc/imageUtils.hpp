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


#include "imageTransforms.hpp"

namespace mx
{
namespace improc
{

/** \ingroup image_utils
  *@{
  */

/// Reflect pixel coordinates across the given center pixel.
/**
  */ 
template<typename realT>
int reflectImageCoords( int & x1, ///< [out] the reflected x coordinate
                        int & y1, ///< [out] the reflected y coordinate
                        int x0,   ///< [in] the input x coordinate
                        int y0,   ///< [in] the input y coordinate
                        realT xc, ///< [in] the center pixel x coordinate
                        realT yc  ///< [in] the center pixel y coordinate
                      )
{
   x1 =  xc - (x0 - xc);
   y1 =  yc - (y0 - yc);

   return 0;
}

/// Zero any NaNs in an image
template<class imageT>
void zeroNaNs( imageT & im /**< [in/out] image which will have any NaN pixels set to zero */)
{
   for(int c=0; c< im.cols(); ++c)
   {
      for(int r=0; r< im.rows(); ++r)
      {
         if( !std::isnormal( im(r,c) ) )
         {
            im(r,c) = 0;
         }
      }
   }
}

/// Zero any NaNs in an image cube
template<class cubeT>
void zeroNaNCube( cubeT & imc /**< [in/out] cube which will have any NaN pixels set to zero */)
{
   for(int p=0; p< imc.planes(); ++p) 
   {
      for(int c=0; c< imc.cols(); ++c)
      {
         for(int r=0; r< imc.rows(); ++r)
         {
            if( !std::isnormal( imc.image(p)(r,c) ) )
            {
               imc.image(p)(r,c) = 0;
            }
         }
      }
   }
}

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
int imageCenterOfLight( typename imageT::Scalar & x,
                        typename imageT::Scalar & y,
                        const imageT & im
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

/// Find the maximum in an image at sub-pixel resolution by interpolation
/** Uses imageMagnify() to expand the image to the desired scale.  Because of the
  * scaling used in imageMagnify, the desired scale may not be exact.  As a result
  * the actual scale is returned in scale_x and scale_y.
  *
  */ 
template<typename floatT, typename imageT, typename magImageT, typename transformT>
int imageMaxInterp( floatT & x,        ///< [out] the x-position of the maximum, in pixels of the input image
                    floatT & y,        ///< [out] the y-position of the maximum, in pixels of the input image
                    floatT & scale_x,  ///< [in/out] the desired scale or resolution, in pixels < 1, in the x direction.  On output contains the actual scale calculated.
                    floatT & scale_y,  ///< [in/out] the desired scale or resolution, in pixels < 1, in the y direction.  On output contains the actual scale calculated.
                    magImageT & magIm, ///< [in] the magnified image.  This is used as working memory, will be resized.
                    const imageT & im, ///< [in] the image to find the maximum of
                    transformT trans   ///< [in] the transform to use for interpolation
                  )
{
   floatT magSize_x = ceil( (im.rows() - 1.0)/ scale_x) + 1;
   floatT magSize_y = ceil( (im.cols() - 1.0)/ scale_y) + 1;
   
   floatT mag_x = ((floatT) magSize_x-1.0)/((floatT) im.rows() - 1.0);
   floatT mag_y = ((floatT) magSize_y-1.0)/((floatT) im.cols() - 1.0);
   
   scale_x = 1.0/mag_x;
   scale_y = 1.0/mag_y;
   
   magIm.resize(magSize_x, magSize_y);
   
   imageMagnify(magIm, im, trans);
   
   int ix, iy;
   magIm.maxCoeff(&ix,&iy);
   x = ix*scale_x;
   y = iy*scale_y;
      
   return 0;
}

/// Find the maximum in an image at sub-pixel resolution by cubic convolution interpolation
/** Uses imageMagnify() to expand the image to the desired scale.  Because of the
  * scaling used in imageMagnify, the desired scale may not be exact.  As a result
  * the actual scale is returned in scale_x and scale_y.
  *
  */ 
template<typename floatT, typename imageT, typename magImageT>
int imageMaxInterp( floatT & x,        ///< [out] the x-position of the maximum, in pixels of the input image
                    floatT & y,        ///< [out] the y-position of the maximum, in pixels of the input image
                    floatT & scale_x,  ///< [in/out] the desired scale or resolution, in pixels < 1, in the x direction.  On output contains the actual scale calculated.
                    floatT & scale_y,  ///< [in/out] the desired scale or resolution, in pixels < 1, in the y direction.  On output contains the actual scale calculated.
                    magImageT & magIm, ///< [in] the magnified image.  This is used as working memory, will be resized.
                    const imageT & im  ///< [in] the image to find the maximum of
                  )
{
   return imageMaxInterp(x,y,scale_x,scale_y, magIm, im, cubicConvolTransform<typename imageT::Scalar>());
}




/// Combine two images, each with their own mask defining good pixels.
/** The combined image is made up of the pixels in im1 where mask1 is 1, and the pixels of im2 where mask2 is 1.
  * If a pixel in both mask1 and mask2 has a value of 1, that pixel in the combo is the average of im1 and im2.
  * All other pixels are set to 0 in the combined image.
  *
  * Separate template types are used for each argument to allow references, etc.
  * 
  * \tparam imageT the eigen-like array type of the combined image 
  * \tparam imageT1 the eigen-like array type of image 1
  * \tparam imageT2 the eigen-like array type of mask 1
  * \tparam imageT3 the eigen-like array type of image 2
  * \tparam imageT4 the eigen-like array type of mask 2
  */
template<typename imageT, typename imageT1, typename imageT2, typename imageT3, typename imageT4>
void combine2ImagesMasked( imageT & combo,        ///< [out] the combined image.  will be resized.
                           const imageT1 & im1,   ///< [in] the first image
                           const imageT2 & mask1, ///< [in] the mask for the first image
                           const imageT3 & im2,   ///< [in] the second image
                           const imageT4 & mask2  ///< [in] the mask for the second image
                         )
{
   combo.resize(im1.rows(), im2.cols());
   
   for(int c=0; c < combo.cols(); ++c)
   {
      for(int r=0;r<combo.rows();++r)
      {
         if(mask1(r,c) == 1 && mask2(r,c) == 0) combo(r,c) = im1(r,c);
         else if(mask2(r,c) == 1 & mask1(r,c) == 0) combo(r,c) = im2(r,c);
         else if(mask1(r,c) == 1 && mask2(r,c) == 1) combo(r,c) = 0.5*(im1(r,c) + im2(r,c));
         else combo(r,c) = 0;
      }
   }
}

template<typename eigenT, typename eigenTin>
void removeRowsAndCols(eigenT & out, const eigenTin & in, int st, int w)
{
   
   out.resize(in.rows() - w, in.cols() - w);
   
   out.topLeftCorner(st,st) = in.topLeftCorner(st,st);
   
   out.bottomLeftCorner(in.rows()-(st+w), st) = in.bottomLeftCorner(in.rows()-(st+w), st);
   
   out.topRightCorner(st, in.cols()-(st+w))  = in.topRightCorner(st, in.cols()-(st+w));
   
   out.bottomRightCorner(in.rows()-(st+w),in.cols()-(st+w)) = in.bottomRightCorner(in.rows()-(st+w),in.cols()-(st+w));
}


template<typename eigenT, typename eigenTin>
void removeRows(eigenT & out,  const eigenTin & in, int st, int w)
{
   
   out.resize(in.rows() - w, in.cols());
   
   out.topLeftCorner(st,in.cols()) = in.topLeftCorner(st,in.cols());
   
   out.bottomLeftCorner(in.rows()-(st+w), in.cols()) = in.bottomLeftCorner(in.rows()-(st+w), in.cols());
   
}

template<typename eigenT, typename eigenTin>
void removeCols(eigenT & out,  const eigenTin & in, int st, int w)
{
   
   out.resize(in.rows(), in.cols() - w);
   
   out.topLeftCorner(in.rows(), st) = in.topLeftCorner(in.rows(), st);
   
   out.topRightCorner(in.rows(),in.cols()-(st+w)) = in.topRightCorner(in.rows(),in.cols()-(st+w));
   
}   

/// Copy one image to another, with no transformation
/** This is merely memcpy 
  * 
  * \returns dest
  */ 
void * imcpy( void * dest,   ///< [out] the address of the first pixel in the destination image
              void * src,    ///< [in] the address of the first pixel in the source image
              size_t width,  ///< [in] the width in pixels of size szof
              size_t height, ///< [in] the height in pixels of size szof
              size_t szof    ///< [in] the size in bytes of a one pixel
            );

/// Copy one image to another, flipping up-down
/** This is a reversed row-by-row memcpy
  * 
  * \returns dest
  */ 
void * imcpy_flipUD( void * dest,   ///< [out] the address of the first pixel in the destination image
                     void * src,    ///< [in] the address of the first pixel in the source image
                     size_t width,  ///< [in] the width in pixels of size szof
                     size_t height, ///< [in] the height in pixels of size szof
                     size_t szof    ///< [in] the size in bytes of a one pixel
                   );

/// Copy one image to another, flipping left-right
/** 
  * 
  * \returns dest
  */ 
void * imcpy_flipLR( void * dest,   ///< [out] the address of the first pixel in the destination image
                     void * src,    ///< [in] the address of the first pixel in the source image
                     size_t width,  ///< [in] the width in pixels of size szof
                     size_t height, ///< [in] the height in pixels of size szof
                     size_t szof    ///< [in] the size in bytes of a one pixel
                   );

/// Copy one image to another, flipping up-down and left-right
/** 
  * 
  * \returns dest
  */ 
void * imcpy_flipUDLR( void * dest,   ///< [out] the address of the first pixel in the destination image
                       void * src,    ///< [in] the address of the first pixel in the source image
                       size_t width,  ///< [in] the width in pixels of size szof
                       size_t height, ///< [in] the height in pixels of size szof
                       size_t szof    ///< [in] the size in bytes of a one pixel
                     );

///@}

} //namespace math
} //namespace mx



#endif // improc_imageUtils_hpp
