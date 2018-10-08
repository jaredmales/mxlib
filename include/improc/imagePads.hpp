/** \file imagePads.hpp
  * \brief Image padding 
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef improc_imagePads_hpp
#define improc_imagePads_hpp

#include <vector>

#include "../mxError.hpp"

namespace mx
{
namespace improc 
{

/** \addtogroup image_padding
  * These functions pad an image by adding rows and columns, and filling in the new pixels with appropriate values.  This can be a constant, or nearest neighbor.  A version
  * is also provided which expands an aribtrary 1/0 mask, which is useful when working with, say, circular regions of an image.
  */

/** \ingroup image_padding
  * @{ 
  */

/// Pad an image with a constant value
/** 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array
  * \tparam imInT is an Eigen-like array
  */ 
template<typename imOutT, typename imInT>
int padImage( imOutT & imOut, ///< [out] On return contains the padded image.  This will be resized.
              imInT & imIn, ///< [in] The image to be padded.
              unsigned int padSz, ///< [in] The size of the pad.  The padded image (imOut) will be 2*padSz rows and cols larger than the input.
              typename imOutT::Scalar value  ///< [in] the value to use for padding.
            )
{
   int nRows = imIn.rows() + 2*padSz;
   int nCols = imIn.cols() + 2*padSz;
   
   imOut = imOutT::Constant(nRows, nCols, value);
   
   imOut.block(padSz, padSz, imIn.rows(),imIn.cols()) = imIn;
   
   return 0;
}

/// Pad an image with a constant value for reference types.
/** This version can be used with Eigen reference types.
  *  
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array reference
  * \tparam imInT is an Eigen-like array reference
  */ 
template<typename imOutT, typename imInT>
int padImageRef( imOutT imOut, ///< [out] On return contains the padded image.  This will be resized.
                 imInT imIn, ///< [in] The image to be padded.
                 unsigned int padSz, ///< [in] The size of the pad.  The padded image (imOut) will be 2*padSz rows and cols larger than the input.
                 typename imOutT::Scalar value ///< [in] the value to use for padding.
               )
{
   int nRows = imIn.rows() + 2*padSz;
   int nCols = imIn.cols() + 2*padSz;

   
   imOut = imOutT::Constant(nRows, nCols, value);
   
   imOut.block(padSz, padSz, imIn.rows(),imIn.cols()) = imIn;
   
   return 0;
}

/// Pad an image by repeating the values in the edge rows and columns
/** Allocates imOut to hold 2*padSz more rows and columns, and copies imIn to the center.  Then fills in the new pixels
  * by copying the edge pixels outward.
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array 
  * \tparam imInT is an Eigen-like array
  */
template<typename imOutT, typename imInT>
int padImage( imOutT & imOut, ///< [out] On return contains the padded image.  This will be resized.
              imInT & imIn, ///< [in] The image to be padded.
              unsigned int padSz ///< [in] The size of the pad.  The padded image (imOut) will be 2*padSz rows and cols larger than the input.
            )
{
   int dim1 = imIn.rows();
   int dim2 = imIn.cols();
   
   imOut.resize(dim1+2*padSz, dim2 + 2*padSz);
   
   imOut.block(padSz, padSz, dim1, dim2) = imIn;
   
   //First do the left and right edge columns and corners
   for(int i = 0; i< padSz; ++i)
   {
      for(int j=0; j < imOut.cols(); ++j)
      {
         if( j < padSz )
         {
            imOut(i,j) = imIn(0, 0);
            imOut( imOut.rows()-1-i, j) = imIn( dim1-1, 0);
         }
         else if( j >= dim2 + padSz)
         {
            imOut(i,j) = imIn(0, dim2-1);
            imOut( imOut.rows()-1-i, j) = imIn( dim1-1, dim2-1);
         }
         else
         {
            imOut( i, j ) = imOut( padSz, j );
            imOut( imOut.rows() - padSz + i,  j )   = imOut( padSz + dim1-1, j );
         }
      }
   }
   
   //Now do the top and bottom rows.
   for(int i=padSz; i < dim1+padSz; ++i)
   {
      for(int j=0; j< padSz; ++j)
      {
         imOut(i,j) = imOut(i, padSz);
         imOut(i, imOut.cols()-1-j) = imOut(i, imOut.cols()-1-padSz);
      }
   }
   
   
   return 0;
}

/// Pad an image by repeating the values at the edge of a 1/0 mask
/** Allocates imOut to match imMask, and copies imIn to the center multiplied by the mask.  Then fills in the new pixels
  * by copying the edge pixels outward. For each pixel on the edge of the mask (touching at least one 1-valued pixel
  * in the mask), the value is chosen as the value of the nearest pixel with a 1 in the mask.  If the closest pixels are more than one equidistant pixels,
  * their mean value is used.    
  * 
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array
  * \tparam imInT is an Eigen-like array 
  */
template<typename imOutT, typename imInT>
int padImage( imOutT & imOut, ///< [out] On return contains the padded image.  This will be resized.
              imInT & imIn, ///< [in] The image to be padded.
              imInT & imMask, ///< [in] The 1/0 mask image.
              unsigned int padSz ///< [in] The number of iterations of the padding loop.
            )
{

   typedef imInT imMaskT;
   
   int dim1 = imMask.rows();
   int dim2 = imMask.cols();
   
   imOut.resize(dim1, dim2);
   
   int stR = 0.5*(dim1 - imIn.rows());
   int stC = 0.5*(dim2 - imIn.cols());


   if(stR < 0)
   {
      mxError("padImage", MXE_INVALIDARG, "imMask must have at least as many rows as imIn.");
      return -1;
   }
   if(stC < 0)
   {
      mxError("padImage", MXE_INVALIDARG, "imMask must have at least as many columns as imIn.");
      return -1;
   }
   
   imOut.block(stR, stC, imIn.rows(), imIn.cols()) = imIn;
   
   imOut *= imMask;
   
   imMaskT tMask = imMask; //Need a modifiable mask, which grows as the image is extended
   
   std::vector<int> nmi, nmj; //Used to track which pixels of the mask need to be updated
   
   for(unsigned int n=0; n < padSz; ++n)
   {
      nmi.clear();
      nmj.clear();
      for(int ii =0; ii< dim1; ++ii)
      {
         for(int jj =0; jj<dim2; ++jj)
         {
            if(tMask(ii,jj) == 1) continue;
            
            typename imOutT::Scalar mv = 0;
            //int r, minr;// = 3;
            int nv = 0;
            
            //Loop over the neighboring pixels
            for(int i=-1; i<=1; ++i)
            {
               for(int j=-1; j<=1; ++j)
               {
                  if(i==0 && j==0) continue;
                  if( ii + i < 0 || ii + i >= dim1) continue;
                  if( jj + j < 0 || jj + j >= dim2) continue;
                  
                  if( tMask( ii+i, jj+j) == 1)
                  {
                     mv += imOut( ii+i, jj+j);
                     ++nv;
  
                  }
               }
            }
            
            //Found at least 1 neighboring pixel in the 1's part of the mask
            if(nv > 0)
            {
               imOut(ii, jj) = mv/nv;
               nmi.push_back(ii);
               nmj.push_back(jj);
            }
         }
      }
      
      for(size_t i=0; i< nmi.size(); ++i)
      {
         tMask( nmi[i], nmj[i]) = 1;
      }
      
   }
      
   return 0;
   
}

///Cut down a padded image
/**
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array
  * \tparam imInT is an Eigen-like array
  */ 
template<typename imOutT, typename imInT>
int cutPaddedImage( imOutT & imOut,  ///< [out] On return contains the cut image.  This will be resized.
                    const imInT & imIn, ///< [in] The image to be cut down.
                    unsigned int padSz ///< [in] The size of the pad.  The output image will be smaller by 2*padSz rows and cols.
                  )
{
   if(2*padSz > imIn.rows())
   {
      mxError("cutPaddedImage", MXE_INVALIDARG, "padSz too large, must be < imIn.rows()/2");
      return -1;
   }
   
   if(2*padSz > imIn.cols())
   {
      mxError("cutPaddedImage", MXE_INVALIDARG, "padSz too large, must be < imIn.cols()/2");
      return -1;
   }
   
   int nRows = imIn.rows() - 2*padSz;
   int nCols = imIn.cols() - 2*padSz;
   
   imOut = imIn.block(padSz, padSz, nRows, nCols);
   
   return 0;
}

///@}

   
}// namespace improc 
}// namespace mx

#endif //improc_imagePads_hpp
