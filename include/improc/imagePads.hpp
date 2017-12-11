/** \file imagePads.hpp
  * \brief Image padding 
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __imagePads_hpp__
#define __imagePads_hpp__

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

///Pad an image with a constant value
/** 
  * \param [out] imOut is resized, and on return contains the padded image.
  * \param [in] imIn the image to be padded.
  * \param [in] padSz the size of the pad to be added.  The image will grow by 2*padSz rows and columns.
  * \param [in] value the value to use for padding.
  *
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array
  * \tparam imInT is an Eigen-like array
  */ 
template<typename imOutT, typename imInT>
int padImage( imOutT & imOut, 
              imInT & imIn, 
              unsigned int padSz, 
              typename imOutT::Scalar value)
{
   int nRows = imIn.rows() + 2*padSz;
   int nCols = imIn.cols() + 2*padSz;
   
   imOut = imOutT::Constant(nRows, nCols, value);
   
   imOut.block(padSz, padSz, imIn.rows(),imIn.cols()) = imIn;
   
   return 0;
}

///Pad an image with a constant value
/** This version can be used with Eigen reference types.
  *  
  * \param [out] imOut is resized, and on return contains the padded image.
  * \param [in] imIn the image to be padded.  
  * \param [in] padSz the size of the pad to be added.  The image will grow by 2*padSz rows and columns.
  * \param [in] value the value to use for padding.
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array reference
  * \tparam imInT is an Eigen-like array reference
  */ 
template<typename imOutT, typename imInT>
int padImageRef( imOutT imOut, 
                 imInT imIn, 
                 unsigned int padSz,
                 typename imOutT::Scalar value)
{
   int nRows = imIn.rows() + 2*padSz;
   int nCols = imIn.cols() + 2*padSz;

   
   imOut = imOutT::Constant(nRows, nCols, value);
   
   imOut.block(padSz, padSz, imIn.rows(),imIn.cols()) = imIn;
   
   return 0;
}

///Pad an image by repeating the values in the edge rows and columns
/** Allocates imOut to hold 2*padSz more rows and columns, and copies imIn to the center.  Then fills in the new pixels
  * by copying the edge pixels outward.
  *
  * \param [out] imOut is resized, and on return contains the padded image.
  * \param [in] imIn the image to be padded.  
  * \param [in] padSz the size of the pad to be added.  The image will grow by 2*padSz rows and columns.
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array 
  * \tparam imInT is an Eigen-like array
  */
template<typename imOutT, typename imInT>
int padImage( imOutT & imOut,
              imInT & imIn,
              unsigned int padSz )
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

///Pad an image by repeating the values at the edge of a 1/0 mask
/** Allocates imOut to match imMask, and copies imIn to the center multiplied by the mask.  Then fills in the new pixels
  * by copying the edge pixels outward. For each pixel on the edge of the mask (touching at least one 1-valued pixel
  * in the mask), the value is chosen as the value of the nearest pixel with a 1 in the mask.  If the closest pixels are more than one equidistant pixels,
  * their mean value is used.    
  * 
  * \param [out] imOut is resized to match imMask, and on return contains the padded image.
  * \param [in] imIn the image to be padded.  
  * \param [in] imMask is the 1/0 mask.
  * \param [in] padSz the number of iterations of the padding loop.
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array
  * \tparam imInT is an Eigen-like array 
  */
template<typename imOutT, typename imInT>
int padImage( imOutT & imOut,
              imInT & imIn,
              imInT & imMask,
              unsigned int padSz )
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
   
   for(int n=0; n < padSz; ++n)
   {
      nmi.clear();
      nmj.clear();
      for(int ii =0; ii< dim1; ++ii)
      {
         for(int jj =0; jj<dim2; ++jj)
         {
            if(tMask(ii,jj) == 1) continue;
            
            typename imOutT::Scalar mv = 0;
            int r, minr = 3;
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
                     //if(nv == 0) nv = 1;
                     
//                     r = i*i + j*j;
  
//                      if(r == minr) //Equidistant pixels
//                      {
//                         mv += imOut( ii+i, jj+j); //So we average
//                         ++nv;
//                      }
//                      else if(r < minr ) //New closest pixel
//                      {
//                         mv = imOut( ii+i, jj+j);
//                         nv = 1;
//                         minr = r;
//                      }
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
      
      for(int i=0; i< nmi.size(); ++i)
      {
         tMask( nmi[i], nmj[i]) = 1;
      }
      
   }
      
   return 0;
   
}

///Cut down a padded image
/**
  * \param imOut is resized, and on return contains the cut image. 
  * \param inIn is the image to be cut.  
  * \param padSz the size of the pad.  The image will shrink by 2*padSz rows and columns.
  * 
  * \retval 0 on success
  * \retval -1 on error
  * 
  * \tparam imOutT is an Eigen-like array
  * \tparam imInT is an Eigen-like array
  */ 
template<typename imOutT, typename imInT>
int cutPaddedImage( imOutT & imOut, 
                    const imInT & imIn, 
                    unsigned int padSz )
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

#endif //__imagePads_hpp__
