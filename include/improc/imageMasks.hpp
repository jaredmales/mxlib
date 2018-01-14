/** \file imageMasks.hpp
  * \brief Declares and defines functions to work with image masks
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
   
#ifndef improc_imageMasks_hpp
#define improc_imageMasks_hpp

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include "imageTransforms.hpp"

namespace mx
{

namespace improc 
{
   

/// Fills in the cells of an Eigen 2D Array with their radius from the center
/** \ingroup image_masks
  *
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class eigenT> 
void radiusImage( eigenT & m, ///< [out] the allocated radius array, will be filled in with radius values. 
                  typename eigenT::Scalar xc, ///< [in] the x center
                  typename eigenT::Scalar yc,  ///< [in] the y center
                  typename eigenT::Scalar scale=1 ///< [in] [optional] a scaling to apply to each value (default = 1)
                )
{
   typedef typename eigenT::Scalar arithT;
   
   arithT f_x, f_y;

   size_t dim1 = m.rows();
   size_t dim2 = m.cols();
   
   for(size_t i=0; i < dim1; i++)
   {
      f_x = (i-xc)*(i-xc);
      
      for(size_t j=0; j < dim2; j++)
      {
         f_y = (j-yc)*(j-yc);
         
         m(i,j) = sqrt( f_x + f_y)*scale;
      }
   }
   
}

/// Fills in the cells of Eigen-like 2D Array with their radius from the canonical center
/** \ingroup image_masks
  *
  * The center is @f$ (x_c, y_c) = (0.5*(dim_1-1), 0.5*(dim_2 -1)) @f$.
  *
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class eigenT> 
void radiusImage( eigenT & m, ///< [out] the allocated radius array, will be filled in with radius values.
                  typename eigenT::Scalar scale=1 ///< [in] [optional] a scaling to apply to each value (default = 1)
                )
{
   typedef typename eigenT::Scalar arithT;
   
   arithT xc, yc;
   
   xc = 0.5*(m.rows()-1);
   yc = 0.5*(m.cols()-1);
   
   radiusImage(m, xc, yc, scale);
}




/// Fills in the cells of an Eigen-like 2D Array with their angle relative to the center
/** \ingroup image_masks
  *
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class eigenT> 
void angleImage( eigenT & m, ///< [out]  the allocated angle array.  Will be filled in with angle values.
                 typename eigenT::Scalar xc, ///< [in] the x center
                 typename eigenT::Scalar yc, ///< [in] the y center
                 typename eigenT::Scalar scale=radian<typename eigenT::Scalar>()  ///< [in] [optional] a scaling to apply to each angle value. Default converts to degrees, set to 1 for radians.
               )
{
   typedef typename eigenT::Scalar arithT;
   arithT f_x, f_y;

   size_t dim1 = m.rows();
   size_t dim2 = m.cols();
   
   for(size_t i=0; i < dim1; i++)
   {
      f_x = (i-xc);
      
      for(size_t j=0; j < dim2; j++)
      {
         f_y = (j-yc);
         
         m(i,j) = fmod(atan2(f_y, f_x) + two_pi<typename eigenT::Scalar>(), two_pi<typename eigenT::Scalar>())  *scale;
      }
   }
   
}

/// Fills in the cells of Eigen 2D Array with their angle relative the canonical center
/** \ingroup image_masks
  *
  * The center is @f$ (x_c, y_c) = (0.5*(dim_1-1), 0.5*(dim_2 -1)) @f$.
  *
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class eigenT> 
void angleImage( eigenT & m, ///< [out] the allocated angle array.  Will be filled in with angle values.
                 typename eigenT::Scalar scale = radian<typename eigenT::Scalar>() ///< [in] [optional] a scaling to apply to each angle value. Default converts to degrees, set to 1 for radians.
               )
{
   typedef typename eigenT::Scalar arithT;
   
   arithT xc, yc;
   
   xc = 0.5*(m.rows()-1);
   yc = 0.5*(m.cols()-1);
   
   angleImage(m, xc, yc, scale);
  
}


/// Fills in the cells of Eigen-like arrays with their radius amd angle relative to the center
/** \ingroup image_masks
  *
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class eigenT> 
void radAngImage( eigenT & rIm, ///< [out] the allocated radius array, will be filled in with radius values.
                  eigenT & qIm, ///< [out] the allocated angle array, same size as rIm.  Will be filled in with angle values.
                  typename eigenT::Scalar xc, ///< [in] the x center
                  typename eigenT::Scalar yc, ///< [in] the y center
                  typename eigenT::Scalar rscale = 1, ///< [in] [optional] a scaling to apply to each radius value. Default is 1.0.
                  typename eigenT::Scalar qscale= radian<typename eigenT::Scalar>()  ///< [in] [optional] a scaling to apply to each angle value. Default converts to degrees, set to 1 for radians.
                )
{
   typedef typename eigenT::Scalar arithT;
   
   arithT f_x, f_y;

   size_t dim1 = rIm.rows();
   size_t dim2 = rIm.cols();
   
   for(size_t i=0; i < dim1; ++i)
   {
      f_x = ( ((arithT)i)-xc);
      
      for(size_t j=0; j < dim2; ++j)
      {
         f_y = (((arithT)j)-yc);
         rIm(i,j) = std::sqrt( f_x*f_x + f_y*f_y)*rscale;
         qIm(i,j) = fmod(atan2(f_y, f_x) + two_pi<arithT>(), two_pi<arithT>()) *qscale;
      }
   }
}



///Get the vector indices of an annular region in an image
/**
  * \ingroup image_masks
  * 
  * \tparam eigenT is an Eigen-like 2D array type
  * 
  * \returns a vector containing the 1D indices of the region defined by the input parameters
  */
template<typename eigenT>
std::vector<size_t> annulusIndices( eigenT &rIm,  ///< [in] a radius image of the type produced by \ref radiusImage
                                    eigenT &qIm,  ///< [in] an angle image of the type produce by \ref angleImage
                                    typename eigenT::Scalar xcen,  ///< [in] the x center of the image
                                    typename eigenT::Scalar ycen,  ///< [in] the y center of the image
                                    typename eigenT::Scalar min_r, ///< [in] the minimum radius of the region
                                    typename eigenT::Scalar max_r, ///< [in] the maximum radius of the region
                                    typename eigenT::Scalar min_q, ///< [in] the minimum angle of the region
                                    typename eigenT::Scalar max_q, ///< [in] the maximum angle of the region
                                    eigenT * mask = 0 ///< [in] [optional] pointer to a mask image, only pixels of value 1 are included in the indices.
                                  )
{

   std::vector<size_t> idx;
   
   int min_x = -max_r, max_x = max_r, min_y = -max_r, max_y = max_r;

   if(max_q == 0) max_q = 360.;
   
   size_t msize = ((pi<double>()*(max_r*max_r - min_r*min_r)) * (max_q-min_q)/360.) *1.01 + 1;
   
   //This was tested, this is slightly faster than resize with an erase.
   idx.reserve(msize);
   
   int x0 = xcen+min_x;
   if(x0 < 0) x0 = 0;
   int x1 = xcen+max_x+1;
   if(x1 > rIm.rows()) x1 = rIm.rows();
   int y0 = ycen+min_y;
   if(y0 < 0) y0 = 0;
   int y1 = ycen+max_y+1;
   if(y1 > rIm.cols()) y1 = rIm.cols();
   
   for(size_t i = x0; i< x1; ++i)
   {
      for(size_t j = y0; j< y1; ++j)
      { 
         if(rIm(i,j) >= min_r && rIm(i,j) < max_r && qIm(i,j) >= min_q && qIm(i,j) < max_q) 
         {
            if( mask )
            {
               if( (*mask)(i,j) == 0) continue;
            }
            
            idx.push_back(j*rIm.cols() + i);
         }
      }
   }
   
   
   return idx;
}



inline void rectangleIndices( std::vector<size_t> & idx, 
                              size_t rows, 
                              size_t cols, 
                              size_t xmin, 
                              size_t xmax, 
                              size_t ymin, 
                              size_t ymax
                            )
{
      
   //if(xmin < 0) xmin = 0;
   if(xmax > rows-1) xmax = rows-1;
   
   //if(ymin < 0) ymin = 0;
   if(ymax > cols-1) ymax = cols-1;
   
   idx.reserve( (xmax-xmin+1)*(ymax-ymin + 1) );
   
   for(int i=xmin; i<=xmax; ++i)
   {
      for(int j=ymin;j<=ymax; ++j)
      {
         idx.push_back( j*rows + i);
      }
   }
}

template<class eigenT>
void rectangleIndices( std::vector<size_t> & idx, 
                       eigenT & mask, 
                       size_t xmin, 
                       size_t xmax, 
                       size_t ymin, 
                       size_t ymax
                     )
{
   rectangleIndices(idx, (size_t) mask.rows(), (size_t) mask.cols(), xmin, xmax, ymin, ymax);  
}

///Apply a mask to an image
/** The pixels indicated by the indices are set to a value.
  * 
  * \ingroup image_masks
  * 
  * \tparam eigenT is an Eigen-like 2D array type
  *  
  */
template<class eigenT> 
void applyMask( eigenT & maskedIm,  ///< [out] the image to mask (will be modified)
                const std::vector<size_t> & idx, ///< [in] the indices of the pixels to mask
                typename eigenT::Scalar maskval ///< [in] the mask value.
              )
{
   for(size_t i=0; i< idx.size(); ++i)
   {
      maskedIm(idx[i]) = maskval;
   }
}


///Mask a circle in an image.
/** The circle is describe by its center coordinates and radius. Any value can be set for the mask,
  * with 0 being the default.  
  * 
  * \tparam arrayT is an Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<class arrayT> 
void maskCircle( arrayT & m, ///< [in/out] the image to be masked, is modified.
                 typename arrayT::Scalar xcen, ///< [in] the x coordinate of the center of the circle
                 typename arrayT::Scalar ycen, ///< [in] the y coordinate of the center of the circle
                 typename arrayT::Scalar rad,  ///< [in] the radius of the circle
                 typename arrayT::Scalar val = 0,  ///< [in] [optional] the mask value.  Default is 0.
                 typename arrayT::Scalar pixbuf = 0.5  ///< [in] [optional] buffer for radius comparison.  Default is 0.5 pixels.
               )
{
   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   typename arrayT::Scalar r;
   
   
   for(size_t i=0; i < l0; i++)
   {
      for(size_t j=0; j < l1; j++)
      {
         r = sqrt( pow(i-xcen, 2) + pow(j-ycen, 2) );
         
         if(r <= rad+pixbuf) m(i,j) = val;
      }
   }
}   

///Populate a mask based on a typical CCD bleeding pattern.
/** Masks a circle for saturated pixels, as well as a horizontal bar for bleeding (in both directions if needed).
  *
  * \returns 0 on success, -1 otherwise.
  * 
  * \tparam imT is an Eigen-like 2D array type.
  */ 
template< typename imT >
int ccdBleedMask( imT & im,  ///< [out] the mask, on output will be 1/0.  Must be allocated prior to call.
                  typename imT::Scalar x0, ///< [in] the x-coordinate of the center of the mask
                  typename imT::Scalar y0, ///< [in] the y-coordinate of the center of the mask
                  typename imT::Scalar rad, ///< [in] the radius of the central circle.
                  typename imT::Scalar height, ///< [in] the half-height of the saturation bar.
                  typename imT::Scalar lwidth, ///< [in] the length of the saturation bar to the left.
                  typename imT::Scalar rwidth ///< [in] the length of the saturation bar to the right.
                )
{
   typename imT::Scalar x, y;
   
   for(int i=0; i<im.rows(); ++i)
   {
      x = i;
      for(int j=0; j<im.cols(); ++j)
      {
         y = j;
         
         if( sqrt( pow(x-x0,2) + pow(y-y0,2)) <= rad)
         {
            im(i,j) = 0;
         }
         else if( fabs(y - y0) <= height && x >= x0-lwidth && x <= x0 + rwidth )
         {
            im(i,j) = 0;
         }
         else
         {
            im(i,j) = 1;
         }
      }
   }
   
   return 0;
}

/// Cut out a region of an image specified by an index-mask.
/** The output sill be a row-image containing the pixel values.
  * 
  * \tparam imageTout is the output image Eigen-like type.
  * \tparam imageTin is the input image Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<typename imageTout, typename imageTin>
void cutImageRegion( imageTout & imout,  ///< [out] a row-image containing the pixel values specified by the indices.
                     const imageTin & imin,  ///< [in] a 2-D image
                     const std::vector<size_t> & idx, ///< [in] the linear indices of the pixel values
                     bool resize = true ///< [in] [optional] flag controls whether imout is resized.  It is if true.
                   )
{
   if(resize)
   {
      imout.resize(idx.size(),1);
   }
   
   //#pragma omp parallel for schedule(static, 1)
   for(int i=0;i<idx.size();++i)
   {
      imout(i) = imin( idx[i] );
   }
   
}

/// Insert a region of an image specified by an index-mask.
/** Both the input and the output are row-images containing the pixel values.
  * 
  * \tparam imageTout is the output image Eigen-like type.
  * \tparam imageTin is the input image Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<typename imageTout, typename imageTin >
void insertImageRegion( imageTout imout, ///< [out] a row-image into which the pixel values specified by the indices are inserted.
                        const imageTin & imin,  ///< [in] a row-image containing the pixel values, same size as idx
                        const std::vector<size_t> & idx ///< [in] the linear indices of the pixel values
                      )
{
   //#pragma omp parallel for schedule(static, 1)
   for(int i=0;i< idx.size();++i)
   {
      imout(idx[i],0) = imin(i,0);
   }
   
} 

template< typename imageT, typename transformT = cubicConvolTransform<typename imageT::Scalar>>
void rotateMask( imageT & rotMask,
                 imageT & mask,
                 typename imageT::Scalar angle
               )
{
   typedef typename imageT::Scalar Scalar;
   
   imageRotate( rotMask, mask, angle, transformT());
   
   for(int ii=0; ii< rotMask.rows(); ++ii)
   {
      for(int jj=0; jj < rotMask.cols(); ++jj)
      {
         if( rotMask(ii,jj) < 0.5) rotMask(ii,jj) = 0;
         else rotMask(ii,jj) = 1;
      }
   }
}
   



} //namespace improc 
} //namespace mx

#endif // improc_imageMasks_hpp
