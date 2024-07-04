/** \file imageMasks.hpp
  * \brief Declares and defines functions to work with image masks
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

#ifndef improc_imageMasks_hpp
#define improc_imageMasks_hpp

#include "../math/constants.hpp"
#include "../math/geo.hpp"

#include "eigenImage.hpp"
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
void radiusImage( eigenT & m,                     ///< [out] the allocated radius array, will be filled in with radius values. 
                  typename eigenT::Scalar xc,     ///< [in] the x center
                  typename eigenT::Scalar yc,     ///< [in] the y center
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
void radiusImage( eigenT & m,                     ///< [out] the allocated radius array, will be filled in with radius values.
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
  * \tparam angleT is the angle type, either radiansT<realT> or degreesT<realT>.  Note that realT sets the type for all arithmetic.
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class angleT, class eigenT> 
void angleImage( eigenT & m,                ///< [out]  the allocated angle array.  Will be filled in with angle values.
                 typename angleT::realT xc, ///< [in] the x center
                 typename angleT::realT yc  ///< [in] the y center
               )
{
   typedef typename angleT::realT realT;
   
   for(size_t j=0; j < m.cols(); j++)
   {
      realT f_y = ( static_cast<realT>(j)-yc);
         
      for(size_t i=0; i < m.rows(); i++)
      {
         realT f_x = ( static_cast<realT>(i)-xc);
      
         m(i,j) = fmod(atan2(f_y, f_x) + math::two_pi<realT>(), math::two_pi<realT>()) * angleT::scale;
      }
   }
   
}

/// Fills in the cells of Eigen 2D Array with their angle relative the canonical center
/** \ingroup image_masks
  *
  * The center is @f$ (x_c, y_c) = (0.5*(dim_1-1), 0.5*(dim_2 -1)) @f$.
  *
  * \tparam angleT is the angle type, either radiansT<realT> or degreesT<realT>.  Note that realT sets the type for all arithmetic.
  * \tparam eigenT is an Eigen-like 2D array type
  */  
template<class angleT, class eigenT> 
void angleImage( eigenT & m /** < [out] the allocated angle array.  Will be filled in with angle values. */)
{
   typedef typename angleT::realT realT;
      
   realT xc = 0.5*(m.rows()-1);
   realT yc = 0.5*(m.cols()-1);
   
   angleImage<angleT>(m, xc, yc);
}


/// Fills in the cells of Eigen-like arrays with their radius amd angle relative to the center
/** \ingroup image_masks
  *
  * \tparam angleT is the angle type, either radiansT<realT> or degreesT<realT>.  Note that realT sets the type for all arithmetic.
  * \tparam eigenT1 is an Eigen-like 2D array type.  Should be resolved by compiler.
  * \tparam eigenT2 is an Eigen-like 2D array type.  Should be resolved by compiler.
  * 
  */  
template<class angleT, class eigenT1, class eigenT2> 
void radAngImage( eigenT1 & rIm,                    ///< [out] the allocated radius array, will be filled in with radius values.
                  eigenT2 & qIm,                    ///< [out] the angle array, will be re-sized to match rIm.  Will be filled in with angle values.
                  typename angleT::realT xc,        ///< [in] the x center
                  typename angleT::realT yc,        ///< [in] the y center
                  typename angleT::realT rscale = 1 ///< [in] [optional] a scaling to apply to each radius value. Default is 1.0.
                )
{
   typedef typename angleT::realT realT;
   
   qIm.resize(rIm.rows(), rIm.cols());

   for(int cc=0; cc < rIm.cols(); ++cc)
   {
      realT f_y = ( (static_cast<realT>(cc))-yc);
      
      for(int rr=0; rr < rIm.rows(); ++rr)
      {
         realT f_x = ((static_cast<realT>(rr))-xc);

         rIm(rr,cc) = std::sqrt( f_x*f_x + f_y*f_y)*rscale;
         qIm(rr,cc) = fmod(atan2(f_y, f_x) + math::two_pi<realT>(), math::two_pi<realT>()) * angleT::scale;
      }
   }
}



///Get the vector indices of an annular region in an image
/**
  * \ingroup image_masks
  * 
  * \tparam angleT is the angle type, either radiansT<realT> or degreesT<realT>.  Note that realT sets the type for all arithmetic.
  * \tparam eigenT1 is an Eigen-like 2D array type.  Should be resolved by compiler.
  * \tparam eigenT2 is an Eigen-like 2D array type.  Should be resolved by compiler.
  * \tparam eigenT3 is an Eigen-like 2D array type.  Should be resolved by compiler.
  * 
  * \returns a vector containing the 1D indices of the region defined by the input parameters
  */
template<typename angleT, typename eigenT1, typename eigenT2, typename eigenT3=eigenT1>
std::vector<size_t> annulusIndices( const eigenT1 & rIm,          ///< [in] a radius image of the type produced by \ref radiusImage
                                    const eigenT2 & qIm,          ///< [in] an angle image of the type produce by \ref angleImage
                                    typename angleT::realT xcen,  ///< [in] the x center of the image
                                    typename angleT::realT ycen,  ///< [in] the y center of the image
                                    typename angleT::realT min_r, ///< [in] the minimum radius of the region
                                    typename angleT::realT max_r, ///< [in] the maximum radius of the region
                                    typename angleT::realT min_q, ///< [in] the minimum angle of the region. 
                                    typename angleT::realT max_q, ///< [in] the maximum angle of the region. 
                                    eigenT3 * mask = 0            /**< [in] [optional] pointer to a mask image, only pixels of value 1 
                                                                                       are included in the indices.*/
                                  )
{

   std::vector<size_t> idx;
   
   int min_x = -max_r;
   int max_x = max_r; 
   int min_y = -max_r;
   int max_y = max_r;
   
   int x0 = xcen+min_x;
   if(x0 < 0) x0 = 0;

   int x1 = xcen+max_x+1;
   if(x1 > rIm.rows()) x1 = rIm.rows();

   int y0 = ycen+min_y;
   if(y0 < 0) y0 = 0;

   int y1 = ycen+max_y+1;
   if(y1 > rIm.cols()) y1 = rIm.cols();
   
   
   //Normalize the angles by mod
   min_q = math::angleMod<angleT>(min_q);
   //If max_q is exactly 360/2pi (to within a tol of 100*eps) we don't mod, other mod.
   if( fabs(max_q - angleT::full) > 100*std::numeric_limits<typename angleT::realT>::epsilon()) max_q = math::angleMod<angleT>(max_q);

   //Find the mid-point to enable sectors wider than 180.
   typename angleT::realT mid_q;
   if(min_q <= max_q) mid_q = 0.5*(min_q + max_q);
   else mid_q = 0.5*((min_q - angleT::full) + max_q);
   
   std::cerr << min_q << " " << mid_q << " " << max_q << "\n";
   for(int j = y0; j < y1; ++j)
   {
      for(int i = x0; i < x1; ++i)
      { 
         if(rIm(i,j) < min_r) continue;
         if(rIm(i,j) >= max_r) continue;

         if( !( (math::angleDiff<angleT>(min_q,qIm(i,j)) >= 0 && math::angleDiff<angleT>(qIm(i,j),mid_q) > 0) ||
              (math::angleDiff<angleT>(mid_q,qIm(i,j)) >= 0 && math::angleDiff<angleT>(qIm(i,j),max_q) > 0)) ) continue;

         if( mask )
         {
            if( (*mask)(i,j) == 0) continue;
         }
            
         idx.push_back(j*rIm.rows() + i);
      }
   }
   
   return idx;
}

///Get the coordinates of the bounding rectangle of an annulus.
/**
  * \ingroup image_masks
  * 
  * \tparam angleT is the angle type, either radiansT<realT> or degreesT<realT>.  Note that realT sets the type for all arithmetic.
  * 
  */
template<typename angleT>
void annulusBoundingRect( int & x0,                     ///< [out] The lower left x-coordinate of the bounding rect.
                          int & y0,                     ///< [out] The lower left y-coordinate of the bounding rect.
                          int & x1,                     ///< [out] The upper right x-coordinate of the bounding rect.
                          int & y1,                     ///< [out] The upper right y-coordinate of the bounding rect.
                          typename angleT::realT xcen,  ///< [in] the x center of the image
                          typename angleT::realT ycen,  ///< [in] the y center of the image
                          typename angleT::realT min_r, ///< [in] the minimum radius of the region
                          typename angleT::realT max_r, ///< [in] the maximum radius of the region
                          typename angleT::realT min_q, ///< [in] the minimum angle of the region. 
                          typename angleT::realT max_q  ///< [in] the maximum angle of the region. 
                        )
{
    typedef typename angleT::realT realT;

    //Get the corners
    realT x00 = xcen + min_r*cos(min_q*angleT::radians);
    realT y00 = ycen + min_r*sin(min_q*angleT::radians);
    realT x01 = xcen + max_r*cos(min_q*angleT::radians);
    realT y01 = ycen + max_r*sin(min_q*angleT::radians);

    realT x10 = xcen + min_r*cos(max_q*angleT::radians);
    realT y10 = ycen + min_r*sin(max_q*angleT::radians);
    realT x11 = xcen + max_r*cos(max_q*angleT::radians);
    realT y11 = ycen + max_r*sin(max_q*angleT::radians);

    //vertex of min_r is probably not necessary, but might as well
    realT x20 = xcen + min_r*cos(math::angleMean<angleT>({min_q, max_q})*angleT::radians);
    realT y20 = ycen + min_r*sin(math::angleMean<angleT>({min_q, max_q})*angleT::radians);

    //vertex of max_r
    realT x21 = xcen + max_r*cos(math::angleMean<angleT>({min_q, max_q})*angleT::radians);
    realT y21 = ycen + max_r*sin(math::angleMean<angleT>({min_q, max_q})*angleT::radians);

    x0 = std::ceil(std::min({x00,x01,x10,x11,x20,x21}));
    y0 = std::ceil(std::min({y00,y01,y10,y11,y20,y21}));
    x1 = std::ceil(std::max({x00,x01,x10,x11,x20,x21}));
    y1 = std::ceil(std::max({y00,y01,y10,y11,y20,y21}));
}

/// Reflect vector indices across the given center pixel.
/** This assumes column major order.
  * \returns the vector of reflected indices
  */
template<typename realT>
std::vector<size_t> reflectImageIndices( const std::vector<size_t> & idxi, ///< [in] the vector indices to reflect
                                         int w,                            ///< [in] the image width
                                         int h,                            ///< [in] the image height
                                         realT xc,                         ///< [in] the image center x coordinate
                                         realT yc                          ///< [in] the image center y coordinate
                                       )
{
   std::vector<size_t> idxr;
   idxr.reserve(idxi.size());

   int x0, y0, x1, y1;
   for(size_t n =0; n < idxi.size(); ++n)
   {
      int y0 = idxi[n] / w;
      int x0 = idxi[n] - y0*w;

      reflectImageCoords(x1, y1, x0, y0, xc, yc);

      if(x1 < w && y1 < h) idxr.push_back( y1*w + x1 );
   }

   return idxr;
}

template<typename sizeT>
void rectangleIndices( std::vector<sizeT> & idx, 
                              sizeT rows, 
                              sizeT cols, 
                              sizeT xmin, 
                              sizeT xmax, 
                              sizeT ymin, 
                              sizeT ymax
                            )
{
      
   //if(xmin < 0) xmin = 0;
   if(xmax > rows-1) xmax = rows-1;
   
   //if(ymin < 0) ymin = 0;
   if(ymax > cols-1) ymax = cols-1;
   
   idx.reserve( (xmax-xmin+1)*(ymax-ymin + 1) );
   
   for(sizeT i=xmin; i<=xmax; ++i)
   {
      for(sizeT j=ymin;j<=ymax; ++j)
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
   rectangleIndices<size_t>(idx, (size_t) mask.rows(), (size_t) mask.cols(), xmin, xmax, ymin, ymax);  
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
/** The circle is describe by its center coordinates and radius. Any value can be set for the mask.
  * Pixels outside the masked circle are not altered.
  * 
  * \tparam arrayT is an Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<class arrayT> 
void maskCircle( arrayT & m,                          ///< [in.out] the image to be masked, is modified.
                 typename arrayT::Scalar xcen,        ///< [in] the x coordinate of the center of the circle
                 typename arrayT::Scalar ycen,        ///< [in] the y coordinate of the center of the circle
                 typename arrayT::Scalar rad,         ///< [in] the radius of the circle
                 typename arrayT::Scalar val,         ///< [in] the mask value. 
                 typename arrayT::Scalar pixbuf = 0.5 ///< [in] [optional] buffer for radius comparison.  Default is 0.5 pixels.
               )
{
   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   typename arrayT::Scalar rr;
   
   
   for(size_t c=0; c < m.cols(); c++)
   {
      for(size_t r=0; r < m.rows(); r++)
      {
         rr = sqrt( pow(r-xcen, 2) + pow(c-ycen, 2) );
         
         if(rr <= rad+pixbuf) m(r,c) = val;
      }
   }
}   

///Mask a circle in an image at the standard center.
/** The circle is centered at 0.5*(rows-1) and 0.5*(cols-1), and described by its radius. Any value can be set for the mask.  
  * Pixels outside the masked circle are not altered.
  * 
  * \tparam arrayT is an Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<class arrayT> 
void maskCircle( arrayT & m,                          ///< [in.out] the image to be masked, is modified.
                 typename arrayT::Scalar rad,         ///< [in] the radius of the circle
                 typename arrayT::Scalar val,         ///< [in] the mask value.
                 typename arrayT::Scalar pixbuf = 0.5 ///< [in] [optional] buffer for radius comparison.  Default is 0.5 pixels.
               )
{
   return maskCircle(m, 0.5*(m.rows()-1.0), 0.5*(m.cols()-1.0), rad, val, pixbuf);
}   

/// Mask an ellipse in an image.
/** The ellipse is describe by its center coordinates and x and y direction radii (the semi-major and -minor axes, 
  * in either order). Any value can be set for the mask,cwith 0 being the default.
  * 
  * \tparam arrayT is an Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<class arrayT> 
void maskEllipse( arrayT & m,                         ///< [in.out] the image to be masked, is modified.
                  typename arrayT::Scalar xcen,        ///< [in] the x coordinate of the center of the ellipse
                  typename arrayT::Scalar ycen,        ///< [in] the y coordinate of the center of the ellipse
                  typename arrayT::Scalar xrad,        ///< [in] the x radius of the ellipse
                  typename arrayT::Scalar yrad,        ///< [in] the y radius of the ellipse
                  typename arrayT::Scalar ang,         ///< [in] the c.c.w. angle to rotate the ellipse by
                  typename arrayT::Scalar val = 0,     ///< [in] [optional] the mask value.  Default is 0.
                  typename arrayT::Scalar pixbuf = 0.5 ///< [in] [optional] buffer for radius comparison.  Default is 0.5 pixels.
                )
{
   typedef typename arrayT::Scalar realT;

   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   realT r;
   realT x;
   realT y;
   realT xe, ye;
   realT rad;
   
   realT cq = cos(ang);
   realT sq = sin(ang);

   realT xr,yr;

   for(size_t i=0; i < l0; i++)
   {
      x = i-xcen;
      for(size_t j=0; j < l1; j++)
      {
         y = j-ycen;
         
         xr = x*cq - y*sq;
         yr = x*sq + y*cq;

         //Coordinate on the ellipse where the line to the point intersects
         xe = (pow(xrad*yrad,2)/ (pow(yrad,2) + pow(xrad*yr/xr,2)));
         ye = (pow(yrad,2) - xe*pow(yrad/xrad,2));
         
         rad = sqrt(xe + ye); 
         
         r = sqrt( pow(xr, 2) + pow(yr, 2) );
         
         if(r <= rad+pixbuf) m(i,j) = val;
      }
   }
}   

/// Mask a wedge in an image.
/** The wedge is describe by its center coordinate, central angle, and angular half-width. Any value can be set for the mask.
  * Pixels outside the masked circle are not altered.
  * 
  * \tparam arrayT is an Eigen-like type.
  * 
  * \ingroup image_masks
  */
template<class arrayT> 
void maskWedge( arrayT & m,                          ///< [in.out] the image to be masked, is modified.
                typename arrayT::Scalar xcen,        ///< [in] the x coordinate of the center of the circle
                typename arrayT::Scalar ycen,        ///< [in] the y coordinate of the center of the circle
                typename arrayT::Scalar angCen,      ///< [in] the central angle of the wedge, in degrees
                typename arrayT::Scalar angHW,       ///< [in] the angular half-wdith of the wedge, in degrees
                typename arrayT::Scalar val = 0      ///< [in] [optional] the mask value.  Default is 0.
              )
{
   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   typedef typename arrayT::Scalar angleT;
   
   angleT dang;
   
   for(size_t c=0; c < m.cols(); c++)
   {
      for(size_t r=0; r < m.rows(); r++)
      {
         
         dang = math::angleDiff<math::degreesT<angleT>>(math::rtod(atan2( c-ycen, r-xcen )), angCen);
         

         if( dang > -angHW && dang <= angHW)
         {
            m(r,c) = val;
         }
      }
   }
}

///Draw a thin (1-pixel) line from one point to another.
/** Calculates the line connecting the two points and sets the pixels along that line to the 
  * supplied value.
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup image_masks
  */  
template<typename realT>
int drawLine( eigenImage<realT> & mask, ///< [in/out] [pre-allocated] The array in which to draw the line.
              realT x0,                 ///< [in] The x coordinate of the first point
              realT y0,                 ///< [in] The y coordinate of the first point
              realT x1,                 ///< [in] The x coordinate of the second point
              realT y1,                 ///< [in] The y coordinate of the second point
              realT val                 ///< [in] The value to set on the pixels in the line
            )
{
   if( fabs(x1-x0) <= 1 && fabs(y1-y0) <= 1)
   {
      //If it is a single pixel, plot at the rounded average point.
      realT xa = 0.5*(x0+x1) + 0.5;
      realT ya = 0.5*(y0+y1) + 0.5;
      if( xa >=0 && xa < mask.rows() && ya >=0 && ya < mask.cols())  mask( (int) xa, (int) ya) = val;
      return 0;
   }
   
   realT _x0, _x1, _y0, _y1;
   
   realT length = sqrt( pow(x1-x0,2) + pow(y1-y0,2));
   
   //If X doesn't change enough, we try to do it as a function of y. 
   if( fabs(x1-x0) <= 1 )
   {
      if( y0 < y1)
      {
         _x0 = x0;
         _x1 = x1;
         _y0 = y0;
         _y1 = y1;
      }
      else
      {
         _x0 = x1;
         _x1 = x0;
         _y0 = y1;
         _y1 = y0;
      }
   
   
      realT m = (_x1-_x0)/(_y1-_y0);
      realT b = _x0 - m*_y0;
      
      realT dy = (_y1-_y0)/(2*length+1);
      for(realT y = _y0; y <= _y1; y += dy)
      {
         realT x = m*y + b;
      
         if( x+0.5 >=0 && x+0.5 < mask.rows() && y+0.5 >=0 && y+0.5 < mask.cols())  mask( (int) (x+0.5), (int) (y+0.5)) = val;
      }
      
      return 0;
   }
   
   //Ordinarily draw in x.
   if( x0 < x1)
   {
      _x0 = x0;
      _x1 = x1;
      _y0 = y0;
      _y1 = y1;
   }
   else
   {
      _x0 = x1;
      _x1 = x0;
      _y0 = y1;
      _y1 = y0;
   }
   
   
   realT m = (_y1-_y0)/(_x1-_x0);
   realT b = _y0 - m*_x0;
   
   realT dx = (_x1-_x0)/(2*length+1);
   for(realT x = _x0; x <= _x1; x += dx)
   {
      realT y = m*x + b;
      
      if( x+0.5 >=0 && x+0.5 < mask.rows() && y+0.5 >=0 && y+0.5 < mask.cols())  mask( (int) (x+0.5), (int) (y+0.5)) = val;
   }
}
      
      
///Draw a thick line from one point to another.
/** Calculates the line connecting the two points and sets the pixels along that line to the 
  * supplied value.  Makes the line thick by drawing lines perpindicular to each point with
  * length equal to the specified width.
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup image_masks
  */  
template<typename realT>
int drawLine( eigenImage<realT> & mask, ///< [in.out] [pre-allocated] The array in which to draw the line.
              realT x0,                 ///< [in] The x coordinate of the first point
              realT y0,                 ///< [in] The y coordinate of the first point
              realT x1,                 ///< [in] The x coordinate of the second point
              realT y1,                 ///< [in] The y coordinate of the second point
              realT width,              ///< [in] The desired full-width of the line
              realT val                 ///< [in] The value to set on the pixels in the line
            )
{
   //You can't do this.
   if( width <= 0 ) return -1;
   
   //If no reason to try this hard, don't.
   if(width < 1) return drawLine(mask, x0, y0, x1, y1, val);
   
   realT _x0, _x1, _y0, _y1;
   
   //If X doesn't change enough, we try to do it as a function of y. 
   if( fabs(x1-x0) <= 1 )
   {
      if( y0 < y1)
      {
         _x0 = x0;
         _x1 = x1;
         _y0 = y0;
         _y1 = y1;
      }
      else
      {
         _x0 = x1;
         _x1 = x0;
         _y0 = y1;
         _y1 = y0;
      }
   
   
      realT m = (_x1-_x0)/(_y1-_y0);
      realT b = _x0 - m*_y0;
      
      realT dy = (_y1-_y0)/(mask.cols()+1);
      for(realT y = _y0; y <= _y1; y += dy)
      {
         realT x = m*y + b;
      
         realT xs, ys, xe, ye; //The start and end point of the perpindicular line.
         
         realT q1 = atan(m)+0.5*math::pi<realT>();
         realT cq = cos(q1);
         realT sq = sin(q1);
         
         xs = x + 0.5*width*cq;
         ys = y + 0.5*width*sq;
         
         xe = x - 0.5*width*cq;
         ye = y - 0.5*width*sq;
         
         drawLine(mask, xs, ys, xe, ye, val);
      }
      
      return 0;
   }
   
   //Ordinarily draw in x.
   if( x0 < x1)
   {
      _x0 = x0;
      _x1 = x1;
      _y0 = y0;
      _y1 = y1;
   }
   else
   {
      _x0 = x1;
      _x1 = x0;
      _y0 = y1;
      _y1 = y0;
   }
   
   
   realT m = (_y1-_y0)/(_x1-_x0);
   realT b = _y0 - m*_x0;
   
   realT dx = (_x1-_x0)/(mask.rows()+1);
   for(realT x = _x0; x <= _x1; x += dx)
   {
      realT y = m*x + b;
      
         realT xs, ys, xe, ye; //The start and end point of the perpindicular line.
         
         realT q1 = atan(m)+0.5*math::pi<realT>();
         realT cq = cos(q1);
         realT sq = sin(q1);
         
         xs = x + 0.5*width*cq;
         ys = y + 0.5*width*sq;
         
         xe = x - 0.5*width*cq;
         ye = y - 0.5*width*sq;
         
         drawLine(mask, xs, ys, xe, ye, val);
         
   }
}








///Populate a mask based on a typical CCD bleeding pattern.
/** Masks a circle for saturated pixels, as well as a horizontal bar for bleeding (in both directions if needed).
  *
  * \returns 0 on success, -1 otherwise.
  * 
  * \tparam imT is an Eigen-like 2D array type.
  * 
  * \ingroup image_masks
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
/** The output will be a row-image containing the pixel values.
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
   for(size_t i=0;i<idx.size();++i)
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
   for(size_t i=0;i< idx.size();++i)
   {
      imout(idx[i],0) = imin(i,0);
   }
   
} 

/// Rotate a binary mask
/** Sets edge pixels to 0 or 1 depending on the interpolation, being below/above 0.5.
  */
template< typename imageT, typename transformT = cubicConvolTransform<typename imageT::Scalar>>
void rotateMask( imageT & rotMask,
                 imageT & mask,
                 typename imageT::Scalar angle
               )
{   
   imageRotate( rotMask, mask, angle, transformT());
   
   for(int jj=0; jj < rotMask.cols(); ++jj)
   {
      for(int ii=0; ii< rotMask.rows(); ++ii)   
      {
         if( rotMask(ii,jj) < 0.5) rotMask(ii,jj) = 0;
         else rotMask(ii,jj) = 1;
      }
   }
}
   



} //namespace improc 
} //namespace mx

#endif // improc_imageMasks_hpp
