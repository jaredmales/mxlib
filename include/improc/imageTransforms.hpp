/** \file imageTransforms.hpp
 * \author Jared R. Males
 * \brief Image interpolation and transformation
 * \ingroup image_processing_files
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

#ifndef improc_imageTransforms_hpp
#define improc_imageTransforms_hpp

#include <cstddef>
#include <cmath>

#include <iostream>
#include <limits>


namespace mx
{
namespace improc
{

///Transformation by bi-linear interpolation
/** \ingroup image_transforms
  */
template<typename _arithT>
struct bilinearTransform
{
   typedef _arithT arithT;

   static const size_t width = 2;
   static const size_t lbuff = 0;

   template<typename arrT, typename arithT>
   void operator()(arrT & kern, arithT x, arithT y)
   {
      kern.resize(width, width);

      kern(0,0) = (1.-x)*(1.-y);
      kern(0,1) = (1.-x)*y;
      kern(1,0) = x*(1.-y);
      kern(1,1) = x*y;
   }
};

///Typedef for bilinearTransform with single precision
/** \ingroup image_transforms
  */
typedef bilinearTransform<float> bilinearTransf;

///Typedef for bilinearTransform with double precision
/** \ingroup image_transforms
  */
typedef bilinearTransform<double> bilinearTransd;


///Transformation by cubic convolution interpolation
/** Uses the cubic convolution interpolation kernel.  See <a href="https://en.wikipedia.org/wiki/Bicubic_interpolation">https://en.wikipedia.org/wiki/Bicubic_interpolation </a>.
  *
  * The parameter \ref cubic should be left as the default -0.5 in most cases, which gives the bicubic spline interpolator.  See
  * <a href="https://en.wikipedia.org/wiki/Cubic_Hermite_spline">https://en.wikipedia.org/wiki/Cubic_Hermite_spline</a>.
  *
  * \tparam _arithT is the type in which to do all calculations.  Should be a floating point type.
  *
  * \ingroup image_transforms
  */
template<typename _arithT>
struct cubicConvolTransform
{
   typedef _arithT arithT; ///< The type in which all calculations are performed.

   static const size_t width = 4;
   static const size_t lbuff = 1;

   arithT cubic {-0.5}; ///< The kernel parameter.  The default value -0.5 gives the bicubic spline interpolator.

   explicit cubicConvolTransform(arithT c)
   {
      cubic = c;
   }

   cubicConvolTransform() {}

   cubicConvolTransform(const cubicConvolTransform & t)
   {
      cubic = t.cubic;
   }

   arithT cubicConvolKernel(arithT d)
   {
      if(d <= 1) return (cubic+2.)*d*d*d - (cubic+3.)*d*d + 1.;

      if(d < 2) return cubic*d*d*d -5.*cubic*d*d + 8.*cubic*d - 4.*cubic;

      return 0;
   }

   template<typename arrT, typename arithT>
   void operator()(arrT & kern, arithT x, arithT y)
   {
      arithT km2x,km1x,kp1x,kp2x;
      arithT km2y,km1y,kp1y,kp2y;

      km2x = cubicConvolKernel((1.+x));
      km1x = cubicConvolKernel(x);
      kp1x = cubicConvolKernel(1.-x);
      kp2x = cubicConvolKernel(2.-x);

      km2y = cubicConvolKernel((1.+y));
      km1y = cubicConvolKernel(y);
      kp1y = cubicConvolKernel(1.-y);
      kp2y = cubicConvolKernel(2.-y);

      kern(0,0) = km2x*km2y;
      kern(0,1) = km2x*km1y;
      kern(0,2) = km2x*kp1y;
      kern(0,3) = km2x*kp2y;

      kern(1,0) = km1x*km2y;
      kern(1,1) = km1x*km1y;
      kern(1,2) = km1x*kp1y;
      kern(1,3) = km1x*kp2y;

      kern(2,0) = kp1x*km2y;
      kern(2,1) = kp1x*km1y;
      kern(2,2) = kp1x*kp1y;
      kern(2,3) = kp1x*kp2y;

      kern(3,0) = kp2x*km2y;
      kern(3,1) = kp2x*km1y;
      kern(3,2) = kp2x*kp1y;
      kern(3,3) = kp2x*kp2y;
   }
};




///Typedef for cubicConvolTransform with single precision
/** \ingroup image_transforms
  */
typedef cubicConvolTransform<float> cubicConvolTransf;

///Typedef for cubicConvolTransform with double precision
/** \ingroup image_transforms
  */
typedef cubicConvolTransform<double> cubicConvolTransd;

/** \ingroup image_transforms
  * @{
  */


/// Rotate an image represented as an eigen array
/** Uses the given transformation type to rotate an image.
  *
  * \tparam transformT specifies the transformation to use [will be resolved by compiler]
  * \tparam arrT is the eigen array type of the output [will be resolved by compiler]
  * \tparam arrT2 is the eigen array type of the input [will be resolved by compiler]
  * \tparam floatT is a floating point type [will be resolved by compiler in most cases]
  *
  */
template<typename transformT, typename arrT, typename arrT2, typename floatT>
void imageRotate( arrT & transim,  ///< [out] The rotated image.  Must be pre-allocated.
                  const arrT2 &im, ///< [in] The image to be rotated.
                  floatT dq,       ///< [in] the angle, in radians, by which to rotate in the c.c.w. direction
                  transformT trans ///< [in] is the transformation to use
                )
{
   typedef typename transformT::arithT arithT;
   arithT cosq, sinq;
   arithT x0, y0, x,y;
   arithT xcen, ycen;

   int Nrows, Ncols;

   int i0, j0;

   const int lbuff = transformT::lbuff;
   const int width = transformT::width;

   cosq = cos(dq);
   sinq = sin(dq);

   Nrows = im.rows();
   Ncols = im.cols();


   transim.resize(Nrows, Ncols);

   //The geometric image center
   xcen = 0.5*(Nrows-1.);
   ycen = 0.5*(Ncols-1.);

   int xulim = Nrows-width+lbuff;// - 1;
   int yulim = Ncols-width+lbuff;// - 1;

   arithT xc_x_cosq = xcen*cosq;
   arithT xc_x_sinq = xcen*sinq;
   arithT yc_x_cosq = ycen*cosq;
   arithT yc_x_sinq = ycen*sinq;

   xc_x_cosq += yc_x_sinq;
   xc_x_sinq -= yc_x_cosq;


   #ifdef MXLIB_USE_OMP
   #pragma omp parallel private(x0,y0,i0,j0,x,y)
   #endif
   {
      arithT i_x_cosq, i_x_sinq;
      arrT kern;
      kern.resize(width,width);

      #ifdef MXLIB_USE_OMP
      #pragma omp for schedule(static, 1)
      #endif
      for(int i=0;i<Nrows; ++i)
      {
         i_x_cosq = i*cosq - xc_x_cosq;// + xcen;
         i_x_sinq = -(i*sinq - xc_x_sinq);// + ycen;

         for(int j=0;j<Ncols; ++j)
         {
            //We are actually doing this rotation matrix:
            //x0 =  (i-xcen)*cosq + (j-ycen)*sinq;
            //y0 = -(i-xcen)*sinq + (j-ycen)*cosq;
            //This is the minimum-op representation of the above rotation matrix:
            x0 =  i_x_cosq + j*sinq;
            y0 =  i_x_sinq + j*cosq;

            //Get lower left index
            i0 = x0 +xcen;
            j0 = y0 +ycen;

            if(i0 <= lbuff || i0 >= xulim || j0 <= lbuff || j0 >= yulim)
            {
               transim(i,j) = 0;
               continue;
            }

            //Get the residual
            x = x0+xcen-i0;
            y = y0+ycen-j0;

            trans(kern, x, y);
            transim(i,j) = (im.block(i0-lbuff,j0-lbuff, width, width) * kern).sum();
         }//for j
      }//for i
   }//#pragma omp parallel

}//void imageRotate(arrT & transim, const arrT2  &im, floatT dq, transformT trans)

/// Shift an image by whole pixels, wrapping around..
/** The output image can be smaller than the input image, in which case the wrapping still occurs for the input image, but only
  * output images worth of pixels are actually shifted.  This is useful, for instance, when propagating large turbulence phase screens
  * where one only needs a small section at a time.
  *
  * \tparam outputArrT is the eigen array type of the output [will be resolved by compiler]
  * \tparam inputArrT is the eigen array type of the input [will be resolved by compiler]
  *
  */
template<typename outputArrT, typename inputArrT>
void imageShiftWP( outputArrT & out,  ///< [out] contains the shifted image.  Must be pre-allocated, but can be smaller than the in array.
                   inputArrT & in,    ///< [in] the image to be shifted.
                   int dx,            ///< [in] the amount to shift in the x direction
                   int dy             ///< [in] the amount to shift in the y direction
                 )
{
   dx %= in.rows();
   dy %= in.cols();

   int outr = out.rows();
   int outc = out.cols();
   int inr = in.rows();
   int inc = in.cols();
   #ifdef MXLIB_USE_OMP
   #pragma omp parallel
   #endif
   {
      int x, y;

      #ifdef MXLIB_USE_OMP
      #pragma omp for
      #endif
      for(int cc=0;cc<outc; ++cc)
      {
         y = cc - dy;

         if (y < 0) y += inc;
         else if (y >= inc) y -= inc;
         
         for(int rr=0; rr<outr; ++rr)
         {
            x = rr - dx;
            
            if(x < 0) x += inr;
            else if (x >= inr) x -= inr;

            out(rr, cc) = in(x,y);
         }
      }
      
   }
}

/// Shift an image by whole pixels, wrapping around, with a scaling image applied to the shifted image.
/** The output image can be smaller than the input image, in which case the wrapping still occurs for the input image, but only
  * output images worth of pixels are actually shifted.  This is useful, for instance, when propagating large turbulence phase screens
  * where one only needs a small section at a time.
  *
  * The scaling is applied to the output image.  The scale image must be the same size as the output image.
  * 
  * \tparam outputArrT is the eigen array type of the output [will be resolved by compiler]
  * \tparam inputArrT is the eigen array type of the input [will be resolved by compiler]
  * \tparam scaleArrT is the eigen array type of the scale image [will be resolved by compiler]
  *
  */
template<typename outputArrT, typename inputArrT, typename scaleArrT>
void imageShiftWP( outputArrT & out,  ///< [out] contains the shifted image.  Must be pre-allocated, but can be smaller than the in array.
                   inputArrT & in,    ///< [in] the image to be shifted.
                   scaleArrT & scale, ///< [in] image of scale values applied per-pixel to the output (shifted) image, same size as out
                   int dx,            ///< [in] the amount to shift in the x direction
                   int dy             ///< [in] the amount to shift in the y direction
                 )
{
   dx %= in.rows();
   dy %= in.cols();

   int outr = out.rows();
   int outc = out.cols();
   int inr = in.rows();
   int inc = in.cols();
   #ifdef MXLIB_USE_OMP
   //#pragma omp parallel
   #endif
   {
      int x, y;

      #ifdef MXLIB_USE_OMP
      //#pragma omp for
      #endif
      for(int cc=0;cc<outc; ++cc)
      {
         y = cc - dy;

         if (y < 0) y += inc;
         else if (y >= inc) y -= inc;
         
         for(int rr=0; rr<outr; ++rr)
         {
            x = rr - dx;
            
            if(x < 0) x += inr;
            else if (x >= inr) x -= inr;

            out(rr, cc) = in(x,y)*scale(rr,cc);
         }
      }
      
   }
}

/// Shift an image.
/** Uses the given transformation type to shift an image.  
  * 
  * Note that this does not treat the edges
  * of the image, determined by the buffer width (lbuff) of the kernel and the size of shift.  If you wish to 
  * treat the edges, you must pad the image by at least lbuff+abs(shift) pixels in each direction, and
  * implement a strategy (zeros, mirror, wrap) prior to calling this function.
  * 
  * \tparam arrOutT is the Eigen-like array type of the output [will be resolved by compiler]
  * \tparam arrInT is the Eigen-like array type of the input [will be resolved by compiler]
  * \tparam floatT1 is a floating point type [will be resolved by compiler]
  * \tparam floatT2 is a floating point type [will be resolved by compiler]
  * \tparam transformT specifies the transformation to use [will be resolved by compiler]
  *
  */
template<typename arrOutT, typename arrInT, typename floatT1, typename floatT2, typename transformT>
void imageShift( arrOutT & transim, ///< [out] Will contain the shifted image.  Will be allocated.
                 const arrInT  &im, ///< [in] the image to be shifted.
                 floatT1 dx,        ///< [in] the amount to shift in the x direction
                 floatT2 dy,        ///< [in] the amount to shift in the y direction
                 transformT trans   ///< [in] trans is the transformation to use
               )
{
   typedef typename transformT::arithT arithT;

   int Nrows, Ncols;

   const int lbuff = transformT::lbuff;
   const int width = transformT::width;

   Nrows = im.rows();
   Ncols = im.cols();

   int xulim = Nrows-width+lbuff;
   int yulim = Ncols-width+lbuff;

   transim.resize(Nrows, Ncols);

   #ifdef MXLIB_USE_OMP
   #pragma omp parallel
   #endif
   {
      int i0, j0;
      // (rx, ry) is fractional residual of shift
      arithT rx = 1-(dx-floor(dx));    
      arithT ry = 1-(dy-floor(dy));
   
      arrOutT kern;
      kern.resize(width,width);
      trans(kern, rx, ry);

      #ifdef MXLIB_USE_OMP
      #pragma omp for
      #endif
      for(int i=0;i<Nrows; ++i)
      {
         // (i,j) is position in new image
         // (i0,j0) is integer position in old image

         i0 = i-dx;

         if(i0 <= lbuff || i0 >= xulim)
         {
            for(int j=0;j<Ncols; ++j)
            {
               transim(i,j) = 0;
            }
            continue;
         }

         for(int j=0;j<Ncols; ++j)
         {
            j0 = j-dy;

            if(j0 <= lbuff || j0 >= yulim)
            {
               transim(i,j) = 0;
               continue;
            }

            transim(i,j) = (im.block(i0-lbuff,j0-lbuff, width, width) * kern).sum();
         }//for j
      }//for i
   }//#pragam omp

} //imageShift


/// Magnify an image.
/** Uses the given transformation type to magnify the input image to the size of the output image.
  *
  * Here we assume that the image center is the mxlib standard:
  * \code
    x_center = 0.5*(im.rows()-1);
    y_center = 0.5*(im.cols()-1);
  * \endcode
  * Some care is necessary to prevent magnification from shifting the image with respect to this center.  The main result is that the
  * magnification factors (which can be different in x and y) are defined thus:
  * \code
    x_mag = (transim.rows()-1.0) / (im.rows()-1.0);
    y_mag = (transim.cols()-1.0) / (im.cols()-1.0);
  * \endcode
  *
  * Example:
  * \code
    im1.resize(512,512);
    //add image to im1
    im2.resize(1024,1024);
    imageMagnify(im2,im1, cubicConvolTransform<double>());
    \endcode
  * In this exmple, the image in im1 will be magnified by `1023.0/511.0 = 2.002x` and placed in im2.
  *
  * This transform function does not handle edges.  If treatment of edges is desired, you must pad the input
  * image using the desired strategy before calling this function.  Note that the padded-size of the input image
  * will affect the magnification factor.
  *
  * \tparam arrOutT is the eigen array type of the output.
  * \tparam arrInT is the eigen array type of the input.
  * \tparam transformT specifies the transformation to use.
  */
template<typename arrOutT, typename arrInT, typename transformT>
void imageMagnify( arrOutT & transim, ///< [out] contains the magnified image.  Must be pre-allocated.
                   const arrInT  &im, ///< [in] is the image to be magnified.
                   transformT trans   ///< [in] is the transformation to use
                 )
{
   typedef typename transformT::arithT arithT;

   arithT x0, y0, x,y;

   int Nrows, Ncols;

   int i0, j0;

   const int lbuff = transformT::lbuff;
   const int width = transformT::width;

   Nrows = transim.rows();
   Ncols = transim.cols();

   int xulim = im.rows()-lbuff - 1;
   int yulim = im.cols()-lbuff - 1;

   arithT x_scale = ( (arithT) im.rows()-1.0)/ (transim.rows()-1.0); //this is 1/x_mag
   arithT y_scale = ( (arithT) im.cols()-1.0)/ (transim.cols()-1.0); //this is 1/y_mag

   arithT xcen = 0.5*( (arithT) transim.rows() - 1.0);
   arithT ycen = 0.5*( (arithT) transim.cols() - 1.0);

   arithT xcen0 = 0.5*( (arithT) im.rows() - 1.0);
   arithT ycen0 = 0.5*( (arithT) im.cols() - 1.0);

//   #pragma omp parallel private(x0,y0,i0,j0,x,y) num_threads(4)
   {
      arrOutT kern;
      kern.resize(width,width);

      for(int j=0;j<Ncols; ++j)
      {
         // (i,j) is position in new image
         // (x0,y0) is true position in old image
         // (i0,j0) is integer position in old image
         // (x, y) is fractional residual of (x0-i0, y0-j0)
        
         y0 = ycen0 + (j-ycen)*y_scale;
         j0 = y0;

         if(j0 < lbuff || j0 >= yulim)
         {
            for(int i=0;i<Nrows; ++i)
            {
               transim(i,j) = 0;
            }
            continue;
         }

//       #pragma omp for
         for(int i=0;i<Nrows; ++i)
         {
            x0 = xcen0 + (i-xcen)*x_scale;
            i0 = x0; //just converting to int
        
            if(i0 < lbuff || i0 >= xulim)
            {
               transim(i,j) = 0;
               continue;
            }


            //Get the residual
            x = x0-i0;
            y = y0-j0;

            trans(kern, x, y);
            transim(i,j) = (im.block(i0-lbuff,j0-lbuff, width, width) * kern).sum();
         }//for j
      }//for i
   }//#pragma omp

}

/// Magnify an image with the cubic convolution interpolator.
/** Uses the cubic convolution interpolator to magnify the input image to the size of the output image.
  * 
  * This is a wrapper for imageMagnify with the transform type specified.
  *
  * \tparam arrOutT is the eigen array type of the output.
  * \tparam arrInT is the eigen array type of the input.
  * 
  * \overload
  */
template<typename arrOutT, typename arrInT>
void imageMagnify( arrOutT & transim, ///< [out] contains the magnified image.  Must be pre-allocated.
                   const arrInT  &im  ///< [in] is the image to be magnified.
                 )
{
   return imageMagnify(transim, im, cubicConvolTransform<typename arrInT::Scalar>());
}

/// Re-bin an image using the sum, reducing its size while conserving the total flux.
/** Optionally this can be the mean instead of the sum filter, in which case total flux is not conserved.
  */
template<typename imageOutT, typename imageInT>
int imageRebinSum( imageOutT & imout,     ///< [out] the re-binned image.  Must be allocated to size which is an integer factor smaller than imin.
                   const imageInT & imin, ///< [in] the image to rebin
                   bool mean = false      ///< [in] if true the output is the mean rather than the sum.
                 )
{
   int rebin = imin.rows()/imout.rows();
   if(imin.cols()/imout.cols() != rebin) return -1;

   int N = 1;
   if(mean) N = rebin*rebin;
   for(int i=0;i<imout.rows(); ++i)
   {
      for(int j=0; j<imout.cols(); ++j)
      {
         imout(i,j) = imin.block( i*rebin, j*rebin, rebin, rebin).sum()/N;
      }
   }
   
   return 0;
}

/// Re-bin an image using the sum, reducing its size while conserving the total flux.  Records the value and position of the re-binned max pixel.
/** Optionally this can be the mean instead of the sum filter, in which case total flux is not conserved.
  * 
  * \overload 
  */
template<typename imageOutT, typename imageInT>
int imageRebinSum( imageOutT & imout,                 ///< [out] the re-binned image.  Must be allocated to size which is an integer factor smaller than imin.
                   int & xMax,                        ///< [out] the x-locatioin of the max pixel
                   int & yMax,                        ///< [out] the y-locatioin of the max pixel
                   typename imageOutT::Scalar & pMax, ///< [out] the value of the max pixel
                   const imageInT & imin,             ///< [in] the image to rebin
                   bool mean = false                  ///< [in] if true the output is the mean rather than the sum.
                 )
{
   int rebin = imin.rows()/imout.rows();
   if(imin.cols()/imout.cols() != rebin) return -1;

   int N = 1;
   if(mean) N = rebin*rebin;
   
   xMax = 0;
   yMax = 0;
   pMax = std::numeric_limits<typename imageOutT::Scalar>::lowest();
   
   for(int i=0;i<imout.rows(); ++i)
   {
      for(int j=0; j<imout.cols(); ++j)
      {
         imout(i,j) = imin.block( i*rebin, j*rebin, rebin, rebin).sum()/N;
         if(imout(i,j) > pMax)
         {
            pMax = imout(i,j);
            xMax = i;
            yMax = j;
         }
      }
      
   }

   return 0;
}

/// Re-bin an image using the mean.
/** This is a wrapper for imageRebinSum with `mean=true`.
  */
template<typename imageOutT, typename imageInT>
int imageRebinMean( imageOutT & imout,    ///< [out] the re-binned image.  Must be allocated to size which is an integer factor smaller than imin.
                    const imageInT & imin ///< [in] the image to rebin
                  )
{
   return imageRebinSum(imout, imin, true);
}

/// Re-bin an image using the mean.  Records the value and position of the re-binned max pixel.
/** This is a wrapper for imageRebinSum with `mean=true`.
  *
  * \overload 
  */
template<typename imageOutT, typename imageInT>
int imageRebinMean( imageOutT & imout,                 ///< [out] the re-binned image.  Must be allocated to size which is an integer factor smaller than imin.
                    int & xMax,                        ///< [out] the x-locatioin of the max pixel
                    int & yMax,                        ///< [out] the y-locatioin of the max pixel
                    typename imageOutT::Scalar & pMax, ///< [out] the value of the max pixel
                    const imageInT & imin,             ///< [in] the image to rebin
                    bool mean = false                  ///< [in] if true the output is the mean rather than the sum.
                  )
{
   return imageRebinSum(imout, xMax, yMax, pMax, imin, true);
}

/// Re-bin an image, takes the mean with a min/max rejection.
/** The mean is calculated after rejecting the minimuma and maximum value.
  */
template<typename imageOutT, typename imageInT>
int imageRebinMeanReject( imageOutT & imout,    ///< [out] the re-binned image.  Must be allocated to size which is an integer factor smaller than imin.
                          const imageInT & imin ///< [in] the image to rebin
                        )
{
   int rebin = imin.rows()/imout.rows();
   if(imin.cols()/imout.cols() != rebin) return -1;

   int N = rebin*rebin - 2;
   
   for(int i=0;i<imout.rows(); ++i)
   {
      for(int j=0; j<imout.cols(); ++j)
      {
         register typename imageOutT::Scalar sum = 0;
         register typename imageOutT::Scalar max = imin(i*rebin, j*rebin);
         register typename imageOutT::Scalar min = imin(i*rebin, j*rebin);
         for(int k=0;k<rebin;++k)
         {
            for(int l=0;l<rebin;++l)
            {
               sum += imin(i*rebin+k, j*rebin+l);
               if(imin(i*rebin+k, j*rebin+l) > max) max = imin(i*rebin+k, j*rebin+l);
               if(imin(i*rebin+k, j*rebin+l) < min) min = imin(i*rebin+k, j*rebin+l);
               
            }
         }
         imout(i,j) = (sum-max-min)/N;/**/
      }
   }
   
   return 0;
}

/// Re-bin an image, takes the mean with a min/max rejection.  Records the value and position of the re-binned max pixel.
/** The mean is calculated after rejecting the minimuma and maximum value.
  * 
  * \overload 
  */
template<typename imageOutT, typename imageInT>
int imageRebinMeanReject( imageOutT & imout,                 ///< [out] the re-binned image.  Must be allocated to size which is an integer factor smaller than imin.
                          int & xMax,                        ///< [out] the x-locatioin of the max pixel
                          int & yMax,                        ///< [out] the y-locatioin of the max pixel
                          typename imageOutT::Scalar & pMax, ///< [out] the value of the max pixel
                          const imageInT & imin              ///< [in] the image to rebin
                        )
{
   int rebin = imin.rows()/imout.rows();
   if(imin.cols()/imout.cols() != rebin) return -1;

   int N = rebin*rebin - 2;
   
   xMax = 0;
   yMax = 0;
   pMax = std::numeric_limits<typename imageOutT::Scalar>::lowest();
   
   for(int i=0;i<imout.rows(); ++i)
   {
      for(int j=0; j<imout.cols(); ++j)
      {
         register typename imageOutT::Scalar sum = 0;
         register typename imageOutT::Scalar max = imin(i*rebin, j*rebin);
         register typename imageOutT::Scalar min = imin(i*rebin, j*rebin);
         for(int k=0;k<rebin;++k)
         {
            for(int l=0;l<rebin;++l)
            {
               sum += imin(i*rebin+k, j*rebin+l);
               if(imin(i*rebin+k, j*rebin+l) > max) max = imin(i*rebin+k, j*rebin+l);
               if(imin(i*rebin+k, j*rebin+l) < min) min = imin(i*rebin+k, j*rebin+l);
               
            }
         }
         imout(i,j) = (sum-max-min)/N;
         
         if(imout(i,j) > pMax)
         {
            pMax = imout(i,j);
            xMax = i;
            yMax = j;
         }
      }
   }
   
   return 0;
}

/// Down-sample an image, reducing its size while conserving the total flux.
/** If the old size is an integer multiple of the new size, this is just a re-bin.  If not an integer multiple,
  * the image is interpolated after performing the closest re-bin, and then re-normalized to conserve flux.
  *
  * \todo Allow selection of interpolator, providing a default version.
  */
template<typename imageOutT, typename imageInT>
void imageDownSample(imageOutT & imout, const imageInT & imin)
{
   typedef typename imageOutT::Scalar Scalar;


   //Record this for normalization later
   Scalar inputTotal = fabs(imin.sum());


   //As a first step, rebin to nearest whole pixel factor which is larger than the desired output size
   int closestRebin = imin.rows()/imout.rows();//, imin.cols()/imout.cols() );

   float sample = ( (float) imin.rows())/ closestRebin;

   while(  sample != floor(sample))
   {
      --closestRebin;
      if(closestRebin == 1) break;
      sample = ( (float) imin.rows())/ closestRebin;
   }

   //Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> temp;
   imageOutT temp;
   temp.resize( imin.rows()/closestRebin, imin.cols()/closestRebin);

   for(int i=0;i<temp.rows(); ++i)
   {
      for(int j=0; j<temp.cols(); ++j)
      {
         temp(i,j) = imin.block( i*closestRebin, j*closestRebin, closestRebin, closestRebin).sum();
      }
   }

   //If the output image is now the requested size return.
   if(temp.rows() == imout.rows() && temp.cols() == imout.cols())
   {
      imout = temp;
      return;
   }
   //Otherwise, re-sample using bilinear interpolation.
   typedef bilinearTransform<Scalar> transformT;

   transformT trans;
   //Eigen::Array<Scalar, -1,-1> kern;
   imageOutT kern;

   const int lbuff = transformT::lbuff;
   const int width = transformT::width;

   for(int i=0;i<imout.rows(); ++i)
   {
      for(int j=0;j<imout.cols(); ++j)
      {
         double x = ( (double) i/ imout.rows())*temp.rows();
         double y = ( (double) j/ imout.cols())*temp.cols();

         trans(kern, x-floor(x), y-floor(y));

         imout(i,j) = (temp.block( floor(x)-lbuff, floor(y)-lbuff, width, width)*kern).sum();

      }
   }

   //Normalize
   Scalar outputTotal = fabs(imout.sum());
   imout *= inputTotal/outputTotal;
}

///@}

} //namespace improc
} //namespace mx




#endif //improc_imageTransforms_hpp
