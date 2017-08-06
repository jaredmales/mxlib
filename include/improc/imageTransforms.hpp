/** \file imageTransforms.hpp
 * \author Jared R. Males
 * \brief Image interpolation and transformation
 * \ingroup image_processing_files
 *
 */

#ifndef __imageTransforms_hpp__
#define __imageTransforms_hpp__

#include <cstddef>
#include <cmath>

#include <iostream>


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
/** \ingroup image_transforms 
  */
template<typename _arithT>
struct cubicConvolTransform
{
   typedef _arithT arithT;
   
   static const size_t width = 4;
   static const size_t lbuff = 1;

   arithT cubic;
   
   cubicConvolTransform()
   {
      cubic = -0.5;
   }
   
   cubicConvolTransform(arithT c)
   {
      cubic = c;
   }
   
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
  * \tparam arrOutT is the eigen array type of the output [will be resolved by compiler]
  * \tparam arrInT is the eigen array type of the input [will be resolved by compiler]
  * \tparam floatT is a floating point type [will be resolved by compiler in most cases]
  * \tparam transformT specifies the transformation to use [will be resolved by compiler]
  *
  * \param [out] transim contains the shifted image.  Must be pre-allocated.
  * \param [in] im is the image to be shifted.
  * \param [in] dq is the amount in radians to rotate in the c.c.w. direction
  * \param [in] trans is the transformation to use
  */
template<typename transformT, typename arrT, typename arrT2, typename floatT>
void imageRotate(arrT & transim, const arrT2 &im, floatT dq, transformT trans)
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
   
         
   #pragma omp parallel private(x0,y0,i0,j0,x,y) 
   {
      arithT i_x_cosq, i_x_sinq;
      arrT kern; 
      kern.resize(width,width);
   
      #pragma omp for schedule(static, 1)
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
  * \param [out] out contains the shifted image.  Must be pre-allocated, but can be smaller than the in array.
  * \param [in] in is the image to be shifted.
  * \param [in] dx is the amount to shift in the x direction
  * \param [in] dy is the amount to shift in the y direction
  *
  * \tparam outputArrT is the eigen array type of the output [will be resolved by compiler]
  * \tparam inputArrT is the eigen array type of the input [will be resolved by compiler]
  *
  */
template<typename outputArrT, typename inputArrT> 
void imageShiftWP( outputArrT & out,
                   inputArrT & in,
                   int dx,
                   int dy )
{
   dx %= in.rows();
   dy %= in.cols();

   #pragma omp parallel num_threads(4)
   {
      int x, y;
   
      #pragma omp for   
      for(int i=0;i<out.rows(); ++i)
      {
         x = i - dx;
         
         if(x < 0) x += in.rows();
         else if (x >= in.rows()) x -= in.rows();

         for(int j=0; j<out.cols(); ++j)
         {
            y = j - dy;
            if (y < 0) y += in.cols();
            else if (y >= in.cols()) y -= in.cols();
         
            out(i, j) = in(x,y);
         }
      }
   }
}

/// Shift an image.
/** Uses the given transformation type to shift an image.
  *
  * \tparam arrOutT is the eigen array type of the output [will be resolved by compiler]
  * \tparam arrInT is the eigen array type of the input [will be resolved by compiler]
  * \tparam floatT is a floating point type [will be resolved by compiler in most cases]
  * \tparam transformT specifies the transformation to use [will be resolved by compiler]
  * 
  * \param [out] transim contains the shifted image.  Must be pre-allocated.
  * \param [in] im is the image to be shifted.
  * \param [in] dx is the amount to shift in the x direction
  * \param [in] dy is the amount to shift in the y direction
  * \param [in] trans is the transformation to use
  * 
  */
template<typename arrOutT, typename arrInT, typename floatT, typename transformT>
void imageShift( arrOutT & transim, 
                 const arrInT  &im, 
                 floatT dx, 
                 floatT dy, 
                 transformT trans )
{
   typedef typename transformT::arithT arithT;
   
   arithT x0, y0, x,y;
   //arithT xcen, ycen;
   
   int Nrows, Ncols;
   
   int i0, j0;
    
   const int lbuff = transformT::lbuff;
   const int width = transformT::width;
   
   Nrows = im.rows();
   Ncols = im.cols();

   
   transim.resize(Nrows, Ncols);
   
   //The geometric image center
   //xcen = 0.5*(Nrows-1.);       
   //ycen = 0.5*(Ncols-1.);
   
   int xulim = Nrows-width+lbuff;// - 1;
   int yulim = Ncols-width+lbuff;// - 1;
  
         
   #pragma omp parallel private(x0,y0,i0,j0,x,y) num_threads(4)
   {
      arrOutT kern; 
      kern.resize(width,width);
      
      #pragma omp for 
      for(int i=0;i<Nrows; ++i)
      {
         // (i,j) is position in new image
         // (x0,y0) is true position in old image
         // (i0,j0) is integer position in old image
         // (x, y) is fractional residual of (x0-i0, y0-j0)

         x0 = i-dx;
         i0 = x0; //just converting to int
            
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

            y0 = j-dy;
            j0 = y0;
            
            if(j0 <= lbuff || j0 >= yulim) 
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
   }//#pragam omp
         
} //void imageShift(arrOutT & transim, const arrInT  &im, floatT dx, floatT dy)


/// Magnify an image.
/** Uses the given transformation type to magnify the input image to the size of the output image.
  *
  * Here we assume that the image center is the mxlib standard: 
  * \code
    x_center = 0.5*(im.rows()-1);
    y_center = 0.5*(im.cols()-1);
  * \endcode
  * Some care is necessary to prevent magnification from shifting with respect to this center.  The main result is that the 
  * magnification factors (which can be different in x and y) are defined thus:
  * \code
    x_mag = (transim.rows()-1.0) / (im.rows()-1.0);
    y_mag = (transim.cols()-1.0) / (im.cols()-1.0);
  * \endcode
  * 
  * \tparam arrOutT is the eigen array type of the output. 
  * \tparam arrInT is the eigen array type of the input.
  * \tparam transformT specifies the transformation to use.
  */
template<typename arrOutT, typename arrInT, typename transformT>
void imageMagnify( arrOutT & transim, ///< [out] contains the magnified image.  Must be pre-allocated.
                   const arrInT  &im, ///< [in] is the image to be magnified.
                   transformT trans ///< [in] is the transformation to use
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
      
//      #pragma omp for 
      for(int i=0;i<Nrows; ++i)
      {
         // (i,j) is position in new image
         // (x0,y0) is true position in old image
         // (i0,j0) is integer position in old image
         // (x, y) is fractional residual of (x0-i0, y0-j0)

         x0 = xcen0 + (i-xcen)*x_scale;
         i0 = x0; //just converting to int
            
         if(i0 < lbuff || i0 >= xulim) 
         {
            for(int j=0;j<Ncols; ++j)
            {
               transim(i,j) = 0;
            }
            continue;
         }
            
         for(int j=0;j<Ncols; ++j)
         {

            y0 = ycen0 + (j-ycen)*y_scale; 
            j0 = y0;
            
            if(j0 < lbuff || j0 >= yulim) 
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
   }//#pragam omp
         
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
      Scalar outputTotal = fabs(imout.sum());
      
      //Normalize
      imout *= 1.0;//(1./    inputTotal/outputTotal;
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




#endif //__imageTransforms_hpp__
