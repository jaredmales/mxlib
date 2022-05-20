/** \file imageFilters.hpp
  * \brief Image filters (smoothing, radial profiles, etc.)
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __imageFilters_hpp__
#define __imageFilters_hpp__

#include <cstdlib>

#include "../math/gslInterpolator.hpp"
#include "../math/vectorUtils.hpp"
#include "../math/geo.hpp"

#include "imageMasks.hpp"

namespace mx
{
namespace improc 
{
   

///Symetric Gaussian smoothing kernel
/** \ingroup image_filters_kernels
  * 
  */
template<typename _arrayT, size_t _kernW=4>
struct gaussKernel
{
   typedef _arrayT arrayT;
   typedef typename _arrayT::Scalar arithT;
   static const int kernW = _kernW;
   
   arrayT kernel;
   
   arithT _fwhm;
   
   explicit gaussKernel(arithT fwhm)
   {
      _fwhm = fwhm;
      
      int w = kernW*_fwhm;
      
      if(w % 2 == 0) w++;
      
      kernel.resize(w, w);

      arithT kcen = 0.5*(w-1.0);
      
      arithT sig2 = _fwhm/2.354820045030949327;
      sig2 *= sig2;
   
      arithT r2;
      for(int i=0; i < w; ++i)
      {
         for(int j=0; j < w; ++j)
         {
            r2 = pow(i-kcen,2) + pow(j-kcen,2) ;
            kernel(i,j) = exp(-r2/(2.0*sig2));
         }
      }
      
      kernel /= kernel.sum();
   }
   
   int maxWidth()
   {
      return _kernW*_fwhm;
   }
   
   void setKernel(arithT x, arithT y, arrayT & kernelArray)
   {
      //Unused parts of interface:
      static_cast<void>(x);
      static_cast<void>(y);
      
      kernelArray = kernel;
   }
   
};

/// Azimuthally variable boxcare kernel.
/** Averages the image in a boxcare defined by a radial and azimuthal extent.
  * 
  * \ingroup image_filters_kernels
  */
template<typename _arrayT, size_t _kernW=2>
struct azBoxKernel
{
   typedef _arrayT arrayT;
   typedef typename _arrayT::Scalar arithT;

   static const int kernW = _kernW;
   
   arithT m_radWidth {0}; ///< the half-width of the averaging box, in the radial direction, in pixels.
   arithT m_azWidth {0};  ///< the half-width of the averaging box, in the azimuthal direction, in pixels.
   arithT m_maxAz {0};    ///< the maximum half-width of the averging box in the azimuthal direction, in degrees. \>= 0. If 0 or \>= 180, then no maximum is enforced.

   int m_maxWidth;
   
   azBoxKernel( arithT radWidth, ///< [in] the half-width of the averaging box, in the radial direction, in pixels.
                arithT azWidth   ///< [in] the half-width of the averaging box, in the azimuthal direction, in pixels.
              ) : m_radWidth(fabs(radWidth)), m_azWidth(fabs(azWidth))
   {      
      setMaxWidth();
   }

   azBoxKernel( arithT radWidth, ///< [in] the half-width of the averaging box, in the radial direction, in pixels.
                arithT azWidth,  ///< [in] the half-width of the averaging box, in the azimuthal direction, in pixels.
                arithT maxAz     ///< [in] the maximum half-width of the averging box in the azimuthal direction, in degrees. \>= 0. If 0 or \>= 180, then no maximum is enforced.
              ) : m_radWidth(fabs(radWidth)), m_azWidth(fabs(azWidth))
   {      
      setMaxWidth();

      maxAz = fabs(maxAz);
      if(maxAz >= 180) maxAz = 0; //Larger than 180 means no limit.

      m_maxAz = math::dtor( maxAz );
      
   }

   /// Sets the max width based on the configured az and rad widths.
   void setMaxWidth()
   {
      // atan is where derivative of width/height is 0.

      arithT qmax = atan(m_azWidth/m_radWidth);
      arithT mx1 = kernW*( (int) ( fabs(m_azWidth*sin(qmax)) + fabs(m_radWidth*cos(qmax)) ) + 1 );
      
      qmax = atan(m_radWidth/m_azWidth);
      arithT mx2 = kernW*( (int) ( fabs(m_azWidth*cos(qmax)) + fabs(m_radWidth*sin(qmax)) ) + 1 );

      m_maxWidth = 0.5*std::max(mx1, mx2);
   }

   int maxWidth()
   {
      return m_maxWidth;
   }
   
   void setKernel(arithT x, arithT y, arrayT & kernel)
   {
      arithT rad0 = sqrt((arithT) (x*x + y*y));
      
      arithT sinq = y/rad0;
      arithT cosq = x/rad0;
      
      //Only calc q if we're going to use it.
      arithT q = 0;
      if(m_maxAz > 0) q = atan2(sinq, cosq);

      int w = kernW*( (int) ( fabs(m_azWidth*sinq) + fabs(m_radWidth*cosq) ) + 1 );
      int h = kernW*( (int) ( fabs(m_azWidth*cosq) + fabs(m_radWidth*sinq) ) + 1 );

      if(w > m_maxWidth*2)
      {
         std::stringstream errs;
         errs << "|" << kernW << "|" << m_radWidth << "|" << m_azWidth << "|" << m_maxWidth << "|" << x << "|" << y << "|" << rad0 << "|" << sinq << "|" << cosq;
         errs << "|" << w << "|" << h << "|";  
         mxThrowException(err::sizeerr, "azBoxKernel::setKernel", "width bigger than 2*maxWidth.  This is a bug.  Details: " + errs.str());
      }

      if(h > m_maxWidth*2)
      {
         std::stringstream errs;
         errs << "|" << kernW << "|" << m_radWidth << "|" << m_azWidth << "|" << m_maxWidth << "|" << x << "|" << y << "|" << rad0 << "|" << sinq << "|" << cosq;
         errs << "|" << w << "|" << h << "|";  
         mxThrowException(err::sizeerr, "azBoxKernel::setKernel", "height bigger than 2*maxWidth.  This is a bug.  Details: " + errs.str());
      }

      kernel.resize(w, h);
      
      arithT xcen = 0.5*(w-1.0);
      arithT ycen = 0.5*(h-1.0);
      
      arithT xP, yP;
      arithT rad, radP;
      arithT sinq2, cosq2, sindq, q2, dq;
      for(int j=0; j < h; ++j)
      {
         yP = j;
         for(int i=0; i < w; ++i)
         {
            xP = i;
            rad = sqrt( pow(xP-xcen,2) + pow(yP-ycen,2) );
            radP = sqrt( pow(x+xP-xcen,2) + pow(y+yP-ycen,2) );
 
            if( fabs(radP-rad0) > m_radWidth) 
            {
               kernel(i,j) = 0;
               continue;
            }

            sinq2 = (yP-ycen)/rad;
            cosq2 = (xP-xcen)/rad;

            sindq = sinq2*cosq - cosq2*sinq;            
            if( fabs(rad*sindq) <= m_azWidth) 
            {
               if(m_maxAz > 0) //Only check this if needed.
               {
                  q2 = atan2(y+yP-ycen, x+xP-xcen);
                  dq = math::angleDiff(q,q2);
                  if(fabs(dq) > m_maxAz)
                  {
                     kernel(i,j) = 0;
                     continue;
                  }
               }
               kernel(i,j) = 1;
            }
            else kernel(i,j) = 0;
         }
      }

      arithT ksum = kernel.sum();
      if(ksum == 0)
      {
         mxThrowException(err::invalidconfig, "azBoxKernel::setKernel", "kernel sum 0 at " + std::to_string(x) + "," + std::to_string(y));
      }
      kernel /= ksum;
   }
   
};



///Filter an image with a mean kernel.
/** Applies the kernel to each pixel in the image and sums, storing the filtered result in the output image.
  * The kernel-type (kernelT) must have the following interface:
  * \code
  * template<typename _arrayT, size_t kernW>
  * struct filterKernel
  * {
  *     typedef _arrayT arrayT; // arrayT must have an eigen-like interface
  *     typedef typename _arrayT::Scalar arithT;
  *   
  *     filterKernel()
  *     {
  *        //constructor
  *     }
  *   
  *     //The maxWidth function returns the maximum possible full-width (in either direction) of the kernel
  *     //Called once
  *     int maxWidth()
  *     {
  *        //returns the maximum half-width given the configuration
  *     }
  * 
  *     //The setKernel function is called for each pixel.
  *     void setKernel(arithT x, arithT y, arrayT & kernel)
  *     {
  *        //This must resize and populate the passed in kernel array each time
  *        //so that it is re-entrant. 
  *        //Note: On output kernel array should be normalized so that sum() = 1.0
  *        //Note: the width and height of the kernel array should always be odd
  *     }
  * };
  * \endcode
  * 
  * \param [out] fim will be allocated with resize, and on output contains the filtered image
  * \param [in] im is the image to be filtered
  * \param [in] kernel a fully configured obect of type kernelT
  * \param [in] maxr is the maximum radius from the image center to apply the kernel.  pixels
  *                  outside this radius are set to 0.
  * 
  * \tparam imageOutT the type of the output image (must have an Eigen like interface)
  * \tparam imageInT the type of the input image (must have an Eigen like interface)
  * \tparam kernelT is the kernel type (see above)
  * 
  * \ingroup image_filters_kernels
  *
  */ 
template<typename imageOutT, typename imageInT, typename kernelT>
void filterImage(imageOutT & fim, imageInT im, kernelT kernel,  int maxr= 0)
{
   fim.resize(im.rows(), im.cols());
  
   float xcen = 0.5*(im.rows()-1);
   float ycen = 0.5*(im.cols()-1);
   
   if(maxr == 0) maxr = 0.5*im.rows() - kernel.maxWidth();
   
   int mini = 0.5*im.rows() - maxr;
   int maxi = 0.5*im.rows() + maxr;
   int minj = 0.5*im.cols() - maxr;
   int maxj = 0.5*im.cols() + maxr;
   

   typename kernelT::arrayT kernelArray;
   
   #pragma omp parallel private(kernelArray)
   {
      int im_i, im_j, im_p,im_q;
      int kern_i, kern_j, kern_p,kern_q;  
      typename imageOutT::Scalar norm;
   
      #pragma omp for
      for(int i=0; i< im.rows(); ++i)
      {
         for(int j=0; j<im.cols(); ++j)
         {
            if((i >= mini && i< maxi) && (j>= minj && j<maxj))
            {
               kernel.setKernel(i-xcen, j-ycen, kernelArray);
               fim(i,j) = (im.block(i-0.5*(kernelArray.rows()-1), j-0.5*(kernelArray.cols()-1), kernelArray.rows(), kernelArray.cols())*kernelArray).sum();
            }
            else
            {
               kernel.setKernel(i-xcen, j-ycen, kernelArray);
         
               im_i = i - 0.5*(kernelArray.rows()-1);
               if(im_i < 0) im_i = 0;
         
               im_j = j - 0.5*(kernelArray.cols()-1);
               if(im_j < 0) im_j = 0;

               im_p = im.rows() - im_i;
               if(im_p > kernelArray.rows()) im_p = kernelArray.rows();
          
               im_q = im.cols() - im_j;
               if(im_q > kernelArray.cols()) im_q = kernelArray.cols();
         
               kern_i = 0.5*(kernelArray.rows()-1) - i;
               if(kern_i < 0) kern_i = 0;
         
               kern_j = 0.5*(kernelArray.cols()-1) - j;
               if(kern_j < 0) kern_j = 0;
      
               kern_p = kernelArray.rows() - kern_i;
               if(kern_p > kernelArray.rows()) kern_p = kernelArray.rows();
         
               kern_q = kernelArray.cols() - kern_j;
               if(kern_q > kernelArray.cols()) kern_q = kernelArray.cols();
       
               //Pick only the smallest widths
               if(im_p < kern_p) kern_p = im_p;
               if(im_q < kern_q) kern_q = im_q;
   
               norm = kernelArray.block(kern_i, kern_j, kern_p, kern_q ).sum();
        
               fim(i,j) = ( im.block(im_i, im_j, kern_p, kern_q) * kernelArray.block(kern_i, kern_j, kern_p, kern_q )).sum()/norm;
               if( !std::isfinite(fim(i,j))) fim(i,j) = 0.0;
            }
         }
      }
   }// pragma omp parallel
}   
   
///@}

   
/// Smooth an image using the mean in a rectangular box, optionally rejecting the highest and lowest values.
/** Calculates the mean value in a rectangular box of imIn, of size meanFullSidth X meanFullWidth and stores it in the corresonding center pixel of imOut.
  * Does not smooth the 0.5*meanFullwidth rows and columns of the input image, and the values of these pixels are not
  * changed in imOut (i.e. you should 0 them before the call). 
  * 
  * imOut is not re-allocated.
  * 
  * If rejectMinMax is `true` then the minimum and maximum value in the box are not included in the mean.  Rejection makes
  * the algorithm somewhat slower, depending on box width.
  *
  * \tparam imageTout is an eigen-like image array
  * \tparam imageTin is an eigen-like image array
  * 
  * \returns 0 on success
  * \returns -1 on error.
  * 
  * \ingroup image_filters_average
  */
template<typename imageTout, typename imageTin>
int meanSmooth( imageTout & imOut,        ///< [out] the smoothed image. Not re-allocated, and the edge pixels are not modified.
                const imageTin & imIn,    ///< [in] the image to smooth
                int meanFullWidth,        ///< [in] the full-width of the smoothing box
                bool rejectMinMax = false ///< [in] whether or not to reject the min and max value.
              )
{
   typedef typename imageTout::Scalar scalarT;
      
   int buff = 0.5*meanFullWidth;
   
   int nPix = meanFullWidth*meanFullWidth;
   
   if(rejectMinMax) //avoid the branch on every pixel
   {
      nPix -= 2;
      for(int jj=buff; jj< imIn.cols()-buff; ++jj)
      {
         for(int ii=buff; ii < imIn.rows()-buff; ++ii)
         {
            scalarT sum = 0;
            scalarT max = sum;
            scalarT min = sum;
            for(int ll = jj-buff; ll < jj+buff+1; ++ll)
            {
               for(int kk = ii-buff; kk < ii+buff+1; ++kk)
               {
                  //scalarT val = imIn(kk,ll);
                  sum += imIn(kk,ll);
                  if(imIn(kk,ll) > max) max = imIn(kk,ll);
                  if(imIn(kk,ll) < min) min = imIn(kk,ll);
               }
            }
            imOut(ii,jj) = (sum-max-min)/nPix;
         }
      }
   }
   else
   {
      for(int jj=buff; jj< imIn.cols()-buff; ++jj)
      {
         for(int ii=buff; ii < imIn.rows()-buff; ++ii)
         {
            scalarT sum = 0;
            for(int ll = jj-buff; ll < jj+buff+1; ++ll)
            {
               for(int kk = ii-buff; kk < ii+buff+1; ++kk)
               {
                  sum += imIn(kk,ll);
               }
            }
            imOut(ii,jj) = sum/nPix;
         }
      }
   }
   
   return 0;
}

/// Smooth an image using the mean in a rectangular box, optionally rejecting the highest and lowest values.  Determines the location and value of the highest pixel.
/** Calculates the mean value in a rectangular box of imIn, of size meanFullSidth X meanFullWidth and stores it in the corresonding center pixel of imOut.
  * Does not smooth the 0.5*meanFullwidth rows and columns of the input image, and the values of these pixels are not
  * changed in imOut (i.e. you should 0 them before the call). 
  * 
  * imOut is not re-allocated.
  * 
  * If rejectMinMax is `true` then the minimum and maximum value in the box are not included in the mean.  Rejection makes
  * the somewhat slower, depending on box width.
  *
  * This version also determines the location and value of the maximum pixel.  This adds some overhead, maybe on the order of 10% slower than without.
  * 
  * \overload
  * 
  * \tparam imageTout is an eigen-like image array
  * \tparam imageTin is an eigen-like image array
  * 
  * \returns 0 on success
  * \returns -1 on error.
  * 
  * \ingroup image_filters_average
  */
template<typename imageTout, typename imageTin>
int meanSmooth( imageTout & imOut,                 ///< [out] the smoothed image. Not re-allocated, and the edge pixels are not modified.
                int & xMax,                        ///< [out] the x-locatioin of the max pixel
                int & yMax,                        ///< [out] the y-locatioin of the max pixel
                typename imageTout::Scalar & pMax, ///< [out] the value of the max pixel
                const imageTin & imIn,             ///< [in] the image to smooth
                int meanFullWidth,                 ///< [in] the full-width of the smoothing box
                bool rejectMinMax = false          ///< [in] whether or not to reject the min and max value.
              )
{
   typedef typename imageTout::Scalar scalarT;
      
   int buff = 0.5*meanFullWidth;
   
   int nPix = meanFullWidth*meanFullWidth;
   
   pMax = std::numeric_limits<scalarT>::lowest();
   
   if(rejectMinMax) //avoid the branch on every pixel.
   {
      nPix -= 2;
      for(int jj=buff; jj< imIn.cols()-buff; ++jj)
      {
         for(int ii=buff; ii < imIn.rows()-buff; ++ii)
         {
            scalarT sum = 0;
            scalarT max = sum;
            scalarT min = sum;
            for(int ll = jj-buff; ll < jj+buff+1; ++ll)
            {
               for(int kk = ii-buff; kk < ii+buff+1; ++kk)
               {
                  //scalarT val = imIn(kk,ll);
                  sum += imIn(kk,ll);
                  if(imIn(kk,ll) > max) max = imIn(kk,ll);
                  if(imIn(kk,ll) < min) min = imIn(kk,ll);
               }
            }
            imOut(ii,jj) = (sum-max-min)/nPix;
            if(imOut(ii,jj) > pMax)
            {
               pMax = imOut(ii,jj);
               xMax = ii;
               yMax = jj;
            }
         }
      }
   }
   else
   {
      for(int jj=buff; jj< imIn.cols()-buff; ++jj)
      {
         for(int ii=buff; ii < imIn.rows()-buff; ++ii)
         {
            scalarT sum = 0;
            for(int ll = jj-buff; ll < jj+buff+1; ++ll)
            {
               for(int kk = ii-buff; kk < ii+buff+1; ++kk)
               {
                  sum += imIn(kk,ll);
               }
            }
            imOut(ii,jj) = sum/nPix;
            if(imOut(ii,jj) > pMax)
            {
               pMax = imOut(ii,jj);
               xMax = ii;
               yMax = jj;
            }
         }
      }
   }
   
   return 0;
}      
         

/// Smooth an image using the median in a rectangular box.  Also Determines the location and value of the highest pixel in the smoothed image.
/** Calculates the median value in a rectangular box of imIn, of size medianFullSidth X medianFullWidth and stores it in the corresonding center pixel of imOut.
  * Does not smooth the 0.5*medianFullwidth rows and columns of the input image, and the values of these pixels are not
  * changed in imOut (i.e. you should 0 them before the call). 
  * 
  * imOut is not re-allocated.
  * 
  * Also determines the location and value of the maximum pixel.  This is a negligble overhead compared to the median operation.
  * 
  * 
  * \tparam imageTout is an eigen-like image array
  * \tparam imageTin is an eigen-like image array
  * 
  * \returns 0 on success
  * \returns -1 on error.
  * 
  * \ingroup  image_filters_average
  */
template<typename imageTout, typename imageTin>
int medianSmooth( imageTout & imOut,                 ///< [out] the smoothed image. Not re-allocated, and the edge pixels are not modified.
                  int & xMax,                        ///< [out] the x-locatioin of the max pixel
                  int & yMax,                        ///< [out] the y-locatioin of the max pixel
                  typename imageTout::Scalar & pMax, ///< [out] the value of the max pixel
                  const imageTin & imIn,             ///< [in] the image to smooth
                  int medianFullWidth                ///< [in] the full-width of the smoothing box
                )
{
   typedef typename imageTout::Scalar scalarT;
      
   int buff = 0.5*medianFullWidth;
   
   std::vector<scalarT> pixs(medianFullWidth*medianFullWidth);
   
   pMax = std::numeric_limits<scalarT>::lowest();
   
   for(int jj=buff; jj< imIn.cols()-buff; ++jj)
   {
      for(int ii=buff; ii < imIn.rows()-buff; ++ii)
      {
         int n=0;
         for(int ll = jj-buff; ll < jj+buff+1; ++ll)
         {
            for(int kk = ii-buff; kk < ii+buff+1; ++kk)
            {
               pixs[n] = imIn(kk,ll);
               ++n;
            }
         }
         
         imOut(ii,jj) = math::vectorMedianInPlace(pixs);
         if(imOut(ii,jj) > pMax)
         {
            pMax = imOut(ii,jj);
            xMax = ii;
            yMax = jj;
         }
      }
   }
   
   return 0;
}      

/// Smooth an image using the median in a rectangular box.  
/** Calculates the median value in a rectangular box of imIn, of size medianFullSidth X medianFullWidth and stores it in the corresponding center pixel of imOut.
  * Does not smooth the outer 0.5*medianFullwidth rows and columns of the input image, and the values of these pixels are not
  * changed in imOut (i.e. you should 0 them before the call). 
  * 
  * imOut is not re-allocated.
  * 
  * \overload 
  * 
  * \tparam imageTout is an eigen-like image array
  * \tparam imageTin is an eigen-like image array
  * 
  * \returns 0 on success
  * \returns -1 on error.
  * 
  * \ingroup image_filters_average
  */
template<typename imageTout, typename imageTin>
int medianSmooth( imageTout & imOut,                 ///< [out] the smoothed image. Not re-allocated, and the edge pixels are not modified.
                  const imageTin & imIn,             ///< [in] the image to smooth
                  int medianFullWidth                ///< [in] the full-width of the smoothing box
                )
{
   int xMax, yMax;
   typename imageTout::Scalar pMax;
   return medianSmooth(imOut, xMax, yMax, pMax, imIn, medianFullWidth);
}

template<typename eigenImT>
void rowEdgeMedSubtract( eigenImT & im, ///< The image to filter
                         int ncols      ///< The number of columns on each side of the image to use as the reference
                       )
{
   typedef typename eigenImT::Scalar realT;

   std::vector<realT> edge(2*ncols);
   for(int rr=0; rr<im.rows(); ++rr)
   {
      for(int cc=0; cc< ncols; ++cc)
      {
         edge[cc] = im(rr,cc);
      }
      
      for(int cc=0; cc< ncols; ++cc)
      {
         edge[ncols+cc] = im(rr,im.cols()-ncols + cc);
      }

      realT med = math::vectorMedian(edge);
      im.row(rr) -= med;
   }

   return;
}

template<typename eigenImT>
void colEdgeMedSubtract( eigenImT & im, ///< The image to filter
                         int nrows      ///< The number of rows on each side of the image to use as the reference
                       )
{
   typedef typename eigenImT::Scalar realT;

   std::vector<realT> edge(2*nrows);
   for(int cc=0; cc<im.cols(); ++cc)
   {
      for(int rr=0; rr< nrows; ++rr)
      {
         edge[rr] = im(rr,cc);
      }
      
      for(int rr=0; rr< nrows; ++rr)
      {
         edge[nrows+rr] = im(im.rows()-nrows + rr, cc);
      }

      realT med = math::vectorMedian(edge);
      im.col(cc) -= med;
   }

   return;
}

//------------ Radial Profile --------------------//

template<typename floatT>
struct radval
{
   floatT r;
   floatT v;
};

template<typename floatT>
struct radvalRadComp
{
   bool operator()(radval<floatT> rv1, radval<floatT> rv2)
   {
      return (rv1.r < rv2.r);
   }
};

template<typename floatT>
struct radvalValComp
{
   bool operator()(radval<floatT> rv1, radval<floatT> rv2)
   {
      return (rv1.v < rv2.v);
   }
};


/// Calculate the the radial profile
/** The median radial profile is calculated by rebinning to a 1 pixel grid.
  * 
  * 
  * \tparam vecT the std::vector-like type to contain the profile 
  * \tparam eigenImT1 the eigen-array-like type of the input image
  * \tparam eigenImT2 the eigen-array-like type of the radius and mask image
  * \tparam eigenImT3 the eigen-array-like type of the mask image
  * 
  * \ingroup rad_prof
  */ 
template<typename vecT, typename eigenImT1, typename eigenImT2, typename eigenImT3>
void radprof( vecT & rad,              ///< [out] the radius points for the profile. Should be empty.
              vecT & prof,             ///< [out] the median image value at the corresponding radius. Should be empty.
              const eigenImT1 & im,    ///< [in] the image of which to calculate the profile
              const eigenImT2 & radim, ///< [in] image of radius values per pixel
              const eigenImT3 * mask,  ///< [in] [optional] 1/0 mask, only pixels with a value of 1 are included in the profile. Set to 0 to not use.
              bool mean = false,       ///< [in] [optional] set to true to use the mean.  If false (default) the median is used.
              typename eigenImT1::Scalar minr = 0
            )
{
   typedef typename eigenImT1::Scalar floatT;
   
   int dim1 = im.rows();
   int dim2 = im.cols();
   
   floatT maxr; 
   //floatT minr;
   
   size_t nPix;
   if(mask)
   {
      nPix = mask->sum();
   }
   else
   {
      nPix = dim1*dim2;
   }
   
   /* A vector of radvals will be sorted, then binned*/
   std::vector<radval<floatT> > rv(nPix);
   
   size_t i=0;
   
   for(int c=0;c<im.cols();++c)
   {
      for(int r=0;r<im.rows();++r)
      {
         if(mask)
         {
            if( (*mask)(r,c) == 0) continue;
         }
         
         rv[i].r = radim(r,c);
         rv[i].v = im(r,c);
         ++i;
      }
   }
   
   sort(rv.begin(), rv.end(), radvalRadComp<floatT>());
   
//    for(auto it=rv.begin(); it != rv.end(); ++it)
//    {
//       std::cout << it->r << " " << it->v << "\n";
//    }
//    
//    exit(0);
   
   
   if(minr == 0) minr = rv[0].r;
   maxr = rv.back().r;
   
   /*Now bin*/
   floatT dr = 1.0;
   floatT r0 = minr;
   floatT r1 = minr + dr;
   int i1=0, i2, n;
   
   floatT med;
  
   while(r1 < maxr)
   {
      while(rv[i1].r < r0) ++i1;
      i2 = i1;
      while(rv[i2].r <= r1) ++i2;
      
      if(mean)
      {
         med = 0;
         for(int in=i1; in<i2; ++in) med += rv[in].v;
         med /= (i2-i1);
      }
      else
      {
         n = 0.5*(i2-i1);

         std::nth_element(rv.begin()+i1, rv.begin()+i1+n, rv.begin()+i2, radvalValComp<floatT>());
      
         med = (rv.begin()+i1+n)->v;
         
         //Average two points if even number of points
         if((i2-i1)%2 == 0)
         {
            med = 0.5*(med + (*std::max_element(rv.begin()+i1, rv.begin()+i1+n, radvalValComp<floatT>())).v);
         }
      }
      
      rad.push_back(.5*(r0+r1));
      prof.push_back(med);
      i1 = i2;
      r0 += dr;
      r1 += dr;
   }
   

   
}

/// Calculate the the radial profile
/** The median radial profile is calculated by rebinning to a 1 pixel grid.
  * This version calculates a centered radius image.
  * 
  * \overload
  * 
  * \tparam vecT the std::vector-like type to contain the profile 
  * \tparam eigenImT1 the eigen-array-like type of the input image
  * \tparam eigenImT2 the eigen-array-like type of the radius and mask image
  * \tparam eigenImT3 the eigen-array-like type of the mask image
  * 
  * \ingroup rad_prof
  */ 
template<typename vecT, typename eigenImT1, typename eigenImT2>
void radprof( vecT & rad,             ///< [out] the radius points for the profile. Should be empty.
              vecT & prof,            ///< [out] the median image value at the corresponding radius. Should be empty.
              const eigenImT1 & im,   ///< [in] the image of which to calculate the profile
              const eigenImT2 & mask, ///< [in] 1/0 mask, only pixels with a value of 1 are included in the profile 
              bool mean = false       ///< [in] [optional] set to true to use the mean.  If false (default) the median is used.
            )
{
   eigenImage<typename eigenImT1::Scalar> radim;
   radim.resize(im.cols(), im.rows());
   
   radiusImage(radim);
   
   radprof(rad, prof, im, radim, &mask, mean);
}

/// Calculate the the radial profile
/** The median radial profile is calculated by rebinning to a 1 pixel grid.
  * This version calculates a centered radius image.
  * 
  * \overload
  * 
  * \tparam vecT the std::vector-like type to contain the profile 
  * \tparam eigenImT1 the eigen-array-like type of the input image
  * 
  * \ingroup rad_prof
  */ 
template<typename vecT, typename eigenImT1>
void radprof( vecT & rad,           ///< [out] the radius points for the profile. Should be empty.
              vecT & prof,          ///< [out] the median image value at the corresponding radius. Should be empty.
              const eigenImT1 & im, ///< [in] the image of which to calculate the profile
              bool mean = false,     ///< [in] [optional] set to true to use the mean.  If false (default) the median is used.
              double dr = 1
            )
{
   eigenImage<typename eigenImT1::Scalar> radim;
   radim.resize(im.cols(), im.rows());
   
   radiusImage(radim);
   
   radprof(rad, prof, im, radim, (eigenImage<typename eigenImT1::Scalar> *) nullptr, mean, dr);
}

///Form a radial profile image, and optionally subtract it from the input
/** The radial profile is calculated using linear interpolation on a 1 pixel grid
  * 
  * 
  * \tparam radprofT the eigen array type of the output
  * \tparam eigenImT1 the eigen array type of the input image
  * \tparam eigenImT2 the eigen array type of the radius image
  * \tparam eigenImT3 the eigen array type of the mask image
  * 
  * \ingroup rad_prof
  */ 
template<typename radprofT, typename eigenImT1, typename eigenImT2, typename eigenImT3>
void radprofim( radprofT & radprofIm,   ///< [out] the radial profile image.  This will be resized.
                eigenImT1 & im,         ///< [in the image to form the profile of. 
                const eigenImT2 & rad,  ///< [in] an array of radius values for each pixel
                const eigenImT3 * mask, ///< [in] [optional 1/0 mask, only pixels with a value of 1 are included in the profile. Can be nullptr.
                bool subtract,          ///< [in] if true, then on ouput im will have had its radial profile subtracted.
                bool mean = false       ///< [in] [optional] set to true to use the mean.  If false (default) the median is used. 
              )
{
  
  
   std::vector<double> med_r, med_v; //Must be double for GSL interpolator
  
   radprof(med_r, med_v, im, rad, mask);
   
   /* And finally, interpolate onto the radius image */
   radprofIm.resize(im.rows(), im.cols() );
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interp(med_r, med_v);
   
   for(int c=0;c<im.cols();++c)
   {
      for(int r=0;r<im.rows();++r)
      {
         if(mask)
         {
            if( (*mask)(r,c) == 0) 
            {
               radprofIm(r,c) = 0;
               continue;
            }
         }
         
         radprofIm(r,c) = interp( ((double) rad(r,c)) );
         if(subtract) im(r,c) -= radprofIm(r,c);
      }
   }
   
}




///Form a radial profile image, and optionally subtract it from the input
/** The radial profile is calculated using linear interpolation on a 1 pixel grid.
  * This version calculates a centered radius image.
  * 
  * \tparam radprofT the eigen array type of the output
  * \tparam eigenImT the eigen array type of the input
  * 
  * \ingroup rad_prof
  */ 
template<typename radprofT, typename eigenImT>
void radprofim( radprofT & radprof,    ///< [out] the radial profile image.  This will be resized.
                eigenImT & im,         ///< [in] the image to form the profile of. 
                bool subtract = false, ///< [in] [optional] if true, then on ouput im will have had its radial profile subtracted.
                bool mean = false      ///< [in] [optional] set to true to use the mean.  If false (default) the median is used. 
              )
{
   eigenImage<typename eigenImT::Scalar> rad;
   rad.resize(im.rows(), im.cols());
   
   radiusImage(rad);
   
   radprofim(radprof, im, rad, (eigenImage<typename eigenImT::Scalar> *) nullptr, subtract);
   
}
  

/** \ingroup std_prof 
  * @{
  */

/// Form a standard deviation image, and optionally divide the input by it to form a S/N map.
/** The standard deviation profile is calculated using linear interpolation on a 1 pixel grid
  * 
  * \tparam eigenImT the eigen array type of the output and non-reference images.  Each image input can be a different type to allow references, etc.
  * 
  */ 
template<typename eigenImT, typename eigenImT1, typename eigenImT2, typename eigenImT3>
void stddevImage( eigenImT & stdIm,                 ///< [out] the standard deviation image.  This will be resized.
                  const eigenImT1 & im,             ///< [in] the image to form the standard deviation profile of, never altered.
                  const eigenImT2 & rad,            ///< [in] array of radius values
                  const eigenImT3 & mask,           ///< [in] a 1/0 mask.  0 pixels are excluded from the std-dev calculations.
                  typename eigenImT::Scalar minRad, ///< [in] the minimum radius to analyze
                  typename eigenImT::Scalar maxRad, ///< [in] the maximum radius to analyze
                  bool divide                ///< [in] if true, the output is the input image is divided by the std-dev profile, i.e. a S/N map.  default is false.
                )
{
   typedef typename eigenImT::Scalar floatT;
   
   int dim1 = im.cols();
   int dim2 = im.rows();
   
   floatT mr = rad.maxCoeff();
   
   /* A vector of radvals will be sorted, then binned*/
   std::vector<radval<floatT> > rv(dim1*dim2);
      
   for(int i=0;i<rv.size();++i)
   {
      if(mask(i) == 0) continue;
      
      rv[i].r = rad(i);
      rv[i].v = im(i);
   }
   
   sort(rv.begin(), rv.end(), radvalRadComp<floatT>());
   
   /*Now bin*/
   floatT dr = 1;
   floatT r0 = 0;
   floatT r1 = dr;
   int i1=0, i2, n;
   
   floatT stdVal;
     
   std::vector<double> std_r, std_v;
   
   while(r1 < mr)
   {
      while(rv[i1].r < r0 && i1 < rv.size())
      {
          ++i1;
      }
      if(i1 == rv.size()) 
      {
          break;
      }
      
      i2 = i1;
      while(rv[i2].r <= r1 && i2 < rv.size()) ++i2;
      n = 0.5*(i2-i1);

      std::vector<double> vals;
      
      for(int i=i1; i< i2; ++i)
      {
         vals.push_back( rv[i].v);
      } 
      
      std_r.push_back(.5*(r0+r1));
            
      std_v.push_back( std::sqrt(math::vectorVariance(vals)) ) ;
      i1 = i2;
      r0 += dr;
      r1 += dr;
   }
      
   /* And finally, interpolate onto the radius image */
   stdIm.resize(dim1, dim2);
   math::gslInterpolator<math::gsl_interp_linear<double>> interp(std_r, std_v);
   
   for(int i=0;i<dim1;++i)
   {
      for(int j=0;j<dim2;++j)
      {
         if(rad(i,j) < minRad || rad(i,j) > maxRad)
         {
            stdIm(i,j) = 0;
         }
         else
         {
            stdIm(i,j) = interp( ((double) rad(i,j)) );
            if(divide) stdIm(i,j) = im(i,j) / stdIm(i,j);
         }
      }
   }
   
}
 
/// Form a standard deviation image, and optionally divide the input by it to form a S/N map.
/** The standard deviation profile is calculated using linear interpolation on a 1 pixel grid
  * 
  * This version creates a radius map on each call, and calls the above version.  This should not 
  * be used for repeated alls, rather create a radius map ahead of time.
  * 
  * \overload 
  * 
  * \tparam eigenImT the eigen array type of the output and non-reference images
  * 
  */ 
template<typename eigenImT, typename eigenImT1, typename eigenImT2>
void stddevImage( eigenImT & stdIm,                 ///< [out] the standard deviation image.  This will be resized.
                  const eigenImT1 & im,             ///< [in] the image to form the standard deviation profile of, never altered.
                  const eigenImT2 & mask,           ///< [in] a 1/0 mask.  0 pixels are excluded from the std-dev calculations.
                  typename eigenImT::Scalar minRad, ///< [in] the minimum radius to analyze
                  typename eigenImT::Scalar maxRad, ///< [in] the maximum radius to analyze
                  bool divide = false               ///< [in] [optional] if true, the output is the input image is divided by the std-dev profile, i.e. a S/N map.  default is false.
                )
{   
   int dim1 = im.cols();
   int dim2 = im.rows();
   
   eigenImage<typename eigenImT::Scalar> rad;
   rad.resize(dim1, dim2);
   
   radiusImage(rad);
   stddevImage(stdIm, im, rad, mask, minRad, maxRad, divide );
   
}

/// Form a standard deviation image for each imamge in a cube, and optionally divide the input by it forming a S/N map cube.
/** The standard deviation profile is calculated using linear interpolation on a 1 pixel grid
  * 
  * 
  * \tparam eigencubeT is the eigen cube type of the input and output cubes.
  * \tparam eigenImT the eigen array type of the output and non-reference images.
  * 
  */ 
template<typename eigenCubeT, typename eigenCubeT1, typename eigenImT>
void stddevImageCube( eigenCubeT & stdImc,              ///< [out]  the standard deviation image cube.  This will be resized.
                      const eigenCubeT1 & imc,          ///< [in] the image cube to form the standard deviation profile of. 
                      const eigenImT & mask,            ///< [in] a 1/0 mask.  0 pixels are excluded from the std-dev calculations.
                      typename eigenImT::Scalar minRad, ///< [in] the minimum radius to analyze
                      typename eigenImT::Scalar maxRad, ///< [in] the maximum radius to analyze
                      bool divide = false               ///< [in] [optional] if true, the output is the input image is divided by the std-dev profile, i.e. a S/N map.  default is false.
                    )
{
   int dim1 = imc.cols();
   int dim2 = imc.rows();
   
   eigenImage<typename eigenCubeT::Scalar> rad;
   rad.resize(dim1, dim2);
   
   radiusImage(rad);
   
   stdImc.resize(imc.rows(), imc.cols(), imc.planes());
   
   //#pragma omp parallel for
   for(int i=0; i< imc.planes(); ++i)
   {
      eigenImage<typename eigenCubeT::Scalar> im, stdIm;
      
      im = imc.image(i);
      
      stddevImage(stdIm, im, rad, mask, minRad, maxRad, divide );

      stdImc.image(i) = stdIm;
      
   }
}

///@}

} //namespace improc 
} //namespace mx

#endif //__imageFilters_hpp__  


