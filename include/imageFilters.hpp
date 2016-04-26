/** \file imageFilters.hpp
  * \brief Image filters (smoothing, radial profiles, etc.)
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __imageFilters_hpp__
#define __imageFilters_hpp__

#include <cstdlib>

#include "gslInterpolation.hpp"

namespace mx
{

///Symetric Gaussian smoothing kernel
/** \ingroup image_processing
  */
template<typename _arrayT, size_t _kernW=4>
struct gaussKernel
{
   typedef _arrayT arrayT;
   typedef typename _arrayT::Scalar arithT;
   static const int kernW = _kernW;
   
   arrayT kernel;
   
   arithT _fwhm;
   
   gaussKernel(arithT fwhm)
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
            kernel(i,j) = exp(-r2/sig2);
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
      kernelArray = kernel;
   }
   
};

///Azimuthally variable boxcare smoothing kernel.
/** \ingroup image_processing
  */
template<typename _arrayT, size_t _kernW=4>
struct azBoxKernel
{
   typedef _arrayT arrayT;
   typedef typename _arrayT::Scalar arithT;


   static const int kernW = _kernW;
   
   //arrayT kernel;
   
   arithT _radWidth;
   arithT _azWidth;
   int _maxWidth;
   
   azBoxKernel(arithT radWidth, arithT azWidth)
   {
      _radWidth = radWidth;
      _azWidth = azWidth;
      
      _maxWidth = std::max(radWidth, azWidth);
   }
   
   int maxWidth()
   {
      return _maxWidth;
   }
   
   void setKernel(arithT x, arithT y, arrayT & kernel)
   {
      arithT rad = sqrt((arithT) (x*x + y*y));
      
      arithT sinq = y/rad;
      arithT cosq = x/rad;
      
      
      int w = 2*( (int) std::max( fabs(_azWidth*sinq), fabs(_radWidth*cosq) ) + 1 );
      int h = 2*( (int) std::max( fabs(_azWidth*cosq), fabs(_radWidth*sinq) ) + 1);
      
      
      kernel.resize(w, h);
      
      arithT xcen = 0.5*(w-1.0);
      arithT ycen = 0.5*(h-1.0);
      
      arithT sinq2, cosq2, sindq;
      for(int i=0; i < w; ++i)
      {
         for(int j=0; j < h; ++j)
         {
            rad = sqrt( pow(i-xcen,2) + pow(j-ycen,2) );
            sinq2 = (j-ycen)/rad;
            cosq2 = (i-xcen)/rad;
            sindq = sinq2*cosq - cosq2*sinq;
            if( rad <= _radWidth && fabs(rad*sindq) <= _azWidth) kernel(i,j) = 1;
            else kernel(i,j) = 0;
         }
      }
      if(kernel.sum() == 0)
      {
         std::cerr << "Kernel sum 0: " << x << " " << y << "\n";
         exit(-1);
      }
      kernel /= kernel.sum();
   }
   
};
            
///Filter an image with a kernel.
/** Applies the kernel to each pixel in the image, storing the filtered result in the output image.
  * The kernel-type (kernelT) must have the following interface:
  * \code
  * template<typename _arrayT, size_t kernW=4>
  * struct filterKernel
  * {
  *     typedef _arrayT arrayT;
  *     typedef typename _arrayT::Scalar arithT;
  *   
  *     filterKernel()
  *     {
  *        //constructor
  *     }
  *   
  *     //The setKernel function is called for each pixel.
  *     void setKernel(arithT x, arithT y, arrayT & kernel)
  *     {
  *        //This may or may not do anything (e.g., kernel could be created in constructor)
  *        //but on output kernel array should be normalized so that sum() = 1.0
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
  * \ingroup image_processing 
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
   
   for(int i=mini; i<maxi; ++i)
   {
      for(int j=minj; j<maxj; ++j)
      {
         kernel.setKernel(i-xcen, j-ycen, kernelArray);
         fim(i,j) = (im.block(i-0.5*kernelArray.rows(), j-0.5*kernelArray.cols(), kernelArray.rows(), kernelArray.cols())*kernelArray).sum();
      }
   }

   //Now handle the edges
   #pragma omp parallel private(kernelArray)
   {
      int im_i, im_j, im_p,im_q;
      int kern_i, kern_j, kern_p,kern_q;  
      typename imageOutT::Scalar norm;
   
      #pragma omp for
      for(size_t i=0; i< im.rows(); i++)
      {
         for(size_t j=0; j<im.cols(); j++)
         {
            //if((i >= maxr && i< im.rows()-maxr) && (j>= maxr && j<im.rows()-maxr)) continue;
            if((i >= mini && i< maxi) && (j>= minj && j<maxj)) continue;
         
            kernel.setKernel(i-xcen, j-ycen, kernelArray);
         
            im_i = i - 0.5*kernelArray.rows();
            if(im_i < 0) im_i = 0;
         
            im_j = j - 0.5*kernelArray.cols();
            if(im_j < 0) im_j = 0;

            im_p = im.rows() - im_i;
            if(im_p > kernelArray.rows()) im_p = kernelArray.rows();
          
            im_q = im.cols() - im_j;
            if(im_q > kernelArray.cols()) im_q = kernelArray.cols();
         
            kern_i = 0.5*kernelArray.rows() - i;
            if(kern_i < 0) kern_i = 0;
         
            kern_j = 0.5*kernelArray.cols() - j;
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
   }// pragma omp parallel
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



///Form a radial profile image, and optionally subtract it from the input
/** The radial profile is calculated using linear interpolation on a 1 pixel grid
  * \ingroup image_processing
  * 
  * \param [out] radprof is the radial profile image.  This will be resized.
  * \param [in] im is the image to form the profile of. 
  * \param [in] rad is an array of radius values
  * \param [in] subtract if true, then on ouput im will have had its radial profile subtracted.
  * 
  * \tparam radprofT the eigen array type of the output
  * \tparam eigenimT the eigen array type of the input
  */ 
template<typename radprofT, typename eigenimT>
void radprofim( radprofT & radprof, 
                eigenimT & im,
                radprofT & rad,
                bool subtract = false)
{
   typedef typename eigenimT::Scalar floatT;
   
   int dim1 = im.cols();
   int dim2 = im.rows();
   
   floatT mr = rad.maxCoeff();
   
   /* A vector of radvals will be sorted, then binned*/
   std::vector<radval<floatT> > rv(dim1*dim2);
   
   for(int i=0;i<rv.size();++i)
   {
      rv[i].r = rad(i);
      rv[i].v = im(i);
   }
   
   sort(rv.begin(), rv.end(), radvalRadComp<floatT>());
   
   /*Now bin*/
   floatT dr = 1;
   floatT r0 = 0;
   floatT r1 = dr;
   int i1=0, i2, n;
   
   floatT med;
  
   std::vector<double> med_r, med_v;
   while(r1 < mr)
   {
      while(rv[i1].r < r0) ++i1;
      i2 = i1;
      while(rv[i2].r <= r1) ++i2;
      
      n = 0.5*(i2-i1);

      std::nth_element(rv.begin()+i1, rv.begin()+i1+n, rv.begin()+i2, radvalValComp<floatT>());
      
      med = (rv.begin()+i1+n)->v;
      
      med_r.push_back(.5*(r0+r1));
      med_v.push_back(med);
      i1 = i2;
      r0 += dr;
      r1 += dr;
   }
   
   /* And finally, interpolate onto the radius image */
   radprof.resize(dim1, dim2);
   gslInterpolator interp(gsl_interp_linear, med_r, med_v);
   
   for(int i=0;i<dim1;++i)
   {
      for(int j=0;j<dim2;++j)
      {
         radprof(i,j) = interp.interpolate( ((double) rad(i,j)) );
         if(subtract) im(i,j) -= radprof(i,j);
      }
   }
   
}




///Form a radial profile image, and optionally subtract it from the input
/** The radial profile is calculated using linear interpolation on a 1 pixel grid.
  * This version calculates a centered radius image.
  * 
  * \ingroup image_processing
  * 
  * \param [out] radprof is the radial profile image.  This will be resized.
  * \param [in] im is the image to form the profile of. 
  * \param [in] subtract if true, then on ouput im will have had its radial profile subtracted.
  * 
  * \tparam radprofT the eigen array type of the output
  * \tparam eigenimT the eigen array type of the input
  */ 
template<typename radprofT, typename eigenimT>
void radprofim( radprofT & radprof, 
                eigenimT & im, 
                bool subtract = false)
{
   typedef typename eigenimT::Scalar floatT;
   
   int dim1 = im.cols();
   int dim2 = im.rows();
   
   radprofT rad;
   rad.resize(dim1, dim2);
   
   radiusImage(rad);
   
   radprofim(radprof, im, rad, subtract);
   
}
  
///Form a standard deviation image, and optionally divide the input by it
/** The standard deviation profile using linear interpolation on a 1 pixel grid
  * \ingroup image_processing
  * 
  * \param [out] stdIm is the standard deviation image.  This will be resized.
  * \param [in] im is the image to form the standard deviation profile of. 
  * \param [in] rad is an array of radius values
  * \param [in] mask is an array which is a 1/0 mask.  0 pixels are excluded from the std-dev calculations.
  * \param [in] minRad is the minimum radius to analyze
  * \param [in] maxRad is the maximum radius to analyze
  * \param [in] divide if true, the input image is divided by the std-def profile
  * 
  * \tparam eigenimT the eigen array type of the output and non-reference images
  */ 
template<typename eigenimT>
void stddevImage( eigenimT & stdIm, 
                  eigenimT & im,
                  eigenimT & rad,
                  eigenimT & mask,
                  typename eigenimT::Scalar minRad,
                  typename eigenimT::Scalar maxRad, 
                  bool divide = false )
{
   typedef typename eigenimT::Scalar floatT;
   
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
      while(rv[i1].r < r0) ++i1;
      i2 = i1;
      while(rv[i2].r <= r1) ++i2;
      
      n = 0.5*(i2-i1);

      std::vector<double> vals;
      
      for(int i=i1; i< i2; ++i)
      {
         vals.push_back( rv[i].v);
      } 
      
      std_r.push_back(.5*(r0+r1));
      
      std_v.push_back( std::sqrt(mx::vectorVariance(vals)) ) ;
      i1 = i2;
      r0 += dr;
      r1 += dr;
   }
   
   /* And finally, interpolate onto the radius image */
   stdIm.resize(dim1, dim2);
   mx::gslInterpolator interp(gsl_interp_linear, std_r, std_v);
   
   for(int i=0;i<dim1;++i)
   {
      for(int j=0;j<dim2;++j)
      {
         if(rad(i,j) < minRad || rad(i,j) > maxRad)
         {
            stdIm(i,j) = 0;
            if(divide) im(i,j) = 0;
         }
         else
         {
            stdIm(i,j) = interp.interpolate( ((double) rad(i,j)) );
            if(divide) im(i,j) /= stdIm(i,j);
         }
      }
   }
   
}
 
///Form a standard deviation image, and optionally divide the input by it
/** The standard deviation profile using linear interpolation on a 1 pixel grid
  * \ingroup image_processing
  * 
  * \param [out] stdIm is the standard deviation image.  This will be resized.
  * \param [in] im is the image to form the standard deviation profile of. 
  * \param [in] mask is an array which is a 1/0 mask.  0 pixels are excluded from the std-dev calculations.
  * \param [in] minRad is the minimum radius to analyze
  * \param [in] maxRad is the maximum radius to analyze
  * \param [in] divide if true, the input image is divided by the std-def profile
  * 
  * \tparam eigenimT the eigen array type of the output and non-reference images
  */ 
template<typename eigenimT>
void stddevImage( eigenimT & stdIm, 
                  eigenimT & im,
                  eigenimT & mask,
                  typename eigenimT::Scalar minRad,
                  typename eigenimT::Scalar maxRad,
                  bool divide = false )
{
   typedef typename eigenimT::Scalar floatT;
   
   int dim1 = im.cols();
   int dim2 = im.rows();
   
   eigenimT rad;
   rad.resize(dim1, dim2);
   
   mx::radiusImage(rad);
   
   stddevImage(stdIm, im, rad, mask, minRad, maxRad, divide );
   
}

///Form a standard deviation image for each imamge in a cube, and optionally divide the input by it
/** The standard deviation profile using linear interpolation on a 1 pixel grid
  * \ingroup image_processing
  * 
  * \param [out] stdImc is the standard deviation image cube.  This will be resized.
  * \param [in] imc is the image cube to form the standard deviation profile of. 
  * \param [in] mask is an array which is a 1/0 mask.  0 pixels are excluded from the std-dev calculations.
  * \param [in] minRad is the minimum radius to analyze
  * \param [in] maxRad is the maximum radius to analyze
  * \param [in] divide if true, the input image is divided by the std-def profile
  * 
  * \tparam eigencubeT is the eigen cube type of the input and output cubes
  * \tparam eigenimT the eigen array type of the output and non-reference images
  */ 
template<typename eigencubeT, typename eigenimT>
void stddevImageCube( eigencubeT & stdImc, 
                      eigencubeT & imc,
                      eigenimT & mask,
                      typename eigenimT::Scalar minRad,
                      typename eigenimT::Scalar maxRad,
                      bool divide = false )
{
   typedef typename eigenimT::Scalar floatT;
   
   int dim1 = imc.cols();
   int dim2 = imc.rows();
   
   eigenimT rad;
   rad.resize(dim1, dim2);
   
   mx::radiusImage(rad);

   eigenimT im, stdIm;
   
   stdImc.resize(imc.rows(), imc.cols(), imc.planes());
   
   for(int i=0; i< imc.planes(); ++i)
   {
      im = imc.image(i);
      
      stddevImage(stdIm, im, rad, mask, minRad, maxRad, divide );

      stdImc.image(i) = stdIm;
      
      if(divide) imc.image(i) = im;
      
   }
}

} //namespace mx

#endif //__imageFilters_hpp__  


