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

template<size_t kernW=4, typename arithT=double, typename arrT, typename fwhmT> 
void gaussKernel(arrT &kernel, fwhmT fwhm)
{
   size_t kernSize = fwhm*kernW; 
   
   if(kernSize % 2 == 0) kernSize++;
   
   arithT kernCen = 0.5*(kernSize-1.);

   kernel.resize(kernSize, kernSize);

   fwhmT r2;

   fwhmT sig2 = fwhm/2.354820045030949327;
   sig2 *= sig2;

   for(size_t i =0; i < kernSize; i++)
   {
      for(size_t j=0; j < kernSize; j++)
      {
         r2 = (i-kernCen)*(i-kernCen) + (j-kernCen)*(j-kernCen);

         kernel(i,j) = exp(-r2/sig2);
      }
   }
   
   kernel /= kernel.sum();
}



template<typename arrT1, typename arrT2, typename arrT3> 
void smoothImage(arrT1 &smim, arrT2 & im, arrT3 & kernel)
{
   
   smim.resize(im.rows(), im.cols());
   smim.setZero();
   
   size_t kernW = kernel.rows();
   size_t fullW = 0.5*kernW;
   
   //First do the main part of the image as fast as possible
   #pragma omp parallel for
   for(int i=fullW; i< im.rows()-fullW; i++)
   {
      for(int j=fullW; j<im.cols()-fullW; j++)
      {
         smim(i,j) = (im.block(i-fullW, j-fullW, kernW, kernW)*kernel).sum();
      }
   }

   //Now handle the edges
   int im_i, im_j, im_p,im_q;
   int kern_i, kern_j, kern_p,kern_q;
   
   typename arrT1::Scalar norm;
   
   for(size_t i=0; i< im.rows(); i++)
   {
      for(size_t j=0; j<im.cols(); j++)
      {
         if((i >= fullW && i< im.rows()-fullW) && (j>= fullW && j<im.rows()-fullW)) continue;
         
         im_i = i - fullW;
         if(im_i < 0) im_i = 0;
         
         im_j = j - fullW;
         if(im_j < 0) im_j = 0;

         im_p = im.rows() - im_i;
         if(im_p > kernW) im_p = kernW;
          
         im_q = im.cols() - im_j;
         if(im_q > kernW) im_q = kernW;
         
         kern_i = fullW - i;
         if(kern_i < 0) kern_i = 0;
         
         kern_j = fullW - j;
         if(kern_j < 0) kern_j = 0;
      
         kern_p = kernW - kern_i;
         if(kern_p > kernW) kern_p = kernW;
         
         kern_q = kernW - kern_j;
         if(kern_q > kernW) kern_q = kernW;
       
         //Pick only the smallest widths
         if(im_p < kern_p) kern_p = im_p;
         if(im_q < kern_q) kern_q = im_q;
   
         norm = kernel.block(kern_i, kern_j, kern_p, kern_q ).sum();
         smim(i,j) = ( im.block(im_i, im_j, kern_p, kern_q) * kernel.block(kern_i, kern_j, kern_p, kern_q )).sum()/norm;
      }
   }
   
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
  
  
} //namespace mx

#endif //__imageFilters_hpp__  


