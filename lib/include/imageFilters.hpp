
#ifndef __imageFilters_hpp__
#define __imageFilters_hpp__

#include <cstdlib>
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



template<typename arrT> void smoothImage(arrT &smim, arrT & im, arrT & kernel)
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
   
   typename arrT::Scalar norm;
   
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
  
} //namespace mx

#endif //__imageFilters_hpp__  


