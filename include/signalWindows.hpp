/** \file signalWindows.hpp
  * \brief Procedures to calculate window functions for signal processing
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing
  *
  */

#ifndef __signalWindows_hpp__
#define __signalWindows_hpp__

#include <boost/math/constants/constants.hpp>

namespace mx
{
   
/** \brief Create a 2-D Tukey window
  * 
  * Function to create a 2-D Tukey window.  
  * 
  * The width of the window is controlled by alpha.  alpha = 0 is a cylinder, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  * \param filt is a pre-allocated array of size dim x dim
  * \param dim is the size of the array
  * \param N is the diameter of the window, nominally N=dim, but this is not required
  * \param alpha controls the window width.
  * \param xc is the desired x center of the window
  * \param yc is the desired y center of the window
  */
template<typename floatT>
void tukey2d(floatT *filt, int dim, floatT N, floatT alpha, floatT xc, floatT yc)
{
   
   int ii, jj;
   floatT x, y, r;

   floatT rad = 0.5*(N-1.0);

   floatT pi = boost::math::constants::pi<floatT>();
   
   for(ii=0; ii<dim; ++ii)
   {
      x = ( (floatT) ii) - xc;
      for(jj=0; jj<dim; ++jj)
      {
         y = ( (floatT) jj) - yc;

         r = sqrt(x*x + y*y);

         //Following mxlib convention of including half pixels
         if(r > rad + 0.5)
         {
            filt[ii*dim + jj] = 0.0;
         }
         else if(rad + r > (N-1)*(1-0.5*alpha) && alpha > 0.)
         {
            //Have to prevent going greater than N-1 due to half pixel inclusion.
            floatT dr = rad+r;
            if(dr > N-1) dr = N-1;
            
            filt[ii*dim + jj] = 0.5*(1.0 + cos(pi * ( 2.*(dr)/(alpha*(N-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[ii*dim + jj] = 1.0;
         }
      }
   }
}

/** \brief Create a 2-D Tukey window on an annulus
  * 
  * Function to create a 2-D Tukey window on an annulus.  
  * 
  * The width of the window is controlled by alpha.  alpha = 0 is a cylinder, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  * \param filt is a pre-allocated array of size dim x dim.
  * \param dim is the size of the array.
  * \param N is the outer diameter of the window, nominally N=dim, but this is not required.
  * \param eps is the ratio of inner diameter to the outer diameter
  * \param alpha controls the window width.
  * \param xc is the desired x center of the window.
  * \param yc is the desired y center of the window.
  *
  * \tparam floatT is a floating point type 
  */
template<typename floatT>
void tukey2dAnnulus(floatT *filt, int dim, floatT N, floatT eps, floatT alpha, floatT xc, floatT yc)
{

   int ii, jj;
   floatT x, y, r, z, Z;

   floatT rad = 0.5*(N-1.0);
   
   Z = floor((1.0-eps)*(rad)) + 1.0;
   
   floatT pi = boost::math::constants::pi<floatT>();
      
   for(ii=0; ii<dim; ++ii)
   {
      x = ( (floatT) ii) - xc;
      for(jj=0; jj<dim; ++jj)
      {
         y = ( (floatT) jj) - yc;

         r = sqrt(x*x + y*y);

         z = (r - eps*(rad));
         
         //Following mxlib convention of including half pixels
         if(r > rad + 0.5 || r < eps*rad)
         {
            filt[ii*dim + jj] = 0.0;
         }
         else if(z <= 0.5*alpha*(Z-1) && alpha > 0)
         {
            filt[ii*dim + jj] = 0.5*(1.0 + cos(pi*(2.*z/(alpha*(Z-1)) -1.0) ));
         }
         else if(z >= (Z-1)*(1.-0.5*alpha) && alpha > 0)
         {
            z = z*((Z-0.5)/Z); //Stretch a little to help with the half pixel
            if(z > Z-1) z = Z-1;
            filt[ii*dim + jj] = 0.5*(1.0 + cos(pi* ( 2.*(z)/(alpha*(Z-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[ii*dim + jj] = 1.0;
         }
      }
   }
}

} //namespace mx
#endif // __signalWindows_hpp__

