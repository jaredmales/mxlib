/** \file signalWindows.hpp
  * \brief Procedures to calculate window functions for signal processing
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef signalWindows_hpp
#define signalWindows_hpp

#include <boost/math/constants/constants.hpp>

namespace mx
{
namespace sigproc 
{
   

/// Create a 1-D Tukey window
/** 
  * Function to create a 1-D Tukey window.  
  * 
  * The width of the window is controlled by alpha.  alpha = 0 is a square wave, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  * \param filt is a pre-allocated array of size N
  * \param N is the size of the array
  * \param alpha controls the window width.
  * 
  * \tparam floatT specifies the floating point type
  * 
  * \ingroup signal_processing
  */
template<typename floatT, typename alphaT>
void tukey1d(floatT *filt, int N, alphaT alpha)
{
   floatT pi = boost::math::constants::pi<floatT>();
   floatT falpha = static_cast<floatT>(alpha); 
   
   floatT lim1 = falpha*(N-1.0)/2.0;
   floatT lim2 = (N-1.0)*(1.0-0.5*falpha);
   
   for(int ii=0; ii<N; ++ii)
   {
    
      if(ii < lim1 && falpha > 0.)
      {
         filt[ii] = 0.5*(1.0 + cos(pi * ( 2.*(ii)/(falpha*(N-1)) - 1.0) ));
      }
      else if(ii > lim2 && falpha > 0.)
      {
         filt[ii] = 0.5*(1.0 + cos(pi * ( 2.*(ii)/(falpha*(N-1)) - 2./falpha + 1.0) ));
      }
      else
      {
         filt[ii] = 1.0;
      }
   }
}
   
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
  * 
  *  \ingroup signal_processing
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
  * 
  *  \ingroup signal_processing
  */
template<typename floatT>
void tukey2dAnnulus(floatT *filt, int dim, floatT N, floatT eps, floatT alpha, floatT xc, floatT yc)
{

   int ii, jj;
   floatT x, y, r, z, Z;

   floatT rad = 0.5*(N-1.0);
   
   Z = (1-eps)*rad+1.0; //floor((1.0-eps)*(rad)) + 1.0;
   
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
         else if(z > (Z-1)*(1.-0.5*alpha) && alpha > 0)
         {
            z = z*((Z-0.5)/Z); //Stretch a little to help with the half pixel
            if(z > Z) z = Z-1;
            filt[ii*dim + jj] = 0.5*(1.0 + cos(pi* ( 2.*(z)/(alpha*(Z-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[ii*dim + jj] = 1.0;
         }
      }
   }
}


} //namespace sigproc 
} //namespace mx

#endif // signalWindows_hpp

