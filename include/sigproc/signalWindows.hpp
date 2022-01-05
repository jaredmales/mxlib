/** \file signalWindows.hpp
  * \brief Procedures to calculate window functions for signal processing
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2015 - 2020 Jared R. Males (jaredmales@gmail.com)
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

#include <cmath>
#include "../math/constants.hpp"

namespace mx
{
namespace sigproc 
{
namespace window
{
   
/// The Tukey Window
/** 
  * The width of the window is controlled by alpha.  alpha = 0 is a square wave, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  * 
  * \tparam realT a floating point type
  * 
  * \ingroup signal_windows1D
  */
template<typename realT>
void tukey( realT *filt, ///< [out] the pre-allocated array to hold the filter
            int N,       ///< [in] the size of the filter
            realT alpha  ///< [in] the width parameter
          )
{
   constexpr realT pi = math::pi<realT>();
   
   realT lim1 = alpha*(N-1.0)/2.0;
   realT lim2 = (N-1.0)*(1.0-0.5*alpha);
   
   for(int ii=0; ii<N; ++ii)
   {
    
      if(ii < lim1 && alpha > 0.)
      {
         filt[ii] = 0.5*(1.0 + cos(pi * ( 2.*(ii)/(alpha*(N-1)) - 1.0) ));
      }
      else if(ii > lim2 && alpha > 0.)
      {
         filt[ii] = 0.5*(1.0 + cos(pi * ( 2.*(ii)/(alpha*(N-1)) - 2./alpha + 1.0) ));
      }
      else
      {
         filt[ii] = 1.0;
      }
   }
}
 
/// The Tukey Window
/** 
  * The width of the window is controlled by alpha.  alpha = 0 is a square wave, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  * 
  * \overload
  * 
  * \tparam realT a floating point type
  * 
  * \ingroup signal_windows1D
  */
template<typename realT>
void tukey( std::vector<realT> & filt, ///< [out] the pre-allocated vector to hold the filter
            realT alpha                ///< [in] the width parameter
          )
{
   tukey(filt.data(), filt.size(), alpha);
}

/// The generalized 2 parameter cosine Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void genCosine( realT * filt, ///< [out] The pre-allocated vector which will store the filter
                size_t N,     ///< [in] the size of the filter vector
                realT a0,     ///< [in] parameter of the generalized cosine window
                realT a1      ///< [in] parameter of the generalized cosine window
              )
{
   constexpr realT pi = math::pi<realT>();

   for( size_t n=0; n<N; ++n)
   {
      filt[n] = a0 - a1*cos(2*pi*n/N);
   }
}

/// The generalized 3 parameter cosine Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void genCosine( realT * filt, ///< [out] The pre-allocated vector which will store the filter
                size_t N,     ///< [in] the size of the filter vector
                realT a0,     ///< [in] parameter of the generalized cosine window
                realT a1,     ///< [in] parameter of the generalized cosine window
                realT a2      ///< [in] parameter of the generalized cosine window
              )
{
   constexpr realT pi = math::pi<realT>();

   for( size_t n=0; n<N; ++n)
   {
      filt[n] = a0 - a1*cos(2*pi*n/N) + a2*cos(4*pi*n/N);
   }
}

/// The generalized 4 parameter cosine Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void genCosine( realT * filt, ///< [out] The pre-allocated vector which will store the filter
                size_t N,     ///< [in] the size of the filter vector
                realT a0,     ///< [in] parameter of the generalized cosine window
                realT a1,     ///< [in] parameter of the generalized cosine window
                realT a2,     ///< [in] parameter of the generalized cosine window
                realT a3      ///< [in] parameter of the generalized cosine window
              )
{
   constexpr realT pi = math::pi<realT>();

   for( size_t n=0; n<N; ++n)
   {
      filt[n] = a0 - a1*cos(2*pi*n/N) + a2*cos(4*pi*n/N) - a3*cos(6*pi*n/N);
   }
}

/// The Hann Window
/** 
  * See https://en.wikipedia.org/wiki/Window_function
  * 
  * \tparam realT a floating point type
  * 
  * \ingroup signal_windows1D
  */
template<typename realT>
void hann( realT *filt, ///< [out] the pre-allocated array to hold the filter
           int N        ///< [in] the size of the filter
         )
{
   genCosine<realT>(filt, N, 0.5, 0.5);
}

/// The Hann Window
/** 
  * See https://en.wikipedia.org/wiki/Window_function
  * 
  * \overload
  * 
  * \tparam realT a floating point type
  * 
  * \ingroup signal_windows1D
  */
template<typename realT>
void hann( std::vector<realT> & filt /**< [out] the pre-allocated vector to hold the filter */)
{
   hann(filt.data(), filt.size());
}

/// The Hamming Window
/** 
  * See https://en.wikipedia.org/wiki/Window_function
  * 
  * \tparam realT a floating point type
  * 
  * \ingroup signal_windows1D
  */
template<typename realT>
void hamming( realT *filt, ///< [out] the pre-allocated array to hold the filter
              int N        ///< [in] the size of the filter
            )
{
   genCosine<realT>(filt, N, 25.0/46.0, 1.0-25.0/46.0);
}

/// The Hamming Window
/** 
  * See https://en.wikipedia.org/wiki/Window_function
  * 
  * \overload
  * 
  * \tparam realT a floating point type
  * 
  * \ingroup signal_windows1D
  */
template<typename realT>
void hamming( std::vector<realT> & filt /**< [out] the pre-allocated vector to hold the filter */)
{
   hamming(filt.data(), filt.size());
}

/// The Blackman Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void blackman( realT * filt, ///< [out] The pre-allocated vector which will store the filter
               size_t N      ///< [in] the size of the filter vector
             )
{
   genCosine<realT>(filt, N, 0.42, 0.5, 0.08);
}

/// The Blackman Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \overload
  * 
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void blackman( std::vector<realT> & filt /**< [out] The pre-allocated vector which will store the filter */)
{
   blackman(filt.data(), filt.size());
}

/// The Exact Blackman Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void exactBlackman( realT * filt, ///< [out] The pre-allocated vector which will store the filter
                    size_t N      ///< [in] the size of the filter vector
                  )
{
   genCosine<realT>(filt, N,  0.42659, 0.49656, 0.076849);
}

/// The Exact Blackman Windwo
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \overload
  * 
  * \tparam realT a real floating point type
  */ 
template<typename realT>
void exactBlackman( std::vector<realT> & filt /**< [out] The pre-allocated vector which will store the filter */)
{
   exactBlackman(filt.data(), filt.size());
}

/// The Nuttal Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void nuttal( realT * filt, ///< [out] The pre-allocated vector which will store the filter
             size_t N      ///< [in] the size of the filter vector
           )
{
   genCosine<realT>(filt, N, 0.355768, 0.487396, 0.144232, 0.012604);
}

/// The Nuttal Window
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \overload
  * 
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void nuttal( std::vector<realT> & filt /**< [out] The pre-allocated vector which will store the filter */)
{
   nuttal(filt.data(), filt.size());
}

/// The Blackman-Nuttal Windwo
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void blackmanNuttal( realT * filt, ///< [out] The pre-allocated vector which will store the filter
                     size_t N      ///< [in] the size of the filter vector
                   )
{
   genCosine<realT>(filt, N, 0.3635819, 0.4891775, 0.1365995, 0.0106411 );
}

/// The Blackman-Nuttal Windwo
/** See https://en.wikipedia.org/wiki/Window_function
  *
  * \overload
  * 
  * \tparam realT a real floating point type
  * 
  * \ingroup signal_windows1D
  */ 
template<typename realT>
void blackmanNuttal( std::vector<realT> & filt /**< [out] The pre-allocated vector which will store the filter */)
{
   blackmanNuttal(filt.data(), filt.size());
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
  *  \ingroup signal_windows2D
  */
template<typename realT>
void tukey2d(realT *filt, int dim, realT N, realT alpha, realT xc, realT yc)
{
   
   int ii, jj;
   realT x, y, r;

   realT rad = 0.5*(N-1.0);

   realT pi = math::pi<realT>();
   
   for(ii=0; ii<dim; ++ii)
   {
      x = ( (realT) ii) - xc;
      for(jj=0; jj<dim; ++jj)
      {
         y = ( (realT) jj) - yc;

         r = sqrt(x*x + y*y);

         //Following mxlib convention of including half pixels
         if(r > rad + 0.5)
         {
            filt[jj*dim + ii] = 0.0;
         }
         else if(rad + r > (N-1)*(1-0.5*alpha) && alpha > 0.)
         {
            //Have to prevent going greater than N-1 due to half pixel inclusion.
            realT dr = rad+r;
            if(dr > N-1) dr = N-1;
            
            filt[jj*dim + ii] = 0.5*(1.0 + cos(pi * ( 2.*(dr)/(alpha*(N-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[jj*dim + ii] = 1.0;
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
  * \tparam realT is a floating point type 
  * 
  *  \ingroup signal_windows2D
  */
template<typename realT>
void tukey2dAnnulus(realT *filt, int dim, realT N, realT eps, realT alpha, realT xc, realT yc)
{

   int ii, jj;
   realT x, y, r, z, Z;

   realT rad = 0.5*(N-1.0);
   
   Z = (1-eps)*rad+1.0; //floor((1.0-eps)*(rad)) + 1.0;
   
   realT pi = math::pi<realT>();
      
   for(ii=0; ii<dim; ++ii)
   {
      x = ( (realT) ii) - xc;
      for(jj=0; jj<dim; ++jj)
      {
         y = ( (realT) jj) - yc;

         r = sqrt(x*x + y*y);

         z = (r - eps*(rad));
         
         //Following mxlib convention of including half pixels
         if(r > rad + 0.5 || r < eps*rad)
         {
            filt[jj*dim + ii] = 0.0;
         }
         else if(z <= 0.5*alpha*(Z-1) && alpha > 0)
         {
            filt[jj*dim + ii] = 0.5*(1.0 + cos(pi*(2.*z/(alpha*(Z-1)) -1.0) ));
         }
         else if(z > (Z-1)*(1.-0.5*alpha) && alpha > 0)
         {
            z = z*((Z-0.5)/Z); //Stretch a little to help with the half pixel
            if(z > Z) z = Z-1;
            filt[jj*dim + ii] = 0.5*(1.0 + cos(pi* ( 2.*(z)/(alpha*(Z-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[jj*dim + ii] = 1.0;
         }
      }
   }
}

} //namespace window
} //namespace sigproc 
} //namespace mx

#endif // signalWindows_hpp

