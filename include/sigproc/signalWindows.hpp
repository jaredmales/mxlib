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
  * The shape of the window is controlled by alpha.  alpha = 0 is a cylinder, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  *  \ingroup signal_windows2D
  */
template<typename realT>
void tukey2d( realT *filt, ///< [out] a pre-allocated array of size \p rows x \p cols (column major)
              int rows,    ///< [in] the number of rows in filt
              int cols,    ///< [in] the number of cols in filt
              realT D,     ///< [in] the diameter of the window
              realT alpha, ///< [in] controls the window shape.  1.0 gives a Hann window, 0.0 gives a cylinder (a.k.a. no window)
              realT xc,    ///< [in] the desired x center of the window.
              realT yc     ///< [in] the desired y center of the window.
            )
{
   
   int ii, jj;

   realT rad = 0.5*(D-1.0);

   realT pi = math::pi<realT>();
   
   for(int cc=0; cc<cols; ++cc)
   {
      realT y = ( (realT) cc) - yc;
      for(int rr=0; rr<rows; ++rr)
      {
         realT x = ( (realT) rr) - xc;

         realT r = sqrt(x*x + y*y);

         //Following mxlib convention of including half pixels
         if(r > rad + 0.5)
         {
            filt[cc*rows + rr] = 0.0;
         }
         else if(rad + r > (D-1)*(1-0.5*alpha) && alpha > 0.)
         {
            //Have to prevent going greater than N-1 due to half pixel inclusion.
            realT dr = rad+r;
            if(dr > D-1) dr = D-1;
            
            filt[cc*rows + rr] = 0.5*(1.0 + cos(pi * ( 2.*(dr)/(alpha*(D-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[cc*rows + rr] = 1.0;
         }
      }
   }
}

/** \brief Create a 2-D Tukey window
  * 
  * Function to create a 2-D Tukey window.  
  * 
  * The shape of the window is controlled by alpha.  alpha = 0 is a cylinder, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam arrT an Eigen-like array type
  * 
  * \overload
  * 
  * \ingroup signal_windows2D
  */
template<typename arrT>
void tukey2d( arrT & filt,                 ///< [in.out] a pre-allocated array
              typename arrT::Scalar D,     ///< [in] the diameter of the window
              typename arrT::Scalar alpha, ///< [in] controls the window shape.  1.0 gives a Hann window, 0.0 gives a cylinder (a.k.a. no window)
              typename arrT::Scalar xc,    ///< [in] the desired x center of the window.
              typename arrT::Scalar yc     ///< [in] the desired y center of the window.
            )
{
   tukey2d(filt.data(), filt.rows(), filt.cols(), D, alpha, xc, yc);
}

/** \brief Create a 2-D Tukey window on an annulus
  * 
  * Function to create a 2-D Tukey window on an annulus.  
  * 
  * The shape of the window is controlled by alpha.  alpha = 0 is a cylinder, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  * \tparam realT is a floating point type 
  * 
  *  \ingroup signal_windows2D
  */
template<typename realT>
void tukey2dAnnulus( realT *filt, ///< [out] a pre-allocated array of size \p rows x \p cols (column major)
                     int rows,    ///< [in] the number of rows in filt
                     int cols,    ///< [in] the number of cols in filt
                     realT D,     ///< [in] the outer diameter of the window
                     realT eps,   ///< [in] the ratio of inner diameter to the outer diameter
                     realT alpha, ///< [in] controls the window shape.  1.0 gives a Hann window, 0.0 gives a cylinder (a.k.a. no window)
                     realT xc,    ///< [in] the desired x center of the window.
                     realT yc     ///< [in] the desired y center of the window.
                   )
{
   realT rad = 0.5*(D-1.0);
   
   int Z = (1-eps)*rad+1.0; //floor((1.0-eps)*(rad)) + 1.0;
   
   realT pi = math::pi<realT>();
      
   for(int cc=0; cc<cols; ++cc)
   {
      realT y = ( (realT) cc) - yc;
      for(int rr=0; rr<rows; ++rr)
      {
         realT x = ( (realT) rr) - xc;

         realT r = sqrt(x*x + y*y);

         realT z = (r - eps*(rad));
         
         //Following mxlib convention of including half pixels
         if(r > rad + 0.5 || r < eps*rad)
         {
            filt[cc*rows + rr] = 0.0;
         }
         else if(z <= 0.5*alpha*(Z-1) && alpha > 0)
         {
            filt[cc*rows + rr] = 0.5*(1.0 + cos(pi*(2.*z/(alpha*(Z-1)) -1.0) ));
         }
         else if(z > (Z-1)*(1.-0.5*alpha) && alpha > 0)
         {
            z = z*((Z-0.5)/Z); //Stretch a little to help with the half pixel
            if(z > Z) z = Z-1;
            filt[cc*rows + rr] = 0.5*(1.0 + cos(pi* ( 2.*(z)/(alpha*(Z-1)) - 2./alpha + 1.0) ));
         }
         else
         {
            filt[cc*rows + rr] = 1.0;
         }
      }
   }
}

/** \brief Create a 2-D Tukey window on an annulus
  * 
  * Function to create a 2-D Tukey window on an annulus.  
  * 
  * The shape of the window is controlled by alpha.  alpha = 0 is a cylinder, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  * \overload
  * 
  * \tparam realT is a floating point type 
  * 
  * \ingroup signal_windows2D
  */
template<typename arrT>
void tukey2dAnnulus( arrT & filt,                 ///< [in,out] a pre-allocated array
                     typename arrT::Scalar D,     ///< [in] the outer diameter of the window
                     typename arrT::Scalar eps,   ///< [in] the ratio of inner diameter to the outer diameter
                     typename arrT::Scalar alpha, ///< [in] controls the window shape.  1.0 gives a Hann window, 0.0 gives a cylinder (a.k.a. no window)
                     typename arrT::Scalar xc,    ///< [in] the desired x center of the window.
                     typename arrT::Scalar yc     ///< [in] the desired y center of the window.
                   )
{
   return tukey2dAnnulus(filt, filt.rows(), filt.cols(), D, eps, alpha, xc, yc);
}

/** \brief Create a 2-D Tukey window on a rectangle
  * 
  * Function to create a 2-D Tukey window on a rectangle.  
  * 
  * The shape of the window is controlled by alpha.  alpha = 0 is a rectangle window, alpha=1.0 is the Hann window.
  * 
  * See https://en.wikipedia.org/wiki/Window_function
  *
  *  \ingroup signal_windows2D
  */
template<typename realT>
void tukey2dSquare( realT *filt, ///< [out] a pre-allocated array of size \p rows x \p cols (column major)
                    int rows,    ///< [in] the number of rows in filt
                    int cols,    ///< [in] the number of cols in filt
                    realT W,     ///< [in] the width of the window, corresponds to rows
                    realT H,     ///< [in] the height of the window, corresponds to cols
                    realT alpha, ///< [in] controls the window shape.  1.0 gives a Hann window, 0.0 gives a cylinder (a.k.a. no window)
                    realT xc,    ///< [in] the desired x center of the window.
                    realT yc     ///< [in] the desired y center of the window.
                  )
{
   
    int ii, jj;

    realT W2 = 0.5*(W-1.0);
    realT H2 = 0.5*(H-1.0);

    realT pi = math::pi<realT>();
   
    for(int cc=0; cc<cols; ++cc)
    {
        realT y = fabs(( (realT) cc) - yc);

        if( y > H2 + 0.5)
        {
            for(int rr=0; rr<rows; ++rr)
            {
                filt[cc*rows + rr] = 0;
            }
            continue;
        }

        for(int rr=0; rr<rows; ++rr)
        {
            realT x = fabs(( (realT) rr) - xc);

            if( x > W2 + 0.5)
            {
                filt[cc*rows + rr] = 0;
                continue;
            }
            else if( (W2 + x > (W-1)*(1-0.5*alpha)) && (H2 + x > (H-1)*(1-0.5*alpha)) && alpha > 0.)
            {
                //Have to prevent going greater than N-1 due to half pixel inclusion.
                realT dx = W2+x;
                if(dx > W-1) dx = W-1;
            
                realT dy = H2+y;
                if(dy > H-1) dy = H-1;

                filt[cc*rows + rr] = 0.5*(1.0 + cos(pi * ( 2.*(dx)/(alpha*(W-1)) - 2./alpha + 1.0) ));
                filt[cc*rows + rr] *= 0.5*(1.0 + cos(pi * ( 2.*(dy)/(alpha*(H-1)) - 2./alpha + 1.0) ));
         }
         else
         {
                filt[cc*rows + rr] = 1.0;
         }
      }
   }
}

} //namespace window
} //namespace sigproc 
} //namespace mx

#endif // signalWindows_hpp

