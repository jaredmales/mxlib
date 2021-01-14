/** \file array2FitGaussian2D.hpp
  * \author Jared R. Males
  * \brief Wrapper for a native array to pass to \ref levmarInterface, with 2D Gaussian details.
  * \ingroup fitting_files
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_fit_array2FitGaussian2D_hpp
#define math_fit_array2FitGaussian2D_hpp

#include "../../mxError.hpp"

namespace mx
{
namespace math
{
namespace fit
{


/// Wrapper for a native array to pass to \ref levmarInterface, with 2D Gaussian details.
/** Supports fixing G0, G, x0, and y0 independently.  The shape and orientation can be fixed, but
  * for the general form, sigma_x, sigma_y, and theta can only be fixed together. 
  * \ingroup gaussian_peak_fit
  */
template<typename realT>
struct array2FitGaussian2D
{
   realT * data {nullptr}; ///< Pointer to the array
   size_t nx {0};          ///< X dimension of the array
   size_t ny {0};          ///< Y dimension of the array
   
   realT * mask {nullptr}; ///< Pointer to the (optional) mask array.  Any 0 pixels are excluded from the fit.
   
   realT m_G0 {0};
   realT m_G {0};
   realT m_x0 {0};
   realT m_y0 {0};
   realT m_sigma_x {0};
   realT m_sigma_y {0};
   realT m_theta {0};
   
   realT m_a {0};
   realT m_b {0};
   realT m_c {0};
   
   realT m_sigma {0};
   
   int m_G0_idx {0};
   int m_G_idx {1};
   int m_x0_idx {2};
   int m_y0_idx {3};
   int m_sigma_x_idx {4}; ///< Index of sigma_x in the parameters.  Re-used for a.
   int m_sigma_y_idx {5}; ///< Index of sigma_y in the parameters.  Re-used for b.
   int m_theta_idx {6};  ///< Index of theta in the parameters.  Re-used for c.
   
   int m_sigma_idx {4};  ///< Index of sigma in the symmetric case
   
   int m_nparams {7};
   int m_maxNparams {7};
   
   void setSymmetric();
   
   void setGeneral();
   
   /// Set whether each parameter is fixed.
   /** Sets the parameter indices appropriately.
     */
   void setFixed( bool G0, 
                  bool G, 
                  bool x0, 
                  bool y0, 
                  bool sigma_x, 
                  bool sigma_y,
                  bool theta
                );
   
   realT G0( realT * p );
   
   void G0( realT * p,
            realT nG0
          );
   
   realT G( realT * p );
   
   void G( realT * p,
           realT nG
         );
   
   realT x0( realT * p );
   
   void x0( realT * p,
            realT nx0
          );
   
   realT y0( realT * p );
   
   void y0( realT * p,
            realT ny0
          );
   
   realT sigma_x( realT * p );
   
   void sigma_x( realT * p,
                 realT nsigma_x
               );
   
   realT sigma_y( realT * p );
   
   void sigma_y( realT * p,
                 realT nsigma_y
               );
   
   realT theta( realT * p );
   
   void theta( realT * p,
               realT ntheta
             );
   
   realT sigma( realT * p );
   
   void sigma( realT * p,
               realT nsigma
             );
   
   realT a( realT * p );
   
   void a( realT * p,
           realT na
         );
   
   realT b( realT * p );
   
   void b( realT * p,
           realT nb
         );
   
   realT c( realT * p );
   
   void c( realT * p,
           realT nc
         );
   
   int nparams();
   
};

template<typename realT>
void array2FitGaussian2D<realT>::setSymmetric()
{
   m_maxNparams = 5;
}

template<typename realT>
void array2FitGaussian2D<realT>::setGeneral()
{
   m_maxNparams = 7;
}

template<typename realT>
void array2FitGaussian2D<realT>::setFixed( bool G0, 
                                           bool G, 
                                           bool x0, 
                                           bool y0, 
                                           bool sigma_x, 
                                           bool sigma_y,
                                           bool theta
                                         )
{
   int idx = 0;
   
   if(G0) m_G0_idx = -1;
   else m_G0_idx = idx++;
   
   if(G) m_G_idx = -1;
   else m_G_idx = idx++;

   if(x0) m_x0_idx = -1;
   else m_x0_idx = idx++;
   
   if(y0) m_y0_idx = -1;
   else m_y0_idx = idx++;
   
   if(m_maxNparams == 5)
   {
      if(sigma_x) m_sigma_idx = -1;
      else m_sigma_idx = idx++;
   }
   else if(sigma_x && sigma_y && theta)
   {
      m_sigma_x_idx = -1;
      m_sigma_y_idx = -1;
      m_theta_idx = -1;
   }
   else
   {
      if(sigma_x || sigma_y || theta)
      {
         mxError("array2FitGaussian2D::setFixed", MXE_NOTIMPL, "cannot fix sigma_x, sigma_y, and theta separately");
      }
      
      m_sigma_x_idx = idx++;
      m_sigma_y_idx = idx++;
      m_theta_idx = idx++;
   }
   
   m_nparams = idx;
}

template<typename realT>
realT array2FitGaussian2D<realT>::G0( realT * p )
{
   if( m_G0_idx < 0 )
   {
      return m_G0;
   }
   else
   {
      return p[m_G0_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::G0( realT * p,
                                     realT nG0
                                   )
{
   if( m_G0_idx < 0 )
   {
      m_G0 = nG0;
   }
   else
   {
      p[m_G0_idx]=nG0;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::G( realT * p )
{
   if( m_G_idx < 0 )
   {
      return m_G;
   }
   else
   {
      return p[m_G_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::G( realT * p,
        realT nG
      )
{
   if( m_G_idx < 0 )
   {
      m_G = nG;
   }
   else
   {
      p[m_G_idx] = nG;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::x0( realT * p )
{
   if( m_x0_idx < 0 )
   {
      return m_x0;
   }
   else
   {
      return p[m_x0_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::x0( realT * p,
         realT nx0
       )
{
   if( m_x0_idx < 0 )
   {
      m_x0 = nx0;
   }
   else
   {
      p[m_x0_idx] = nx0;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::y0( realT * p )
{
   if( m_y0_idx < 0 )
   {
      return m_y0;
   }
   else
   {
      return p[m_y0_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::y0( realT * p,
                                     realT ny0
                                   )
{
   if( m_y0_idx < 0 )
   {
      m_y0 = ny0;
   }
   else
   {
      p[m_y0_idx] = ny0;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::sigma_x( realT * p )
{
   if( m_sigma_x_idx < 0 )
   {
      return m_sigma_x;
   }
   else
   {
      return p[m_sigma_x_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::sigma_x( realT * p,
                                          realT nsigma_x
                                        )
{
   if( m_sigma_x_idx < 0 )
   {
      m_sigma_x = nsigma_x;
   }
   else
   {
      p[m_sigma_x_idx] = nsigma_x;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::sigma_y( realT * p )
{
   if( m_sigma_y_idx < 0 )
   {
      return m_sigma_y;
   }
   else
   {
      return p[m_sigma_y_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::sigma_y( realT * p,
                                          realT nsigma_y
                                        )
{
   if( m_sigma_y_idx < 0 )
   {
      m_sigma_y = nsigma_y;
   }
   else
   {
      p[m_sigma_y_idx] = nsigma_y;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::theta( realT * p )
{
   if( m_theta_idx < 0 )
   {
      return m_theta;
   }
   else
   {
      return p[m_theta_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::theta( realT * p,
                                        realT ntheta
                                      )
{
   if( m_theta_idx < 0 )
   {
      m_theta = ntheta;
   }
   else
   {
      p[m_theta_idx] = ntheta;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::sigma( realT * p )
{
   if( m_sigma_idx < 0 )
   {
      return m_sigma;
   }
   else
   {
      return p[m_sigma_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::sigma( realT * p,
                                        realT nsigma
                                      )
{
   if( m_sigma_idx < 0 )
   {
      m_sigma = nsigma;
   }
   else
   {
      p[m_sigma_idx] = nsigma;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::a( realT * p )
{
   //This aliases sigma_x after param normalization
   if( m_sigma_x_idx < 0 )
   {
      return m_a;
   }
   else
   {
      return p[m_sigma_x_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::a( realT * p,
        realT na
      )
{
   //This aliases sigma_x after param normalization
   if( m_sigma_x_idx < 0 )
   {
      m_a = na;
   }
   else
   {
      p[m_sigma_x_idx] = na;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::b( realT * p )
{
   //This aliases sigma_y after param normalization
   if( m_sigma_y_idx < 0 )
   {
      return m_b;
   }
   else
   {
      return p[m_sigma_y_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::b( realT * p,
                                    realT nb
                                  )
{
   //This aliases sigma_y after param normalization
   if( m_sigma_y_idx < 0 )
   {
      m_b = nb;
   }
   else
   {
      p[m_sigma_y_idx] = nb;
   }
}

template<typename realT>
realT array2FitGaussian2D<realT>::c( realT * p )
{
   //This aliases theta after param normalization
   if( m_theta_idx < 0 )
   {
      return m_c;
   }
   else
   {
      return p[m_theta_idx];
   }
}

template<typename realT>
void array2FitGaussian2D<realT>::c( realT * p,
                                    realT nc
                                  )
{
   //This aliases theta after param normalization
   if( m_theta_idx < 0 )
   {
      m_c = nc;
   }
   else
   {
      p[m_theta_idx] = nc;
   }
}

template<typename realT>
int array2FitGaussian2D<realT>::nparams()
{
   return m_nparams;
}

} //namespace fit
} //namespace math

} //namespace mx

#endif //math_fit_array2FitGaussian2D_hpp

