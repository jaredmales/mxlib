/** \file gslInterpolation.hpp
  * \brief Wrappers for using the GNU Scientific Library 1-D interpolation functions
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup interpolation
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

#ifndef gslInterpolation_hpp
#define gslInterpolation_hpp

#include <vector>
#include <gsl/gsl_interp.h>

namespace mx
{
namespace math
{
   
template<typename _realT>
struct gslInterpolator
{
   typedef _realT realT;
   
   static_assert( std::is_same<double, typename std::remove_cv<realT>::type>::value, "GSL Interpolation only works with double");
   
   gsl_interp * interp {0};
   gsl_interp_accel * acc {0};
   
   realT *_xin;
   realT *_yin;
   
   void setup(const gsl_interp_type * interpT, realT *xin, realT *yin, size_t Nin)
   {
      if(interp) gsl_interp_free(interp);
      if(acc) gsl_interp_accel_free(acc);
      
      interp =  gsl_interp_alloc(interpT, Nin);
      acc = gsl_interp_accel_alloc ();

      gsl_interp_init(interp, xin, yin, Nin);
   
      gsl_interp_accel_reset(acc);
      
      _xin = xin;
      _yin = yin;
   }
   
   void setup(const gsl_interp_type * interpT, std::vector<realT> & xin, std::vector<realT> & yin)
   {
      setup(interpT, xin.data(), yin.data(), xin.size());
   }
   
   gslInterpolator()
   {
      
   }
   
   gslInterpolator(const gsl_interp_type * interpT, realT *xin, realT *yin, size_t Nin)
   {
      setup(interpT, xin,yin,Nin);
   }
   
   gslInterpolator(const gsl_interp_type * interpT, std::vector<realT> & xin, std::vector<realT> & yin)
   {
      setup(interpT, xin.data(), yin.data(), xin.size());
   }
   
   ~gslInterpolator()
   {
      if(interp) gsl_interp_free(interp);
      if(acc) gsl_interp_accel_free (acc);
   }
   
   realT interpolate(const realT & x)
   {
      realT y;
      gsl_interp_eval_e (interp, _xin, _yin, x, acc, &y);
      return y;
   }
   
   realT operator()(const realT & x)
   {
      return interpolate(x);
   }
   
};

///Interpolate a 1-D data X vs Y discrete function onto a new X axis
/**
  * \param interpT one of the <a href="https://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html#Interpolation-Types">gsl interpolation types</a>.
  * \param [in]  xin the input x-axis
  * \param [in]  yin the input y-values
  * \param [in]  Nin the size of the input x-y data series
  * \param [in]  xout the desired x-axis
  * \param [out] yout the output interpolated y-values, pre-allocated
  * \param [in]  Nout the size of the output x-y axis
  * 
  * \retval
  * 
  * \ingroup interpolation
  */
template<typename realT>
int gsl_interpolate( const gsl_interp_type * interpT, 
                     realT *xin, 
                     realT *yin, 
                     size_t Nin, 
                     realT *xout, 
                     realT *yout, 
                     size_t Nout )
{
   static_assert( std::is_same<double, typename std::remove_cv<realT>::type>::value, "GSL Interpolation only works with double");
   
   gsl_interp * interp =  gsl_interp_alloc(interpT, Nin);
   gsl_interp_accel * acc = gsl_interp_accel_alloc ();

   gsl_interp_init(interp, xin, yin, Nin);
   
   gsl_interp_accel_reset(acc);
   
   for(size_t i=0;i<Nout; ++i)
   {
      gsl_interp_eval_e (interp, xin, yin, xout[i], acc, &yout[i]);
   }
   
   gsl_interp_free(interp);
   gsl_interp_accel_free (acc);
   
   return 0;
}

///Interpolate a 1-D data X vs Y discrete function onto a new X axis (vector version)
/**
  * \param interpT of the <a href="https://www.gnu.org/software/gsl/doc/html/interp.html#d-interpolation-types">gsl interpolation types</a>.
  * \param [in] xin the input x-axis
  * \param [in] yin the input y-values
  * \param [in] xout the desired x-axis
  * \param [out] yout the output interpolated y-values, does not need to be allocated
  * 
  * \retval
  * 
  * \ingroup interpolation
  */
template<typename realT>
int gsl_interpolate( const gsl_interp_type * interpT,
                     std::vector<realT> & xin,
                     std::vector<realT> & yin,
                     std::vector<realT> & xout,
                     std::vector<realT> & yout )
{
   yout.resize(xout.size());
   return gsl_interpolate(interpT, xin.data(), yin.data(), xin.size(), xout.data(), yout.data(), xout.size());
}

} //namespace math
} //namespace mx

#endif //gslInterpolation_hpp
