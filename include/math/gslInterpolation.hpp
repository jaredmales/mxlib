/** \file gslInterpolation.hpp
  * \brief Wrappers for using the GNU Scientific Library 1-D interpolation functions
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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
#include <gsl/gsl_errno.h>

#include "../mxException.hpp"

namespace mx
{
namespace math
{
   
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
  * \returns 0 on success
  * \returns -1 on a gsl error.
  * 
  * \todo report errors iaw mxlib standard in gsl_interpolate
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
   
   gsl_set_error_handler_off();

   gsl_interp * interp =  gsl_interp_alloc(interpT, Nin);
   if(interp == nullptr)
   {
      std::cerr << "gsl_interpolate: gsl_interp_alloc failed\n";
      return -1;
   }

   gsl_interp_accel * acc = gsl_interp_accel_alloc ();
   if(acc == nullptr)
   {
      std::cerr << "gsl_interpolate: gsl_interp_accel_alloc failed\n";
      return -1;

   }
   int gsl_errno = gsl_interp_init(interp, xin, yin, Nin);

   if(gsl_errno != 0)
   {
      std::cerr << "gsl_interpolate: error from gsl_interp_init [" << gsl_strerror(gsl_errno) << "]\n";
      gsl_interp_free(interp);
      gsl_interp_accel_free (acc);
      return -1;
   }
   
   gsl_errno = 0;
   for(size_t i=0;i<Nout; ++i)
   {
      //Don't error check, let it set NAN
      gsl_errno += gsl_interp_eval_e (interp, xin, yin, xout[i], acc, &yout[i]);
   }
   
   gsl_interp_free(interp);
   gsl_interp_accel_free (acc);
   
   if(gsl_errno)
   {
      std::cerr << "gsl_interpolate: error(s) reported by gsl_interp_eval_e\n";
      return -1;
   }

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
