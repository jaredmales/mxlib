/** \file logInterpolator.hpp
  * \brief Interpolation in log space.
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2022 Jared R. Males (jaredmales@gmail.com)
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

#include "../mxException.hpp"
#include "gslInterpolator.hpp"

namespace mx
{
namespace math 
{

/// Interpolate a function in log space
/** Given a discrete function, conduct linear interpolation in log space.
  * The input vectors are converted to their log10 values.  Linear interpolation
  * using gslInterpolator is conducted on the log10 of the input value.  The output
  * is converted back. So the result is
  * \code 
  *   y = pow(10, interp(log10(x))). 
  * \endcode
  *  
  * \ingroup interpolation
  */  
template<typename realT>
class logInterpolator
{

protected:
   gslInterpolator<gsl_interp_linear<realT>> m_interp; ///< The interpolator

   std::vector<realT>  m_logx; ///< Internal storage of the log10 values of the x values
   std::vector<realT>  m_logy; ///< Internal storage of the lgo10 values of the y values

public:
   /// Default constructor
   logInterpolator(){}

   /// Convert the inputs to their log10 values, and construct the interpolator.
   /**
     * \throws mxException if vectors are not the same size.
     */
   logInterpolator( const std::vector<realT> & x, /// [in] the input x-axis
                    const std::vector<realT> & y  /// [in] the input y-axis
                  )
   {
      setup(x,y);
   }

   /// Convert the inputs to their log10 values, and construct the interpolator.
   /**
     * \throws mx::err::sizeerr if vectors are not the same size.
     */ 
   void setup( const std::vector<realT> & x, /// [in] the input x-axis
               const std::vector<realT> & y  /// [in] the input y-axis
             )
   {
      if(x.size() != y.size())
      {
         mxThrowException(err::sizeerr, "radprofIntegral", "vectors must have same size");
      }

      m_logx.resize(x.size());
      m_logy.resize(y.size());

      for(size_t n = 0; n < x.size(); ++n)
      {
         m_logx[n] = log10(x[n]);
         m_logy[n] = log10(y[n]);
      }

      m_interp.setup(m_logx, m_logy);
   }

   /// Calculate the interpolated value at the input \par x.
   /**
     * \returns the interpolated value
     */ 
   realT operator()(const realT & x)
   {
      return pow(static_cast<realT>(10), m_interp(log10(x)));
   }

};

} //namespace math
} //namespace mx
