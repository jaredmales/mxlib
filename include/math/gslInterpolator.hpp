/** \file gslInterpolator.hpp
  * \brief Class for managing 1-D interpolation using the GNU Scientific Library
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

#ifndef gslInterpolator_hpp
#define gslInterpolator_hpp

#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "../mxException.hpp"

namespace mx
{
namespace math
{
   
/// GSL Linear Interpolation
/** \ingroup interpolation
  */ 
template<typename _realT>
struct gsl_interp_linear
{
   typedef _realT realT;

   static const gsl_interp_type * interpolator()
   {
      return ::gsl_interp_linear;
   }
};

#ifndef MX_OLD_GSL
/// GSL Steffen Interpolation
/** \ingroup interpolation
  */
template<typename _realT>
struct gsl_interp_steffen
{
   typedef _realT realT;

   static const gsl_interp_type * interpolator()
   {
      return ::gsl_interp_steffen;
   }
};
#endif

/// Class to manage interpolation using the GSL interpolation library.
/**
  *
  * \tparam interpT is the interpolation type, which also specifies the floating point precision.
  * 
  * \ingroup interpolation
  */ 
template<typename interpT>
class gslInterpolator
{
public:
   typedef typename interpT::realT realT;
   
   static_assert( std::is_same<double, typename std::remove_cv<realT>::type>::value, "GSL Interpolation only works with double");

protected:

   gsl_interp * m_interp {nullptr}; ///< the gsl interpolator structure
   gsl_interp_accel * m_acc {nullptr}; ///< the gsl interpolation accelerator structure
   
   realT *m_xin; ///< the input x data
   realT *m_yin; ///< the input y data

public:

   /// Default constructor
   gslInterpolator();

   /// Raw pointer constructor   
   gslInterpolator( realT *xin, ///< [in] the input x data, this pointer is stored and must remain valid
                    realT *yin, ///< [in] the input y data, this pointer is stored and must remain valid
                    size_t Nin  ///< [in] the size of data vectors
                  );
   
   /// Vector constructor
   gslInterpolator( std::vector<realT> & xin, ///< [in] the input x data, the pointer to xin.data() is stored and xin must remain unchanged
                    std::vector<realT> & yin  ///< [in] the input y data, the pointer to yin.data() is stored and yin must remain unchanged
                  );
   
   /// Destructor
   /** Deallocates working memory.
     */  
   ~gslInterpolator();
   
   /// Setup the interpolator for the supplied data pointers
   /**
     * \throws mx::err::allocerr if allocation fails
     * \throws mx::err::liberr if GSL initialization fails 
     */ 
   void setup( realT *xin, ///< [in] the input x data, this pointer is stored and must remain valid
               realT *yin, ///< [in] the input y data, this pointer is stored and must remain valid
               size_t Nin  ///< [in] the size of data vectors
             );
   
   /// Setup the interpolator for the supplied data vectors
   /**
     * \throws mx::err::sizeerr if the vectors are not the same size
     * \throws mx::err::allocerr if allocation fails
     * \throws mx::err::liberr if GSL initialization fails 
     */
   void setup( std::vector<realT> & xin, ///< [in] the input x data, the pointer to xin.data() is stored and xin must remain unchanged
               std::vector<realT> & yin  ///< [in] the input y data, the pointer to yin.data() is stored and yin must remain unchanged
             );   

   /// Calculate the interpolated function value at a point
   /**
     * \returns the interpolated value at \p x 
     */ 
   realT operator()( const realT & x /**< [in] the point at which to interpolate */ );
   
};

template<typename interpT>
gslInterpolator<interpT>::gslInterpolator() 
{   
}

template<typename interpT>
gslInterpolator<interpT>::gslInterpolator( realT *xin, 
                                           realT *yin, 
                                           size_t Nin  
                                         )
{
   setup(xin,yin,Nin);
}

template<typename interpT>
gslInterpolator<interpT>::gslInterpolator( std::vector<realT> & xin, 
                                           std::vector<realT> & yin  
                                        )
{
   setup(xin.data(), yin.data(), xin.size());
}

template<typename interpT>
gslInterpolator<interpT>::~gslInterpolator()
{
   if(m_interp != nullptr) gsl_interp_free(m_interp);
   if(m_acc != nullptr) gsl_interp_accel_free (m_acc);
}

template<typename interpT>
void gslInterpolator<interpT>::setup( realT *xin, 
                                      realT *yin, 
                                      size_t Nin  
                                    )
{
   if(m_interp != nullptr) gsl_interp_free(m_interp);
   if(m_acc != nullptr) gsl_interp_accel_free(m_acc);
   
   m_interp =  gsl_interp_alloc(interpT::interpolator(), Nin);
   if(!m_interp) mxThrowException(err::allocerr, "gslInterpolation::setup", "gsl_interp_alloc failed");
   m_acc = gsl_interp_accel_alloc ();
   if(!m_acc) mxThrowException(err::allocerr, "gslInterpolation::setup", "gsl_interp_accel_alloc failed");
   int errv = gsl_interp_init(m_interp, xin, yin, Nin);

   if(errv != 0)
   {
      mxThrowException(err::liberr, "gslInterpolation::setup", "gsl_interp_init failed");
   }
   errv = gsl_interp_accel_reset(m_acc);
   
   if(errv != 0)
   {
      mxThrowException(err::liberr, "gslInterpolation::setup", "gsl_interp_init failed");
   }
   m_xin = xin;
   m_yin = yin;
}

template<typename interpT>
void gslInterpolator<interpT>::setup( std::vector<realT> & xin, 
                                      std::vector<realT> & yin  
                                    )
{
   if(xin.size() != yin.size())
   {
      mxThrowException(err::sizeerr, "gslInterpolator<interpT>::setup", "input vectors must be the same size");
   }
   setup(xin.data(), yin.data(), xin.size());
}

template<typename interpT>
typename interpT::realT gslInterpolator<interpT>::operator()(const realT & x)
{
   realT y;
   int errv = gsl_interp_eval_e (m_interp, m_xin, m_yin, x, m_acc, &y);
   if(errv != 0 && errv != GSL_EDOM)
   {
      mxThrowException(err::liberr, "gslInterpolation::setup", "gsl_interp_eval_e failed");
   }
   return y;
}
   


} //namespace math
} //namespace mx

#endif //gslInterpolator_hpp
