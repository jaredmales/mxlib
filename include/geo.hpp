/** \file geo.hpp
  * \author Jared R. Males
  * \brief Utilities for working with angles
  * \ingroup gen_math_files
  *
  */

#ifndef __geo_hpp__
#define __geo_hpp__

#include <boost/math/constants/constants.hpp>

#include <sofa.h>

namespace mx
{
   
/** \ingroup geo
  * @{
  */

///Convert from degrees to radians
/**
  * Uses boost constants:
  * \code
  * return q*boost::math::constants::degree<realT>();
  * \endcode
  * 
  * \param q is the angle to convert
  * 
  * \return the angle q converted to radians
  */
template<typename realT>
realT dtor(realT q)
{
   return q * boost::math::constants::degree<realT>();
}

///Convert from radians to degrees
/**
  * Uses boost constants:
  * \code
  * return q*boost::math::constants::radian<realT>();
  * \endcode
  * 
  * \param q is the angle to convert
  * 
  * \return the angle q converted to degrees
  */
template<typename realT>
realT rtod(realT q)
{
   return q * boost::math::constants::radian<realT>();
}

///Calculate the difference between two angles, correctly across 0/360.
/** Calculates \f$ dq = q2- q1 \f$, but accounts for crossing 0/360.  This implies
  * that \f$ dq \le 180 \f$.
  * 
  * \param q1 angle to subtract from q2, in degrees.
  * \param q2 angle to subtract q1 from, in degrees.
  *
  * \returns the difference of q2 and q1
  *
  * \tparam realT is the type in which to do arithmetic
  * 
  */
template<int degrad = 0, typename realT>
realT angleDiff(realT q1, realT q2)
{ 
   static_assert(std::is_floating_point<realT>::value);
   
   realT dq = q2-q1;

   realT half;
   
   if(degrad) half = boost::math::constants::pi<realT>();
   else half = static_cast<realT>(180);
   
   if (abs(dq) > half)
   {
      if(dq < 0)
      {
         dq = dq + static_cast<realT>(2)*half;
      }
      else
      {
         dq = dq - static_cast<realT>(2)*half;
      }
   }
   
   return dq;
}


/// @}

} //namespace mx

#endif //__geo_hpp__
