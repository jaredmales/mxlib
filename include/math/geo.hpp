/** \file geo.hpp
  * \author Jared R. Males
  * \brief Utilities for working with angles
  * \ingroup gen_math_files
  *
  */

#ifndef __geo_hpp__
#define __geo_hpp__

#include <vector>


#include <boost/math/constants/constants.hpp>



namespace mx
{
   
/** \ingroup geo
  * @{
  */

///Calculate the semi-latus rectum of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>
#define semilatrect(a,e) (e == 0.0 ? a : (e == 1.0 ? 2.*a : (e < 1. ? a*(1-e*e) : a*(e*e-1))))

///Calculate the focal parameter of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>
#define focus(a,e) (e == 0.0 ? 1e34 : (e == 1.0 ? 2.*a : (e < 1. ? a*(1-e*e)/e : a*(e*e-1)/e)))

///Calculate the semi-major axis of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>, given the focal parameter and the eccentricity
#define semimaj(p,e) (e == 1.0 ? 1e34 : (e < 1 ? p*e/(1-e*e) : p*e/(e*e-1) ))

///Calculate the eccentricity of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a> given the semi-major axis and the focal parameter
#define eccent(a, p) (a == 0.0 ? 1e34 : (p >= 1e9 ? 0.0 : (p>0 ? (-p/(2*a)+0.5*std::sqrt(p*p/(a*a) + 4)) : (p/(2*a)+0.5*std::sqrt(p*p/(a*a) + 4)) ) ))


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
   static_assert(std::is_floating_point<realT>::value, "realT must be floating point");
   
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


///Calculate the mean of a set of angles, correctly across 0/360.
/** Calculates the mean by decomposing into vectors and averaging the components. This accounts for crossing 0/360.  
  * 
  * \param q vector of angles to average.
  *
  * \returns the mean angle
  *
  * \tparam degrad controls whether angles are degrees (false) or radians (true)
  * \tparam realT is the type in which to do arithmetic
  * 
  */
template<typename realT, int degrad = 0>
realT angleMean(std::vector<realT> & q)
{ 
   static_assert(std::is_floating_point<realT>::value, "realT must be floating point");
   
   realT s = 0;
   realT c = 0;

   realT d2r;
   
   if( degrad ) d2r = 1;
   else d2r  = boost::math::constants::pi<realT>()/static_cast<realT>(180);

   for(int i=0; i< q.size(); ++i)
   {
      s += sin( q[i]*d2r );
      c += cos( q[i]*d2r );
   }
   
   s /= q.size();
   c /= q.size();
   
   return atan(s/c)/d2r ;

}


/// @}

} //namespace mx

#endif //__geo_hpp__
