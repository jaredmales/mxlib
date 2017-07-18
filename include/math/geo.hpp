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
namespace math
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

///Calculate the angle modulo full-circle, normalizing to a positive value.
/** The output will be betweeen 0 and 360 (or 0 and 2pi).
  * 
  * \returns the value of q normalized to 0 <= q < 360[2pi]
  *
  * \tparam degrad controls whether this is in degrees (0, default) or radians (1)
  * \tparam realT is the type in which to do arithmetic
  * 
  */
template<int degrad = 0, typename realT>
realT angleMod(realT q /**< [in] the angle */)
{ 
   static_assert(std::is_floating_point<realT>::value, "angleMod: realT must be floating point");
   
   realT full;
   
   if(degrad) full = boost::math::constants::two_pi<realT>();
   else full = static_cast<realT>(360);
   
   q = fmod(q,full);
   
   if(q < 0) q += full;
   
   return q;
}


///Calculate the difference between two angles, correctly across 0/360.
/** Calculates \f$ dq = q2- q1 \f$, but accounts for crossing 0/360.  This implies
  * that \f$ dq \le 180 \f$.
  * 
  *
  * \returns the difference of q2 and q1
  * 
  * \tparam degrad controls whether this is in degrees (0, default) or radians (1)
  * \tparam realT is the type in which to do arithmetic
  * 
  */
template<int degrad = 0, typename realT>
realT angleDiff( realT q1, ///< [in] angle to subtract from q2, in degrees.
                 realT q2 ///< [in] angle to subtract q1 from, in degrees.
               )
{ 
   static_assert(std::is_floating_point<realT>::value, "angleDiff: realT must be floating point");
   
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
template<int degrad = 0, typename realT>
realT angleMean(std::vector<realT> & q)
{ 
   static_assert(std::is_floating_point<realT>::value, "angleMean: realT must be floating point");
   
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

///Make a vector of angles continuous, fixing the 0/360 crossing.
/** The vector is modified so it is continuous.
  *
  * 
  * \tparam degrad controls whether angles are degrees (false) or radians (true)
  * \tparam realT is the type in which to do arithmetic
  */
template<int degrad = 0, typename realT>
int continueAngles( std::vector<realT> & angles, ///< [in] the vector of angles
                    realT threshold=0.75 ///< [in] [optional] the fraction of a full circle at which to consider a difference in angle discontinuous.
                  )
{
   realT full;
   
   if(degrad) full = boost::math::constants::two_pi<realT>();
   else full = static_cast<realT>(360);

   threshold*=full;
   
   realT adj = 0;
   
   realT last = angles[0];
   
   if( fabs(angles[1] - angles[0]) > threshold)
   {
      if( angles[1] > angles[0]) adj = -full;
      else adj = full;
      
      angles[1] += adj;
   }
   
   if(angles.size() == 2) return 0;
   
   for(int i=2; i< angles.size(); ++i)
   {
      angles[i] += adj;
    
      if( fabs(angles[i] - angles[i-1]) > threshold)
      {
         if( angles[i] > angles[i-1]) adj += -full;
         else adj += full;
      
         angles[i] += adj;
      }
   }
   
   return 0;
}


/// @}

} //namespaace math
} //namespace mx

#endif //__geo_hpp__
