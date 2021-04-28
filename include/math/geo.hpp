/** \file geo.hpp
  * \author Jared R. Males
  * \brief Utilities for working with angles
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_geo_hpp
#define math_geo_hpp

#include <vector>

#include <cmath>

#include "constants.hpp"


namespace mx
{
namespace math
{
   

///Calculate the semi-latus rectum of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>
/**
  * \ingroup geo
  */  
#define semilatrect(a,e) (e == 0.0 ? a : (e == 1.0 ? 2.*a : (e < 1. ? a*(1-e*e) : a*(e*e-1))))

///Calculate the focal parameter of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>
/**
  * \ingroup geo
  */
#define focus(a,e) (e == 0.0 ? 1e34 : (e == 1.0 ? 2.*a : (e < 1. ? a*(1-e*e)/e : a*(e*e-1)/e)))

///Calculate the semi-major axis of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>, given the focal parameter and the eccentricity
/**
  * \ingroup geo
  */
#define semimaj(p,e) (e == 1.0 ? 1e34 : (e < 1 ? p*e/(1-e*e) : p*e/(e*e-1) ))

///Calculate the eccentricity of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a> given the semi-major axis and the focal parameter
/**
  * \ingroup geo
  */
#define eccent(a, p) (a == 0.0 ? 1e34 : (p >= 1e9 ? 0.0 : (p>0 ? (-p/(2*a)+0.5*std::sqrt(p*p/(a*a) + 4)) : (p/(2*a)+0.5*std::sqrt(p*p/(a*a) + 4)) ) ))


///Convert from degrees to radians
/**
  * 
  * \param q is the angle to convert
  * 
  * \return the angle q converted to radians
  * 
  * \ingroup geo
  */
template<typename realT>
realT dtor(realT q)
{
   return q * pi<realT>()/static_cast<realT>(180);
}

///Convert from radians to degrees
/**
  * 
  * \param q is the angle to convert
  * 
  * \return the angle q converted to degrees
  * 
  * \ingroup geo
  */
template<typename realT>
realT rtod(realT q)
{
   return q * static_cast<realT>(180)/pi<realT>();
}

///Calculate the angle modulo full-circle, normalizing to a positive value.
/** The output will be betweeen 0 and 360 (or 0 and 2pi).
  * 
  * \returns the value of q normalized to 0 <= q < 360[2pi]
  *
  * \tparam degrad controls whether this is in degrees (0, default) or radians (1)
  * \tparam realT is the type in which to do arithmetic
  * 
  * \ingroup geo
  */
template<int degrad = 0, typename realT>
realT angleMod(realT q /**< [in] the angle */)
{ 
   static_assert(std::is_floating_point<realT>::value, "angleMod: realT must be floating point");
   
   realT full;
   
   if(degrad) full = two_pi<realT>();
   else full = static_cast<realT>(360);
   
   q = fmod(q,full);
   
   if(q < 0) q += full;
   
   return q;
}

struct radians;

struct degrees;

template<typename degrad, typename realT>
struct degradT;

template<typename realT>
struct degradT<degrees,realT>
{
   static constexpr realT half = static_cast<realT>(180.0);
};

template<typename realT>
struct degradT<radians,realT>
{
   static constexpr realT half = pi<realT>();
};

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
  *  \ingroup geo
  */
template<typename degrad = radians, typename realT>
realT angleDiff( realT q1, ///< [in] angle to subtract from q2, in degrees.
                                   realT q2 ///< [in] angle to subtract q1 from, in degrees.
                                 )
{ 
   //typedef typename degradT::realT realT;
   
   static_assert(std::is_floating_point<realT>::value, "angleDiff: realT must be floating point");
   
   realT dq = q2-q1;

   realT half = degradT<degrad,realT>::half;//   degradT::half;
  
   if(std::abs(dq) > half)
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
  *  \ingroup geo
  */
template<int degrad = 0, typename realT>
realT angleMean(std::vector<realT> & q)
{ 
   static_assert(std::is_floating_point<realT>::value, "angleMean: realT must be floating point");
   
   realT s = 0;
   realT c = 0;

   realT d2r;
   
   if( degrad ) d2r = 1;
   else d2r  = pi<realT>()/static_cast<realT>(180);

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
  * 
  *  \ingroup geo
  */
template<int degrad = 0, typename realT>
int continueAngles( std::vector<realT> & angles, ///< [in] the vector of angles
                    realT threshold=0.75 ///< [in] [optional] the fraction of a full circle at which to consider a difference in angle discontinuous.
                  )
{
   realT full;
   
   if(degrad) full = two_pi<realT>();
   else full = static_cast<realT>(360);

   threshold*=full;
   
   realT adj = 0;
   
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

/// Rotate a point about the origin.
/** The rotation is counter-clockwise for positive angles.
  *
  * \tparam realT a real floating point type
  * 
  *  \ingroup geo
  */ 
template<typename realT>
void rotatePoint( realT & x0,  ///< [in/out] the x-coordinate of the point to rotate.  On exit contains the rotated value.
                  realT & y0,  ///< [in/out] the y-coordinate of the point to rotate.  On exit contains the rotated value.
                  realT angle  ///< [in] the angle by which to rotate [radians]
                )
{
   realT x1;
   realT y1;
   
   realT cq = cos(angle);
   realT sq = sin(angle);
   
   x1 = x0*cq - y0*sq;
   y1 = x0*sq + y0*cq;
   
   x0 = x1;
   y0 = y1;
}

} //namespaace math
} //namespace mx

#endif //math_geo_hpp
