/** \file orbitUtils.hpp
 * \author Jared R. Males
 * \brief Utilities related to orbits.
 * \ingroup astrofiles
 *
 */

#ifndef mx_astro_orbitUtils_hpp
#define mx_astro_orbitUtils_hpp

#include "../math/constants.hpp"
#include "kepler.hpp"

namespace mx
{
namespace astro
{

///Calculate the period of an orbit, given the masses and semi-major axis.
/** 
  * \tparam units specifies the unit system and precision, see \ref astrounits. 
  * 
  * \returns the period in the specified units
  * 
  * \ingroup orbits
  */ 
template<typename units>
typename units::realT orbitPeriod( typename units::realT m1, ///< [in] mass of one object. 
                                   typename units::realT m2, ///< [in] mass of the other object (can be 0).
                                   typename units::realT a ///< [in] the semi-major axis of the orbit.
                                 )
{
   return math::two_pi<typename units::realT>()* sqrt( (a*a*a) / (constants::G<units>()*(m1+m2)));
}


///Calculate the semi-major axis of an orbit, given the masses and Period, using SI units.
/** 
  * \tparam units specifies the unit system and precision, see \ref astrounits.
  * 
  * \returns the semi-major axis in the specified units
  * 
  * \ingroup orbits
  */ 
template<typename units>
typename units::realT orbitSemiMaj( typename units::realT m1, ///< [in] mass of one object.
                                    typename units::realT m2, ///< [in] mass of the other object (can be 0).
                                    typename units::realT P ///< [in] the period of the orbit.
                                  )
{
   typedef typename units::realT realT;
   
   constexpr realT two = static_cast<realT>(2); //Just for convenience.
   
   return  pow( pow(P/(math::two_pi<realT>()),two)  * constants::G<units>() *(m1+m2), math::third<realT>());
}

///Calculate the mean anomaly at time t, given the time of pericenter passage t0 and period P.  
/**
  * \ingroup orbits
  */
template<typename realT>
realT orbitMeanAnol( realT t,
                     realT t0,
                     realT P
                   )
{
   return math::two_pi<realT>()*(t-t0)/P;
}

///Calculate the separation and true anomaly for an orbit at the specified mean anomaly.
/** 
  * \tparam realT is the type used for arithmetic
  * 
  * \returns the numer of iterations in solving Kepler's equation, if successful.
  * \returns -1 if itmax is exceeded.
  * 
  * \ingroup orbits
  */
template<typename realT>
long orbitElements( realT & r,             ///< [out] is the orbital separation
                    realT & f,             ///< [out] is the true anomaly
                    realT E,               ///< [out] is the eccentric anomaly, provided since it is free
                    realT D,               ///< [out] is the error after the last iteration of the Kepler solution.
                    realT e,               ///< [in] is the eccentricity
                    realT M,               ///< [in] is the mean anomaly
                    realT a,               ///< [in] is the semi-major axis
                    realT tol=KEPLER_TOL,  ///< [in] is the desired tolerance
                    long itmax=KEPLER_ITMAX ///< [in] is the maximum number of iterations
                  )
{
   long its;
   
   if(e < 0.0)
   {
      std::cerr << "e < 0 in rf_elements\n";
      return -1;
   }

   if(e == 1.0)
   {
      std::cerr << "e = 1 in rf_elements, which is not currently handled.\n";
      return -1;
   }

   its = solve_kepler(E, D,  e, M, tol, itmax);
   
   if(e > 1.0)
   {
      f = 2.0*std::atan(std::sqrt((e+1.0)/(e-1.0))*std::tanh(E/2.0));
      r = a*(1.0-e*std::cosh(E));
   } 
   else if(e<1.0)
   {
      if(e == 0.0) 
      {
         f = M;
         r = a;
      }
      else 
      {
         f = 2.0*std::atan(std::sqrt((1.0+e)/(1.0-e))*std::tan(E/2.0));
         r = a * (1.0-e*std::cos(E));
      }
   }

   return its;
}

///Get the orbital phase at true anomaly f. Calculates the cos(alfa) where alfa is the orbital phase.
/** Only calculates cos(alfa) to save an operation when possible.
  *    
  * \tparam realT is the type used for arithmetic
  * 
  * \returns the cosine of the phase angle
  * 
  * \ingroup orbits
  */ 
template<typename realT>
realT orbitPhaseCosine( realT f,   ///< [in] is the true anomoaly at which to calculate alfa
                        realT w,   ///< [in] the argument of pericenter of the orbit
                        realT inc  ///< [in] is the inclinatin of the orbit
                      )
{
   return sin(f+w)*sin(inc);
}

///Get the lambertian phase function at orbital phase specified as cos(alfa), where alfa is the phase angle.
/** Uses cos(alfa) to save an operation after calculating orbit phase.
  * 
  * \tparam realT is the type used for arithmetic
  * 
  * \returns the lambert phase function at alpha
  * 
  * \ingroup orbits
  */
template<typename realT>
realT orbitLambertPhase( realT cos_alf /**< [in] the cosine of the phase angle*/)
{
   realT alf = std::acos(cos_alf);
   return (sin(alf) + (math::pi<realT>()-alf)*cos_alf)/math::pi<realT>();   
}

///Calculate the radial velocity of an orbiting body inclined relative to an observer.
/** \returns the radial velocity of the body with mass m1 as it orbits the barycenter of its orbit w.r.t. m2.
  * 
  * \ingroup orbits 
  */
template<typename units>
typename units::realT orbitRV( typename units::realT m1,  ///< [in] the mass of the body for which the RV is computed
                               typename units::realT m2,  ///< [in] the mass of the other body in the orbiting pair
                               typename units::realT inc, ///< [in] the inclination of the orbit of m2.
                               typename units::realT a,   ///< [in] the semi-major axis of the orbit
                               typename units::realT e,   ///< [in] the eccentricity of the orbit
                               typename units::realT w,   ///< [in] the argument of pericenter
                               typename units::realT f    ///< [in] the true anomaly at which the RV is to be calculated.
                             )
{
   typedef typename units::realT realT;
   
   realT msini = m2*sin(inc); 

   realT rv = (msini)*std::sqrt((constants::G<units>()/(m1+m2))/(a*(1-e*e)));
   
   rv *= (std::cos(w+f) + e*std::cos(w));
   
   return rv;
}


///Calculate various quantities of an orbit given the keplerian elements and a vector of times.
/** Only those quantities with non-null pointers are actually calculated.
  * 
  * The orbital parameters with dimensions, (a, P, and t0) must be in consistent units.
  * 
  * \tparam vectorT a type whose elements are accessed with the [] operator.  No other requirements are placed on this type, so it could be a native pointer (i.e. double *).
  * \tparam realT the type in which to perform calculations.  Does not need to match the storage type of vectorT.
  *  
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  * 
  * \ingroup orbits
  */
template<typename vectorT, typename realT>
int orbitCartesianWork( vectorT  * x,       ///< [out] [optional] If not NULL, will be the projected x positions of the orbit.  Must be at least as long as t.
                        vectorT  * y,       ///< [out] [optional] If not NULL, will be the projected y positions of the orbit.  Must be at least as long as t.
                        vectorT  * z,       ///< [out] [optional] If not NULL, will be the projected z positions of the orbit.  Must be at least as long as t. 
                        vectorT  * r,       ///< [out] [optional] If not NULL, will be the separation of the orbit.  Must be at least as long as t. 
                        vectorT  * rProj,   ///< [out] [optional] If not NULL, will be the projected separation of the orbit.  Must be at least as long as t. 
                        vectorT  * f,       ///< [out] [optional] If not NULL, will be the true anomaly of the orbit.  Must be at least as long as t. 
                        vectorT  * cos_alf, ///< [out] [optional] If not NULL, will be the phase angle of the orbit.  Must be at least as long as t.  
                        vectorT  * phi,     ///< [out] [optional] If not NULL, will be the Lambertian phase function.  Must be at least as long as t. 
                        const vectorT & t,  ///< [in] the times at which to calculate the projected positions
                        const size_t N,     ///< [in] the number of points contained in t, and allocated in nx, ny, nz.  Avoids requiring a size() member (i.e. if vectorT is native pointers).
                        realT a,            ///< [in] the semi-major axis of the orbit.
                        realT P,            ///< [in] the orbital period.
                        realT e,            ///< [in] the eccentricity of the orbit.
                        realT t0,           ///< [in] the time of pericenter passage of the orbit.
                        realT i,            ///< [in] the inclination of the orbit [radians].
                        realT w,            ///< [in] the argument of pericenter of the orbit [radians].
                        realT W             ///< [in] the longitude of the ascending node of the orbit [radians].
                      )
{
   long rv;
   realT _M, _r, _f, _E, _D;
   
   //Just do these calcs once
   realT cos_i = cos(i);
   realT sin_i = sin(i);
   realT cos_W = cos(W);
   realT sin_W = sin(W);
   
   realT cos_wf, sin_wf;
   
   for(size_t j=0; j < N; j++)
   {
      _M = orbitMeanAnol(t[j], t0, P);
      
      rv = orbitElements(_r, _f, _E, _D, e, _M, a);
      
      if(rv < 0) return -1;
      
      //one calc
      cos_wf = cos(w+_f);
      sin_wf = sin(w+_f);
      
      //Assumption is that these will be optimized away if NULL
      if(x) (*x)[j] = _r*(cos_W*cos_wf-sin_W*sin_wf*cos_i);
      if(y) (*y)[j] = _r*(sin_W*cos_wf+cos_W*sin_wf*cos_i);
      if(z) (*z)[j] = _r*sin_wf*sin_i;
      if(r) (*r)[j] = _r;
      if(rProj) 
      {
         if(x && y)
         {
            (*rProj)[j] = sqrt((*x)[j]*(*x)[j] + (*y)[j]*(*y)[j]);
         }
         else if(x && !y)
         {
            realT _y = _r*(sin_W*cos_wf+cos_W*sin_wf*cos_i);
            (*rProj)[j] = sqrt((*x)[j]*(*x)[j] + _y*_y);
         }
         else if(!x && y)
         {
            realT _x = _r*(cos_W*cos_wf-sin_W*sin_wf*cos_i);
            (*rProj)[j] = sqrt( _x*_x + (*y)[j]*(*y)[j]);
         }
         else
         {
            realT _x = _r*(cos_W*cos_wf-sin_W*sin_wf*cos_i);
            realT _y = _r*(sin_W*cos_wf+cos_W*sin_wf*cos_i);
            (*rProj)[j] = sqrt( _x*_x + _y*_y);
         }
      }
         
      if(f) (*f)[j] = _f;
      
      if(cos_alf) (*cos_alf)[j] = sin_wf*sin_i;
     
      if(phi) (*phi)[j] = orbitLambertPhase<realT>(sin_wf*sin_i);
      
   }
   return 0;
}


///Calculate the cartesian x-y position and the Lambert phase function of an orbit given Keplerian elements and a vector of times.
/**
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
  * 
  * \retval -1 on error calculating the orbit (from orbitCartesianWork)
  * \retval 0 on success.
  * 
  * \ingroup orbits
  */
template<typename vectorT, typename arithT>
int orbitCartesian2DPhi( vectorT &x, ///< [out] the projected x positions of the orbit.  Must be at least as long as t.
                         vectorT &y, ///< [out] the projected y positions of the orbit.  Must be at least as long as t.
                         vectorT &r, ///< [out] the separation of the orbit.  Must be at least as long as t.
                         vectorT &rProj, ///< [out] the projected separation of the orbit.  Must be at least as long as t.
                         vectorT &phi, ///< [out] the Lambert phase function.  Must be at least as long as t.
                         const vectorT & t, ///< [in] the times at which to calculate the projected positions
                         const size_t N,   ///< [in] the number of points contained in N, and allocated in nx, ny, nz
                         arithT a,  ///< [in] the semi-major axis of the orbit
                         arithT P,  ///< [in] the orbital period
                         arithT e,  ///< [in] the eccentricity of the orbit
                         arithT t0,  ///< [in] the time of pericenter passage of the orbit
                         arithT i,  ///< [in] the inclination of the orbit
                         arithT w,  ///< [in] the argument of pericenter of the orbit
                         arithT W  ///< [in] the longitude of the ascending node of the orbit.
                       )
{
   return orbitCartesianWork(&x, &y, (vectorT *)0, &r, &rProj, (vectorT *)0, (vectorT *)0, &phi, t, N, a, P, e, t0, i, w, W);
}



} //namespace astro
} //namespace mx

#endif //mx_astro_orbitUtils_hpp
