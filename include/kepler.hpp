/** \file kepler.hpp
 * \author Jared R. Males
 * \brief Declarations for the utilities related to the Kepler problem.
 *
 */

#ifndef __kepler_hpp__
#define __kepler_hpp__

#include <iostream>
#include <cmath>
#include <cstdlib>

//#include "../vmop/MMatrix1"

#include "astroconstants.h"

//typedef long floatT floatT;

//#include "astrotypes.h"


/** \addtogroup kepler 
  * @{
  */
   
///The default tolerance for solutions to the Kepler problem
#define KEPLER_TOL   1e-8

///The default maximum number of iterations for solutions to the Kepler problem
#define KEPLER_ITMAX 1000

//#define CUBE_ROOT( X)  (exp( log( X) / 3.))

///Calculate the semi-latus rectum of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>
#define semilatrect(a,e) (e == 0.0 ? a : (e == 1.0 ? 2.*a : (e < 1. ? a*(1-e*e) : a*(e*e-1))))

///Calculate the focal parameter of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>
#define focus(a,e) (e == 0.0 ? 1e34 : (e == 1.0 ? 2.*a : (e < 1. ? a*(1-e*e)/e : a*(e*e-1)/e)))

///Calculate the semi-major axis of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a>, given the focal parameter and the eccentricity
#define semimaj(p,e) (e == 1.0 ? 1e34 : (e < 1 ? p*e/(1-e*e) : p*e/(e*e-1) ))

///Calculate the eccentricity of a <a href="http://en.wikipedia.org/wiki/Conic_section">conic section</a> given the semi-major axis and the focal parameter
#define eccent(a, p) (a == 0.0 ? 1e34 : (p >= 1e9 ? 0.0 : (p>0 ? (-p/(2*a)+0.5*std::sqrt(p*p/(a*a) + 4)) : (p/(2*a)+0.5*std::sqrt(p*p/(a*a) + 4)) ) ))

///Calculate the period of an orbit, given the masses and semi-major axis, using solar units.
/** 
  * \param m1 is the first mass, in solar masses.
  * \param m2 is the second mass, in solar masses.
  * \param a is the semi-major axis, in AU
  * 
  * \returns the period in days
  */ 
template<typename arithT>
arithT get_period_solar(arithT m1, arithT m2, arithT a)
{
   return 2*DPI*std::sqrt((a*a*a)/(GM_SOL_AUD*(m1+m2)));
}

///Calculate the period of an orbit, given the masses and semi-major axis, using Earth units.
/** 
  * \param m1 is the first mass, in Earth masses.
  * \param m2 is the second mass, in Earth masses.
  * 
  * \param a is the semi-major axis, in m
  * 
  * \returns the period in secs
  */ 
template<typename arithT>
arithT get_period_earth(arithT m1, arithT m2, arithT a)
{
   return 2*DPI*std::sqrt((a*a*a)/(GRAV_G*(m1+m2)*MASS_EARTH));
}

///Calculate the semi-major axis of an orbit, given the masses and Period, using SI units.
/** 
  * \param m1 is the first mass, in kg.
  * \param m2 is the second mass, in kg.
  * 
  * \param P is the period, in secs
  * 
  * \returns the semi-major axis, in m
  *
  * \tparam arithT is the type used for arithmetic 
  */ 
template<typename arithT>
arithT get_semimaj(arithT m1, arithT m2, arithT P)
{
   return  pow( (P*P/(4.*DPI*DPI)) * GRAV_G*(m1+m2), 1./3.);
}


///Calculate the mean anomaly at time t, given the time of pericenter passage t0 and period P.  
#define MEANANOL(t, t0, P) (2.*DPI*((t)-(t0))/(P))

///Solve the hyperbolic kepler equation (for e> 1).
/** NOTE: this is not tested very well (as of 2014.08.09)
  * 
  * \param E [output] is the eccentric anomaly
  * \param D [output] is the error after the last iteration
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * 
  * \returns -1 on exceeding itmax, otherwise reurns the number of iterations.
  * 
  * \tparam arithT is the type used for arithmetic
  */
template<typename arithT>
long hyperbolic_kepler(arithT *E, arithT *D, arithT e, arithT M, arithT tol, long itmax)
{
   arithT curr, err, thresh;
   int is_negative = 0;
   long n_iter = 0;

   if(e == 1) e = 1.000001;

   if(M == 0) return( 0.);

   err = e * sinh( curr) - curr - M;
   while( fabs( err) > tol && n_iter < 10000)
   {
      n_iter++;
      curr -= err / (e * cosh( curr) - 1.);
      err = e * sinh( curr) - curr - M;
   }
   
   *E = ( is_negative ? -curr : curr);

   if(n_iter > itmax)
   {
      std::cerr << "hyperbolic kepler failed to converge e=" << e << "\n";
      return -1;
   }

   return n_iter;
}

///Calculate the next iteration of Danby's quartic Newton-Raphson method.
//floatT kepler_1(floatT e, floatT M, floatT Ei);
template<typename arithT>
arithT kepler_danby_1(arithT e, arithT M, arithT Ei)
{
   arithT cosE, sinE, fi, fi1, fi2, fi3, di1, di2, di3;
   
   //These are expensive, do them just once.
   cosE = cos(Ei);
   sinE = sin(Ei);
   
   fi = Ei - e*sinE - M;
   fi1 = 1.0 - e*cosE;
   fi2 = e*sinE;
   fi3 = e*cosE;
   
   di1 = -fi / fi1;
   di2 = -fi/(fi1 + 0.5*di1*fi2);
   di3 = -fi/(fi1 + 0.5*di2*fi2 + di2*di2*fi3/6.0);
   
   return Ei + di3;
}

///Solve Kepler's equation using Danby's quartic Newton-Raphson method.
/** \param E [output] is the eccentric anomaly
  * \param D [output] is the error after the last iteration
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval -1 on exceeding itmax, otherwise reurns the number of iterations.
  */
template<typename arithT>
long solve_kepler_danby(arithT *E, arithT *D, arithT e, arithT M, arithT tol, long itmax)
{
   long i;
   arithT lastE, sinE, sign;
   
   sinE = sin(M);
   sign = 1.0;
   if(sinE < 0.0) sign = -1.0;
   
   (*E) = M + 0.85*e*sign;
   
   for(i=0; i < itmax; i++)
   {
      lastE = (*E);
      (*E) = kepler_danby_1(e, M, (*E));
      
      //Test for convergence to within tol
      //Make sure we have iterated at least twice to prevent early convergence
      if(i > 0. && fabs((*E)-lastE) < tol)
      {
         (*D) = M - (*E-e*sin(*E));
         return i;
      }
   }
   
   return -1;  //This means itmax exceeded.
}

///Solve Kepler's equation for any e.  Uses solve_kepler_danby if e < 1.0, hyperbolic_kepler otherwise.
/** \param E [output] is the eccentric anomaly
  * \param D [output] is the error after the last iteration
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval -1 on exceeding itmax, otherwise reurns the number of iterations.
  */
template<typename arithT>
long solve_kepler(arithT *E, arithT *D, arithT e, arithT M, arithT tol=KEPLER_TOL, long itmax=KEPLER_ITMAX)
{   
   if (e < 1.0)
   {
      return solve_kepler_danby(E, D, e, M, tol, itmax);
   }
   else return hyperbolic_kepler(E, D, e, M, tol, itmax);
   
}

///Calculate the next iteration of Danby's quartic Newton-Raphson method (difference form).
//floatT keplerdiff_1(floatT EC, floatT ES, floatT dM, floatT dEi);
template<typename arithT>
arithT kepler_danby_diff_1(arithT EC, arithT ES, arithT dM, arithT dEi)
{
   arithT cosE0, sinE0, cosdE, sindE, fi, fi1, fi2, fi3, di1, di2, di3;
   
   //These are expensive, do them just once.
   //cosE0 = cos(E0);
   //sinE0 = sin(E0);
   cosdE = cos(dEi);
   sindE = sin(dEi);
   
   fi = dEi - EC*sindE + ES*(1-cosdE) - dM;
   fi1 = 1.0 - EC*cosdE + ES*sindE;
   fi2 = EC*sindE + ES*cosdE;
   fi3 = EC*cosdE - ES*sindE;
   
   di1 = -fi / fi1;
   di2 = -fi/(fi1 + 0.5*di1*fi2);
   di3 = -fi/(fi1 + 0.5*di2*fi2 + di2*di2*fi3/6.0);
   
   return dEi + di3;
}

///Solve Kepler's equation in difference form using Danby's quartic Newton-Raphson method.
/** \param dE [output] is the eccentric anomaly
 * \param D [output] is the error after the last iteration
 * \param e [input] is the eccentricity
 * \param EC [input]
 * \param ES [input]
 * \param dM [input]
 * \param tol [input] is the desired tolerance
 * \param itmax [input] is the maximum number of iterations
 * \retval -1 on exceeding itmax, otherwise reurns the number of iterations.
 */
template<typename arithT>
long solve_keplerdiff_danby(arithT *dE, arithT *D, arithT e, arithT EC, arithT ES, arithT dM, arithT tol, long itmax)
{
   long i;
   arithT lastdE, sindE, sign;
   
   sindE = ES*cos(dM-ES) + EC*sin(dM-ES);
   sign = 1.0;
   if(sindE < 0.0) sign = -1.0;
   
   (*dE) = dM + 0.85*sign*e-ES;
   
   for(i=0; i < itmax; i++)
   {
      lastdE = (*dE);
      (*dE) = kepler_danby_diff_1(EC, ES, dM, (*dE));
      
      //Test for convergence to within tol
      //Make sure we have iterated at least twice to prevent early convergence
      if(i > 0. && fabs((*dE)-lastdE) < tol)
      {
         (*D) = dM - (*dE - EC*sin(*dE) + ES*(1-cos(*dE)));
         return i;
      }
   }
   
   return -1;  //This means itmax exceeded.
}

///Calculate the radius and true anomaly for an orbit at the specified mean anomaly.
/** \param r [output] is the radius
  * \param f [output] is the true anomaly
  * \param E [output] is the eccentric anomaly, provided since it is free
  * \param D [output] is the error after the last iteration of the Kepler solution.
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param a [input] is the semi-major axis
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval -1 on exceeding itmax, otherwise reurns the number of iterations in solving Kepler's equation.
  */
template<typename arithT>
long rf_elements(arithT *r, arithT *f, arithT *E, arithT *D, arithT e, arithT M, arithT a, arithT tol=KEPLER_TOL, long itmax=KEPLER_ITMAX)
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
      *f = 2.0*std::atan(std::sqrt((e+1.0)/(e-1.0))*std::tanh(*E/2.0));
      *r = a*(1.0-e*std::cosh(*E));
   } 
   else if(e<1.0)
   {
      if(e == 0.0) 
      {
         *f = M;
         *r = a;
      }
      else 
      {
         *f = 2.0*std::atan(std::sqrt((1.0+e)/(1.0-e))*std::tan(*E/2.0));
         *r = a * (1.0-e*std::cos(*E));
      }
   }

   return its;
}

///Calculate the radial velocity at a point in an orbit specfied by time t.
/** \param rv [output] the radial velocity in meters per second.
  * \param r [output] is the radius
  * \param f [output] is the true anomaly
  * \param t [input] the time
  * \param mstar [input] is the mass of the central star
  * \param msini [input] is the mass of the orbiting body, times the sin of inclination.  Set to 0 if unknown.
  * \param a [input] is the semi-major axis.
  * \param e [input] is the eccentricity
  * \param t0 [input] is the time of pericenter passage
  * \param w [input] is the argument of pericenter
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval 0 on success, -1 otherwise.
  * \todo get_rv needs to be put into the same format as other astro functions
  * \todo specify units for get_rv
  * \todo error handling in astrodynamics!
  */
template<typename arithT>
long get_rv(arithT *rv, arithT * r, arithT *f, arithT *E, arithT *D, arithT t, arithT mstar, arithT msini, arithT a, arithT e, arithT t0, arithT w, arithT tol=KEPLER_TOL, long itmax=KEPLER_ITMAX)
//int get_rv(arithT *rv, arithT * r, arithT *f, arithT *E, arithT *D, arithT t, arithT mstar, arithT msini, arithT a, arithT e, arithT t0, arithT w, arithT tol, long itmax)
{
   arithT fpa, M, m2;
   long its;
  
   if(e < 0.0)
   {
      std::cerr << "e < 0.0 in get_rv\n";
      exit(0);
   }
   if(e == 1.0)
   {
      std::cerr << "e == 1.0 in get_rv\n";
      exit(0);
   }

   m2 = 0.0;//msini/MR_JUDPITER;
   if(e < 1.0) M = std::sqrt(GM_SOL_AUD*(mstar+m2) / POW_F(a,3))*(t-t0);
   else if(a < 0) M = std::sqrt(GM_SOL_AUD*(mstar+m2) / POW_F(-a,3))*(t-t0);
   else
   {
      std::cerr << "Bad conic parameters in get_rv\n";
      std::cerr << "a=" << a << " e=" << e << "\n";
      exit(0);
   }

   its = rf_elements(r, f, E, D, e, M, a, tol, itmax);

   if(*r <= 0 || (*r > 2.0*a && a > 0.0))
   {
      std::cerr << "r <= 0 || r < 0.5a in get_rv\n";
      std::cerr << "a=" << a << " e=" << e << " r=" << *r <<"\n";
      return -1;
   }
   else 
   {
      //*rv = (msini/MR_JUDPITER)*std::sqrt((GM_SOL_AUD/(mstar+m2))*(2.0/(*r) - 1.0/a)) * AU_M /DAY_SEC;
      *rv = (msini/MR_JUPITER)*std::sqrt((GM_SOL_AUD/(mstar+m2))/(a*(1-e*e))) * AU_M /DAYSEC;
   }
   //fpa = std::atan(e*SIN_F(*f)/(1+e*std::cos(*f)));
   //*rv *= std::cos(*f + w - fpa);
   *rv *= (std::cos(w+*f) + e*std::cos(w));
   //*rv = (msini/MR_JUDPITER)*std::sqrt(GM_SOL_AUD/(mstar*a*(1-e*e)))*(cos(*f+w)+e*cos(w))*AU_M /DAY_SEC;
   return its;
}

///Get the lambertian phase function at orbital phase specified as cos(alfa), where alfa is the phase angle.
/** Uses cos(alfa) to save an operation after calculating orbit phase.
  * 
  * \tparam arithT is the type used for arithmetic
  * 
  * \param cos_alf [input] specifies the phase angle, as its cosine
  * 
  * \retval the lambert phase function
  */
template<typename arithT>
arithT get_lambert_phi(const arithT cos_alf)
{
   arithT alf = std::acos(cos_alf);
   return (sin(alf) + (DPI-alf)*cos_alf)/DPI;   
}

//cartesian_orbit_work
///Calculate various quantities of an orbit given the keplerian elements and a vector of times.
/** Only those quantities with non-null pointers are actually calculated.
  * 
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
  *  
  * \param x       [output] the projected x positions of the orbit.  Must be at least as long as t.
  * \param y       [output] the projected y positions of the orbit.  Must be at least as long as t.
  * \param z       [output] the projected z positions of the orbit.  Must be at least as long as t.
  * \param f       [output] the true anomaly of the orbit.  Must be at least as long as t.
  * \param cos_alf [output] the phase angle of the orbit.  Must be at least as long as t.
  * \param phi  [output] the Lambertian phase function.  Must be at least as long as t.
  * \param t       [input] the times at which to calculate the projected positions
  * \param N       [input] the number of points contained in N, and allocated in nx, ny, nz
  * \param a       [input] the semi-major axis of the orbit
  * \param P       [input] the orbital period
  * \param e       [input] the eccentricity of the orbit
  * \param t0      [input] the time of pericenter passage of the orbit
  * \param i       [input] the inclination of the orbit
  * \param w       [input] the argument of pericenter of the orbit
  * \param W       [input] the longitude of the ascending node of the orbit.
  * 
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
template<typename vectorT, typename arithT>
int cartesian_orbit_work( vectorT  * x, 
                          vectorT  * y, 
                          vectorT  * z,
                          vectorT  * f,
                          vectorT  * cos_alf,
                          vectorT  * phi,
                          const vectorT & t, 
                          const size_t N,  
                          arithT a, 
                          arithT P, 
                          arithT e, 
                          arithT t0, 
                          arithT i, 
                          arithT w, 
                          arithT W)
{
   long rv;
   arithT M;
   arithT r, _f, E, D;
   
   //Just do these calcs once
   arithT cos_i = cos(i);
   arithT sin_i = sin(i);
   arithT cos_W = cos(W);
   arithT sin_W = sin(W);
   
   arithT cos_wf, sin_wf;
   
   for(size_t j=0; j < N; j++)
   {
      M = MEANANOL(t[j], t0, P);
      
      rv = rf_elements(&r, &_f, &E, &D, e, M, a);
      
      if(rv < 0) return -1;
      
      //one calc
      cos_wf = cos(w+_f);
      sin_wf = sin(w+_f);
      
      //Assumption is that these will be optimized away if NULL
      if(x) (*x)[j] = r*(cos_W*cos_wf-sin_W*sin_wf*cos_i);
      if(y) (*y)[j] = r*(sin_W*cos_wf+cos_W*sin_wf*cos_i);
      if(z) (*z)[j] = r*sin_wf*sin_i;
      
      if(f) (*f)[j] = _f;
      
      if(cos_alf) (*cos_alf)[j] = sin_wf*sin_i;
     
      if(phi) (*phi)[j] = get_lambert_phi<arithT>(sin_wf*sin_i);
      
   }
   return 0;
}


///Calculate the cartesian x-y position of an orbit given keplerian elements and a vector of times.
/**
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
  * 
  * \param nx [output] the projected x positions of the orbit.  Must be at least as long as t.
  * \param ny [output] the projected y positions of the orbit.  Must be at least as long as t.
  * \param t  [input] the times at which to calculate the projected positions
  * \param N  [input] the number of points contained in N, and allocated in nx, ny, nz
  * \param a  [input] the semi-major axis of the orbit
  * \param P  [input] the orbital period
  * \param e  [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i  [input] the inclination of the orbit
  * \param w  [input] the argument of pericenter of the orbit
  * \param W  [input] the longitude of the ascending node of the orbit.
  * 
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
template<typename vectorT, typename arithT>
int cartesian_orbit_2D( vectorT &nx, 
                        vectorT &ny, 
                        const vectorT & t, 
                        const size_t N,  
                        arithT a, 
                        arithT P, 
                        arithT e, 
                        arithT t0, 
                        arithT i, 
                        arithT w, 
                        arithT W)
{
   return cartesian_orbit_work(&nx, &ny, 0, 0, 0, 0, t, N, a, P, e, t0, i, w, W);
}


///Calculate the cartesian x-y-z position of an orbit given keplerian elements and a vector of times.
/**
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
  * 
  * \param nx [output] the projected x positions of the orbit.  Must be at least as long as t.
  * \param ny [output] the projected y positions of the orbit.  Must be at least as long as t.
  * \param nz [output] the projected z positions of the orbit.  Must be at least as long as t.
  * \param t  [input] the times at which to calculate the projected positions
  * \param N  [input] the number of points contained in N, and allocated in nx, ny, nz
  * \param a  [input] the semi-major axis of the orbit
  * \param P  [input] the orbital period
  * \param e  [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i  [input] the inclination of the orbit
  * \param w  [input] the argument of pericenter of the orbit
  * \param W  [input] the longitude of the ascending node of the orbit.
  * 
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
template<typename vectorT, typename arithT>
int cartesian_orbit3D( vectorT &nx, 
                       vectorT &ny, 
                       vectorT &nz, 
                       const vectorT & t, 
                       const size_t N,  
                       arithT a, 
                       arithT P, 
                       arithT e, 
                       arithT t0, 
                       arithT i, 
                       arithT w, 
                       arithT W)
{
   return cartesian_orbit_work(&nx, &ny, &nz, 0, 0, 0, t, N, a, P, e, t0, i, w, W);
}

///Calculate the cartesian x-y position and true anomaly of an orbit given keplerian elements and a vector of times.
/**
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
  * 
  * \param nx [output] the projected x positions of the orbit.  Must be at least as long as t.
  * \param ny [output] the projected y positions of the orbit.  Must be at least as long as t.
  * \param f  [output] the true anomaly of the orbit.  Must be at least as long as t.
  * \param t  [input] the times at which to calculate the projected positions
  * \param N  [input] the number of points contained in N, and allocated in nx, ny, nz
  * \param a  [input] the semi-major axis of the orbit
  * \param P  [input] the orbital period
  * \param e  [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i  [input] the inclination of the orbit
  * \param w  [input] the argument of pericenter of the orbit
  * \param W  [input] the longitude of the ascending node of the orbit.
  * 
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
template<typename vectorT, typename arithT>
int cartesian_orbit2D_f( vectorT &nx, 
                         vectorT &ny, 
                         vectorT &f, 
                         const vectorT & t, 
                         const size_t N,  
                         arithT a, 
                         arithT P, 
                         arithT e, 
                         arithT t0, 
                         arithT i, 
                         arithT w, 
                         arithT W )
{
   return cartesian_orbit_work<vectorT, arithT>(&nx, &ny, 0, &f, 0, 0, t, N, a, P, e, t0, i, w, W);
}


///Calculate the cartesian x-y-z position of an orbit and the true anomaly given keplerian elements and a vector of times.
/**
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
* 
  * \param nx [output] the projected x positions of the orbit.  Must be at least as long as t.
  * \param ny [output] the projected y positions of the orbit.  Must be at least as long as t.
  * \param nz [output] the projected z positions of the orbit.  Must be at least as long as t.
  * \param f  [output] the true anomaly of the orbit.  Must be at least as long as t.
  * \param t  [input] the times at which to calculate the projected positions
  * \param N  [input] the number of points contained in N, and allocated in nx, ny, nz
  * \param a  [input] the semi-major axis of the orbit
  * \param P  [input] the orbital period
  * \param e  [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i  [input] the inclination of the orbit
  * \param w  [input] the argument of pericenter of the orbit
  * \param W  [input] the longitude of the ascending node of the orbit.
  * 
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
template<typename vectorT, typename arithT>
int cartesian_orbit3D_f( vectorT &nx, 
                         vectorT &ny, 
                         vectorT &nz, 
                         vectorT &f, 
                         const vectorT & t, 
                         const size_t N,  
                         arithT a, 
                         arithT P, 
                         arithT e, 
                         arithT t0, 
                         arithT i, 
                         arithT w, 
                         arithT W )
{
   return cartesian_orbit_work(&nx, &ny, &nz, &f, 0, 0, t, N, a, P, e, t0, i, w, W);
}

///Calculate the cartesian x-y position and the Lambert phase function of an orbit given keplerian elements and a vector of times.
/**
  * \tparam vectorT a type whose elements are accessed with the [] operator.
  * \tparam arithT the type in which to perform arithmetic
  * 
  * \param nx [output] the projected x positions of the orbit.  Must be at least as long as t.
  * \param ny [output] the projected y positions of the orbit.  Must be at least as long as t.
  * \param phi  [output] the true anomaly of the orbit.  Must be at least as long as t.
  * \param t  [input] the times at which to calculate the projected positions
  * \param N  [input] the number of points contained in N, and allocated in nx, ny, nz
  * \param a  [input] the semi-major axis of the orbit
  * \param P  [input] the orbital period
  * \param e  [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i  [input] the inclination of the orbit
  * \param w  [input] the argument of pericenter of the orbit
  * \param W  [input] the longitude of the ascending node of the orbit.
  * 
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
template<typename vectorT, typename arithT>
int cartesian_orbit2D_phi( vectorT &nx, 
                         vectorT &ny, 
                         vectorT &phi, 
                         const vectorT & t, 
                         const size_t N,  
                         arithT a, 
                         arithT P, 
                         arithT e, 
                         arithT t0, 
                         arithT i, 
                         arithT w, 
                         arithT W )
{
   return cartesian_orbit_work(&nx, &ny, (vectorT *)0, (vectorT *)0, (vectorT *)0, &phi, t, N, a, P, e, t0, i, w, W);
}



///Calculate the projection of an orbit onto the reference plane for a set of times.
/** NOTE: this is depreciated.  Use cartesian_orbit instead.
  * 
  * \param nx [output] the projected x positions of the orbit
  * \param ny [output] the projected y positions of the orbit
  * \param t [input] the times at which to calculate the projected positions
  * \param a [input] the semi-major axis of the orbit
  * \param P [input] the orbital period
  * \param e [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i [input] the inclination of the orbit
  * \param w [input] the argument of pericenter of the orbit
  * \param W [input] the longitude of the ascending node of the orbit.
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
//int project_orbit( mx::Vectord &nx, mx::Vectord &ny, mx::Vectord &f, const mx::Vectord & t, double a, double P, double e, double t0, double i, double w, double W);

///Get the orbital phase at true anomaly f. Calculates the cos(alfa) where alfa is the orbital phase.
/** Only calculates cos(alfa) to save an operation when possible.
  * \param cos_alf [output] is the result
  * \param f [input] is the true anomoaly at which to calculate alfa
  * \param w [input] the argument of pericenter of the orbit
  * \param inc [input] is the inclinatin of the orbit
  * \retval -1 on error
  * \retval 0 on success
  */ 
int get_orbit_phase(double &cos_alf, double f, double w, double inc)
{
   cos_alf = sin(f+w)*sin(inc);
}

///Get the orbital phase at true anomaly f. Calculates the cos(alfa) where alfa is the orbital phase.
/** Only calculates cos(alfa) to save an operation when possible.
  * \param cos_alf [output] is the result
  * \param f [input] is the true anomoaly at which to calculate alfa
  * \param w [input] the argument of pericenter of the orbit
  * \param inc [input] is the inclinatin of the orbit
  * \retval -1 on error
  * \retval 0 on success
  */ 
template<typename vectorT>
int get_orbit_phase(vectorT &cos_alf, const vectorT &f, const double w, const double inc)
{
   cos_alf.allocate(f.size());
   
   for(size_t i=0; i<f.size(); i++)
   {
      cos_alf[i] = sin(f[i]+w)*sin(inc);
   }
   return 0;
}





///Get the lambertian phase function at orbital phase specified as cos(alfa), where alfa is the phase angle.
/** Uses cos(alfa) to save an operation after calculating orbit phase.
  * \param phi [output] is the result
  * \param cos_alf [input] specifies the phase angle, as its cosine
  * \retval -1 on error
  * \retval 0 on success
  */
//int get_lambert_phi(mx::Vectord &phi, const mx::Vectord &cos_alf);

///Calculate the inclination and the longitude of the ascending node given a point on the orbit projected onto the plane of the sky, the argument of pericenter, the current true radius, and the current true anomaly.
/** Solutions are not always possible.  You must check the return value to know if the results are valid.
 * \param i [output] the inclination (there will be 2)
 * \param W [output] the longitude of the ascending node (there will be 2)
 * \param x [input] the x coordinate of the projection onto the reference plane
 * \param y [input] the y coordinate of the projection onto the reference plane
 * \param rho [input] the projected separation (yes this is redundant to x and y, but it is passed independently for efficiency)
 * \param r [input] the physical (3D) separation
 * \param f [input] the true anomaly
 * \param w [input] the argument of pericenter
 * \retval 0 on success
 * \retval -1 if no solution is possible due to ...
 * \retval -2 if no solution is possible due to z = 0
 * \retval -3 if no solution is possible due to ...
 */
//int get_iW_1pt(mx::Vectord &i, mx::Vectord &W, double x, double y, double rho, double r, double f, double w);

///Calculate the longitude of the ascending node given a point on the orbit projected onto the plane of the sky, the argument of pericenter, the current radius, and the current true anomaly.
/** Used for the z = 0, i != 0/PI case.  i.e. sin(w+f) = 0, cos(w+f) = +/- 1.
 * w must be either -f or PI - f.
 * \param W [output] the longitude of the ascending node.
 * \param x [input] the x coordinate of the starting point
 * \param y [input] the y coordinate of the starting point
 * \param r [input] the current 3D separation from the star
 * \param f [input] the current true anomaly
 * \param w [input] the argument of pericenter
 * \retval 0 on success
 * \retval -1 if the conditions aren't  met by w an df
 */
template<typename arithT>
int get_W_1pt_z0(arithT &W, arithT x, arithT y, arithT r, arithT f, arithT w)
{
   arithT sinW, cosW;
   
   if(w == -f)
   {
      sinW = y/r;
      cosW = x/r;
   }
   else if(w == DPI - f)
   {
      sinW = -y/r;
      cosW = -x/r;
   }
   else return -1;
   
   W = atan2(sinW, cosW);
   return 0;
}
/// @}

#endif //__kepler_hpp__
