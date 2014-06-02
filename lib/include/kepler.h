/** \file kepler.h
 * \author Jared R. Males
 * \brief Declarations for the utilities related to the Kepler problem.
 *
 */

#ifndef __kepler_h__
#define __kepler_h__

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "../vmop/MMatrix1"

#include "astroconstants.h"

//typedef long floatT floatT;

#include "astrotypes.h"


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
#define eccent(a, p) (a == 0.0 ? 1e34 : (p >= 1e9 ? 0.0 : (p>0 ? (-p/(2*a)+0.5*SQRT_F(p*p/(a*a) + 4)) : (p/(2*a)+0.5*SQRT_F(p*p/(a*a) + 4)) ) ))

///Calculate the period of an orbit, given the masses and semi-major axis, using solar units.
/** \param m1 is the first mass, in solar masses.
  * \param m2 is the second mass, in solar masses.
  * \param a is the semi-major axis, in AU
  * \retval T the period in days
  */ 
floatT get_period_solar(floatT m1, floatT m2, floatT a);//Msol, AU, days

///Calculate the period of an orbit, given the masses and semi-major axis, using Earth units.
/** \param m1 is the first mass, in Earth masses.
  * \param m2 is the second mass, in Earth masses.
  * \param a is the semi-major axis, in m
  * \retval T the period in secs
  */ 
floatT get_period_earth(floatT m1, floatT m2, floatT a);//Mearth, m, secs

floatT get_semimaj(floatT m1, floatT m2, floatT P);


///Calculate the mean anomaly at time t, given the time of pericenter passage t0 and period T.  
#define MEANANOL(t, t0, T) (2*PI*(t-t0)/T)

///Solve the hyperbolic kepler equation (for e> 1).
/** \param E [output] is the eccentric anomaly
  * \param D [output] is the error after the last iteration
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval -1 on exceeding itmax, otherwise reurns the number of iterations.
  */
long hyperbolic_kepler(floatT *E, floatT *D, floatT e, floatT M, floatT tol, long itmax);

//Calculate the next iteration of Danby's quartic Newton-Raphson method.
floatT kepler_1(floatT e, floatT M, floatT Ei);

///Solve Kepler's equation using Danby's quartic Newton-Raphson method.
/** \param E [output] is the eccentric anomaly
  * \param D [output] is the error after the last iteration
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval -1 on exceeding itmax, otherwise reurns the number of iterations.
  */
long solve_kepler_danby(floatT *E, floatT *D, floatT e, floatT M, floatT tol, long itmax);

///Solve Kepler's equation for any e.  Uses solve_kepler_danby if e < 1.0, hyperbolic_kepler otherwise.
/** \param E [output] is the eccentric anomaly
  * \param D [output] is the error after the last iteration
  * \param e [input] is the eccentricity
  * \param M [input] is the mean anomaly
  * \param tol [input] is the desired tolerance
  * \param itmax [input] is the maximum number of iterations
  * \retval -1 on exceeding itmax, otherwise reurns the number of iterations.
  */
long solve_kepler(floatT *E, floatT *D, floatT e, floatT M, floatT tol=KEPLER_TOL, long itmax=KEPLER_ITMAX);

//Calculate the next iteration of Danby's quartic Newton-Raphson method.
floatT keplerdiff_1(floatT EC, floatT ES, floatT dM, floatT dEi);

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
int solve_keplerdiff_danby(floatT *dE, floatT *D, floatT e, floatT EC, floatT ES, floatT dM, floatT tol, int itmax);

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
long rf_elements(floatT *r, floatT *f, floatT *E, floatT *D, floatT e, floatT M, floatT a, floatT tol=KEPLER_TOL, long itmax=KEPLER_ITMAX);

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
int get_rv(floatT *rv, floatT * r, floatT *f, floatT *E, floatT *D, floatT t, floatT mstar, floatT msini, floatT a, floatT e, floatT t0, floatT w, floatT tol=KEPLER_TOL, long itmax=KEPLER_ITMAX);

///Calculate the cartesian x-y-z position of an orbit given keplerian elements and a vector of times.
/** \param nx [output] the projected x positions of the orbit.  Should be same length as t.
  * \param ny [output] the projected y positions of the orbit.  Should be same length as t.
  * \param nz [output] the projected z positions of the orbit.  Should be same length as t.
  * \param t  [input] the times at which to calculate the projected positions
  * \param a  [input] the semi-major axis of the orbit
  * \param P  [input] the orbital period
  * \param e  [input] the eccentricity of the orbit
  * \param t0 [input] the time of pericenter passage of the orbit
  * \param i  [input] the inclination of the orbit
  * \param w  [input] the argument of pericenter of the orbit
  * \param W  [input] the longitude of the ascending node of the orbit.
  * \retval -1 on error calculating the orbit (from rf_elements)
  * \retval 0 on success.
  */
int cartesian_orbit( mx::Vectord &nx, mx::Vectord &ny, mx::Vectord &nz, const mx::Vectord & t, double a, double P, double e, double t0, double i, double w, double W);

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
int project_orbit( mx::Vectord &nx, mx::Vectord &ny, mx::Vectord &f, const mx::Vectord & t, double a, double P, double e, double t0, double i, double w, double W);

///Get the orbital phase at true anomaly f. Calculates the cos(alfa) where alfa is the orbital phase.
/** Only calculates cos(alfa) to save an operation when possible.
  * \param cos_alf [output] is the result
  * \param f [input] is the true anomoaly at which to calculate alfa
  * \param w [input] the argument of pericenter of the orbit
  * \param inc [input] is the inclinatin of the orbit
  * \retval -1 on error
  * \retval 0 on success
  */ 
int get_orbit_phase(double &cos_alf, double f, double w, double inc);

///Get the orbital phase at true anomaly f. Calculates the cos(alfa) where alfa is the orbital phase.
/** Only calculates cos(alfa) to save an operation when possible.
  * \param cos_alf [output] is the result
  * \param f [input] is the true anomoaly at which to calculate alfa
  * \param w [input] the argument of pericenter of the orbit
  * \param inc [input] is the inclinatin of the orbit
  * \retval -1 on error
  * \retval 0 on success
  */ 
int get_orbit_phase(mx::Vectord &cos_alf, const mx::Vectord &f, double w, double inc);

///Get the lambertian phase function at orbital phase specified as cos(alfa), where alfa is the phase angle.
/** Uses cos(alfa) to save an operation after calculating orbit phase.
  * \param phi [output] is the result
  * \param cos_alf [input] specifies the phase angle, as its cosine
  * \retval -1 on error
  * \retval 0 on success
  */
int get_lambert_phasef(double &phi, double cos_alf);

///Get the lambertian phase function at orbital phase specified as cos(alfa), where alfa is the phase angle.
/** Uses cos(alfa) to save an operation after calculating orbit phase.
  * \param phi [output] is the result
  * \param cos_alf [input] specifies the phase angle, as its cosine
  * \retval -1 on error
  * \retval 0 on success
  */
int get_lambert_phasef(mx::Vectord &phi, const mx::Vectord &cos_alf);

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
int get_iW_1pt(mx::Vectord &i, mx::Vectord &W, double x, double y, double rho, double r, double f, double w);

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
int get_W_1pt_z0(double &W, double x, double y, double r, double f, double w);

/// @}

#endif //__KEPLER_H__
