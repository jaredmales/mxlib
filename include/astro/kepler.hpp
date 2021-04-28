/** \file kepler.hpp
 * \author Jared R. Males
 * \brief Declarations for the utilities related to the Kepler problem.
 * \ingroup astrofiles
 *
 */

#ifndef __mx_astro_kepler_hpp__
#define __mx_astro_kepler_hpp__

#include <iostream>
#include <cmath>
#include <cstdlib>


#include "units.hpp"

namespace mx
{
   
namespace astro
{


   
#ifndef KEPLER_TOL
///The default tolerance for solutions to the Kepler problem
#define KEPLER_TOL (1e-8)
#endif

#ifndef KEPLER_ITMAX
///The default maximum number of iterations for solutions to the Kepler problem
#define KEPLER_ITMAX (1000)
#endif




///Solve the hyperbolic kepler equation (for e> 1).
/** \note this is not tested very well (as of 2014.08.09)
  * 
  * 
  * \returns -1 on exceeding itmax, otherwise reurns the number of iterations.
  * 
  * \tparam realT is the type used for arithmetic
  * 
  * \ingroup kepler
  */
template<typename realT>
long hyperbolic_kepler( realT & E,   ///< [out] is the eccentric anomaly
                        realT & err, ///< [out] is the error after the last iteration
                        realT e,     ///< [in] is the eccentricity
                        realT M,     ///< [in] is the mean anomaly
                        realT tol,   ///< [in] is the desired tolerance
                        long itmax   ///< [in] is the maximum number of iterations
                      )
{
   realT curr, thresh;
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
   
   E = ( is_negative ? -curr : curr);

   if(n_iter > itmax)
   {
      std::cerr << "hyperbolic kepler failed to converge e=" << e << "\n";
      return -1;
   }

   return n_iter;
}

///Calculate the next iteration of Danby's quartic Newton-Raphson method.
/** \ingroup kepler
  */
template<typename realT>
realT kepler_danby_1( realT e, 
                      realT M, 
                      realT Ei)
{
   realT cosE, sinE, fi, fi1, fi2, fi3, di1, di2, di3;
   
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
/** 
  * \returns -1 on exceeding itmax.
  * \returns the number of iterations on success.
  * 
  * \ingroup kepler
  */
template<typename realT>
long solve_kepler_danby( realT & E, ///< [out] is the eccentric anomaly
                         realT & D, ///< [out] is the error after the last iteration
                         realT e,   ///< [in] is the eccentricity
                         realT M,   ///< [in] is the mean anomaly
                         realT tol, ///< [in] is the desired tolerance
                         long itmax  ///< [in] is the maximum number of iterations
                       )
{
   long i;
   realT lastE, sinE, sign;
   
   sinE = sin(M);
   sign = 1.0;
   if(sinE < 0.0) sign = -1.0;
   
   E = M + 0.85*e*sign;
   
   for(i=0; i < itmax; i++)
   {
      lastE = E;
      E = kepler_danby_1(e, M, E);
      
      //Test for convergence to within tol
      //Make sure we have iterated at least twice to prevent early convergence
      if(i > 1 && fabs(E-lastE) < tol)
      {
         D = M - (E-e*sin(E));
         return i;
      }
   }
   
   return -1;  //This means itmax exceeded.
}

///Solve Kepler's equation for any e.  Uses solve_kepler_danby if e < 1.0, hyperbolic_kepler otherwise.
/** 
  * \returns the number of iterations on success.
  * \returns -1 if itmax is exceeded.
  * 
  * \ingroup kepler
  */
template<typename realT>
long solve_kepler( realT & E, ///< [out] is the eccentric anomaly
                   realT & D, ///< [out] is the error after the last iteration
                   realT e,   ///< [in] is the eccentricity
                   realT M,   ///< [in] is the mean anomaly
                   realT tol=KEPLER_TOL,  ///< [in] is the desired tolerance
                   long itmax=KEPLER_ITMAX  ///< [in] is the maximum number of iterations
                 )
{   
   if (e < 1.0)
   {
      return solve_kepler_danby(E, D, e, M, tol, itmax);
   }
   else return hyperbolic_kepler(E, D, e, M, tol, itmax);
   
}





#if 0

#endif


#if 0

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

#endif

} //namespace astro
} //namespace mx
#endif //__mx_astro_kepler_hpp__
