#include "include/keplerdanby.h"


//Calculate the next iteration of Danby's quartic Newton-Raphson method.
double kepler_1(double e, double M, double Ei)
{
   double cosE, sinE, fi, fi1, fi2, fi3, di1, di2, di3;

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

//Solve Kepler's equation using Danby's quartic Newton-Raphson method.
int solve_kepler_danby(double *E, double *D, double e, double M, double tol, int itmax)
{
   int i;
   double lastE, sinE, sign;

   sinE = sin(M);
   sign = 1.0;
   if(sinE < 0.0) sign = -1.0;

   (*E) = M + 0.85*e*sign;

   for(i=0; i < itmax; i++)
   {
      lastE = (*E);
      (*E) = kepler_1(e, M, (*E));

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

/************ Difference Form *************************/
//Calculate the next iteration of Danby's quartic Newton-Raphson method.
double keplerdiff_1(double EC, double ES, double dM, double dEi)
{
   double cosE0, sinE0, cosdE, sindE, fi, fi1, fi2, fi3, di1, di2, di3;
   
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

//Solve Kepler's equation using Danby's quartic Newton-Raphson method.
int solve_keplerdiff_danby(double *dE, double *D, double e, double EC, double ES, double dM, double tol, int itmax)
{
   int i;
   double lastdE, sindE, sign;
   
   sindE = ES*cos(dM-ES) + EC*sin(dM-ES);
   sign = 1.0;
   if(sindE < 0.0) sign = -1.0;
   
   (*dE) = dM + 0.85*sign*e-ES;
   
   for(i=0; i < itmax; i++)
   {
      lastdE = (*dE);
      (*dE) = keplerdiff_1(EC, ES, dM, (*dE));
      
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

