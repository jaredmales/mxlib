/** \file kepler.cpp
 * \author Jared R. Males
 * \brief Definitions of some utilities for working with the Kepler problem
 *
 */

//#include "kepler.hpp"






//long hyperbolic_kepler(floatT *E, floatT *D, floatT ecc, floatT mean_anom, floatT tol, long itmax)


// int solve_barker(floatT *f, floatT M)
// {
//    floatT u, y;
// 
//    u = 3*M;
//    if(M > 0) u+=SQRT_F(8*M*M+1);
//    else u-=SQRT_F(8*M*M+1);
// 
//    y = POW_F(fabs(u), .3333333);
//    if(y == 0) std::cerr << "0\n";
//    *f = 2*ATAN_F(1/y - y);
// 
//    if(M > 0) *f *=-1;
//    return 0;
// }

//****** Solution to Kepler's equation **********//
//Calculate the next iteration of Danby's quartic Newton-Raphson method.
// floatT kepler_danby_1(floatT e, floatT M, floatT Ei)
// {
//    floatT cosE, sinE, fi, fi1, fi2, fi3, di1, di2, di3;
//    
//    //These are expensive, do them just once.
//    cosE = cos(Ei);
//    sinE = sin(Ei);
//    
//    fi = Ei - e*sinE - M;
//    fi1 = 1.0 - e*cosE;
//    fi2 = e*sinE;
//    fi3 = e*cosE;
//    
//    di1 = -fi / fi1;
//    di2 = -fi/(fi1 + 0.5*di1*fi2);
//    di3 = -fi/(fi1 + 0.5*di2*fi2 + di2*di2*fi3/6.0);
//    
//    return Ei + di3;
// }

//Solve Kepler's equation using Danby's quartic Newton-Raphson method.
// long solve_kepler_danby(floatT *E, floatT *D, floatT e, floatT M, floatT tol, long itmax)
// {
//    long i;
//    floatT lastE, sinE, sign;
//    
//    sinE = sin(M);
//    sign = 1.0;
//    if(sinE < 0.0) sign = -1.0;
//    
//    (*E) = M + 0.85*e*sign;
//    
//    for(i=0; i < itmax; i++)
//    {
//       lastE = (*E);
//       (*E) = kepler_danby_1(e, M, (*E));
//       
//       //Test for convergence to within tol
//       //Make sure we have iterated at least twice to prevent early convergence
//       if(i > 0. && fabs((*E)-lastE) < tol)
//       {
//          (*D) = M - (*E-e*sin(*E));
//          return i;
//       }
//    }
//    
//    return -1;  //This means itmax exceeded.
// }


// long solve_kepler(floatT *E, floatT *D, floatT e, floatT M, floatT tol, long itmax)
// {   
//    if (e < 1.0)
//    {
//       return solve_kepler_danby(E, D, e, M, tol, itmax);
//    }
//    else return hyperbolic_kepler(E, D, e, M, tol, itmax);
//    
// }

/************ Difference Form *************************/
//Calculate the next iteration of Danby's quartic Newton-Raphson method.
// floatT kepler_danby_diff_1(floatT EC, floatT ES, floatT dM, floatT dEi)
// {
//    floatT cosE0, sinE0, cosdE, sindE, fi, fi1, fi2, fi3, di1, di2, di3;
//    
//    //These are expensive, do them just once.
//    //cosE0 = cos(E0);
//    //sinE0 = sin(E0);
//    cosdE = cos(dEi);
//    sindE = sin(dEi);
//    
//    fi = dEi - EC*sindE + ES*(1-cosdE) - dM;
//    fi1 = 1.0 - EC*cosdE + ES*sindE;
//    fi2 = EC*sindE + ES*cosdE;
//    fi3 = EC*cosdE - ES*sindE;
//    
//    di1 = -fi / fi1;
//    di2 = -fi/(fi1 + 0.5*di1*fi2);
//    di3 = -fi/(fi1 + 0.5*di2*fi2 + di2*di2*fi3/6.0);
//    
//    return dEi + di3;
// }

//Solve Kepler's equation using Danby's quartic Newton-Raphson method.
// long solve_keplerdiff(floatT *dE, floatT *D, floatT e, floatT EC, floatT ES, floatT dM, floatT tol, long itmax)
// {
//    long i;
//    floatT lastdE, sindE, sign;
//    
//    sindE = ES*cos(dM-ES) + EC*sin(dM-ES);
//    sign = 1.0;
//    if(sindE < 0.0) sign = -1.0;
//    
//    (*dE) = dM + 0.85*sign*e-ES;
//    
//    for(i=0; i < itmax; i++)
//    {
//       lastdE = (*dE);
//       (*dE) = kepler_danby_diff_1(EC, ES, dM, (*dE));
//       
//       //Test for convergence to within tol
//       //Make sure we have iterated at least twice to prevent early convergence
//       if(i > 0. && fabs((*dE)-lastdE) < tol)
//       {
//          (*D) = dM - (*dE - EC*sin(*dE) + ES*(1-cos(*dE)));
//          return i;
//       }
//    }
//    
//    return -1;  //This means itmax exceeded.
// }

// long rf_elements(floatT *r, floatT *f, floatT *E, floatT *D, floatT e, floatT M, floatT a, floatT tol, long itmax)
// {
//    long its;
//    
//    if(e < 0.0)
//    {
//       std::cerr << "e < 0 in rf_elements\n";
//       return -1;
//    }
// 
//    if(e == 1.0)
//    {
//       std::cerr << "e = 1 in rf_elements, which is not currently handled.\n";
//       return -1;
//    }
// 
//    its = solve_kepler(E, D,  e, M, tol, itmax);
//    
//    if(e > 1.0)
//    {
//       *f = 2.0*ATAN_F(SQRT_F((e+1.0)/(e-1.0))*TANH_F(*E/2.0));
//       *r = a*(1.0-e*COSH_F(*E));
//    } 
//    else if(e<1.0)
//    {
//       if(e == 0.0) 
//       {
//          *f = M;
//          *r = a;
//       }
//       else 
//       {
//          *f = 2.0*ATAN_F(SQRT_F((1.0+e)/(1.0-e))*TAN_F(*E/2.0));
//          *r = a * (1.0-e*COS_F(*E));
//       }
//    }
// 
//    return its;
// }




// int get_rv(floatT *rv, floatT * r, floatT *f, floatT *E, floatT *D, floatT t, floatT mstar, floatT msini, floatT a, floatT e, floatT t0, floatT w, floatT tol, long itmax)
// {
//    floatT fpa, M, m2;
//    long its;
//   
//    if(e < 0.0)
//    {
//       std::cerr << "e < 0.0 in get_rv\n";
//       exit(0);
//    }
//    if(e == 1.0)
//    {
//       std::cerr << "e == 1.0 in get_rv\n";
//       exit(0);
//    }
// 
//    m2 = 0.0;//msini/MR_JUDPITER;
//    if(e < 1.0) M = SQRT_F(GM_SOL_AUD*(mstar+m2) / POW_F(a,3))*(t-t0);
//    else if(a < 0) M = SQRT_F(GM_SOL_AUD*(mstar+m2) / POW_F(-a,3))*(t-t0);
//    else
//    {
//       std::cerr << "Bad conic parameters in get_rv\n";
//       std::cerr << "a=" << a << " e=" << e << "\n";
//       exit(0);
//    }
// 
//    its = rf_elements(r, f, E, D, e, M, a, tol, itmax);
// 
//    if(*r <= 0 || (*r > 2.0*a && a > 0.0))
//    {
//       std::cerr << "r <= 0 || r < 0.5a in get_rv\n";
//       std::cerr << "a=" << a << " e=" << e << " r=" << *r <<"\n";
//       return -1;
//    }
//    else 
//    {
//       //*rv = (msini/MR_JUDPITER)*SQRT_F((GM_SOL_AUD/(mstar+m2))*(2.0/(*r) - 1.0/a)) * AU_M /DAY_SEC;
//       *rv = (msini/MR_JUPITER)*SQRT_F((GM_SOL_AUD/(mstar+m2))/(a*(1-e*e))) * AU_M /DAYSEC;
//    }
//    //fpa = ATAN_F(e*SIN_F(*f)/(1+e*COS_F(*f)));
//    //*rv *= COS_F(*f + w - fpa);
//    *rv *= (COS_F(w+*f) + e*COS_F(w));
//    //*rv = (msini/MR_JUDPITER)*SQRT_F(GM_SOL_AUD/(mstar*a*(1-e*e)))*(cos(*f+w)+e*cos(w))*AU_M /DAY_SEC;
//    return its;
// }

// template<typename vectorT>
// int cartesian_orbit( vectorT &nx, vectorT &ny, vectorT &nz, const vectorT & t, const size_t N, double a, double P, double e, double t0, double i, double w, double W)
// {
//    long rv;
//    double M;
//    floatT r, _f, E, D;
//    
//    //Just do these calcs once
//    double cos_i = cos(i);
//    double sin_i = sin(i);
//    double cos_W = cos(W);
//    double sin_W = sin(W);
//    
//    double cos_wf, sin_wf;
//    
//    for(size_t j=0; j < N; j++)
//    {
//       M = MEANANOL(t[j], t0, P);
//       
//       rv = rf_elements(&r, &_f, &E, &D, e, M, a);
//       
//       if(rv < 0) return -1;
//       
//       //one calc
//       cos_wf = cos(w+_f);
//       sin_wf = sin(w+_f);
//       
//       nx[j] = r*(cos_W*cos_wf-sin_W*sin_wf*cos_i);
//       ny[j] = r*(sin_W*cos_wf+cos_W*sin_wf*cos_i);
//       nz[j] = r*sin_wf*sin_i;
//       
//    }
//    return 0;
// }

// int project_orbit( mx::Vectord &nx, mx::Vectord &ny, mx::Vectord &f, const mx::Vectord & t, double a, double P, double e, double tau, double i, double w, double W)
// {
//    long rv;
//    double M;
//    floatT r, _f, E, D;
//    for(int j=0; j<t.length(0); j++)
//    {
//       M = MEANANOL(t(j), tau, P);
//       
//       rv = rf_elements(&r, &_f, &E, &D, e, M, a);
//       
//       if(rv < 0) return -1;
//       
//       nx(j) = r*(cos(W)*cos(w+_f)-sin(W)*sin(w+_f)*cos(i));
//       ny(j) = r*(sin(W)*cos(w+_f)+cos(W)*sin(w+_f)*cos(i));
//       f(j) = _f;
//       
//    }
//    return 0;
// }
// 
// int get_orbit_phase(double &cos_alf, double f, double w, double inc)
// {
//    cos_alf = sin(w+f)*sin(inc);
//    return 0;
// }

// int get_orbit_phase(mx::Vectord &cos_alf, const mx::Vectord &f, double w, double inc)
// {
//    cos_alf.allocate(f.length(0));
//    
//    for(size_t i=0; i<f.length(0); i++)
//    {
//       cos_alf(i) = sin(f(i)+w)*sin(inc);
//    }
//    return 0;
// }

// int get_lambert_phasef(double &phi, const double cos_alf)
// {
//    double alf = acos(cos_alf);
//    phi = (sin(alf) + (DPI-alf)*cos_alf)/DPI;
//    
//    return 0;
// }

// int get_lambert_phasef(mx::Vectord &phi, const mx::Vectord &cos_alf)
// {
//    mx::Vectord alf = cos_alf.acos(); //acos(cos_alf);
//    
//    phi.allocate(cos_alf.length());
//    
//    for(size_t i=0;i<cos_alf.length(0); i++)
//    {
//       phi(i) = (sin(alf(i)) - (alf(i)-DPI)*cos_alf(i))/DPI;
//    }
//    
//    return 0;
// }



// int get_iW_1pt(mx::Vectord &i, mx::Vectord &W, double x, double y, double rho, double r, double f, double w)
// {
//    double coswf, coswf2, sinwf, sinwf2, r_rho, r_rho2, cosi, sini, sinW, cosW, bot, projx, projy;
// 
//    //std::cerr << x << " " << y << "\n";
//    
//    coswf = cos(w+f);
//    sinwf = sin(w+f);
//    //std::cout << "sinwf = " << sinwf << "\n";
//    
//    r_rho = rho/r;
//    
//    r_rho2 = r_rho*r_rho;
//    
//    coswf2 = coswf*coswf;
//    
//    if(coswf2 > r_rho2)
//    {
//       i(0) = 0;
//       W(0) = 0;
//       return -1;
//    }
//    
//    if(r_rho2 > 1.)
//    {
//       i(0) = 0;
//       W(0) = 0;
//       return -1;
//    }
//    
//    sinwf2 = sinwf*sinwf;
//    
//    if(sinwf2 == 0)
//    {
//       i(0) = 0;
//       W(0) = 0;
//       return -2;
//    }
//    
//    if(sinwf2 < (r_rho2 - coswf2))
//    {
//       i(0) = 0;
//       W(0) = 0;
//       return -3;
//    }
//    
//    
//    //sini = sqrt((1.-r_rho2)/sinwf2);
//    cosi = sqrt((r_rho2-coswf2)/sinwf2);
//    //std::cout << "sini = " << sini << "\n";
//    
//    //i[0] = asin(sini);//acos(cosi);
//    i(0) = acos(cosi);
//    i(1) = DPI - i(0); //save the second acos
//    //i[1] = asin(-1.*sini);//acos(-1.*cosi);
//    //cosi = cos(i[0]);
//    
//    bot = r*(coswf*coswf+sinwf*sinwf*cosi*cosi);
//    
//    cosW = (y*sinwf*cosi + x*coswf)/bot;
//    sinW = (y*coswf - x*sinwf*cosi)/bot;
//    
//    
//    W(0) = atan2(sinW, cosW);
//    
//    cosW = (-y*sinwf*cosi + x*coswf)/bot;
//    sinW = (y*coswf + x*sinwf*cosi)/bot;
//    
//    W(1) = atan2(sinW, cosW);
//    
//    return 0;
//    
// }

// int get_W_1pt_z0(double &W, double x, double y, double r, double f, double w)
// // {
//    double sinW, cosW;
//    
//    if(w == -f)
//    {
//       sinW = y/r;
//       cosW = x/r;
//    }
//    else if(w == DPI - f)
//    {
//       sinW = -y/r;
//       cosW = -x/r;
//    }
//    else return -1;
//    
//    W = atan2(sinW, cosW);
//    return 0;
// }

/*
int main()
{

   floatT E, r, f, x[3], v[3], xa[3], va[3], sun[3], sunv[3];
   kepler_elements k;

   k.To = 2453979.5;
   k.a  = 9.222655036162030E-01;
   k.e  = 1.910572784113192E-01;
   k.w  = 1.263964368281515E+02 * DPI / 180;
   k.W  = 2.044599682923624E+02 * DPI / 180;
   k.i  = 3.331323146322601E+00 * DPI / 180;
   k.M  = 6.141677776140474E+01 * DPI / 180;
   //k.f  = 8.260735077607306E+01 * DPI / 180;
   k.Tp = 2453924.309173170011;

   cartesian_elements(x, v, &k);
   std::cout.precision(16);

   xa[0] = 5.195787228952771E-01;
   xa[1] = 6.997694784948483E-01;
   xa[2] = -2.454767512223948E-02;

   sun[0] = 2.965806106647883E-03;
   sun[1] = 3.573956163006967E-03;
   sun[2] = -1.115886122514115E-04;

   sunv[0] = -4.600447924254690E-06;
   sunv[1] =  4.841253361277791E-06;
   sunv[2] =  4.585020044498540E-08;

   va[0] = -1.295640093206572E-02;
   va[1] =  1.388616930086464E-02;
   va[2] = -1.047600880742424E-03;

   MX::vector_add(x,x, sun,3);
   MX::vector_add(v,v, sunv,3);

   std::cout << xa[0] << "\t" << xa[1] << "\t" << xa[2] << "\n";
   std::cout << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
   std::cout << v[0] << "\t" << v[1] << "\t" << v[2] << "\n";
   std::cout << va[0] << "\t" << va[1] << "\t" << va[2] << "\n";

   return 0;
}
*/
/*
int main()
{
   floatT M, E, e,f,r, rv,w;
   std::cout << HUGE_VAL << " " << 1/HUGE_VAL << std::endl;
   std::cout << focus(1, 0.5) << " " << focus(0,1) << " " << focus(1,0.0) << std::endl;
   std::cout << semimaj(1, 0.5) << " " << semimaj(0,1) << std::endl;
   std::cout << eccent(1, 0.5) << " " << eccent(0,1) << " " << eccent(HUGE_VAL, 1) << std::endl;
   std::cout << eccent(semimaj(focus(1,.5),.5), focus(1,.5)) << std::endl;
}
   /*  for(int i=0; i< 1000; i++)
   {
      //M = -25*DPI + 25*DPI*((floatT)i/500.0);
      e = .5;//2.0*(floatT)i/1000;
      w = -2*DPI + 2*DPI*((floatT)i/500.0);
      //rf_elements(&r, &f, e, M, 1, -DPI);
      //solve_kepler(&E, e, M, 1e-8);
      get_rv(&rv, &r, &f, 1000, 1.0, 2.0, 2.0, e, 0.0, w);
      std::cout << i << " " << w << " " << r << " " << f << " " << rv << "\n";
      //if(e < 1) std::cout << M+e*sin(E) << "\n";
      //else std::cout << e*SINH_F(E) - M << "\n";
   }
}
/**/
