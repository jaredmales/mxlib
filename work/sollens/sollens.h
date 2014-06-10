

#include <mx/astroconstants.h>


///Minimum focal distance from Sun (m)
#define Z_MIN ((RAD_SUN*RAD_SUN)/(2.*SRS_M))

///The coefficient for projection distance in focal plane (unitless)
#define X_FP_COEFF (Z_MIN/PC_M)

///The coefficient for the spot diameter
#define SPOTD_COEFF (2./(DPI*DPI)/SQRT_2*sqrt(1./SRS_M)*sqrt(Z_MIN)*1e-6)

///The coefficient for resolution
#define RES_COEFF (SPOTD_COEFF/X_FP_COEFF)

///The coefficient for v_foc in m/sec for a 1AU 1Msol orbit
#define V_FOC_COEFF (D2PI/(DJY*DAYSEC)*AU_M)


template<typename arithT>
arithT sollens_x_sp_rel( arithT d, arithT z)
{
   return X_FP_COEFF *z/d;
}

template<typename arithT>
arithT sollens_spotdiam(arithT lambda, arithT  z)
{
   return SPOTD_COEFF*lambda*sqrt(z);
}

template<typename arithT>
arithT sollens_resolution(arithT lambda, arithT d, arithT z)
{
   return RES_COEFF*d*lambda/sqrt(z);
}

///V_foc for a 1AU planet around a 1Msun star, in meters per sec
template<typename arithT>
arithT sollens_vfoc_1au1Msol(arithT d, arithT z)
{
   return V_FOC_COEFF*sollens_x_sp_rel(d, z);
}

///V_foc for a 1AU planet around a 1Msun star, in resolution elements per sec
template<typename arithT>
arithT sollens_vfoc_1au1Msol_resel(arithT lambda, arithT d, arithT z)
{
   return V_FOC_COEFF*X_FP_COEFF/SPOTD_COEFF/(d*lambda)*sqrt(z);
}

void makevect(double *x, double *y, int n)
{
   for(int i=0;i<n;i++)
   {
      x[i] = i;
      y[i] = exp(i);
   }
}