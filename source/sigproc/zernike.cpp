/** \file zernike.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Working with the Zernike polynomials.
  * 
  * \ingroup signal_processing_files
  * 
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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




#include "sigproc/zernike.hpp"

namespace mx
{
namespace sigproc 
{
   
int noll_nm( int & n,
             int & m,
             int j   
           )
{
   if( j < 1)
   {
      mxError("noll_nm", MXE_INVALIDARG, "The Noll index j cannot be less than 1 in the Zernike polynomials");
      return -1;
   }
   
   
   n = ceil(-1.5 + sqrt(0.25 + 2*j) - 1e-10); // 1e-10 is to avoid wrong rounding due to numerical precision

   int  jrem = j - (n * (n+1)/2+1);
   m = (jrem + (jrem % 2) * abs( (n % 2)-1) + fabs( (jrem % 2)-1) * (n %  2)) * (-math::func::sign( (j % 2)-0.5));
   
   return 0;
}

int noll_j( unsigned int n,
            int m
          )
{
    if( ((n-m) % 2) ) return -1; //Check if odd

    int mn = n % 4;

    int dm = 0;
    if( m >= 0 && (mn == 2 || mn == 3)) dm = 1;
    else if( m <= 0 && (mn == 0 || mn == 1)) dm = 1;

    int j = (n*(n+1)) / 2 + abs(m) + dm;

    return j;
}

int nZernRadOrd(unsigned int n)
{
    if(n % 2) //odd
    {
        return noll_j(n, -n);
    }
    else
    {
        return noll_j(n, n);
    }
}

//Explicit instantiations:
template
int zernikeRCoeffs<float>( std::vector<float> & c, int n, int m);

template
int zernikeRCoeffs<double>( std::vector<double> & c, int n, int m);

template
int zernikeRCoeffs<long double>( std::vector<long double> & c, int n, int m);

#ifdef HASQUAD
template
int zernikeRCoeffs<__float128>( std::vector<__float128> & c, int n, int m);
#endif

template
float zernikeR<float, double>(float rho, int n, int m, std::vector<double> & c);

template
double zernikeR<double, double>(double rho, int n, int m, std::vector<double> & c);

template
long double zernikeR<long double, long double>(long double rho, int n, int m, std::vector<long double> & c);

#ifdef HASQUAD
template
__float128 zernikeR<__float128,__float128>(__float128 rho, int n, int m, std::vector<__float128> & c);
#endif


template
float zernikeR<float, double>( float rho, int n, int m);

template
double zernikeR<double,double>( double rho, int n, int m);

template
long double zernikeR<long double, long double>( long double rho, int n, int m);

#ifdef HASQUAD
template
__float128 zernikeR<__float128, __float128>( __float128 rho, int n, int m);
#endif

template
float zernike<float,double>(float rho, float phi, int n, int m, std::vector<double> & c);

template
double zernike<double,double>(double rho, double phi, int n, int m, std::vector<double> & c);

template
long double zernike<long double, long double>(long double rho, long double phi, int n, int m, std::vector<long double> & c);

#ifdef HASQUAD
template
__float128 zernike<__float128,__float128>(__float128 rho, __float128 phi, int n, int m, std::vector<__float128> & c);
#endif

template
float zernike<float,double>( float rho, float phi, int n, int m);

template
double zernike<double,double>( double rho, double phi, int n, int m);

template
long double zernike<long double, long double>( long double rho, long double phi, int n, int m);

#ifdef HASQUAD
template
__float128 zernike<__float128,__float128>( __float128 rho, __float128 phi, int n, int m);
#endif


template
float zernike<float, double>(float rho, float phi, int j);

template
double zernike<double, double>(double rho, double phi, int j);

template
long double zernike<long double, long double>(long double rho, long double phi, int j);

#ifdef HASQUAD
template
__float128 zernike<__float128, __float128>(__float128 rho, __float128 phi, int j);
#endif

template
float zernikeQNorm<float>(float k, float phi, int n, int m);

template
double zernikeQNorm<double>(double k, double phi, int n, int m);

template
long double zernikeQNorm<long double>(long double k, long double phi, int n, int m);

#ifdef HASQUAD
template
__float128 zernikeQNorm<__float128>(__float128 k, __float128 phi, int n, int m);
#endif

} //namespace sigproc 
} //namespace mx


