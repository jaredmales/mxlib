
//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef levinsonRecursion_hpp
#define levinsonRecursion_hpp


/*Solves the Toeplitz system PN
j=1 R(N+i−j)xj = yi (i = 0,...,N-1). The Toeplitz matrix need
not be symmetric. y[0..n-1] and r[0..2*n-2] are input arrays; x[0..n-1] is the output array.*/


template<typename floatT>
void levinsonRecursion(floatT * _r, floatT * _x, floatT * _y, int n)
{
   
   int j,k,m,m1,m2;
   floatT pp,pt1,pt2,qq,qt1,qt2,sd,sgd,sgn,shn,sxn;
   
   /*Have to adjust for the stupid NR convention of 1-counting*/
   floatT *r = _r-1;
   floatT *x = _x-1;
   floatT *y = _y-1;
   
   floatT *g,*h, *_g, *_h;

   if (r[n] == 0.0) 
   {
      std::cerr << "toeplz-1 singular principal minor" << "\n";
      return;
   }
   

   /*Have to adjust for the stupid NR convention of 1-counting*/
   _g = new floatT[n]; 
   g = _g-1;
   _h = new floatT[n]; 
   h = _h-1;
   
   x[1]=y[1]/r[n]; //Initialize for the recursion.

   if(n == 1) 
   {
      delete[] _g;
      delete[] _h;
      return;
   }
   

   g[1]=r[n-1]/r[n];

   h[1]=r[n+1]/r[n];

   for (m = 1; m <= n; ++m) 
   { //Main loop over the recursion.

      m1=m+1;
      sxn = -y[m1]; //Compute numerator and denominator for x,
      sd = -r[n];
      
      for (j=1;j <= m; ++j) 
      {
         sxn += r[n+m1-j]*x[j];
         sd += r[n+m1-j]*g[m-j+1];
      }

      if (sd == 0.0) 
      {
         std::cerr<< "toeplz-2 singular principal minor" << "\n";
         
         delete[] _g;
         delete[] _h;
         return;
      }
      
      x[m1]=sxn/sd; //whence x.

      for (j=1; j <= m; ++j) x[j] -= x[m1]*g[m-j+1];

      if (m1 == n) 
      {
         delete[] _g;
         delete[] _h;
         return;
      }

      sgn = -r[n-m1]; //Compute numerator and denominator for G and H,

      shn = -r[n+m1];

      sgd = -r[n];

      for (j=1;j <= m; ++j) 
      {
         sgn += r[n+j-m1]*g[j];
         shn += r[n+m1-j]*h[j];
         sgd += r[n+j-m1]*h[m-j+1];
      }
      
      if (sd == 0.0 || sgd == 0.0) 
      {
         std::cerr << "toeplz-3 singular principal minor" << "\n";
         
         delete[] _g;
         delete[] _h;
         return;
      }
      
      g[m1]=sgn/sgd; //whence G and H.
      h[m1]=shn/sd;

      k=m;
      m2=(m+1) >> 1;
      pp=g[m1];
      qq=h[m1];
      for (j = 1; j <= m2; ++j) 
      {
         pt1=g[j];
         pt2=g[k];
         qt1=h[j];
         qt2=h[k];
         g[j]=pt1-pp*qt2;
         g[k]=pt2-pp*qt1;
         h[j]=qt1-qq*pt2;  
         h[k--]=qt2-qq*pt1;
      }
   } //Back for another recurrence.

   std::cerr << "toeplz - should not arrive here!" << "\n";
   
   delete[] _g;
   delete[] _h;
   return;
   
}

#endif //levinsonRecursion_hpp
