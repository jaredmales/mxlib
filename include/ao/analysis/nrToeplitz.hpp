#ifndef __nrToeplitz_hpp__
#define __nrToeplitz_hpp__

#define FREERETURN { delete _g; delete _h; return;}

/*Solves the Toeplitz system PN
j=1 R(N+i−j)xj = yi (i = 0,...,N-1). The Toeplitz matrix need
not be symmetric. y[0..n-1] and r[0..2*n-2] are input arrays; x[0..n-1] is the output array.*/


template<typename floatT>
void nrToeplitz(floatT * _r, floatT * _x, floatT * _y, int n)
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

   if (n == 1) FREERETURN;
   

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
         return;
      }
      
      x[m1]=sxn/sd; //whence x.

      for (j=1; j <= m; ++j) x[j] -= x[m1]*g[m-j+1];

      if (m1 == n) FREERETURN;

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
   FREERETURN;
   
}

#endif //__nrToeplitz_hpp__
