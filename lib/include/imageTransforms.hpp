
#ifndef __imageTransforms_hpp__
#define __imageTransforms_hpp__

namespace mx
{

template<typename _arithT=double>
struct bilinearTransform
{
   typedef _arithT arithT;
   
   static const size_t width = 2;
   static const size_t lbuff = 0;
   
   template<typename arrT, typename arithT>
   void operator()(arrT & kern, arithT x, arithT y)
   {
      kern.resize(width, width);
      
      kern(0,0) = (1.-x)*(1.-y);
      kern(0,1) = (1.-x)*y;
      kern(1,0) = x*(1.-y);
      kern(1,1) = x*y;
   }
};
  
template<typename arithT>
arithT cubicConvolKernel(arithT cubic, arithT d)
{   
   if(d <= 1) return (cubic+2.)*d*d*d - (cubic+3.)*d*d + 1.;
   
   if(d < 2) return cubic*d*d*d -5.*cubic*d*d + 8.*cubic*d - 4.*cubic;
   
   return 0;
}
   
template<typename _arithT=double>
struct cubicConvolTransform
{
   typedef _arithT arithT;
   
   static const size_t width = 4;
   static const size_t lbuff = 1;

 
   arithT cubic;
   
   cubicConvolTransform()
   {
      cubic = -0.5;
   }
   
   cubicConvolTransform(arithT c)
   {
      cubic = c;
   }
   
   template<typename arrT, typename arithT>
   void operator()(arrT & kern, arithT x, arithT y)
   {
      kern.resize(width, width);
      
      arithT km2x,km1x,kp1x,kp2x;
      arithT km2y,km1y,kp1y,kp2y;
      
      km2x = cubicConvolKernel<arithT>(cubic, (1.+x));
      km1x = cubicConvolKernel<arithT>(cubic, x);
      kp1x = cubicConvolKernel<arithT>(cubic, 1.-x);
      kp2x = cubicConvolKernel<arithT>(cubic, 2.-x);
      
      km2y = cubicConvolKernel<arithT>(cubic, (1.+y));
      km1y = cubicConvolKernel<arithT>(cubic, y);
      kp1y = cubicConvolKernel<arithT>(cubic, 1.-y);
      kp2y = cubicConvolKernel<arithT>(cubic, 2.-y);
      
      
      kern(0,0) = km2x*km2y;
      kern(0,1) = km2x*km1y;
      kern(0,2) = km2x*kp1y;
      kern(0,3) = km2x*kp2y;
      
      kern(1,0) = km1x*km2y;
      kern(1,1) = km1x*km1y;
      kern(1,2) = km1x*kp1y;
      kern(1,3) = km1x*kp2y;
      
      kern(2,0) = kp1x*km2y;
      kern(2,1) = kp1x*km1y;
      kern(2,2) = kp1x*kp1y;
      kern(2,3) = kp1x*kp2y;
      
      kern(3,0) = kp2x*km2y;
      kern(3,1) = kp2x*km1y;
      kern(3,2) = kp2x*kp1y;
      kern(3,3) = kp2x*kp2y;
      
      
   }
};

template<typename arrT, typename floatT, typename transformT>
void imageRotate(arrT & transim, arrT & im, floatT dq, transformT trans)
{
   typedef typename transformT::arithT arithT;
   arithT cosq, sinq;
   arithT x0, y0, x,y;
   arithT xcen, ycen;
   
   int Nrows, Ncols;
   
   int i0, j0;

   //The kernel
   arrT kern;
   
   int lbuff = transformT::lbuff;
   int width = transformT::width;
   
   cosq = cos(dq);
   sinq = sin(dq);
   
   Nrows = im.rows();
   Ncols = im.cols();

   xcen = 0.5*(Nrows-1.);
   ycen = 0.5*(Ncols-1.);
   
   transim.resize(Nrows, Ncols);
   
   //#pragma omp parallel for private(x0,y0,i0,j0,x,y, trans, kern)
   for(int i=0;i<Nrows;i++)
   {
      for(int j=0;j<Ncols;j++)
      {
         x0 =  (i-xcen)*cosq + (j-ycen)*sinq;
         y0 = -(i-xcen)*sinq + (j-ycen)*cosq;
         
         i0 = floor(x0+xcen);
         j0 = floor(y0+ycen);
         
         x = x0+xcen-i0;
         y = y0+ycen-j0;
         
         if(i0 > lbuff && i0 < Ncols -width +lbuff && j0 > lbuff && j0 < Ncols - width+lbuff)
         {
            trans(kern, x, y);
            transim(i,j) = (im.block(i0-lbuff,j0-lbuff,width, width) * kern).sum();
            //im(i0,j0)*(1-x)*(1-y)+im(i0+1,j0)*x*(1-y)+im(i0,j0+1)*(1-x)*y + im(i0+1,j0+1)*x*y;
         }
         else transim(i,j) = 0;
      }
   }
         
}


} //namespace mx




#endif //__imageTransforms_hpp__

