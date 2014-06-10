
#include "geo.h"


template<typename arithT>
std::vector<size_t> imageRegionIndices(arithT min_r, arithT max_r, arithT min_q, arithT max_q, arithT xcen, arithT ycen)
{
   min_q = DTOR(min_q);
   max_q = DTOR(max_q);

   int min_x = -max_r, max_x = max_r, min_y = -max_r, max_y = max_r;
   
   if(abs(max_q - min_q) < DTOR(90))
   {
      min_x = min_r*cos(max_q) - 1;
      max_x = max_r*cos(min_q) + 1;

      if(min_x > max_x) std::swap(min_x, max_x);
   
      min_y = min_r*sin(min_q) - 1;
      max_y = max_r*sin(max_q) + 1;

      if(min_y > max_y) std::swap(min_y, max_y);
   }

   size_t msize = DPI*(max_r*max_r - min_r*min_r) + 1.;
   
   std::vector<size_t> v;
   
   v.resize(msize);
   
   size_t npts = 0;
   
   for(size_t i = xcen+min_x; i< xcen+max_x; i++)
   {
      for(size_t j = ycen+min_y; j< ycen+max_y; j++)
      {  
         r2 = (i-xcen)*(i-xcen) + (j-ycen)*(j-ycen);
         q = atan2( j-ycen, i-xcen );
         
         if(r2 >= min_r*min_r && r2 <= max_r*max_r && 
      }
   }
   
   return v;
}