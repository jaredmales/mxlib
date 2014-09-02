#ifndef __crossCorrelation_hpp__
#define __crossCorrelation_hpp__

namespace mx
{

template<class eigenT, class eigenTin1, class eigenTin2> 
void discreteCrossCorrelation( eigenT & cc, 
                               eigenTin1 & m1, 
                               eigenTin2 & m2,
                               int maxlag = 0)
{

   size_t dim1 = m1.rows();
   size_t dim2 = m2.cols();

   if( dim1 != m2.rows() || dim2 != m2.cols())
   {
      //error
      return;
   }
   
   if( maxlag <= 0 ) maxlag = std::min(dim1,dim2); 
   
   cc.resize(2*maxlag+1, 2*maxlag+1);

   for(int i=(-maxlag); i<maxlag+1; i++)
   {
      for(int j=(-maxlag); j<maxlag+1; j++)
      {         
         cc(i+maxlag, j+maxlag) = 0;

         for(int k=0; k<dim1; k++)
         {
            if(k + i < 0) continue;
            if(k + i >= dim1) continue;
            for(int l=0; l < dim2; l++)
            {
               if(l + j < 0) continue;
               if(l + j >= dim2) continue;
               cc(i+maxlag,j+maxlag) += m1(k,l) * m2(k + i, l + j);
            }
         }
      }
   }
}

template<class eigenT, class eigenTin1, class eigenTin2, class eigenTmask> 
void discreteCrossCorrelation( eigenT & cc, 
                               eigenTin1 & m1, 
                               eigenTin2 & m2,
                               eigenTmask & mask,
                               int maxlag  = 0)
{

   size_t dim1 = m1.rows();
   size_t dim2 = m2.cols();

   if( dim1 != m2.rows() || dim2 != m2.cols() || dim1 != mask.rows() || dim2 != mask.cols() )
   {
      //error
      return;
   }
   
   if( maxlag <= 0 ) maxlag = std::min(dim1,dim2); 
   
   cc.resize(2*maxlag+1, 2*maxlag+1);

   for(int i=(-maxlag); i<maxlag+1; i++)
   {
      for(int j=(-maxlag); j<maxlag+1; j++)
      {         
         cc(i+maxlag, j+maxlag) = 0;

         for(int k=0; k<dim1; k++)
         {
            if(k + i < 0) continue;
            if(k + i >= dim1) continue;
            for(int l=0; l < dim2; l++)
            {
               if(l + j < 0) continue;
               if(l + j >= dim2) continue;
               cc(i+maxlag,j+maxlag) += m1(k,l)*mask(k,l) * m2(k + i, l + j)*mask(k+i, l+j);
            }
         }
      }
   }
}

} //namespace mx

#endif //__crossCorrelation_hpp__
