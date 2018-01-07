#ifndef __crossCorrelation_hpp__
#define __crossCorrelation_hpp__

namespace mx
{

//make these _worker
//change to pre-allocate native arrays
//add x and y lag differences  
   
//template<class eigenT, class eigenTin1, class eigenTin2>
template<typename floatT, typename sizeT=size_t>
void calcDiscreteCrossCorrelation( floatT * cc, 
                                   floatT * m1, 
                                   floatT * m2,
                                   sizeT dim1,
                                   sizeT dim2,
                                   sizeT min_x_lag,
                                   sizeT max_x_lag,
                                   sizeT min_y_lag,
                                   sizeT max_y_lag)
{

   sizeT cc_dim1, cc_dim2;
//   size_t dim1 = m1.rows();
//   size_t dim2 = m2.cols();

//   if( dim1 != m2.rows() || dim2 != m2.cols())
 //  {
      //error
 //     return;
 //  }
   
//    if( min_x_lag <= 0 ) 
//    {
//       max_x_lag = dim1; //std::min(dim1,dim2); 
//       min_x_lag = -max_x_lag;
//    }
//    
//    if( min_y_lag <= 0 ) 
//    {
//       max_y_lag = dim2; //std::min(dim1,dim2); 
//       min_y_lag = -max_y_lag;
//    }
   
   //cc.resize(maxlag-minlag+1, maxlag-minlag+1);

   cc_dim1 = max_x_lag-min_x_lag + 1;
   cc_dim2 = max_y_lag-min_y_lag + 1;
   
   //pout("cc_dims ", cc_dim1, cc_dim2);
   for(sizeT i=min_x_lag; i<max_x_lag+1; ++i)
   {
      for(sizeT j=min_y_lag; j<max_y_lag+1; ++j)
      {         
         //cc(i-minlag, j-minlag) = 0;
         cc[(j-min_y_lag)*cc_dim1 + (i-min_x_lag)] = 0;
         
         //pout((j-min_y_lag)*cc_dim1 + (i-min_x_lag));
         
         for(sizeT k=0; k<dim1; ++k)
         {
            if(k + i < 0) continue;
            if(k + i >= dim1) continue;
            
            for(sizeT l=0; l < dim2; ++l)
            {
               if(l + j < 0) continue;
               if(l + j >= dim2) continue;
               cc[(j-min_y_lag)*cc_dim1 + (i-min_x_lag)] += m1[l*dim1 + k] * m2[(l+j)*dim1 +  (k + i)];
            }
         }
         
      }
   }
}

template<class eigenT, class eigenTin1, class eigenTin2, class eigenTmask> 
void xdiscreteCrossCorrelation( eigenT & cc, 
                               const eigenTin1 & m1, 
                               const eigenTin2 & m2,
                               eigenTmask & mask,
                               int minlag = 0,
                               int maxlag = 0)
{

   size_t dim1 = m1.rows();
   size_t dim2 = m2.cols();

   if( dim1 != m2.rows() || dim2 != m2.cols() || dim1 != mask.rows() || dim2 != mask.cols() )
   {
      //error
      return;
   }
   
   if( minlag <= 0 ) 
   {
      maxlag = std::min(dim1,dim2); 
      minlag = -maxlag;
   }
   
   cc.resize(maxlag-minlag+1, maxlag-minlag+1);

   for(int i=minlag; i<maxlag+1; i++)
   {
      for(int j=minlag; j<maxlag+1; j++)
      {         
         cc(i-minlag, j-minlag) = 0;

         for(int k=0; k<dim1; k++)
         {
            if(k + i < 0) continue;
            if(k + i >= dim1) continue;
            
            for(int l=0; l < dim2; l++)
            {
               if(l + j < 0) continue;
               if(l + j >= dim2) continue;
               cc(i-minlag,j-minlag) += m1(k,l)*mask(k,l) * m2(k + i, l + j)*mask(k+i, l+j);
            }
         }
      }
   }
}


template<class eigenT, class eigenTin1, class eigenTin2> 
void discreteCrossCorrelation( typename eigenT::Scalar & xlag,
                               typename eigenT::Scalar & ylag,
                               eigenT & cc, 
                               const eigenTin1 & m1, 
                               const eigenTin2 & m2,
                               int min_x_lag = 0,
                               int max_x_lag = 0,
                               int min_y_lag = 0,
                               int max_y_lag = 0)
{
   fitGaussian2D<gaussian2D_gen_fitter<float> > fitG;
   
   
   size_t dim1 = m1.rows();
   size_t dim2 = m1.cols();

   if( dim1 != m2.rows() || dim2 != m2.cols())
   {
      //error
      return;
   }
   
//    if( min_x_lag == 0 ) 
//    {
//       max_x_lag = dim1;//std::min(dim1,dim2); 
//       min_x_lag = -max_x_lag;
//    }
//    
//    if( min_y_lag <= 0 ) 
//    {
//       max_y_lag = dim2;//std::min(dim1,dim2); 
//       min_y_lag = -max_y_lag;
//    }
   
   cc.resize(max_x_lag-min_x_lag+1, max_y_lag-min_y_lag+1);
   
   calcDiscreteCrossCorrelation<float, int>((float*)cc.data(),(float*)m1.data(),(float*)m2.data(),dim1, dim2, min_x_lag,max_x_lag, min_y_lag, max_y_lag);
   
   fitG.setArray(cc.data(), cc.rows(), cc.cols());
   
   typename eigenT::Index row, col;
   typename eigenT::Scalar maxval;
   maxval = cc.maxCoeff(&row, &col);
   cc /= maxval; //Normalize to avoid crazy nans in fit, etc.
   //pout(maxval, row, col);
   
   fitG.setGuess(0., 1., row, col, 4., 4., 0.);
   fitG.fit();
   
   //fitG.dump_report();
   
   xlag = min_x_lag+fitG.get_params()[2];
   ylag = min_y_lag+fitG.get_params()[3];
}



} //namespace mx

#endif //__crossCorrelation_hpp__
