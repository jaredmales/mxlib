/** \file vectorUtils.hpp
  * \author Jared R. Males
  * \brief  Header for the std::vector utilities
  * \ingroup gen_math_files
  *
  */

#ifndef __vectorUtils_hpp__
#define __vectorUtils_hpp__ 


#include <vector>
#include <algorithm>

namespace mx
{

/** \ingroup vectorutils
  *@{
  */

///Fill in a vector with a regularly spaced scale
/** Fills in the vector with a 0....N-1 scale.  The spacing of the points
  * can be changed with the scale parameter, and the starting point can be 
  * changed with the offset.
  *
  * \tparam vectorT is a std::vector type.
  */
template<typename vectorT>
void vectorScale( vectorT & vec,  ///< [out] the vector to fill in, can be pre-allocated or not
                  size_t N = 0,  ///< [in] [optional] if specified > 0, then vec is resize()-ed. Default is 0, and vec is not resize()-ed. 
                  typename vectorT::value_type scale = 0, ///< [in] [optional] if specified !=0, then the points are spaced by this value.  Default spacing is 1.
                  typename vectorT::value_type offset = 0 ///< [in] [optional] if specified !=0, then the starting point of the scale is this value.
                )
{
   if(scale == 0) scale = 1.0;
   
   if(N > 0) vec.resize(N);
   
   for(int i=0;i<vec.size(); ++i) vec[i] = i*scale + offset;
}
   
   
///Calculate the mean of a vector.
/** 
  *
  * \returns the mean of vec 
  *
  * \tparam vectorT the std::vector type of vec 
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorMean(const vectorT & vec /**< [in] the vector */)
{
   typename vectorT::value_type mean = 0;
   
   for(size_t i=0; i< vec.size(); ++i) mean += vec[i];
   
   mean/=vec.size();
   
   return mean;
}

///Calculate the weighted mean of a vector.
/** 
  * \returns the weighted mean of vec 
  *
  * \tparam vectorT the std::vector type of vec and w
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorMean( const vectorT & vec, ///< [in] the vector */
                                         const vectorT & w ///< [in] the weights
                                       )
{
   typename vectorT::value_type mean = 0, wsum=0;
   
   for(size_t i=0; i< vec.size(); ++i) 
   {
      mean += w[i]*vec[i];
      wsum += w[i];
   }
   
   mean /= wsum;
   
   return mean;
}

///Calculate median of a vector in-place, altering the vector.
/** Returns the center element if vec has an odd number of elements.  Returns the mean of the center 2 elements if vec has
  * an even number of elements.
  * 
  * \returns the median of vec 
  *
  * \tparam vectorT the std::vector type of vec 
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorMedianInPlace(vectorT & vec /**< [in] the vector, is altered by std::nth_element*/)
{
   typename vectorT::value_type med;
   
   int n = 0.5*vec.size();
   
   std::nth_element(vec.begin(), vec.begin()+n, vec.end());
   
   med = vec[n];
   
   //Average two points if even number of points
   if(vec.size()%2 == 0)
   {
      med = 0.5*(med + *std::max_element(vec.begin(), vec.begin()+n)); 
   }
            
   return med;
} 

///Calculate median of a vector, leaving the vector unaltered.
/** Returns the center element if vec has an odd number of elements.  Returns the mean of the center 2 elements if vec has
  * an even number of elements.
  * 
  * \returns the median of vec 
  *
  * \tparam vectorT the std::vector type of vec 
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorMedian( const vectorT & vec, ///< [in] the vector for which the median is desired
                                           vectorT * work =0 ///< [in] [optional] an optional vector to use as workspace, use to avoid re-allocation
                                         )
{
   typename vectorT::value_type med;
   
   bool localWork = false;
   if(work == 0) 
   {
      work = new vectorT;
      localWork = true;
   }
   
   work->resize(vec.size());
   
   for(int i=0;i<vec.size();++i)
   {
      (*work)[i] = vec[i];
   }

   med = vectorMedianInPlace(*work);
         
   if(localWork) delete work;
   
   return med;
} 

///Calculate the variance of a vector relative to a supplied mean value.
/** 
  * \returns the variance of vec w.r.t. mean 
  *
  * \tparam vectorT the std::vector type of vec 
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorVariance( const vectorT & vec,  ///< [in] the vector
                                             typename vectorT::value_type & mean ///< [in] the mean value with which to calculate the variance
                                           )
{
   typename vectorT::value_type var;
   
   var = 0;
   for(size_t i=0; i< vec.size(); ++i) var += (vec[i]-mean)*(vec[i]-mean);
   
   var /= vec.size();
   
   return var;
}

///Calculate the variance of a vector.
/** 
  * \returns the variance of vec 
  *
  * \tparam vectorT the std::vector type of vec 
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorVariance( const vectorT & vec /**< [in] the vector */)
{
   typename vectorT::value_type mean;
   mean = vectorMean(vec);
   
   return vectorVariance(vec, mean);
}


///Calculate the sigma-clipped mean of a vector
/** Performas sigma-clipping relative to the median, removing any values with deviation from the median > sigma.
  * Continues until either no values are removed, or maxPasses iterations.  If maxPasses == 0, then it is ignored.
  * 
  * \returns the sigma clipped mean of vec 
  *
  * \tparam vectorT the std::vector type of vec 
  * 
  */
template<typename vectorT>
typename vectorT::value_type vectorSigmaMean( const vectorT & vec,  ///<  [in] the vector (unaltered)
                                              const vectorT * weights,  ///<  [in] [optional] the weights (unaltered)
                                              typename vectorT::value_type sigma, ///< [in] the standard deviation threshold to apply.
                                              int & maxPasses  ///< [in/out] [optional] the maximum number of sigma-clipping passes.  Set to actual number of passes on return.
                                            )
{
   vectorT work, wwork;
   
   typename vectorT::value_type med, var, Vsig, dev;

   
   bool doWeight = false;
   if( weights )
   {
      if( weights->size() == vec.size()) doWeight = true;
   }
   
   Vsig = sigma*sigma;
   
   med = vectorMedian(vec, &work);
   var = vectorVariance( work, med );
  
   int nclip;
   int passes = 0;
  
   
   //If weighting, have to synchronize work with weights since work will be
   //partially sorted by median.
   if(doWeight)
   {
      wwork.resize(vec.size());
      for(int i=0; i< vec.size(); ++i)
      {
         work[i] = vec[i];
         wwork[i] = (*weights)[i];
      }
   }
      
   
   while(passes < maxPasses || maxPasses == 0)
   {
      ++passes;
      
      nclip = 0;
      
      for(size_t i=0; i<work.size(); ++i)
      {
         dev = pow(work[i] - med,2)/var;
         if( dev > Vsig )
         {
            work.erase(work.begin()+i);
            
            if(doWeight) wwork.erase(wwork.begin() + i);
            
            --i;
            ++nclip;
         }
      }
      
      if(nclip == 0) break;
      med = vectorMedian(work);
      var = vectorVariance( work, med );
   }
      
   maxPasses = passes;
   
   if(doWeight)
   {
      return vectorMean(work, wwork);
   }
   else
   {
      return vectorMean(work);
   }
}


/** 
  * \overload
  */
template<typename vectorT>
typename vectorT::value_type vectorSigmaMean( const vectorT & vec,  ///<  [in] the vector (unaltered)
                                              typename vectorT::value_type sigma ///< [in] the standard deviation threshold to apply.
                                            )
{
   int maxPasses = 0;
   
   return vectorSigmaMean(vec, (vectorT *) 0, sigma, maxPasses);
}

/** 
  * \overload
  */
template<typename vectorT>
typename vectorT::value_type vectorSigmaMean( const vectorT & vec,  ///<  [in] the vector (unaltered)
                                              const vectorT & weights, ///<  [in] [optional] the weights (unaltered)
                                              typename vectorT::value_type sigma, ///< [in] the standard deviation threshold to apply.
                                              int & maxPasses  ///< [in/out] [optional] the maximum number of sigma-clipping passes.  Set to actual number of passes on return.
                                            )
{
   
   return vectorSigmaMean(vec, &weights, sigma, maxPasses);
}

/** 
  * \overload
  */
template<typename vectorT>
typename vectorT::value_type vectorSigmaMean( const vectorT & vec,  ///<  [in] the vector (unaltered)
                                              const vectorT & weights, ///<  [in] [optional] the weights (unaltered)
                                              typename vectorT::value_type sigma ///< [in] the standard deviation threshold to apply.
                                            )
{
   int maxPasses = 0;
   
   return vectorSigmaMean(vec, &weights, sigma, maxPasses);
}




///@}

} //namespace mx



#endif // __vectorUtils_hpp__

