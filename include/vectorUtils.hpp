/** \file vectorUtils.hpp
  * \author Jared R. Males
  * \brief Some utilities for std::vector
  * \ingroup utils
  *
  */

#ifndef __vectorUtils_hpp__
#define __vectorUtils_hpp__ 


#include <vector>
#include <algorithm>

namespace mx
{

///Calculate the mean of a vector.
/** 
  * \param vec [in] the vector
  *
  * \returns the mean of vec 
  *
  * \tparam vectorT the std::vector<> type of vec 
  * 
  * \ingroup utils
  */
template<typename vectorT>
typename vectorT::value_type vectorMean(const vectorT & vec)
{
   typename vectorT::value_type mean = 0;
   
   for(size_t i=0; i< vec.size(); ++i) mean += vec[i];
   
   mean/=vec.size();
   
   return mean;
}


///Calculate median of a vector in-place, altering the vector.
/** 
  * \param vec [in] the vector, is altered by std::nth_element
  *
  * \returns the median of vec 
  *
  * \tparam vectorT the std::vector<> type of vec 
  * 
  * \ingroup utils
  */
template<typename vectorT>
typename vectorT::value_type vectorMedianInPlace(vectorT & vec)
{
   typename vectorT::value_type med;
   
   int n = 0.5*vec.size();
   
   std::nth_element(vec.begin(), vec.begin()+n, vec.end());
   
   med = vec[n];
   
   if(vec.size()%2 == 0)
   {
      med = 0.5*(med + *std::max_element(vec.begin(), vec.begin()+n)); 
   }
            
   return med;
} 

///Calculate median of a vector, leaving the vector unaltered.
/** 
  * \param vec [in] the vector
  * \param work [in] and optional vector to use as workspace, use to avoid re-allocation
  * 
  * \returns the median of vec 
  *
  * \tparam vectorT the std::vector<> type of vec 
  * 
  * \ingroup utils
  */
template<typename vectorT>
typename vectorT::value_type vectorMedian(const vectorT & vec, vectorT * work =0)
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
  * \param vec [in] the vector
  * \param mean [in] the mean value with which to calculate the variance
  * 
  * \returns the variance of vec w.r.t. mean 
  *
  * \tparam vectorT the std::vector<> type of vec 
  * 
  * \ingroup utils
  */
template<typename vectorT>
typename vectorT::value_type vectorVariance(const vectorT & vec, typename vectorT::value_type & mean)
{
   typename vectorT::value_type var;
   
   var = 0;
   for(size_t i=0; i< vec.size(); ++i) var += (vec[i]-mean)*(vec[i]-mean);
   
   var /= vec.size();
   
   return var;
}

///Calculate the variance of a vector.
/** 
  * \param vec [in] the vector
  *
  * \returns the variance of vec 
  *
  * \tparam vectorT the std::vector<> type of vec 
  * 
  * \ingroup utils
  */
template<typename vectorT>
typename vectorT::value_type vectorVariance(const vectorT & vec)
{
   typename vectorT::value_type mean;
   mean = vectorMean(vec);
   
   return vectorVariance(vec, mean);
}


///Calculate the sigma-clipped mean of a vector
/** Performas sigma-clipping relative to the median, removing any values with deviation from the median > sigma.
  * Continues until either no values are removed, or maxPasses iterations.  If maxPasses == 0, then it is ignored.
  * 
  * \param vec [in] the vector (unaltered)
  * \param sigma [in] the standard deviation threshold to apply.
  * \param var [out] the resulting variance of the vector.
  * \param maxPasses [in] (optional) the maximum number of sigma-clipping passes.  
  * 
  * \returns the mean of vec 
  *
  * \tparam vectorT the std::vector<> type of vec 
  * 
  * \ingroup utils
  */
template<typename vectorT>
typename vectorT::value_type vectorSigmaMean(const vectorT & vec, 
                                             typename vectorT::value_type sigma,
                                             typename vectorT::value_type & var,
                                             int maxPasses = 0 )
{
   vectorT work;
   typename vectorT::value_type med, Vsig, dev;

   Vsig = sigma*sigma;
   
   med = vectorMedian(vec, &work);
  
   var = vectorVariance(vec, med);
   
   int nclip;
   int passes = 0;
   
   while(passes < maxPasses || maxPasses == 0)
   {
      ++passes;
      
      nclip = 0;
      for(size_t i=0; i<work.size(); ++i)
      {
         dev = pow(work[i] - med,2);
         if( dev > Vsig )
         {
            work.erase(work.begin()+i);
            --i;
            ++nclip;
         }
      }
      if(nclip == 0) break;
      var = vectorVariance(vec, med);
   }
      
   return vectorMean(work);
}


} //namespace mx



#endif // __vectorUtils_hpp__

