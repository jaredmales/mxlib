/** \file vectorUtils.hpp
  * \author Jared R. Males
  * \brief  Header for the std::vector utilities
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_vectorUtils_hpp
#define math_vectorUtils_hpp


#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

#include "func/gaussian.hpp"

namespace mx
{
namespace math
{

/** \ingroup vectorutils
  *@{
  */

///Fill in a vector with a regularly spaced scale
/** Fills in the vector with a 0....N-1 scale.  The spacing of the points
  * can be changed with the scale parameter, and the starting point can be
  * changed with the offset.
  *
  * Example:
   \code
   std::vector vec;

   mx::vectorScale(vec, 1000, 0.001, 0.001); // Produces a vector with values 0.001,0.002,.... 1.000
   \endcode
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


///Return the indices of the vector in sorted order, without altering the vector itself
/** Example:
  \code
  std::vector<double> x;
  // -> fill in x with values

  std::vector<size_t> idx;
  idx = mx::sortOrder(x);

  //Print x to stdout in sorted order
  for(int i=0; i< x.size(); ++i) std::cout << x[idx[i]] << "\n";
  \endcode

  * \tparam memberT is the member type of the vector.  Must have < comparison defined.
  */
template <typename memberT>
std::vector<size_t> vectorSortOrder( std::vector<memberT> const& values /**< [in] the vector to sort */)
{
    std::vector<size_t> indices(values.size());

    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort( begin(indices), end(indices), [&](size_t a, size_t b) { return values[a] < values[b]; } );

    return indices; /// \returns the indices of the vector in sorted order.
}

///Calculate the sum of a vector.
/**
  *
  * \returns the sum of vec
  *
  * \tparam vectorT the std::vector type of vec
  *
  */
template<typename valueT>
valueT vectorSum( const valueT * vec, ///< [in] the vector 
                  size_t sz           ///< [in] the size of the vector
                )
{
   valueT sum = 0;

   for(size_t i=0; i< sz; ++i) sum += vec[i];

   return sum;
}

///Calculate the sum of a vector.
/**
  *
  * \returns the sum of vec
  *
  * \tparam vectorT the std::vector type of vec
  *
  */
template<typename vectorT>
typename vectorT::value_type vectorSum(const vectorT & vec /**< [in] the vector */)
{
   typename vectorT::value_type sum = 0;

   for(size_t i=0; i< vec.size(); ++i) sum += vec[i];

   return sum;
}

///Calculate the mean of a vector.
/**
  *
  * \returns the mean of vec
  *
  * \tparam vectorT the std::vector type of vec
  *
  */
template<typename valueT>
valueT vectorMean( const valueT * vec, ///< [in] the vector 
                   size_t sz           ///< [in] the size of the vector
                 )
{
   valueT mean = 0;

   for(size_t i=0; i< sz; ++i) mean += vec[i];

   mean /= sz;

   return mean;
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
  * \tparam valueT the data type
  *
  */
template<typename valueT>
valueT vectorVariance( const valueT * vec, ///< [in] the vector
                       size_t sz,          ///< [in] the size of the vector
                       valueT mean         ///< [in] the mean value with which to calculate the variance
                     )
{
   valueT var;

   var = 0;
   for(size_t i=0; i<sz; ++i) var += pow(vec[i]-mean,2);

   var /= (sz-1);

   return var;
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
                                             const typename vectorT::value_type & mean ///< [in] the mean value with which to calculate the variance
                                           )
{
   typename vectorT::value_type var;

   var = 0;
   for(size_t i=0; i< vec.size(); ++i) var += (vec[i]-mean)*(vec[i]-mean);

   var /= (vec.size()-1);

   return var;
}

///Calculate the variance of a vector.
/**
  * \returns the variance of vec
  *
  * \tparam vectorT the std::vector type of vec
  *
  */
template<typename valueT>
valueT vectorVariance( const valueT * vec, ///< [in] the vector
                       size_t sz           ///< [in] the size of the vector
                     )
{
   valueT mean;
   mean = vectorMean(vec, sz);

   return vectorVariance(vec, sz, mean);
}

///Calculate the variance of a vector.
/**
  * \returns the variance of vec
  *
  * \tparam vectorT the std::vector type of vec
  *
  * \overload
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
typename vectorT::value_type vectorSigmaMean( const vectorT & vec,               ///< [in] the vector (unaltered)
                                              const vectorT & weights,           ///< [in] [optional] the weights (unaltered)
                                              typename vectorT::value_type sigma ///< [in] the standard deviation threshold to apply.
                                            )
{
   int maxPasses = 0;

   return vectorSigmaMean(vec, &weights, sigma, maxPasses);
}


/// Subtract a constant value from a vector
template<typename valueT, typename constT>
void vectorSub( valueT *vec,     ///< [in/out] the vector, each element will have the constant subtracted from it
                size_t sz,       ///< [in] the size of the vector
                const constT & c ///< [in] the constant to subtract from each element
              )
{
   for(size_t n=0; n< sz;++n) vec[n] -= c;
}

/// Subtract a constant value from a vector
template<typename vecT, typename constT>
void vectorSub( vecT &vec,       ///< [in/out] the vector, each element will have the constant subtracted from it
                const constT & c ///< [in] the constant to subtract from each element
              )
{
   vectorSub(vec.data(), vec.size(), c);
}

/// Subtract the mean from a vector
template<typename valueT>
void vectorMeanSub( valueT *vec, ///< [in/out] the vector, each element will have the mean subtracted from it
                    size_t sz    ///< [in] the vector size
                  )
{
   valueT m = vectorMean(vec,sz);
   vectorSub(vec,sz,m);
}

/// Subtract the mean from a vector
template<typename vecT>
void vectorMeanSub( vecT &vec /**< [in/out] the vector, each element will have the mean subtracted from it*/)
{
   vectorMeanSub( vec.data(), vec.size());
}

/// Subtract the median from a vector
template<typename vecT>
void vectorMedianSub( vecT &vec /**<  [in/out] the vector, each element will have the median subtracted from it*/)
{
   typename vecT::value_type m = vectorMedian(vec);
   vectorSub(vec,m);
}

///Smooth a vector using the mean in a window specified by its full-width
template<typename realT>
int vectorSmoothMean( realT * smVec,  ///< [out] the smoothed version of the vector.  At least as large as \p vec.
                      realT * vec,    ///< [in] the input vector, unaltered.
                      size_t vecSize, ///< [in] the size of \p vec
                      int win         ///< [in] the full-width of the smoothing window.  Should be even.  0 results in a slow memcpy.
                    )
{
   realT sum;
   int n;
   for(int i=0; i < vecSize; ++i)
   {
      int j = i - 0.5*win;
      if(j < 0) j = 0;

      sum = 0;
      n = 0;
      while( j <= i+0.5*win && j < vecSize)
      {
         sum += vec[j];
         ++n;
         ++j;
      }

      smVec[i] = sum/n;
   }

   return 0;
}

///Smooth a vector using the mean in a window specified by its full-width
/** \overload
  */
template<typename realT>
int vectorSmoothMean( std::vector<realT> &smVec, ///< [out] the smoothed version of the vector.  Will be resize()-ed.
                      std::vector<realT> & vec,  ///< [in] the input vector, unaltered.
                      int win                    ///< [in] the full-width of the smoothing window.  Should be even.  0 results in a slow memcpy.
                    )
{
   smVec.resize(vec.size());
   return vectorSmoothMean(smVec.data(), vec.data(), vec.size(), win);
}

/// Smooth a vector using the mean in windows specified by their full-widths
/** You supply a window width for each point.  This is useful for, say, logarithmically growing bin sizes in a 
  * PSD.
  * 
  * \returns 0 on success
  * \returns -1 but who are we kidding it we don't check for errors
  * 
  */ 
template<typename realT>
int vectorSmoothMean( std::vector<realT> &smVec, ///< [out] the smoothed version of the vector
                      std::vector<realT> & vec,  ///< [in] the input vector, unaltered.
                      std::vector<int> & wins,   ///< [in] the full-widths of the smoothing windows, same size as vec.
                      bool norm = false          ///< [in] if true the output will normalized to have the same integral as the input. 
                    )
{
   smVec = vec;

   realT sum;
   int n;
   for(int i=0; i < vec.size(); ++i)
   {
      int j = i - 0.5*wins[i];
      if(j < 0) j = 0;

      sum = 0;
      n = 0;
      while( j <= i+0.5*wins[i] && j < vec.size())
      {
         sum += vec[j];
         ++n;
         ++j;
      }

      smVec[i] = sum/n;

   }

   if(norm)
   {
      realT sumin = 0, sumout = 0;
      
      for(int i=0; i < vec.size(); ++i)
      {
         sumin += vec[i];
         sumout += smVec[i];
      }
      
      for(int i=0; i < vec.size(); ++i)
      {
         smVec[i] *= sumin/sumout;
      }
   }
   
   return 0;
}

///Smooth a vector using the median in a window specified by its full-width
template<typename realT>
int vectorSmoothMedian( std::vector<realT> &smVec, ///< [out] the smoothed version of the vector
                        std::vector<realT> & vec, ///< [in] the input vector, unaltered.
                        int win                   ///< [in] the full-width of the smoothing window
                      )
{
   smVec = vec;

   std::vector<realT> tvec;
   int n;
   for(int i=0; i < vec.size(); ++i)
   {
      int j = i - 0.5*win;
      if(j < 0) j = 0;

      tvec.clear();
      while( j <= i+0.5*win && j < vec.size())
      {
         tvec.push_back(vec[j]);
         ++j;
      }

      smVec[i] = vectorMedianInPlace(tvec);

   }

   return 0;
}

///Smooth a vector using the max in a window specified by its full-width
template<typename realT>
int vectorSmoothMax( std::vector<realT> &smVec, ///< [out] the smoothed version of the vector
                     std::vector<realT> & vec, ///< [in] the input vector, unaltered.
                     int win                   ///< [in] the full-width of the smoothing window
                   )
{
   smVec = vec;

   for(int i=0; i < vec.size(); ++i)
   {
      int j = i - 0.5*win;
      if(j < 0) j = 0;

      smVec[i] = vec[j];
      ++j;
      while( j <= i+0.5*win && j < vec.size())
      {
         if(vec[j] > smVec[i]) smVec[i] = vec[j];
         ++j;
      }
   }

   return 0;
}

/// Re-bin a vector by summing (or averaging) in bins of size n points.
/**
  * \returns 0 on success
  * \returns -1 on error
  *
  * \tparam vectorT is any vector-like type with resize(), size(), and the operator()[].
  */
template<typename vectorT>
int vectorRebin( vectorT & binv,           ///< [out] the re-binned vector.  will be resized.
                 const vectorT & v,        ///< [in] the vector to bin.
                 unsigned n,               ///< [in] the size of the bins, in points
                 bool binMean = false      ///< [in] [optional] flag controlling whether sums (false) or means (true) are calculated.
               )
{
   if( n==0 ) return -1;

   binv.resize( v.size()/n );

   for(size_t i=0; i< binv.size(); ++i)
   {
      binv[i] = 0;

      unsigned j;
      for(j=0; j < n; ++j)
      {
         if(i*n + j >= v.size()) break;

         binv[i] += v[ i*n + j];
      }

      if(binMean)
      {
         if(j > 0) binv[i] /= j;
      }
   }

   return 0;
}

/// Calculate and accumulate the means of a timeseries in bins of various sizes.
/** Useful mainly to calculate the variance of the mean as a function of sample size.
  * The output is a vector of vectors, where each element is a vector which contains the means in the 
  * unique bins of size corresponding to the same index in the in put binSzs vector.
  * 
  * \returns 0 on success.
  */ 
template<typename vectorT, typename binVectorT>
int vectorBinMeans( std::vector<vectorT> & means, ///< [out] the means in each distinct bin.  Not cleared, but Will be resized with new means appended.
                    binVectorT & binSzs,          ///< [in] the bin sizes in which to calculate the means
                    const vectorT & v             ///< [in] the input vector to bin .
                  )
{
   vectorT binv;

   means.resize(binSzs.size());
   
   for(size_t i=0; i< binSzs.size(); ++i)
   {
      vectorRebin(binv, v, binSzs[i], true);

      means[i].resize(means[i].size() + binv.size());
      for(size_t j=0; j< binv.size(); ++j)  means[i][ means[i].size() - binv.size() + j] = binv[j];
   }
   
   return 0;
}

/// Convolve (smooth) a vector with a Gaussian.
/**
  * \returns 0 on success.
  * \returns -1 on error.
  *
  * \tparam realT is the real floating point type of the data.
  */
template<typename realT>
int vectorGaussConvolve( std::vector<realT> & dataOut,      ///< [out] The smoothed data vector.  Resized.
                         const std::vector<realT> & dataIn, ///< [in] The data vector to smooth.
                         const std::vector<realT> & scale,  ///< [in] The scale vector used to calculate the kernel.
                         const realT fwhm,                  ///< [in] The FWHM of the Gaussian kernel, same units as scale.
                         const int winhw                    ///< [in] The window half-width in pixels.
                       )
{

   if(dataIn.size() != scale.size()) return -1;

   dataOut.resize(dataIn.size());

   realT sigma = func::fwhm2sigma<realT>(fwhm);

   for(int i=0; i< dataIn.size(); ++i)
   {
      realT G;
      realT sum = 0;
      realT norm = 0;
      for(int j=i-winhw; j<i+winhw; ++j)
      {
         if(j < 0) continue;
         if(j > dataIn.size()-1) continue;

         G = func::gaussian<realT>( scale[j], 0.0, 1.0, scale[i], sigma);

         sum += G*dataIn[j];
         norm += G;
      }

      dataOut[i] = sum/norm;
   }

   return 0;
}




/// Calculate a cumulative histogram of a vector.
/** Sorts the vector and sums.
  *
  * \retval 0 on success, -1 otherwise.
  *
  * \tparam floatT the floating point type of the vector contens.
  */
template<typename floatT>
int vectorCumHist( std::vector<floatT> & svec, ///< [out] Contains the sorted vector.
             std::vector<floatT> & sum,  ///< [out] Contains the cumulative or running sum of the sorted vector
             std::vector<floatT> & vec   ///< [in]  The vector to sort and sum.
           )
{
   svec = vec;

   std::sort(svec.begin(), svec.end());

   sum.resize(svec.size());

   sum[0] = svec[0];

   for(int i=1; i< svec.size(); ++i)
   {
      sum[i] = sum[i-1] + svec[i];
   }

   return 0;
}

/// Calculate a reverse cumulative histogram of a vector.
/** Reverse-sorts the vector and sums.
  *
  * \retval 0 on success, -1 otherwise.
  *
  * \tparam floatT the floating point type of the vector contens.
  */
template<typename floatT>
int vectorCumHistReverse( std::vector<floatT> & svec, ///< [out] Contains the reverse-sorted vector.
                    std::vector<floatT> & sum,  ///< [out] Contains the cumulative or running sum of the reverse-sorted vector
                    std::vector<floatT> & vec   ///< [in]  The vector to reverse-sort and sum.
                  )
{
   svec = vec;

   std::sort(svec.begin(), svec.end(), std::greater<floatT>());

   sum.resize(svec.size());

   sum[0] = svec[0];

   for(int i=1; i< svec.size(); ++i)
   {
      sum[i] = sum[i-1] + svec[i];

   }

   return 0;
}

///@}

} //namespace math
} //namespace mx



#endif // math_vectorUtils_hpp
