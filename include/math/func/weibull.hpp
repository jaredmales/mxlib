/** \file weibull.hpp
  * \brief The Weibull distribution.
  * \ingroup gen_math_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

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

#ifndef weibull_hpp
#define weibull_hpp

///The MLE of the Weibull distribution lambda parameter.
/**
  * \tparam realT a real floating point type 
  * 
  * \returns the MLE of lambda given the data.
  * 
  * \ingroup functions
  */
template<typename realT>
realT weibull_lambda( std::vector<realT> & x, ///<[in] the data points
                      realT k, ///< [in] the shape parameter
                      realT x0 = 0 ///< [in] [optional] the location parameter
                    )
{
   int n = x.size();

   realT lambda = 0;

   for(int i=0; i< n; ++i) lambda += pow( x[i] - x0, k);

   lambda /= n;

   return pow(lambda, 1.0/k);
}

///The general shifted Weibull distribution at a point.
/** Calculates the value of the Weibull distribution at a location specified by x.
  *
  * \f[
  * 
  * f(x; x_o, \lambda, k) = 
  * \begin{cases}
  * \frac{k}{\lambda}\left( \frac{x-x_o}{\lambda}\right)^{k-1} \exp \left( -\left(\frac{x-x_o}{\lambda}\right)^k \right), & x \ge 0 \\
  * 0, & x < 0
  * \end{cases}
  * \f]
  * 
  * \tparam realT a real floating point type 
  *
  * \returns the value of the Weibull distribution at x.
  * 
  * \ingroup functions
  */  
template<typename realT>
realT weibull( realT x, ///< [in] the location at which to calculate the distribution
               realT x0,  ///< [in] the location parameter
               realT k,  ///< [in] the shape parameter
               realT lambda  ///< [in] the scale parameter
             )
{
   if(x - x0 < 0) return 0.0;
   
   realT rat = (x-x0)/lambda;
   
   return (k/lambda) * pow(rat, k-1) * exp( -pow(rat, k) );
   
}

///The Weibull distribution at a point.
/** Calculates the value of the Weibull distribution at a location specified by x.
  *
  * \f[
  * 
  * f(x; \lambda, k) = 
  * \begin{cases}
  * \frac{k}{\lambda}\left( \frac{x}{\lambda}\right)^{k-1} \exp \left( -\left(\frac{x}{\lambda}\right)^k \right), & x \ge 0 \\
  * 0, & x < 0
  * \end{cases}
  * \f]
  * 
  * \tparam realT a real floating point type 
  *
  * \returns the value of the Weibull distribution at x.
  * 
  * \overload 
  * 
  * \ingroup functions
  */  
template<typename realT>
realT weibull( realT x, ///< [in] the location at which to calculate the distribution
               realT k,  ///< [in] the shape parameter
               realT lambda  ///< [in] the scale parameter
             )
{
   return weibull(x, static_cast<realT>(0), k, lambda);
}

   
#endif //weibull_hpp
