/** \file randomT.hpp
  * \author Jared R. Males
  * \brief Defines a random number type
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef mx_math_randomT_hpp
#define mx_math_randomT_hpp

#include <random>

#include "randomSeed.hpp"

namespace mx
{
namespace math
{

///A random number type, which functions like any other arithmetic type.  
/** 
  * Combines a random engine, and a random distribution.  Using the type conversion operator
  * randomT returns the next random deviate whenever it is referenced. 
  *
  * Example:
  * 
    \code
    //This can also be done using the alias definition mx::math::normDistT.  
    randomT<double, std::mt19937_64, std::normal_distribution<double> > norm_distd;
    norm_distd.seed(); // 
  
    double d1 = norm_distd; //get a normally distributed value
    double d2 = norm_distd; //get the next normally distributed value
   \endcode
  * 
  * \test Verify compilation and basic operation of randomT with std::distributions \ref tests_math_randomT_basic "[test doc]" 
  * \test Verify compilation and basic operation of randomT with the Laplace distribution \ref tests_math_randomT_basic_laplace "[test doc]" 
  *
  * \ingroup random
  */
template<class typeT, class _ranengT, class _randistT> class randomT
{

public:

   ///Typedef for the distribution type
   typedef _randistT randistT;
      
   ///Typedef for the engine type
   typedef _ranengT ranengT;

   ///Constructor
   /** By default this calls the seed method, which will use /dev/random to seed the generator on linux, and time(0) on other systems.
     * Set to false to suppress seeding, and/or set a seed with seed(x).
     */ 
   randomT(bool doSeed = true /**< [in] [optional] if true then the seed method is called upon construction.*/ )
   {
      if(doSeed) seed();
   }
   
   ///The random distribution
   _randistT distribution;

   ///The random engine
   ranengT engine;
   
   ///The conversion operator, returns the next value in the sequence, according to the distribution.
   operator typeT()  
   {
      return distribution(engine);
   }
   
   ///Set the seed of the random engine.
   /** Calls the engines seed member function.
     * 
     */
   void seed( typename ranengT::result_type seedval /**< [in] the argument to pass to ranengT::seed() */)
   {
      engine.seed(seedval);
   }
   
   ///Seed the random engine with a good value
   /** Calls \ref mx::randomSeed to get the value.  On linux this uses /dev/urandom.  On other sytems, this uses time(0).
    */
   void seed()
   {
      typename ranengT::result_type seedval;
      randomSeed(seedval);
      seed(seedval);
   }
   
}; 


/**
   * @brief The Laplace (double exponential) continuous distribution for random numbers.
   *
   * The formula for the exponential probability density function is
   * @f$p(x|\lambda) = \frac{\lambda}{2} e^{-\lambda x}@f$.
   *
   * <table border=1 cellpadding=10 cellspacing=0>
   * <caption align=top>Distribution Statistics</caption>
   * <tr><td>Mean</td><td>@f$0@f$</td></tr>
   * <tr><td>Median</td><td>@f$0@f$</td></tr>
   * <tr><td>Mode</td><td>@f$0@f$</td></tr>
   * <tr><td>Range</td><td>@f$[-\infty, \infty]@f$</td></tr>
   * <tr><td>Standard Deviation</td><td>@f$\frac{\sqrt{2}}{\lambda}@f$</td></tr>
   * </table>
   * 
   * This is based on the implementation of the exponential distribution in the GNU ISO C++ Library version 4.6.
   * 
   * \test Verify compilation and basic operation of randomT with the Laplace distribution \ref tests_math_randomT_basic_laplace "[test doc]" 
   * 
   * \ingroup random
   */
template<typename _RealType = double>
class laplace_distribution
{
   static_assert(std::is_floating_point<_RealType>::value, "template argument not a floating point type");

public:
   /** The type of the range of the distribution. */
   typedef _RealType result_type;

   /** Parameter type. */
   struct param_type
   {
      typedef laplace_distribution<_RealType> distribution_type;

      explicit param_type(_RealType __lambda = _RealType(1)) : _M_lambda(__lambda)
      {}

      _RealType lambda() const
      { 
         return _M_lambda;          
      }
      
      friend bool operator==(const param_type& __p1, const param_type& __p2)
      { 
         return __p1._M_lambda == __p2._M_lambda; 
      }

      private:
        _RealType _M_lambda;
   };

public:
   /**
     * @brief Constructs a Laplace (double exponential) distribution with inverse scale
     *        parameter @f$\lambda@f$.
     */
   explicit laplace_distribution(const result_type& __lambda = result_type(1)) : _M_param(__lambda){ }

   explicit laplace_distribution(const param_type& __p) : _M_param(__p)  { }

   /**
     * @brief Resets the distribution state.
     *
     * Has no effect on Laplace distributions.
     */
   void reset() { }

   /**
     * @brief Returns the inverse scale parameter of the distribution.
     */
   _RealType lambda() const
   { 
      return _M_param.lambda(); 
   }

   /**
     * @brief Returns the parameter set of the distribution.
     */
   param_type param() const
   { 
      return _M_param; 
   }

   /**
     * @brief Sets the parameter set of the distribution.
     * @param __param The new parameter set of the distribution.
     */
   void param(const param_type& __param)
   { 
      _M_param = __param;
   }

   /**
     * @brief Returns the greatest lower bound value of the distribution.
    */
   result_type min() const
   { 
      return std::numeric_limits<result_type>::min();  
   }

   /**
     * @brief Returns the least upper bound value of the distribution.
     */
   result_type max() const
   { 
      return std::numeric_limits<result_type>::max(); 
   }

   /**
     * @brief Generating functions.
     */
   template<typename _UniformRandomNumberGenerator> result_type operator()(_UniformRandomNumberGenerator& __urng)
   { 
      return this->operator()(__urng, this->param());
   }

   template<typename _UniformRandomNumberGenerator> result_type operator()
                                  (_UniformRandomNumberGenerator& __urng, const param_type& __p)
   {
      int sgnx = 1;
      
      result_type __aurng = 0.5 - std::generate_canonical<result_type, std::numeric_limits<result_type>::digits,  _UniformRandomNumberGenerator>(__urng);
      
      if(__aurng < 0) 
      {
         sgnx = -1;
         __aurng *= -1;
      }
      
      
      return 0. -__p.lambda()*sgnx*std::log(1.-2.*__aurng);
   }

private:
   param_type _M_param;
};

/**
  * @brief Return true if two exponential distributions have the same
  *        parameters.
  */
template<typename _RealType> 
bool operator== ( const laplace_distribution<_RealType>& __d1,
                  const laplace_distribution<_RealType>& __d2
                )
{ 
   return __d1.param() == __d2.param(); 
}

/**
  * @brief Return true if two exponential distributions have different
  *        parameters.
  */
template<typename _RealType> bool operator!=( const laplace_distribution<_RealType>& __d1,
                                              const laplace_distribution<_RealType>& __d2
                                            )
{ 
   return !(__d1 == __d2); 
}

/**
  * @brief Inserts a %laplace_distribution random number distribution
  * @p __x into the output stream @p __os.
  *
  * @param __os An output stream.
  * @param __x  A %laplace_distribution random number distribution.
  *
  * @returns The output stream with the state of @p __x inserted or in
  * an error state.
  */
template<typename _RealType, typename _CharT, typename _Traits> std::basic_ostream<_CharT, _Traits>&
   operator<<(std::basic_ostream<_CharT, _Traits>&, const laplace_distribution<_RealType>&);

/**
  * @brief Extracts a %laplace_distribution random number distribution
  * @p __x from the input stream @p __is.
  *
  * @param __is An input stream.
  * @param __x A %laplace_distribution random number
  *            generator engine.
  *
  * @returns The input stream with @p __x extracted or in an error state.
  */
template<typename _RealType, typename _CharT, typename _Traits> std::basic_istream<_CharT, _Traits>&
   operator>>(std::basic_istream<_CharT, _Traits>&, laplace_distribution<_RealType>&);




/** \ingroup random
 * @{
 */

///Alias for a uniform random variate
template<typename realT>
using uniDistT = randomT<realT, std::mt19937_64, std::uniform_real_distribution<realT>>;

///Alias for a standard normal random variate
template<typename realT>
using normDistT = randomT<realT, std::mt19937_64, std::normal_distribution<realT>>;


///Alias for an exponential random variate
template<typename realT>
using expDistT = randomT<realT, std::mt19937_64, std::exponential_distribution<realT>>;

///Alias for a laplace random variate
template<typename realT>
using lapDistT = randomT<realT, std::mt19937_64, laplace_distribution<realT>>;

///Alias for a poisson random variate
template<typename intT>
using poissonDistT = randomT<intT, std::mt19937_64, std::poisson_distribution<intT>>;

///Alias for a log normal variate
template<typename realT>
using lognormDistT = randomT<realT, std::mt19937_64, std::lognormal_distribution<realT>>;

/// @}

} //namespace math
} //namespace mx

#endif //mx_math_randomT_hpp




/*5/22/06: added 64 bit mersenne twister support
*/
