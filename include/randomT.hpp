/** \file randomT.hpp
  * \author Jared R. Males
  * \brief Defines a random number type
  * \ingroup random
  *
  */

#ifndef __randomT__
#define __randomT__

#include <random>

namespace mx
{

///Access to a global random number engine  
/** 
  * @returns Returns a reference to the global engine.
  * 
  * \ingroup random
  */
// template<class rangenT> rangenT & random_engine()
// {
//    static rangenT engine;
//    
//    return engine;
// }

///A random number type, which functions like any other arithmetic type.  
/** 
  * Combines a random engine, and a random distribution.  Using the type conversion operator
  * randomT returns the next random deviate whenever it is referenced. The random engine is 
  * the global engine of that type accessed by random_engine().
  *
  * Example:
  * 
  * \code
  * #include <mx/randomSeed.hpp>
  * 
  * 
  * randomT<double, std::mt19937_64, std::normal_distribution<double> > norm_distd;
  * 
  * unsigned long seedval;
  * randomSeed(seedval); //Get seed from urandom
  * norm_distd.engine().seed(seedval);
  * 
  * double d1 = norm_distd; //get a normally distributed value
  * double d2 = norm_distd; //get the next normally distributed value
  * \endcode
  * 
  * \ingroup random
  */
template<class typeT, class _ranengT, class _randistT> class randomT
{

public:

   ///Constructor
   randomT(){}

   ///Typedef for the distribution type
   typedef _randistT randistT;
      
   ///Typedef for the engine type
   typedef _ranengT ranengT;

//    ///Access to the global engine
//    _ranengT & engine()
//    {
//       return random_engine<_ranengT>();
//    }
   
   ///The random distribution
   _randistT distribution;

   ///The random engine
   ranengT engine;
   
   ///The conversion operator, returns the next value in the sequence, according to the distribution.
   operator typeT()  
   {
      //return distribution(random_engine<_ranengT>());
      return distribution(engine);
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
   * <tr><td>Mean</td><td>@f$zero@f$</td></tr>
   * <tr><td>Median</td><td>@f$zero@f$</td></tr>
   * <tr><td>Mode</td><td>@f$zero@f$</td></tr>
   * <tr><td>Range</td><td>@f$[-\infty, \infty]@f$</td></tr>
   * <tr><td>Standard Deviation</td><td>@f$\frac{\sqrt{2}}{\lambda}@f$</td></tr>
   * </table>
   * 
   * This is based on the implementation of the exponential distribution in the GNU ISO C++ Library version 4.6.
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
template<typename _RealType> inline bool operator==(const mx::laplace_distribution<_RealType>& __d1,
                                                    const mx::laplace_distribution<_RealType>& __d2)
{ 
   return __d1.param() == __d2.param(); 
}

/**
  * @brief Return true if two exponential distributions have different
  *        parameters.
  */
template<typename _RealType> inline bool operator!=(const mx::laplace_distribution<_RealType>& __d1,
                                                    const mx::laplace_distribution<_RealType>& __d2)
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
   operator<<(std::basic_ostream<_CharT, _Traits>&, const mx::laplace_distribution<_RealType>&);

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
   operator>>(std::basic_istream<_CharT, _Traits>&, mx::laplace_distribution<_RealType>&);




/** \ingroup random
 * @{
 */

///Typedef for a standard normal distribution of floats
typedef randomT<float, std::mt19937_64, std::normal_distribution<float> > norm_distf;

///Typedef for an exponential distribution of floats
typedef randomT<float, std::mt19937_64, std::exponential_distribution<float> > exp_distf;

///Typedef for a Laplace or double exponential distribution of floats
typedef randomT<float, std::mt19937_64, mx::laplace_distribution<float> > lap_distf;

///Typedef for a standard uniform distribution of doubles
typedef randomT<double, std::mt19937_64, std::uniform_real_distribution<double> > uni_distd;

///Typedef for a standard normal distribution of doubles
typedef randomT<double, std::mt19937_64, std::normal_distribution<double> > norm_distd;

///Typedef for an exponential distribution of doubles
typedef randomT<double, std::mt19937_64, std::exponential_distribution<double> > exp_distd;

///Typedef for a Laplace or double exponential distribution of doubles
typedef randomT<double, std::mt19937_64, mx::laplace_distribution<double> > lap_distd;


///Typedef for a standard normal distribution of long doubles
typedef randomT<long double, std::mt19937_64, std::normal_distribution<long double> > norm_distld;

///Typedef for an exponential distribution of long doubles
typedef randomT<long double, std::mt19937_64, std::exponential_distribution<long double> > exp_distld;

///Typedef for a Laplace or long double exponential distribution of long doubles
typedef randomT<long double, std::mt19937_64, mx::laplace_distribution<long double> > lap_distld;

/// @}

} //namespace MX

#endif //__randomT__




/*5/22/06: added 64 bit mersenne twister support
*/
