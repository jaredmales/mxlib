/** \file trueFalseT.hpp
  * \brief Declares and defines a true-false virtual type used for boolean template overrides
  * \ingroup utils_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __trueFalseT_hpp__
#define __trueFalseT_hpp__

namespace mx
{

namespace meta
{
   
///Template declaration of a trueFalseT type.
/** The specialized types can be used for tag dispatching, that is achieving SFINAE-like behavior but with function overloading when
  * identical function signatures are desired.  This works because trueFalseT\<true\> and trueFalseT\<false\> are
  * dfferent types.
  * 
  * This is a slightly simpler construct than std::true_type and std::false_type, which are derived from std::integral_constant.
  * 
  * This template type is not defined, only the specializations are. The specializations have a typedef of True or False, accordingly,
  * and a static const member value which evaluates to true or false accordingly.
  * 
  * \tparam trueFalse bool to choose which of the specializations to invoke.
  * 
  * \ingroup meta
  */   
template<bool trueFalse>
struct trueFalseT;

///The true specialization of trueFalseT.
/** \ingroup meta
  */
template<>
struct trueFalseT<true>
{
   typedef bool True; ///< Typedef which can be used for SFINAE
   static const bool value = true; ///< bool member value = true
};

///The false specialization of trueFalseT.
/** \ingroup meta
  */
template<>
struct trueFalseT<false>
{
   typedef bool False; ///< Typedef which can be used for SFINAE
   static const bool value = false; ///< bool member value = false
};

} //namespace meta
} //namespace mx


#endif //__trueFalseT_hpp__
