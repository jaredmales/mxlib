/** \file typeTraits.hpp
  * \brief A collection of type trait evaluations.
  * \ingroup utils_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __typeTraits_hpp__
#define __typeTraits_hpp__

namespace mx
{

namespace meta
{
   
template <typename... Ts> using void_t = void;

template <typename T, typename = void>
struct has_value_type : trueFalseT<false> {};

///Test whether a type has a typedef of "value_type"
/** Used for compile-time determination of type.
  * 
  * Example usage:
  * \code
  * bool has_vt = has_value_type<std::vector<float> >; //Evaluates to true
  * bool has_not_vt = has_value_type<float>; //Evaluates to false
  * \endcode
  * 
  * This was taken directly from the example at http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
  * 
  * \ingroup meta
  */
template <typename T>
struct has_value_type<T, void_t<typename T::value_type>> : trueFalseT<true> {};



///Check whether a type is std::vector or not.
/** First resolves whether the type as a typedef value_type.  If not, then member value will be false.
  * If it does have a member type value_type, then it uses std::is_same to compare the type to std::vector. 
  * Then value is equal to std::is_same<T, std::vector<typenameT::value_type>>::value.
  * 
  * \tparam T is the type to evaulate for vector-ness.
  * \tparam has_value_type is a boolean type based on has_value_type SFINAE test.
  * 
  * \ingroup meta
  */ 
template <typename T, bool if_value_type = has_value_type<T>::value >
struct is_std_vector; 

///Partial specialization for the case with value_type member type in T, invoking std::is_same.
/** \ingroup meta
  */
template <typename T>
struct is_std_vector<T, true> 
{
   static const bool value = std::is_same< T, std::vector< typename T::value_type> >::value;
};

///Partial specialization for the case with no value_type member type in T.
/** \ingroup meta
  */
template <typename T>
struct is_std_vector<T,false> 
{
   static const bool value = false;
};

} //namespace meta
} //namespace mx


#endif //__typeTraits_hpp__
