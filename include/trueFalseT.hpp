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

template<bool trueFalse>
struct trueFalseT;

template<>
struct trueFalseT<true>
{
   typedef bool True;
   static const bool value = true;
};

template<>
struct trueFalseT<false>
{
   typedef bool False;
   static const bool value = false;
};


} //namespace mx


#endif //__trueFalseT_hpp__
