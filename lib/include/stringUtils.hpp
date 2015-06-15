/** \file stringUtils.hpp
  * \brief utilities for working with strings
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup stringutils
  *
  */

#include <string>
#include <sstream>

#ifndef __stringUtils_hpp__
#define __stringUtils_hpp__

/// The mxlib c++ namespace
namespace mx
{

/** \addtogroup stringutils
  * @{
  */

///Convert a numerical value to a string
/** \details The default version uses the stream library to convert. 
  *
  * Example:
  * \code
  * double d = 2.898434;
  * std::string str;
  * str = convertToString(d);  //note that you will not normally need to specify <typeT>
  * \endcode
  *
  * \tparam typeT is the type of value to convert
  * 
  * \param value the value of type typeT to be converted
  * 
  * \returns a string representation of value
  */
template<typename typeT> 
std::string convertToString(const typeT & value)
{
   std::ostringstream convert;
   
   convert << value;
   
   return convert.str();
}

/// Specialization of convertToString to avoid converting a string to a string
template<> 
std::string convertToString<std::string>(const std::string & value);


/// Convert a string to a numerical value.
/** The default version attempts to do the conversion with a simple c style cast.  The various template specializations
  * handle conversions to the basic types.
  * 
  * Example:
  * \code
  * std::string str = "2.34567";
  * double d;
  * d = convertFromString<double>(str);
  * \endcode
  * 
  * \tparam typeT is the type of the numerical value desired
  * 
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<typename typeT> 
typeT convertFromString(const std::string & str)
{
   //If no specialization exists, we try to cast
   return (typeT) str;
}

/// Template specialization of convertFromString for char
template<>
char convertFromString<char>(const std::string & str);

/// Template specialization of convertFromString for unsigned char
template<>  
unsigned char convertFromString<unsigned char>(const std::string & str);

/// Template specialization of convertFromString for short
template<> 
short convertFromString<short>(const std::string & str);

/// Template specialization of convertFromString for unsigned short
template<>  
unsigned short convertFromString<unsigned short>(const std::string & str);

/// Template specialization of convertFromString for int
template<>  
int convertFromString<int>(const std::string & str);

/// Template specialization of convertFromString for unsigned int
template<>  
unsigned int convertFromString<unsigned int>(const std::string & str);

/// Template specialization of convertFromString for long
template<> 
long convertFromString<long>(const std::string & str);

/// Template specialization of convertFromString for unsigned long
template<>  
unsigned long convertFromString<unsigned long>(const std::string & str);

/// Template specialization of convertFromString for long long
template<>  
long long convertFromString<long long>(const std::string & str);

/// Template specialization of convertFromString for unsigned long long
template<>  
unsigned long long convertFromString<unsigned long long>(const std::string & str);

/// Template specialization of convertFromString for float
template<>  
float convertFromString<float>(const std::string & str);

/// Template specialization of convertFromString for double
template<>  
double convertFromString<double>(const std::string & str);

/// Template specialization of convertFromString for long double
template<>  
long double convertFromString<long double>(const std::string & str);


/// @}

} //namespace mx

#endif //__stringUtils_hpp__
