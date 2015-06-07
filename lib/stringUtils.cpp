/** \file stringUtils.cpp
  * \brief utilities for working with strings
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup stringutils
  *
  */



#include "stringUtils.hpp"

namespace mx
{

/// Specialization of convertToString to avoid converting a string to a string
template<> 
  std::string convertToString<std::string>(const std::string & value)
{
   return value;
}

/// Template specialization of convertFromString for char
template<>   
char convertFromString<char>(const std::string & str)
{
   return (char) atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned char
template<>   
unsigned char convertFromString<unsigned char>(const std::string & str)
{
   return (unsigned char) atoi(str.c_str());
}

/// Template specialization of convertFromString for short
template<>   
short convertFromString<short>(const std::string & str)
{
   return (short) atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned short
template<>   
unsigned short convertFromString<unsigned short>(const std::string & str)
{
   return (unsigned short) atoi(str.c_str());
}

/// Template specialization of convertFromString for int
template<>   
int convertFromString<int>(const std::string & str)
{
   return atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned int
template<>   
unsigned int convertFromString<unsigned int>(const std::string & str)
{
   return (unsigned int) strtoul(str.c_str(),0,0);
}

/// Template specialization of convertFromString for long
template<>   
long convertFromString<long>(const std::string & str)
{
   return strtol(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for unsigned long
template<>   
unsigned long convertFromString<unsigned long>(const std::string & str)
{
   return strtoul(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for long long
template<>   
long long convertFromString<long long>(const std::string & str)
{
   return strtoll(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for unsigned long long
template<>   
unsigned long long convertFromString<unsigned long long>(const std::string & str)
{
   return strtoull(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for float
template<>   
float convertFromString<float>(const std::string & str)
{
   return strtof(str.c_str(), 0);
}

/// Template specialization of convertFromString for double
template<>   
double convertFromString<double>(const std::string & str)
{
   return strtod(str.c_str(), 0);
}

/// Template specialization of convertFromString for long double
template<>   
long double convertFromString<long double>(const std::string & str)
{
   return strtold(str.c_str(), 0);
}

} //namespace mx

