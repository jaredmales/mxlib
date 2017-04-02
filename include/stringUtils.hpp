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

#include <limits>
/// The mxlib c++ namespace
namespace mx
{

/** \addtogroup stringutils
  * @{
  */

///Convert a numerical value to a string
/** The default version uses the stream library to convert.  A specialization is provided to prevent conversion of std::string.
  *
  * The precision parameter is only used for floating point types.  If the default precision==0 is passed, then the maximum useful precision is 
  * used, from the value of std::numeric_limits\<typeT\>::max_digits10.
  * 
  * The width and pad template parameters can be used to set a fixed maximum width and a pad character.
  * 
  * Examples:
  * 
  * To convert a double precision number to string:
  * \code
    double d = 2.898434;
    std::string str;
    str = convertToString(d);  //note that you will not normally need to specify <typeT>
  * \endcode
  *
  * To produce a fixed with 0 padded string from an integer:
  * \code
    int i = 23;
    std::string str;
    str = convertToString<int, 5, '0'>(i);  //result is "00023".
  * \endcode
  *
  * \tparam typeT is the type of value to convert
  * \tparam width specifies the maximum width, not including the '\0'
  * \tparam pad specifies the pad character
  *
  * \returns a string representation of value
  *
  */
template< typename typeT, unsigned width=0, char pad=' '> 
std::string convertToString( const typeT & value, ///< [in] the value of type typeT to be converted
                             int precision = 0   ///< [in] [optional] the precision (see http://www.cplusplus.com/reference/ios/ios_base/precision/) to use for floating point types.
                           )
{
   
   std::ostringstream convert;
   
   if( std::is_floating_point<typeT>::value )
   {
      if(precision == 0 )
      {
         convert.precision(  std::numeric_limits<typeT>::max_digits10);
      }
      else
      {
         convert.precision(precision);
      }
   }
   
   convert << value;
      
   if(width == 0)
   {
      return convert.str();
   }
   else
   {
      std::string c = convert.str();
      
      if( c.length() > width )
      {
         c.erase( width, c.length()-width);
         return c;
      }
      
      return std::string( width-c.length(), pad) + c;
   }
      
      
      
      
}

/// Specialization of convertToString to avoid converting a string to a string
template<> inline
std::string convertToString<std::string>( const std::string & value,
                                          int precision
                                        )
{
   return value;
}



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
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>  inline
char convertFromString<char>(const std::string & str)
{
   return (char) atoi(str.c_str());
}


/// Template specialization of convertFromString for unsigned char
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>   
unsigned char convertFromString<unsigned char>(const std::string & str)
{
   return (unsigned char) atoi(str.c_str());
}


/// Template specialization of convertFromString for short
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>  inline
short convertFromString<short>(const std::string & str)
{
   return (short) atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned short
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
unsigned short convertFromString<unsigned short>(const std::string & str)
{
   return (unsigned short) atoi(str.c_str());
}

/// Template specialization of convertFromString for int
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
int convertFromString<int>(const std::string & str)
{
   return atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned int
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
unsigned int convertFromString<unsigned int>(const std::string & str)
{
   return (unsigned int) strtoul(str.c_str(),0,0);
}

/// Template specialization of convertFromString for long
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
long convertFromString<long>(const std::string & str)
{
   return strtol(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for unsigned long
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>  inline
unsigned long convertFromString<unsigned long>(const std::string & str)
{
   return strtoul(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for long long
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
long long convertFromString<long long>(const std::string & str)
{
   return strtoll(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for unsigned long long
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
unsigned long long convertFromString<unsigned long long>(const std::string & str)
{
   return strtoull(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for float
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>   
float convertFromString<float>(const std::string & str)
{
   return strtof(str.c_str(), 0);
}


/// Template specialization of convertFromString for double
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>   
double convertFromString<double>(const std::string & str)
{
   return strtod(str.c_str(), 0);
}

/// Template specialization of convertFromString for long double
/**  
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<>   
long double convertFromString<long double>(const std::string & str)
{
   return strtold(str.c_str(), 0);
}


/// Template specialization of convertFromString for bool
/** First looks for 0/1,f/t, or F/T in the first non-space character of str.
  * otherwise uses convertFromString<int>.
  * 
  * \param str is the std::string object to convert.
  * 
  * \returns the converted numerical value.
  */
template<> inline
bool convertFromString<bool>(const std::string & str)
{
   char c = str[0];
   int i=0;
   while(isspace(c) && i < str.length()) c = str[i++];
   
   if(c == '0' || c == 'f' || c == 'F') return false;
   if(c == '1' || c == 't' || c == 'T') return true;
   
   return (bool) convertFromString<int>(str);
}




/// Convert a string to all lower case.
/** Calls the c tolower function for each character in instr.
  * 
  * \param outstr will be resized and populated with the lower case characters
  * \param instr the string to convert.
  */
inline
void toLower(std::string &outstr, const std::string & instr)
{
   outstr.resize(instr.size());
   
   for(int i=0; i < instr.size(); ++i) outstr[i] = tolower(instr[i]);
   
}




/// Convert a string to all lower case.
/** Calls the c tolower function for each character in instr.
  * 
  * \param instr the string to convert.
  * 
  * \return the all lower case string
  */
inline
std::string toLower(const std::string & instr)
{
   std::string outstr;
   
   toLower(outstr, instr);
   
   return outstr;
}



/// @}

} //namespace mx

#endif //__stringUtils_hpp__
