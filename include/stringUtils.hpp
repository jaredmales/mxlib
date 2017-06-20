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
#include <vector>
#include <limits>
#include <algorithm>

#ifndef __stringUtils_hpp__
#define __stringUtils_hpp__

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
  * To produce a fixed width 0 padded string from an integer:
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
  *
  * 
  * \returns the converted numerical value.
  */
template<typename typeT> 
typeT convertFromString( const std::string & str /**< [in] is the string to convert*/ )
{
   //If no specialization exists, we try to cast
   return (typeT) str;
}



/// Template specialization of convertFromString for char
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>  inline
char convertFromString<char>(const std::string & str /**< [in] is the string to convert*/ )
{
   return (char) atoi(str.c_str());
}


/// Template specialization of convertFromString for unsigned char
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>   
unsigned char convertFromString<unsigned char>(const std::string & str /**< [in] is the string to convert*/ )
{
   return (unsigned char) atoi(str.c_str());
}


/// Template specialization of convertFromString for short
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>  inline
short convertFromString<short>(const std::string & str /**< [in] is the string to convert*/ )
{
   return (short) atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned short
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<> inline
unsigned short convertFromString<unsigned short>(const std::string & str /**< [in] is the string to convert*/ )
{
   return (unsigned short) atoi(str.c_str());
}

/// Template specialization of convertFromString for int
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<> inline
int convertFromString<int>(const std::string & str /**< [in] is the string to convert*/ )
{
   return atoi(str.c_str());
}

/// Template specialization of convertFromString for unsigned int
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<> inline
unsigned int convertFromString<unsigned int>(const std::string & str /**< [in] is the string to convert*/ )
{
   return (unsigned int) strtoul(str.c_str(),0,0);
}

/// Template specialization of convertFromString for long
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<> inline
long convertFromString<long>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtol(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for unsigned long
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>  inline
unsigned long convertFromString<unsigned long>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtoul(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for long long
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<> inline
long long convertFromString<long long>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtoll(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for unsigned long long
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<> inline
unsigned long long convertFromString<unsigned long long>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtoull(str.c_str(), 0, 0);
}

/// Template specialization of convertFromString for float
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>   
float convertFromString<float>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtof(str.c_str(), 0);
}


/// Template specialization of convertFromString for double
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>   
double convertFromString<double>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtod(str.c_str(), 0);
}

/// Template specialization of convertFromString for long double
/**  
  *
  * 
  * \returns the converted numerical value.
  */
template<>   
long double convertFromString<long double>(const std::string & str /**< [in] is the string to convert*/ )
{
   return strtold(str.c_str(), 0);
}


/// Template specialization of convertFromString for bool
/** First looks for 0/1,f/t, or F/T in the first non-space character of str.
  * otherwise uses convertFromString\<int\>.
  * 
  * \returns the converted numerical value.
  */
template<> inline
bool convertFromString<bool>(const std::string & str /**< [in] is the string to convert*/ )
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
  */
inline
void toLower( std::string &outstr, ///< [out]  will be resized and populated with the lower case characters
              const std::string & instr ///< [in] is the string to convert
            )
{
   outstr.resize(instr.size());
   
   for(int i=0; i < instr.size(); ++i) outstr[i] = tolower(instr[i]);
   
}




/// Convert a string to all lower case.
/** Calls the c tolower function for each character in instr.
  * 
  * 
  * \return the all lower case string
  */
inline
std::string toLower(const std::string & instr /**< [in] is the string to convert*/ )
{
   std::string outstr;
   
   toLower(outstr, instr);
   
   return outstr;
}

/// Convert a string to all upper case.
/** Calls the c toupper function for each character in instr.
  * 
  */
inline
void toUpper( std::string &outstr, ///< [out]  will be resized and populated with the lower case characters
              const std::string & instr ///< [in] is the string to convert
            )
{
   outstr.resize(instr.size());
   
   for(int i=0; i < instr.size(); ++i) outstr[i] = toupper(instr[i]);
   
}




/// Convert a string to all upper case.
/** Calls the c toupper function for each character in instr.
  * 
  * 
  * \return the all lower case string
  */
inline
std::string toUpper(const std::string & instr /**< [in] is the string to convert*/ )
{
   std::string outstr;
   
   toUpper(outstr, instr);
   
   return outstr;
}

///Remove all white space from a string.
/**
  * Uses std::remove_if.
  */ 
inline
void removeWhiteSpace( std::string & outstr, /// [out] will contain the new string with no whitespace.
                       const std::string & instr ///< [in] is the string to remove whitespace from
                     )
{
   outstr = instr;
   
   outstr.erase(std::remove_if(outstr.begin(), outstr.end(), ::isspace), outstr.end());
}

///Remove all white space from a string.
/**
  * Uses std::remove_if.
  * 
  * \returns the modified string.
  */ 
inline
std::string removeWhiteSpace( const std::string & instr /**< [in] is the string to remove whitespace from*/ )
{
   std::string outstr;
   
   removeWhiteSpace(outstr, instr);
   
   return outstr;
}


/// Wrap a string by breaking it into smaller sized portions of a desired width
/** Whenever possible breaks at spaces.  A single space is discarded at the break.
  */
inline
int stringWrap( std::vector<std::string> & lines, ///< [out] each new entry contains a wrapped portion of the string.  Not cleared, so can accumulate.
                const std::string & str, ///< [in] the string to wrap
                int width ///< [in] the maximum width of the output strings 
              )
{
   int L = str.length();
   
   if(L == 0) return 0;
   int startPos, tmpPos, endPos;
   
   bool done = false;
   
   startPos = 0;

   while( !done )
   {
      if(startPos == L) --startPos; //Just to make sure
      
      endPos = startPos + width;
      
      if(endPos >= L)
      {
         lines.push_back( str.substr( startPos, L-startPos ));
         done = true;
         break;
      }
      
      //Backup to nearest space
      tmpPos = endPos;      
      while( !isspace(str[tmpPos]) && tmpPos >= startPos ) --tmpPos; 
   
      //If we aren't at the beginning (i.e. splitting consecutive characters) we use new end position
      if(tmpPos > startPos) endPos = tmpPos;
      
      
      lines.push_back( str.substr( startPos, endPos-startPos) );
      
            
      startPos = endPos;
      
      //Clear 1 space
      if( str[startPos] == ' ') ++startPos;
   }
   
   return 0;
}

/// @}

} //namespace mx

#endif //__stringUtils_hpp__
