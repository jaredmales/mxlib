/** \file stringUtils.hpp
  * \brief Utilities for working with strings
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup stringutils
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

#ifndef ioutils_stringUtils_hpp
#define ioutils_stringUtils_hpp


#include <string>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>

#include "../mxlib.hpp"

namespace mx
{
namespace ioutils
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
                             int precision = 0    ///< [in] [optional] the precision (see http://www.cplusplus.com/reference/ios/ios_base/precision/) to use for floating point types.
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

// Specialization of convertToString to avoid converting a string to a string
template<> 
std::string convertToString<std::string>( const std::string & value,
                                          int precision
                                        );


/// Convert a string to a numerical value.
/** The default version attempts to do the conversion with a simple c style cast.  Template specializations
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
  * \returns the converted numerical value.
  */
template<typename typeT>
typeT convertFromString( const std::string & str /**< [in] is the string to convert*/ )
{
   //If no specialization exists, we try to cast
   return (typeT) str;
}

template<> 
char convertFromString<char>(const std::string & str /* [in] is the string to convert*/ );

template<> 
char16_t convertFromString<char16_t>(const std::string & str /* [in] is the string to convert*/ );

template<> 
char32_t convertFromString<char32_t>(const std::string & str /* [in] is the string to convert*/ );

template<> 
wchar_t convertFromString<wchar_t>(const std::string & str /* [in] is the string to convert*/ );

template<>
signed char convertFromString<signed char>(const std::string & str /* [in] is the string to convert*/ );

template<>
unsigned char convertFromString<unsigned char>(const std::string & str /* [in] is the string to convert*/ );

template<> 
short convertFromString<short>(const std::string & str /* [in] is the string to convert*/ );

template<> 
unsigned short convertFromString<unsigned short>(const std::string & str /* [in] is the string to convert*/ );

template<>
int convertFromString<int>(const std::string & str /* [in] is the string to convert*/ );

template<>
unsigned int convertFromString<unsigned int>(const std::string & str /* [in] is the string to convert*/ );

template<>
long convertFromString<long>(const std::string & str /* [in] is the string to convert*/ );

template<>
unsigned long convertFromString<unsigned long>(const std::string & str /* [in] is the string to convert*/ );

template<>
long long convertFromString<long long>(const std::string & str /* [in] is the string to convert*/ );

template<>
unsigned long long convertFromString<unsigned long long>(const std::string & str /* [in] is the string to convert*/ );

template<>
float convertFromString<float>(const std::string & str /* [in] is the string to convert*/ );

template<>
double convertFromString<double>(const std::string & str /* [in] is the string to convert*/ );

template<>
long double convertFromString<long double>(const std::string & str /* [in] is the string to convert*/ );


/// Template specialization of convertFromString for bool
/** First looks for 0/1, f/t, or F/T in the first non-space character of str.
  * Otherwise, we use convertFromString\<int\>.
  *
  * \returns the converted numerical value.
  */
template<>
bool convertFromString<bool>(const std::string & str /**< [in] is the string to convert*/ );


/// Convert a string to all lower case.
/** Calls the c tolower function for each character in instr.
  *
  */
void toLower( std::string &outstr, ///< [out]  will be resized and populated with the lower case characters
              const std::string & instr ///< [in] is the string to convert
            );

/// Convert a string to all lower case.
/** Calls the c tolower function for each character in instr.
  *
  *
  * \return the all lower case string
  */
std::string toLower(const std::string & instr /**< [in] is the string to convert*/ );

/// Convert a string to all upper case.
/** Calls the c toupper function for each character in instr.
  *
  */
void toUpper( std::string &outstr, ///< [out]  will be resized and populated with the lower case characters
              const std::string & instr ///< [in] is the string to convert
            );

/// Convert a string to all upper case.
/** Calls the c toupper function for each character in instr.
  *
  *
  * \return the all lower case string
  */
std::string toUpper(const std::string & instr /**< [in] is the string to convert*/ );

///Remove all white space from a string.
/**
  * Uses std::remove_if.
  */
void removeWhiteSpace( std::string & outstr, ///< [out] will contain the new string with no whitespace.
                       const std::string & instr ///< [in] is the string to remove whitespace from
                     );

///Remove all white space from a string.
/**
  * Uses std::remove_if.
  *
  * \returns the modified string.
  */
std::string removeWhiteSpace( const std::string & instr /**< [in] is the string to remove whitespace from*/ );

/// Wrap a string by breaking it into smaller sized portions of a desired width
/** Whenever possible breaks at spaces.  A single space is discarded at the break.
  */
int stringWrap( std::vector<std::string> & lines, ///< [out] each new entry contains a wrapped portion of the string.  Not cleared, so can accumulate.
                const std::string & str, ///< [in] the string to wrap
                int width ///< [in] the maximum width of the output strings
              );

/// Parses a string into a vector of tokens delimited by a character
/** E.g., the string 
  * \code
    std::string s={"0,1,2,3,4"};
    std::vector<int> v;
    parseStringVector(v,s);
    \endcode
  * is parsed to a vector as if it was initialized with 
    \code
    std::vector<int> v = {0,1,2,3,4};
    \endcode 
  *
  * \tparam typeT the type to convert the tokens too.
  */
template<typename typeT>
void parseStringVector( std::vector<typeT> & v, ///< [out] the vector holding the parsed and converted tokens.  Is cleared.
                        const std::string & s,  ///< [in] the string to parse
                        char delim = ','        ///< [in] [optional] the delimiter.  Default is comma \p ','.
                      )
{
   size_t st;
   size_t com;

   st = 0;
   com = s.find(delim, st);

   v.clear();

   while(com != std::string::npos)
   {
      v.push_back( convertFromString<typeT>(s.substr(st, com-st)) );
      st = com + 1;
      com = s.find(delim, st);
   }
   v.push_back( convertFromString<typeT>(s.substr(st, s.size()-st)));

}

/// Parses a string into a vector of tokens delimited by a set of characters
/** E.g., the string 
  * \code
    std::string s={"0,1:2 3,4"};
    std::vector<int> v;
    parseStringVector(v, s, ",: ");
    \endcode
  * is parsed to a vector as if it was initialized with 
    \code
    std::vector<int> v = {0,1,2,3,4};
    \endcode 
  *
  * \tparam typeT the type to convert the tokens too.
  */
template<typename typeT>
void parseStringVector( std::vector<typeT> & v,    ///< [out] the vector holding the parsed and converted tokens.  Is cleared.
                        const std::string & s,     ///< [in] the string to parse
                        const std::string & delims ///< [in] the delimiters. 
                      )
{
   size_t st;
   size_t com;

   st = 0;
   com = s.find_first_of(delims, st);

   v.clear();

   while(com != std::string::npos)
   {
      v.push_back( convertFromString<typeT>(s.substr(st, com-st)) );
      st = com + 1;
      com = s.find_first_of(delims, st);
   }
   v.push_back( convertFromString<typeT>(s.substr(st, s.size()-st)));

}
/// @}

} //namespace ioutils
} //namespace mx

#endif //ioutils_stringUtils_hpp
