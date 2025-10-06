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

#include "../meta/tagT.hpp"
#include "../mxlib.hpp"

namespace mx
{
namespace ioutils
{

/** \addtogroup stringutils
 * @{
 */

/// Convert a numerical value to a string
/** The default version uses the stream library to convert.  A specialization is provided to prevent conversion of
  std::string.
  *
  * The precision parameter is only used for floating point types.  If the default precision==0 is passed, then the
  maximum useful precision is
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
template <typename typeT, unsigned width = 0, char pad = ' '>
[[deprecated( "Use std::format instead" )]]
std::string convertToString(
    const typeT &value, ///< [in] the value of type typeT to be converted
    int precision = 0   ///< [in] [optional] the precision (see
                      ///< http://www.cplusplus.com/reference/ios/ios_base/precision/) to use for floating point types.
)
{

    std::ostringstream convert;

    if( std::is_floating_point<typeT>::value )
    {
        if( precision == 0 )
        {
            convert.precision( std::numeric_limits<typeT>::max_digits10 );
        }
        else
        {
            convert.precision( precision );
        }
    }

    convert << value;

    if( width == 0 )
    {
        return convert.str();
    }
    else
    {
        std::string c = convert.str();

        if( c.length() > width )
        {
            c.erase( width, c.length() - width );
            return c;
        }

        return std::string( width - c.length(), pad ) + c;
    }
}

// Specialization of convertToString to avoid converting a string to a string
template <>
std::string convertToString<std::string>( const std::string &value, int precision );

namespace stoTImpl
{
// tag overloads for stoT

char stoT( const std::string &str, error_t *errc, meta::tagT<char> );
signed char stoT( const std::string &str, error_t *errc, meta::tagT<signed char> );
unsigned char stoT( const std::string &str, error_t *errc, meta::tagT<unsigned char> );
short stoT( const std::string &str, error_t *errc, meta::tagT<short> );
unsigned short stoT( const std::string &str, error_t *errc, meta::tagT<unsigned short> );
int stoT( const std::string &str, error_t *errc, meta::tagT<int> );
unsigned int stoT( const std::string &str, error_t *errc, meta::tagT<unsigned int> );
long stoT( const std::string &str, error_t *errc, meta::tagT<long> );
unsigned long stoT( const std::string &str, error_t *errc, meta::tagT<unsigned long> );
long long stoT( const std::string &str, error_t *errc, meta::tagT<long long> );
unsigned long long stoT( const std::string &str, error_t *errc, meta::tagT<unsigned long long> );
bool stoT( const std::string &str, error_t *errc, meta::tagT<bool> );
float stoT( const std::string &str, error_t *errc, meta::tagT<float> );
double stoT( const std::string &str, error_t *errc, meta::tagT<double> );
long double stoT( const std::string &str, error_t *errc, meta::tagT<long double> );

std::string stoT( const std::string &str, error_t *errc, meta::tagT<std::string> );

#ifdef HAS_QUAD
__float128 stoT( const std::string &str, error_t *errc, meta::tagT<__float128> );
#endif

} // namespace stoTImpl

/// Convert a string to a numerical value.
/** Provides exception-less string conversion.
 *
 * Example:
 * \code
 * std::string str = "2.34567";
 * mx::error_t errc;
 * double d;
 * d = stoT<double>(str, &errc);
 * if(errc != mx::error_t::noerror)
 * {
 * //do something
 * }
 * \endcode
 *
 * Values of `typeT=bool` are converted from strings that start with 0/f/F an 1/t/T as false and
 * true respectively.  If that fails the string is converted to long long and then to bool, so if it is
 * a valid number that fits in long long, it will evaluate to true or false based on whether or not it is 0.
 *
 * Complex types (`std::complex<realT>`) are not supported.
 *
 * \tparam typeT is the type of the numerical value desired
 *
 * \returns the converted numerical value.
 *
 * \b Error \b Codes
 * - mx::error_t::noerror on success
 * - mx::error_t::erange for an out of range value for typeT
 * - mx::error_t::invalidarg for a string that can't be converted
 */
template <typename typeT>
typeT stoT( const std::string &str, /**< [in] the string to convert*/
            error_t *errc = nullptr /**< [out] [optional] pointer to an mxlib error code set during the conversion */
)
{
    return stoTImpl::stoT( str, errc, meta::tagT<typeT>() );
}

/// Convert a string to a numerical value.
/** Provides exception-less string conversion.
 *
 * \overload
 *
 * \tparam typeT is the type of the numerical value desired
 *
 * \returns the converted numerical value.
 *
 * \b Error \b Codes
 * - mx::error_t::noerror on success
 * - mx::error_t::erange for an out of range value for typeT
 * - mx::error_t::invalidarg for a string that can't be converted
 */
template <typename typeT>
typeT stoT( const std::string &str, /**< [in] the string to convert*/
            error_t &errc           /**< [out] [optional] mxlib error code set during the conversion */
)
{
    return stoTImpl::stoT( str, &errc, meta::tagT<typeT>() );
}

/// Convert a string to a numerical value.
/** see \ref stoT
 *
 * \deprecated
 */
template <typename typeT>
[[deprecated( "Use mx::stoT<typeT> instead" )]]
typeT convertFromString( const std::string &str, /**< [in] the string to convert*/
                         error_t *errc = nullptr /**< [out] [optional] mxlib error code */ )
{
    // If no specialization exists, we try to cast
    return stoT<typeT>( str, errc );
}

/// Convert a string to all lower case.
/** Calls the c tolower function for each character in \p instr.
 *
 * \returns mx::error_t::noerror on success
 * \returns mx::error_t::std_length_error if string::resize throws a length error
 * \returns mx::error_t::std_bad_alloc if string::resize throws a bad alloc
 * \returns mx::error_t::std_exception if string::resize throws a std::exception
 * \returns mx::error_t::exception if string::resize throws any other exception
 */
mx::error_t toLower( std::string &outstr,     ///< [out]  will be resized and populated with the lower case characters
                     const std::string &instr ///< [in] is the string to convert
);

/// Convert a string to all lower case.
/** Calls the c tolower function for each character in \p instr.
 *
 * \return the all lower case string
 *
 * \b Errors
 * - mx::error_t::noerror on success
 * - mx::error_t::std_length_error if string::resize throws a length error
 * - mx::error_t::std_bad_alloc if string::resize throws a bad alloc
 * - mx::error_t::std_exception if string::resize throws a std::exception
 * - mx::error_t::exception if string::resize throws any other exception
 *
 */
std::string toLower( const std::string &instr, /**< [in] is the string to convert*/
                     mx::error_t *errc = nullptr /**< [out] [optional] the error code indicating success or error */ );

/// Convert a string to all upper case.
/** Calls the c toupper function for each character in \p instr.
 *
 * \returns mx::error_t::noerror on success
 * \returns mx::error_t::std_length_error if string::resize throws a length error
 * \returns mx::error_t::std_bad_alloc if string::resize throws a bad alloc
 * \returns mx::error_t::std_exception if string::resize throws a std::exception
 * \returns mx::error_t::exception if string::resize throws any other exception
 */
mx::error_t toUpper( std::string &outstr,     ///< [out]  will be resized and populated with the lower case characters
                     const std::string &instr ///< [in] is the string to convert
);

/// Convert a string to all upper case.
/** Calls the c toupper function for each character in \p instr.
 *
 * \returns the all upper case string
 *
 * \b Errors
 * - mx::error_t::noerror on success
 * - mx::error_t::std_length_error if string::resize throws a length error
 * - mx::error_t::std_bad_alloc if string::resize throws a bad alloc
 * - mx::error_t::std_exception if string::resize throws a std::exception
 * - mx::error_t::exception if string::resize throws any other exception
 */
std::string toUpper( const std::string &instr, /**< [in] is the string to convert*/
                     mx::error_t *errc = nullptr /**< [out] [optional] the error code indicating success or error */ );

/// Remove all white space from a string.
/**
 * Uses std::remove_if.
 */
void removeWhiteSpace( std::string &outstr,     ///< [out] will contain the new string with no whitespace.
                       const std::string &instr ///< [in] is the string to remove whitespace from
);

/// Remove all white space from a string.
/**
 * Uses std::remove_if.
 *
 * \returns the modified string.
 */
std::string removeWhiteSpace( const std::string &instr /**< [in] is the string to remove whitespace from*/ );

/// Wrap a string by breaking it into smaller sized portions of a desired width
/** Whenever possible breaks at spaces.  A single space is discarded at the break.
 */
int stringWrap( std::vector<std::string> &lines, /**< [out] each new entry contains a wrapped portion of the string.
                                                            Not cleared, so can accumulate. */
                const std::string &str,          ///< [in] the string to wrap
                int width                        ///< [in] the maximum width of the output strings
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
template <typename typeT>
void parseStringVector(
    std::vector<typeT> &v, ///< [out] the vector holding the parsed and converted tokens.  Is cleared.
    const std::string &s,  ///< [in] the string to parse
    char delim = ','       ///< [in] [optional] the delimiter.  Default is comma \p ','.
)
{
    size_t st;
    size_t com;

    st = 0;
    com = s.find( delim, st );

    v.clear();

    while( com != std::string::npos )
    {
        v.push_back( convertFromString<typeT>( s.substr( st, com - st ) ) );
        st = com + 1;
        com = s.find( delim, st );
    }
    v.push_back( convertFromString<typeT>( s.substr( st, s.size() - st ) ) );
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
template <typename typeT>
void parseStringVector(
    std::vector<typeT> &v,    ///< [out] the vector holding the parsed and converted tokens.  Is cleared.
    const std::string &s,     ///< [in] the string to parse
    const std::string &delims ///< [in] the delimiters.
)
{
    size_t st;
    size_t com;

    st = 0;
    com = s.find_first_of( delims, st );

    v.clear();

    while( com != std::string::npos )
    {
        v.push_back( convertFromString<typeT>( s.substr( st, com - st ) ) );
        st = com + 1;
        com = s.find_first_of( delims, st );
    }
    v.push_back( convertFromString<typeT>( s.substr( st, s.size() - st ) ) );
}
/// @}

} // namespace ioutils
} // namespace mx

#endif // ioutils_stringUtils_hpp
