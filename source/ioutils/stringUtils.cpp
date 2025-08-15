/** \file stringUtils.cpp
 * \brief Implementation of utilities for working with strings
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup stringutils
 *
 */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "../../include/ioutils/stringUtils.hpp"

namespace mx
{
namespace ioutils
{

// Specialization of convertToString to avoid converting a string to a string
template <>
std::string convertToString<std::string>( const std::string &value, int precision )
{
    static_cast<void>( precision );
    return value;
}

namespace stoTImpl
{

template <typename typeT>
typeT stoInt32s( const std::string &str, error_t *errc )
{
    error_t _errc;

    long val = stoT( str, &_errc, long() );

    if( _errc == error_t::noerror )
    {
        if( val < std::numeric_limits<typeT>::lowest() )
        {
            _errc = error_t::erange;
            val = std::numeric_limits<typeT>::lowest();
        }
        else if( val > std::numeric_limits<typeT>::max() )
        {
            _errc = error_t::erange;
            val = std::numeric_limits<typeT>::max();
        }
    }

    if( errc )
    {
        *errc = _errc;
    }

    return static_cast<typeT>( val );
}

char stoT( const std::string &str, error_t *errc, const char & )
{
    return stoInt32s<char>( str, errc );
}

unsigned char stoT( const std::string &str, error_t *errc, const unsigned char & )
{
    return stoInt32s<unsigned char>( str, errc );
}

short stoT( const std::string &str, error_t *errc, const short & )
{
    return stoInt32s<short>( str, errc );
}

unsigned short stoT( const std::string &str, error_t *errc, const unsigned short & )
{
    return stoInt32s<unsigned short>( str, errc );
}

int stoT( const std::string &str, error_t *errc, const int & )
{
    return stoInt32s<int>( str, errc );
}

unsigned int stoT( const std::string &str, error_t *errc, const unsigned int & )
{
    error_t _errc;

    typedef unsigned long ulong;
    unsigned long val = stoT( str, &_errc, ulong() );

    if( _errc == error_t::noerror )
    {
        if( val > std::numeric_limits<unsigned int>::max() )
        {
            _errc = error_t::erange;
            val = std::numeric_limits<unsigned int>::max();
        }
    }

    if( errc )
    {
        *errc = _errc;
    }

    return static_cast<unsigned int>( val );
}

long stoT( const std::string &str, error_t *errc, const long & )
{
    error_t _errc;

    char *end;

    errno = 0;

    long val = std::strtol( str.c_str(), &end, 10 );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

unsigned long stoT( const std::string &str, error_t *errc, const unsigned long & )
{
    error_t _errc;

    char *end;

    errno = 0;

    unsigned long val = std::strtoul( str.c_str(), &end, 10 );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

long long stoT( const std::string &str, error_t *errc, const long long & )
{
    error_t _errc;

    char *end;

    errno = 0;

    long long val = std::strtoll( str.c_str(), &end, 10 );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

unsigned long long stoT( const std::string &str, error_t *errc, const unsigned long long & )
{
    error_t _errc;

    char *end;

    errno = 0;

    unsigned long long val = std::strtoull( str.c_str(), &end, 10 );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

/* bool specialization
 * First looks for 0/1, f/t, or F/T in the first non-space character of str.
 * Otherwise, uses long long and casts to bool.
 *
 * returns the converted numerical value.
 */
bool stoT( const std::string &str, error_t *errc, const bool & )
{
    char c = str[0];
    size_t i = 0;
    while( i < str.length() && isspace( c ) )
    {
        c = str[i++];
    }

    if( c == '0' || c == 'f' || c == 'F' )
    {
        return false;
    }
    else if( c == '1' || c == 't' || c == 'T' )
    {
        return true;
    }

    typedef long long llong;

    return static_cast<bool>( stoT( str, errc, llong() ) );
}

float stoT( const std::string &str, error_t *errc, const float & )
{
    error_t _errc;

    char *end;

    errno = 0;

    float val = std::strtof( str.c_str(), &end );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

double stoT( const std::string &str, error_t *errc, const double & )
{
    error_t _errc;

    char *end;

    errno = 0;

    double val = std::strtod( str.c_str(), &end );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

long double stoT( const std::string &str, error_t *errc, const long double & )
{
    error_t _errc;

    char *end;

    errno = 0;

    long double val = std::strtold( str.c_str(), &end );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

#ifdef HAS_QUAD
__float128 stoT( const std::string &str, error_t *errc, const __float128 & )
{
    error_t _errc;

    char *end;

    errno = 0;

    __float128 val = ::strtof128( str.c_str(), &end );

    if( errno != 0 )
    {
        _errc = errno2error_t( errno );
    }
    else if( end == str.c_str() )
    {
        _errc = error_t::invalidarg;
    }
    else
    {
        _errc = error_t::noerror;
    }

    if( errc )
    {
        *errc = _errc;
    }

    return val;
}

#endif // has_quad

std::string stoT( const std::string &str, error_t *errc, const std::string & )
{
    if(errc)
    {
        *errc = mx::error_t::noerror;
    }

    return str;
}

} // namespace stoTImpl

// Convert a string to all lower case.
mx::error_t toLower( std::string &outstr, const std::string &instr )
{
    try
    {
        outstr.resize( instr.size() );
    }
    catch( const std::length_error &e )
    {
        return mx::error_t::std_length_error;
    }
    catch( const std::bad_alloc &e )
    {
        return mx::error_t::std_bad_alloc;
    }
    catch( const std::exception &e )
    {
        return mx::error_t::std_exception;
    }
    catch( ... )
    {
        return mx::error_t::exception;
    }

    for( size_t i = 0; i < instr.size(); ++i )
    {
        outstr[i] = tolower( instr[i] );
    }

    return mx::error_t::noerror;
}

// Convert a string to all lower case.
std::string toLower( const std::string &instr, mx::error_t *errc )
{
    std::string outstr;

    mx::error_t _errc = toLower( outstr, instr );

    if( errc )
    {
        *errc = _errc;
    }

    return outstr;
}

// Convert a string to all upper case.
mx::error_t toUpper( std::string &outstr, const std::string &instr )
{
    try
    {
        outstr.resize( instr.size() );
    }
    catch( const std::length_error &e )
    {
        return mx::error_t::std_length_error;
    }
    catch( const std::bad_alloc &e )
    {
        return mx::error_t::std_bad_alloc;
    }
    catch( const std::exception &e )
    {
        return mx::error_t::std_exception;
    }
    catch( ... )
    {
        return mx::error_t::exception;
    }

    for( size_t i = 0; i < instr.size(); ++i )
    {
        outstr[i] = toupper( instr[i] );
    }

    return mx::error_t::noerror;
}

// Convert a string to all upper case.
std::string toUpper( const std::string &instr, mx::error_t *errc )
{
    std::string outstr;

    mx::error_t _errc = toUpper( outstr, instr );

    if( errc )
    {
        *errc = _errc;
    }

    return outstr;
}

// Remove all white space from a string.
void removeWhiteSpace( std::string &outstr, const std::string &instr )
{
    outstr = instr;

    outstr.erase( std::remove_if( outstr.begin(), outstr.end(), ::isspace ), outstr.end() );
}

// Remove all white space from a string.
std::string removeWhiteSpace( const std::string &instr )
{
    std::string outstr;

    removeWhiteSpace( outstr, instr );

    return outstr;
}

// Wrap a string by breaking it into smaller sized portions of a desired width
int stringWrap( std::vector<std::string> &lines, const std::string &str, int width )
{
    int L = str.length();

    if( L == 0 )
        return 0;
    int startPos, tmpPos, endPos;

    bool done = false;

    startPos = 0;

    while( !done )
    {
        if( startPos == L )
            --startPos; // Just to make sure

        endPos = startPos + width;

        if( endPos >= L )
        {
            lines.push_back( str.substr( startPos, L - startPos ) );
            done = true;
        }
        else
        {
            // Backup to nearest space
            tmpPos = endPos;
            while( !isspace( str[tmpPos] ) && tmpPos >= startPos )
                --tmpPos;

            // If we aren't at the beginning (i.e. splitting consecutive characters) we use new end position
            if( tmpPos > startPos )
                endPos = tmpPos;

            lines.push_back( str.substr( startPos, endPos - startPos ) );

            startPos = endPos;

            // Clear 1 space
            if( str[startPos] == ' ' )
                ++startPos;
        }
    }

    return 0;
}

} // namespace ioutils
} // namespace mx
