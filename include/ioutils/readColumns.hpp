/** \file readColumns.hpp
 * \author Jared R. Males
 * \brief A utility to read in columns from a text file.
 * \ingroup asciiutils
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

#ifndef __readColumns_hpp__
#define __readColumns_hpp__

#include <cstring>
#include <string>
#include <fstream>
#include <format>

#include "../mxlib.hpp"

#include "stringUtils.hpp"

#define MX_READCOL_MISSINGVALSTR "-99"

namespace mx
{
namespace ioutils
{

struct readColSpaceDelim
{
    static constexpr char delim = ' ';
    static constexpr char strDelim = '"';
    static constexpr char eol = '\n';
    static constexpr char comment = '#';
    static constexpr const char *missingValStr = MX_READCOL_MISSINGVALSTR;
};

struct readColCommaDelim
{
    static constexpr char delim = ',';
    static constexpr char strdelim = '"';
    static constexpr char eol = '\n';
    static constexpr char comment = '#';
    static constexpr const char *missingValStr = MX_READCOL_MISSINGVALSTR;
};

template <class delimT, class verboseT>
error_t readcol( [[maybe_unused]] const char *sin, [[maybe_unused]] int sz, [[maybe_unused]] int & colno )
{
    return error_t::noerror;
}

template <class delimT, class verboseT, typename arrT, typename... arrTs>
error_t readcol( const char *sin, int sz, int & colno, arrT &array, arrTs &...arrays )
{
    try
    {

        // static const unsigned short int nargs = sizeof...(arrTs);
        std::string str;

        int i = 0;
        int l = strlen( sin );

        if( l < 1 )
        {
            return error_t::noerror;
        }

        // Eat white space
        while( isspace( sin[i] ) && sin[i] != delimT::eol && i < l )
        {
            ++i;
        }
        sin = sin + i;
        sz = sz - i;

        // If there's nothing here, we still need to populate the vector
        if( sz <= 1 )
        {
            array.push_back( convertFromString<typename arrT::value_type>( "" ) );
            return error_t::noerror;
        }

        std::stringstream sinstr( sin );

        std::getline( sinstr, str, delimT::delim );

        // Last entry in line might contain eol
        if( str[str.size() - 1] == delimT::eol )
        {
            str.erase( str.size() - 1 );
        }

        if( str.size() == 0 )
        {
            array.push_back( convertFromString<typename arrT::value_type>( MX_READCOL_MISSINGVALSTR ) );
        }
        else
        {
            array.push_back( convertFromString<typename arrT::value_type>( str ) );
        }

        sin += ( str.size() + 1 ) * sizeof( char );
        sz -= ( str.size() + 1 ) * sizeof( char );
    }
    catch( const std::invalid_argument &e )
    {
        return internal::mxlib_error_report<verboseT>( error_t::std_invalid_argument,
                                             std::format( "processing column {}: {}", colno, e.what() ) );
    }
    catch( const std::out_of_range &e )
    {
        return internal::mxlib_error_report<verboseT>( error_t::std_out_of_range,
                                             std::format( "processing column {}: {}", colno, e.what() ) );
    }
    catch( const std::exception &e )
    {
        return internal::mxlib_error_report<verboseT>( error_t::exception,
                                             std::format( "processing column {}: {}", colno, e.what() ) );
    }
    catch( ... )
    {
        return internal::mxlib_error_report<verboseT>( error_t::exception, std::format( "processing column {}.", colno ) );
    }

    ++colno;
    return readcol<delimT, verboseT>( sin, sz, colno, arrays... );
}

/// Read in columns from a text file
/** This function opens a file containing data formatted in columns and reads in the data row by row.
 * The data are stored in std::vectors, which should not be pre-allocated (though they could be reserve()-ed).
 *
 * Example:
 * \code
 * std::vector<int> i1;
 * std::vector<float> f1;
 * std::vector<double> d1;
 *
 * readColumns("data_file.txt", i1, f1, d1);
 * \endcode
 *
 * Note that the types of the vectors do not need to be specified as template arguments.
 *
 * The format of the file can be specified with template arguments like
 * \code
 * readColumns<',', ';', '\r'>("data_file.csv", i1, f1, d1);
 * \endcode
 * which sets the delimmiter to comma, the comment character to ;, and the end-of-line to \\r.
 *
 * Columns can be skipped using mx::ioutils::skipCol.
 *
 * \tparam delim is the character separating columns,  by default this is space.
 * \tparam comment is the character starting a comment.  by default this is #
 * \tparam eol is the end of line character.  by default this is \n
 * \tparam arrTs a variadic list of array types. this is not specified by the user.
 *
 * \todo lineSize should be configurable
 *
 * \ingroup asciiutils
 */
template <class delimT = readColSpaceDelim, class verboseT = verbose::vvv, typename... arrTs>
error_t readColumns( const std::string &fname, ///< [in] is the file name to read from
                     arrTs &...arrays          /**< [out] a variadic list of std::vectors. Any number with mixed
                                                          value_type can be specified. Neither allocated nor cleared,
                                                          so repeated calls will append data.*/
)
{
    // open file
    errno = 0;
    std::ifstream fin;
    fin.open( fname );

    if( !fin.good() )
    {
        error_t errc;
        if( errno != 0 )
        {
            errc = errno2error_t( errno );
        }
        else
        {
            errc = error_t::fileoerr;
        }

        return internal::mxlib_error_report<verboseT>( errc, "Opening " + fname + " for reading" );
    }

    std::string line;

    int64_t lineno = -1;

    while( fin.good() )
    {
        ++lineno;
        try
        {
            std::getline( fin, line, delimT::eol );
        }
        catch( const std::exception &e )
        {
            return internal::mxlib_error_report<verboseT>(
                error_t::exception,
                std::format( "Reading from {} at line {}. {}.", fname, lineno, e.what() ) );
        }
        catch( ... )
        {
            return internal::mxlib_error_report<verboseT>(
                error_t::exception,
                std::format( "Reading from {} at line {}.", fname, lineno ) );
        }

        if( line.size() == 0 )
        {
            continue;
        }

        // Find start of comment and end line at that point.
        size_t i = 0;
        bool nonspace = false; // record if we find a non-space character before the comment
        while( i < line.size() && line[i] != delimT::comment )
        {
            if( !nonspace && !isspace( line[i] ) )
            {
                nonspace = true;
            }
            ++i;
        }

        // Check if line is all comment
        if( i == 0 || !nonspace )
        {
            continue;
        }

        if( i < line.size() )                           // i is > 0 if we're here
        {
            line.erase( line.begin() + i, line.end() ); // does not throw
        }

        int colno = 0;
        error_t errc = readcol<delimT, verboseT>( line.c_str(), line.size(), colno, arrays... );

        if(errc != error_t::noerror)
        {
            return internal::mxlib_error_report<verboseT>( errc, std::format("Reading from {} at line {} column {}",fname, lineno+1, colno+1) );
        }
    }

    // getline will have set fail if there was no new line on the last line.
    if( fin.bad() && !fin.fail() )
    {
        error_t errc;
        if( errno != 0 )
        {
            errc = errno2error_t( errno );
        }
        else
        {
            errc = error_t::filererr;
        }

        return internal::mxlib_error_report<verboseT>( errc, "Reading from " + fname );
    }

    fin.clear(); // Clear the fail bit which may have been set by getline
    errno = 0;
    fin.close();

    if( fin.fail() )
    {
        error_t errc;
        if( errno != 0 )
        {
            errc = errno2error_t( errno );
        }
        else
        {
            errc = error_t::filecerr;
        }

        return internal::mxlib_error_report<verboseT>( errc, "Closing" + fname );
    }

    return error_t::noerror;
}

/// A dummy class to allow mx::readColumns to skip a column(s) in a file without requiring memory allocation.
/** The alternative is to use dummy vectors, which result in excess memory allocations and deallocations.
  * Usage:
  \code
  std::vector<T> col1, col5;
  skipCol sk;
  readColumns("filename.txt", col1, sk, sk, sk, col5); //This results in only columns 1 and 5 being stored.
  \endcode
  *
  * \ingroup asciiutils
  */
struct skipCol
{
    typedef std::string value_type; ///< value_type is defined as std::string so that no conversions take place.

    template <typename T>
    void push_back( const T &arg )
    {
        return;
    }
};

} // namespace ioutils
} // namespace mx

#endif //__readColumns_hpp__
