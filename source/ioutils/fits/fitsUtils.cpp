/** \file fitsUtils.cpp
 * \brief Implementation of utilities to work with FITS files
 * \ingroup fits_processing_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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

#include "ioutils/fits/fitsUtils.hpp"

namespace mx
{
namespace fits
{

int fitsStripApost( std::string &s )
{
    int stripped = 0;
    int p = s.find_first_not_of( " \t\n\t" );
    if( s[p] == '\'' )
    {
        s.erase( 0, p + 1 );
        ++stripped;
    }

    p = s.find_last_not_of( " \t\n\t" );
    if( s[p] == '\'' )
    {
        s.erase( p );
        ++stripped;
    }

    --p;

    while( s[p] == ' ' && p >= 0 )
    {
        s.erase( p );
        --p;
        ++stripped;
    }

    return stripped;
}

void fitsPopulateCard( char headStr[81], char *keyword, char *value, char *comment )
{
    memset( headStr, ' ', 80 );
    headStr[80] = '\0';

    int rpos = 0;

    if( strlen( keyword ) > 8 )
    {
        rpos += snprintf( headStr, 81, "HIERARCH " );

        rpos += snprintf( headStr + rpos, 81 - rpos, "%s =", keyword );
    }
    else
    {
        rpos += snprintf( headStr, 81, "%-8s=", keyword );
    }

    if( strlen( value ) < stdValWidth )
    {
        char fstr[10];
        snprintf( fstr, 10, "%%+%ds", stdValWidth );
        rpos += snprintf( headStr + rpos, 81 - rpos, fstr, value );
    }
    else
    {
        rpos += snprintf( headStr + rpos, 81 - rpos, "%s", value );
    }

    headStr[rpos] = ' ';
    ++rpos;

    headStr[rpos] = '/';
    ++rpos;

    headStr[rpos] = ' ';
    ++rpos;

    snprintf( headStr + rpos, 81 - rpos, "%s", comment );
}

template <>
int fits_write_key<char *>( fitsfile *fptr, char *keyword, void *value, char *comment )
{
    int fstatus = 0;

    fits_write_key_longwarn( fptr, &fstatus );

    fits_write_key_longstr( fptr, keyword, (const char *)value, comment, &fstatus );

    return fstatus;
}

template <>
int fits_write_key<std::string>( fitsfile *fptr, char *keyword, void *value, char *comment )
{
    return fits_write_key<char *>( fptr, keyword, value, comment );
}

template <>
int fits_write_key<bool>( fitsfile *fptr, char *keyword, void *value, char *comment )
{
    unsigned char bc = *( (bool *)value );

    int fstatus = 0;

    fits_write_key( fptr, fitsType<unsigned char>(), keyword, &bc, comment, &fstatus );

    return fstatus;
}

template <>
int fits_write_key<fitsUnknownType>( fitsfile *fptr, char *keyword, void *value, char *comment )
{

    int fstatus = 0;

    char *vstr = (char *)value;

    // If type is unkown, do it verbatim
    char headStr[81];

    fitsPopulateCard( headStr, keyword, vstr, comment );

    fits_write_record( fptr, headStr, &fstatus );
    return fstatus;
}

int fits_write_comment( fitsfile *fptr, char *comment )
{
    int fstatus = 0;

    fits_write_comment( fptr, comment, &fstatus );

    return fstatus;
}

int fits_write_history( fitsfile *fptr, char *history )
{
    int fstatus = 0;

    fits_write_history( fptr, history, &fstatus );

    return fstatus;
}

void fitsErrText( std::string &explan, const std::string &filename, int fstatus )
{
    char emnem[31];

    fits_get_errstatus( fstatus, emnem );

    explan += ": ";
    explan += filename;
    explan += ". CFITSIO: ";

    explan += emnem;
    explan += " (";
    explan += ioutils::convertToString( fstatus );
    explan += ")";
}

} // namespace fits
} // namespace mx
