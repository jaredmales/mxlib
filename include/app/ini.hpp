/** \file ini.hpp
 * \author Jared R. Males
 * \brief The inih  ini-style, file parser modified for mxlib.
 *
 * \ingroup mxApp_files
 *
 */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

// This file was originally ini.c
// Modifications made by jaredmales@gmail.com to support inclusion
// in a c++ library and other features:
//   #ifdef protection
//   made comment character configurable at compile time by adding INI_COMMENT_CHAR
//   copied defines from ini.h
//   copied documentation from ini.h, and updated for doxygen
//   commented out ini.h
// NOTE: original ini.h and ini.c are in the inih directory


#ifndef ini_hpp
#define ini_hpp


/* inih -- simple .INI file parser

inih is released under the New BSD license (see LICENSE.txt). Go to the project
home page for more info:

http://code.google.com/p/inih/

*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>

//#include "ini.h"

/// Parse ini file
/** Same as ini_parse(), but takes a FILE* instead of filename. This doesn't
   close the file when it's finished -- the caller must do that. */
int ini_parse_file(FILE* file,
                   int (*handler)(void*, const char*, const char*,
                                  const char*),
                   void* user
                  );

/// Parse given INI-style file. 
/** May have [section]s, name=value pairs
    (whitespace stripped), and comments.  By default comments start with ';' (semicolon)
    but this can be configured with the INI_COMMENT_CHAR macro. Section
    is "" if name=value pair parsed before any section heading. name:value
    pairs are also supported as a concession to Python's ConfigParser.

    For each name=value pair parsed, call handler function with given user
    pointer as well as section, name, and value (data only valid for duration
    of handler call). Handler should return nonzero on success, zero on error.

    Returns 0 on success, line number of first error on parse error (doesn't
    stop on first error), -1 on file open error, or -2 on memory allocation
    error (only when INI_USE_STACK is zero).
*/
int ini_parse(const char* filename,
              int (*handler)(void*, const char*, const char*, const char*),
              void* user
             );

#endif //ini_hpp
