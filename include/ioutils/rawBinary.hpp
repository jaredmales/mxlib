/** \file rawBinary.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Provides functions for working with raw binary files.
 * \ingroup utils_files
 *
 */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef rawBinary_hpp
#define rawBinary_hpp

#include <cstdio>
#include <string>

#include "../mxlib.hpp"

namespace mx
{
namespace ioutils
{

/// Read an array of data from a file as raw binary.
/** Here raw binary means no formatting or metadata.
 *
 * \ingroup ioutils
 *
 * \returns 0 on success
 * \returns -1 on error
 */
template <typename T>
int readRawBinary( T *data,                    ///< [out] the data pointer
                   size_t szData,              ///< [in] number of elements of sizeof(T) to read
                   const std::string &fileName ///< [in] the file to read from
)
{
    FILE *fout;

    fout = fopen( fileName.c_str(), "rb" );
    if( fout == 0 )
    {
        throw( mx::exception( errno2error_t( errno ), "Error from fopen [" + fileName + "]" ) );
    }

    int nrd = fread( data, sizeof( T ), szData, fout );

    if( nrd != szData )
    {
        int en = errno; // get this before fclose
        fclose( fout );
        // Have to handle case where EOF reached but no error.
        if( en != 0 )
        {
            throw( mx::exception( error_t::filererr, "Error from file [" + fileName + "]" ) );
        }
        else
        {
            throw( mx::exception( error_t::filererr,
                                  "Error reading from file, did not read all elements. [" + fileName + "]" ) );
        }
    }

    int res = fclose( fout );

    if( res != 0 )
    {
        throw(mx::exception(errno2error_t(errno), "Error closing file [" + fileName + "]" );
    }

    return 0;
}

/// Write an array of data to file as raw binary.
/** Here raw binary means no formatting or metadata,
 * just the bytes pointed to by the array are written to
 * disk.
 *
 * \ingroup ioutils
 *
 * \returns 0 on success
 * \returns -1 on error
 */
template <typename T>
int writeRawBinary( const std::string &fileName, ///< [in] the file to write to
                    T *data,                     ///< [in] the data pointer
                    size_t szData )              ///< [in] number of elements of sizeof(T) to write
{
    FILE *fout;

    fout = fopen( fileName.c_str(), "wb" );
    if( fout == 0 )
    {
        throw( mx::exception( errno2error_t( errno ), "Error from fopen [" + fileName + "]" ) );
    }

    int nwr = fwrite( data, sizeof( T ), szData, fout );

    if( nwr != szData )
    {
        int en = errno; // save before close call
        fclose( fout );

        // Have to handle case where EOF reached but no error.
        if( en != 0 )
        {
            throw(mx::exception(error_t::filewerr, "Error writing to file [" + fileName + "]" );
        }
        else
        {
            throw(mx::exception(error_t::filewerr,
                              "Error writing to file, did not write all elements. [" + fileName + "]" ));
        }
    }

    int res = fclose( fout );

    if( res != 0 )
    {
        throw(mx::exception(error_t::filecerr,"Error closing file [" + fileName + "]" ));
    }

    return 0;
}

} // namespace ioutils
} // namespace mx

#endif // rawBinary_hpp
