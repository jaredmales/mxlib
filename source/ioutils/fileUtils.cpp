/** \file fileUtils.cpp
 * \brief Definitions of utilities for working with files
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup fileutils
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

#include "ioutils/fileUtils.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <libgen.h>
#include <cmath>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <filesystem>

#include "../../include/mxException.hpp"

namespace mx
{
namespace ioutils
{

template bool exists<verbose::d>( const std::string &, error_t &);

template bool dir_exists_is<verbose::d>( const std::string &, error_t &);

error_t createDirectories( const std::string &path )
{
    // Use the non throwing version and silently ignore EEXIST errors
    std::error_code ec;
    std::filesystem::create_directories( path, ec );
    if( ec.value() != 0 && ec.value() != EEXIST )
    {
        mx::error_t errc = errno2error_t(ec.value());
        if(errc == error_t::error)
        {
            errc = error_t::filesystem;
        }
        return errc;
    }

    return error_t::noerror;
}

std::string pathStem( const std::string &fname )
{
    std::filesystem::path p( fname );
    return p.stem().string();
}

std::string pathFilename( const std::string &fname )
{
    std::filesystem::path p( fname );
    return p.filename().string();
}

std::string parentPath( const std::string &fname )
{
    std::filesystem::path p( fname );
    return p.parent_path().string();
}




template
error_t getFileNames<verbose::d>( std::vector<std::string> &fileNames,
                                  const std::string &directory,
                                  const std::string &prefix,
                                  const std::string &substr,
                                  const std::string &extension );



std::string fileNamePrependAppend( const std::string &fname, const std::string &prepend, const std::string &append )
{
    std::string dir, base, ext;

    std::filesystem::path p = fname;
    dir = p.parent_path().string();
    base = p.stem().string();
    ext = p.extension().string();

    return dir + '/' + prepend + base + append + ext;
}

std::string fileNameAppend( const std::string &fname, const std::string &append )
{
    return fileNamePrependAppend( fname, "", append );
}

std::string fileNamePrepend( const std::string &fname, const std::string &prepend )
{
    return fileNamePrependAppend( fname, prepend, "" );
}

std::string
getSequentialFilename( const std::string &basename, const std::string &extension, const int startat, const int ndigit )
{
    // int maxdig = 1;
    // for(int j=0;j<ndigit;++j) maxdig *= 10;
    int maxdig = pow( 10, ndigit );

    char formstr[64];
    snprintf( formstr, sizeof( formstr ), "%%0%dd", ndigit );

    char digstr[64];
    int i = startat;

    std::stringstream outn;

    snprintf( digstr, sizeof( digstr ), formstr, i );

    outn << basename;
    outn << digstr;
    outn << extension;

    while( std::filesystem::exists( outn.str() ) && i < maxdig )
    {
        ++i;
        outn.str( "" );

        snprintf( digstr, sizeof( digstr ), formstr, i );

        outn << basename;
        outn << digstr;

        outn << extension;
    }

    return outn.str();
}

off_t fileSize( int fd )
{
    if( fd == -1 )
    {
        return -1;
    }

    struct stat stbuf;

    if( ( fstat( fd, &stbuf ) != 0 ) || ( !S_ISREG( stbuf.st_mode ) ) )
    {
        return -1;
    }

    return stbuf.st_size;
}

off_t fileSize( FILE *f )
{
    return fileSize( fileno( f ) );
}

} // namespace ioutils
} // namespace mx
