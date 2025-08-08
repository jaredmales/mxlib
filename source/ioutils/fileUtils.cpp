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
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <filesystem>

#include "../../include/mxException.hpp"

namespace mx
{
namespace ioutils
{

bool exists( const std::string &path )
{
    return std::filesystem::exists( std::filesystem::path( path ) );
}

int createDirectories( const std::string &path )
{
    // Use the non throwing version and silently ignore EEXIST errors
    std::error_code ec;
    std::filesystem::create_directories( path, ec );
    if( ec.value() != 0 && ec.value() != EEXIST )
    {
        return -1;
    }

    return 0;
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

std::vector<std::string> getFileNames( const std::string &directory,
                                       const std::string &prefix,
                                       const std::string &substr,
                                       const std::string &extension )
{
    std::vector<std::string> vect;
    if( std::filesystem::exists( directory ) )
    {
        if( std::filesystem::is_directory( directory ) )
        {
            try
            {
                bool hasext = false;
                std::string _ext;
                if( extension.size() > 0 )
                {
                    if( extension[0] != '.' )
                    {
                        _ext = '.';
                    }

                    _ext += extension;

                    hasext = true;
                }

                bool hasprefix = false;
                if( prefix != "" )
                {
                    hasprefix = true;
                }

                bool hassub = false;
                if( substr != "" )
                {
                    hassub = true;
                }

                std::filesystem::directory_iterator it{ directory };
                auto it_end = std::filesystem::directory_iterator{};
                for( it; it != it_end; ++it )
                {
                    if( hasext )
                    {
                        if( it->path().extension() != _ext )
                        {
                            continue;
                        }
                    }

                    std::string p = it->path().filename().generic_string();

                    if( hasprefix )
                    {
                        if( p.size() < prefix.size() )
                        {
                            continue;
                        }
                        else
                        {
                            if( p.compare( 0, prefix.size(), prefix ) != 0 )
                            {
                                continue;
                            }
                        }
                    }

                    if( hassub )
                    {
                        if( p.size() < 2 )
                        {
                            continue;
                        }

                        size_t sspos;

                        try
                        {
                            sspos = p.find( substr, 1 ); // only match if not prefix
                        }
                        catch( ... )
                        {
                            std::throw_with_nested( err::exceptthrown( "fileUtils::getFileNames",
                                                                       __FILE__,
                                                                       __LINE__,
                                                                       "from std::string::find" ) );
                        }

                        if( sspos == std::string::npos )
                        {
                            continue;
                        }
                    }

                    // If here then it passed all checks
                    vect.push_back( it->path().native() );
                }

                sort( vect.begin(), vect.end() );
            }
            catch( const err::exceptthrown & e) //this is thrown from within the loop so just pass it on
            {
                throw; //don't add anything
            }
            catch( const std::filesystem::filesystem_error & e)
            {
                std::throw_with_nested( err::exceptthrown( "fileUtils::getFileNames",
                                                           __FILE__,
                                                           __LINE__,
                                                           std::string("from std::filesystem: ") + e.what() ) );
            }
            catch( const std::exception & e )
            {
                std::throw_with_nested( err::exceptthrown( "fileUtils::getFileNames",
                                                           __FILE__,
                                                           __LINE__,
                                                           std::string("(possibly from std::string) ") + e.what() ) );
            }
            catch( ... )
            {
                std::throw_with_nested( err::exceptthrown( "fileUtils::getFileNames",
                                                           __FILE__,
                                                           __LINE__,
                                                           "unknown exception" ) );
            }
        }
        else
        {
            mxThrowException( err::invalidarg, "fileUtils::getFileNames", directory + " is not a directory" );
        }
    }
    else
    {
        mxThrowException( err::notfound, "fileUtils::getFileNames", directory + " does not exist" );
    }

    return vect;
}

std::vector<std::string> getFileNames( const std::string &directory, const std::string &extension )
{
    return getFileNames( directory, "", "", extension );
}

std::vector<std::string> getFileNames( const std::string &directory )
{
    return getFileNames( directory, "", "", "" );
}

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
