
/** \file fileUtils.hpp
 * \brief Declarations of utilities for working with files
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup fileutils
 *
 */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef ioutils_fileUtils_hpp
#define ioutils_fileUtils_hpp

#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>

#include "../mxlib.hpp"

namespace mx
{
namespace ioutils
{

#ifndef ioutils_fileUtils_detail_hpp
#define ioutils_fileUtils_detail_hpp
namespace fileUtilsDetail
{

/** \cond */
/// File-utility stages exposed to deterministic exception tests.
enum class operation
{
    stringToPath,
    getFileNames,
    createDirectories
};

/// Signature of the resettable file-utility operation hook.
using operationHookT = void ( * )( operation );

/// Access the process-wide file-utility operation hook.
operationHookT &operationHook();

/// Signature of the resettable file-utility error-code hook.
using errorCodeHookT = void ( * )( operation, std::error_code & );

/// Access the process-wide file-utility error-code hook.
errorCodeHookT &errorCodeHook();
/** \endcond */

} // namespace fileUtilsDetail
#endif // ioutils_fileUtils_detail_hpp

#ifdef MXLIBTEST_NAMESPACE
namespace MXLIBTEST_NAMESPACE
{
#endif

/** \addtogroup fileutils
 * @{
 */

/// Convert a string to a path, handling exceptions.
/** Wrapper for `path = str` assignment that handles exceptions
 * and implements mxlib standard error handling.
 *
 * \returns error_t::noerror if no exceptions
 * \returns error_t::std_bad_alloc if std::bad_alloc is caught
 * \returns error_t::std_filesystem_error if std::filesystem::filesystem_error is caught
 * \returns error_t::std_exception if any other exceptions are caught
 *
 * \throws a nested mx::exception for any uncaught exceptions.
 */
template <class verboseT>
error_t string2path( std::filesystem::path &path, const std::string &str );

/// Check if a path exists
/**
 * \returns true if the path exists and no errors occur
 * \returns false otherwise
 */
template <class verboseT = verbose::d>
bool exists( const std::string &path, /**< [in] the path to check for existence */
             mx::error_t &errc        /**< [out] error code. Typically convereted as errno from std::filesystem*/
);

/// Check if a path exists and is a directory
/**
 * \returns true only if \p dir both exists and is a directory, and no errors occur
 * \returns false otherwise
 */
template <class verboseT = verbose::d>
bool dir_exists_is( const std::string &dir, /**< [in] the path to check */
                    mx::error_t &errc       /**< [out] error code. Typically convereted as errno from std::filesystem*/
);

/// Create a directory or directories
/** This will create any directories in path that don't exist.  It silently ignores already existing directories.
 *
 * \returns error_t::noerror on success, indicating the directories were created or already existed.
 * \returns other codes, error_t::exxxx (from errno) or error_t::filesystem, on errors.
 */
error_t createDirectories( const std::string &path /**< [in] the path of the directory(ies)to create */ );

/// Get the stem of the filename
/**
 * \returns the stem for the filename, that is without the path or extension
 */
std::string pathStem( const std::string &fname );

/// Get the base filename
/**
 * \returns the filename, including the extension but without the path
 */
std::string pathFilename( const std::string &fname );

/// Get the parent path from a filename
/**
 * \returns the parent path of the file
 */
std::string parentPath( const std::string &fname );

/// Get a list of file names from the specified directory, specifying a prefix, a substring to match, and an extension
/**
 * \returns mx::error_t::success on success
 * \returns mx::error_t::invalidarg if \p directory is not a directory
 * \returns mx::error_t::dirnotfound if \p directory does not exist
 * \returns mx::error_t::exception if an exception is thrown from the standard library
 *
 * \tparam verbose if true then error messages are printed as they occur
 *
 *
 */
template <class verboseT = verbose::d>
error_t getFileNames( std::vector<std::string> &fileNames, /** [out] The populated list of file names.*/
                      const std::string &directory,        /**< [in] The path to the directory to search.
                                                                     Can not be empty.*/
                      const std::string &prefix,           /**< [in] The file name prefix (the beginning
                                                                     characters of the file name) to search
                                                                     for. If "" then not used.*/
                      const std::string &substr,           /**< [in] A substring of the filename to search
                                                                     for. If "" then not used. Only matches
                                                                     after the first character.*/
                      const std::string &extension         /**< [in] The file name extension to search for.
                                                                     If "" then not used. This does not need
                                                                     to include the ".", as in".ext".*/
);

/// Prepend and/or append strings to a file name, leaving the directory and extension unaltered.
/**
 * \returns the new file name
 */
std::string fileNamePrependAppend( const std::string &fname,   /**< [in] the original file name, possibly including a
                                                                         directory and extension*/
                                   const std::string &prepend, /**< [in] is the string to insert at the beginning of the
                                                                         file name after the path*/
                                   const std::string &append   /**< [in] is the string to insert at the end of the file
                                                                         name, before the extension*/
);

/// Append a string to a file name, leaving the directory and extension unaltered.
/**
 * \returns the new file name
 */
std::string fileNameAppend( const std::string &fname, /**< [in] the original file name, possibly including
                                                                a directory and extension*/
                            const std::string &append /**< [in] is the string to insert at the end
                                                                of the file name, before the extension*/
);

/// Prepend strings to a file name, leaving the directory and extension unaltered.
/**
 * \returns the new file name
 */
std::string fileNamePrepend( const std::string &fname,  /**< [in] the original file name, possibly including
                                                                  a directory and extension*/
                             const std::string &prepend /**< [in] is the string to insert at the beginning of
                                                                  the file name after the path*/
);

/// Get the next file in a numbered sequence
/** Searches for files in the path designated by basename of the form basenameXXXXextension
 * where the number of digits in XXXX is set by the \a ndigit parameter.
 *
 * \warning this does not currently detect missing files in the sequence, e.g. if you have 0,1,3 in the directory this
 * will start with 2!
 *
 * \todo switch to using a regex or something so we can detect the missing file.
 *
 * \retval std::string containing the next filename.
 *
 */
std::string getSequentialFilename( const std::string &basename,       ///< [in] path and initial name of the file*/
                                   const std::string &extension = "", /**< [in] [optional] extension to append after the
                                                                                           number. Default is empty.*/
                                   const int startat = 0,             /**< [in] [optional] number to start the
                                                                                           search from.
                                                                                           Default is 0.*/
                                   int ndigit = 4                     /**< [in] [optional] number of digits in string
                                                                                           representation
                                                                                           of the number.Default is 4. */
);

/// Get the size in bytes of a file
/** Uses fstat.
 *
 * \returns the file size if fd is valid and no errors occur
 * \returns -1 on an error
 */
off_t fileSize( int fd /**< [in] an open file descriptor */ );

/// Get the size in bytes of a file pointed to by a FILE pointer
/** Uses fileno to get the associated descriptor, then uses fstat.
 *
 * \returns the file size if fd is valid and no errors occur
 * \returns -1 on an error
 *
 * \overload
 */
off_t fileSize( FILE *f /**< [in] an open file */ );

///@} -fileutils

/* ===================================================================== */
/*                       implementations                                 */

template <class verboseT>
error_t string2path( std::filesystem::path &path, const std::string &str )
{
    try
    {
        fileUtilsDetail::operationHook()( fileUtilsDetail::operation::stringToPath );
        path = str;

        return error_t::noerror;
    }
    catch( const std::bad_alloc &e )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc, e.what() );
        #else
            std::throw_with_nested(mx::exception<verboseT>(error_t::std_bad_alloc));
        #endif
        // clang-format on
    }
    catch( const std::filesystem::filesystem_error &e )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined(MXLIB_CATCH_NONALLOC_EXCEPTIONS)
            return internal::mxlib_error_report<verboseT>( error_t::std_filesystem_error, e.what() );
        #else
            std::throw_with_nested(mx::exception<verboseT>(error_t::std_filesystem_error));
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined(MXLIB_CATCH_NONALLOC_EXCEPTIONS)
            return internal::mxlib_error_report<verboseT>( error_t::std_exception, e.what() );
        #else
            std::throw_with_nested(mx::exception<verboseT>());
        #endif
        // clang-format on
    }
    catch( ... )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined(MXLIB_CATCH_NONALLOC_EXCEPTIONS)
            return internal::mxlib_error_report<verboseT>( error_t::exception, "unknown exception");
        #else
            std::throw_with_nested(mx::exception<verboseT>());
        #endif
        // clang-format on
    }
}

template <class verboseT>
bool exists( const std::string &strpath, mx::error_t &errc )
{
    std::error_code ec;

    std::filesystem::path path;

    try
    {
        errc = string2path<verboseT>( path, strpath );

        if( !!errc )
        {
            internal::mxlib_error_report<verboseT>( errc, "converting path" );
            return false;
        }
    }
    catch( const mx::exception<verboseT> &e )
    {
        std::throw_with_nested( mx::exception<verboseT>( e.code() ) );
    }

    bool ex = std::filesystem::exists( path, ec );

    if( ec.value() != 0 )
    {
        errc = mx::errno2error_t( ec.value() );
        if( errc == error_t::error )
        {
            errc = error_t::filesystem;
        }

        internal::mxlib_error_report<verboseT>( errc, ec.message() );

        return false;
    }

    errc = error_t::noerror;
    return ex;
}

template <class verboseT>
bool dir_exists_is( const std::string &dir, mx::error_t &errc )
{
    std::error_code ec;

    std::filesystem::path path;

    try
    {
        errc = string2path<verboseT>( path, dir );

        if( !!errc )
        {
            internal::mxlib_error_report<verboseT>( errc, "converting path" );
            return false;
        }
    }
    catch( const mx::exception<verboseT> &e )
    {
        std::throw_with_nested( mx::exception<verboseT>( e.code() ) );
    }

    bool exists = std::filesystem::exists( path, ec );

    // clang-format off
    #ifdef MXLIBTEST_DIREXISTSIS_ISEXISTSERR
        ec = std::error_code( EEXIST, std::system_category() ); // LCOV_EXCL_LINE
    #endif
    // clang-format on

    if( ec.value() != 0 )
    {
        errc = mx::errno2error_t( ec.value() );
        if( errc == error_t::error )
        {
            errc = error_t::filesystem;
        }

        internal::mxlib_error_report<verboseT>( errc, ec.message() );

        return false;
    }

    if( !exists )
    {
        return false;
    }

    bool isdir = std::filesystem::is_directory( path, ec );

    // clang-format off
    #ifdef MXLIBTEST_DIREXISTSIS_ISDIRERR
        ec = std::error_code( EACCES, std::system_category() ); // LCOV_EXCL_LINE
    #endif
    // clang-format on

    if( ec.value() != 0 )
    {
        errc = errno2error_t( ec.value() );
        if( errc == mx::error_t::error )
        {
            errc = mx::error_t::filesystem;
        }

        internal::mxlib_error_report<verboseT>( errc, ec.message() );

        return false;
    }

    errc = error_t::noerror;
    return isdir;
}

template <class verboseT>
error_t getFileNames( std::vector<std::string> &fileNames,
                      const std::string &directory,
                      const std::string &prefix,
                      const std::string &substr,
                      const std::string &extension )
{
    try // there are several things that can throw here
    {
        fileUtilsDetail::operationHook()( fileUtilsDetail::operation::getFileNames );
        fileNames.clear();

        if( std::filesystem::exists( directory ) )
        {
            if( std::filesystem::is_directory( directory ) )
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

                bool hasprefix = ( prefix.size() > 0 );

                bool hassub = ( substr.size() > 0 );

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
                            // This won't throw because:
                            //  - prefix has size > 0
                            //  - p.size() >= prefix.size()
                            //  - therefore prefix.size() > 0
                            //  - so pos1 = 0 will not throw.
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

                        size_t sspos = p.find( substr, 1 ); // only match if not prefix

                        if( sspos == std::string::npos )
                        {
                            continue;
                        }
                    }

                    // If here then it passed all checks
                    // this could throw
                    fileNames.push_back( it->path().native() );
                }

                std::sort( fileNames.begin(), fileNames.end() );
            }
            else
            {
                return internal::mxlib_error_report<verboseT>( error_t::invalidarg, directory + " is not a directory" );
            }
        }
        else
        {
            return internal::mxlib_error_report<verboseT>( error_t::dirnotfound, directory + " was not found" );
        }

        return error_t::noerror;
    }
    catch( const std::bad_alloc &e )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc, e.what() );;
        #else
            std::throw_with_nested(mx::exception<verboseT>(error_t::std_bad_alloc));
        #endif
        // clang-format on
    }
    catch( const std::filesystem::filesystem_error &e )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return internal::mxlib_error_report<verboseT>( error_t::std_filesystem_error, e.what() );
        #else
            std::throw_with_nested(mx::exception(error_t::std_filesystem_error));
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return internal::mxlib_error_report<verboseT>( error_t::std_exception, e.what() );
        #else
            std::throw_with_nested(mx::exception<verboseT>(error_t::std_exception));
        #endif
        // clang-format on
    }
    catch( ... )
    {
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return internal::mxlib_error_report<verboseT>( error_t::exception );
        #else
            std::throw_with_nested(mx::exception<verboseT>());
        #endif
        // clang-format on
    }
}

#ifdef MXLIBTEST_NAMESPACE
} // namespace MXLIBTEST_NAMESPACE
#endif

} // namespace ioutils
} // namespace mx

#endif // fileUtils_hpp
