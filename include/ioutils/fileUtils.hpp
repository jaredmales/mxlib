
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

#ifndef filtUtils_hpp
#define filtUtils_hpp

#include <string>
#include <vector>

#include "../mxlib.hpp"

namespace mx
{
namespace ioutils
{

/** \addtogroup fileutils
 * @{
 */

/// Check if a path exists
/**
 * \returns true if the path exists
 * \returns false if it doesn't
 */
bool exists( const std::string &path /**< [in] the path to check for existence */ );

/// Create a directory or directories
/** This will create any directories in path that don't exist.  It silently ignores already existing directories.
 *
 * \returns 0 on success, indicating the directories were created or already existed.
 * \returns -1 on error
 */
int createDirectories( const std::string &path /**< [in] the path of the directory(ies)to create */ );

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
 * \b Tests:
 * - Getting a list of files \ref tests_ioutils_fileUtils_getFileNames "[test doc]"
 *
 */
template<class verboseT=verbose::vvv>
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

template<>
error_t getFileNames<verbose::o>( std::vector<std::string> &fileNames,
                        const std::string &directory,
                        const std::string &prefix,
                        const std::string &substr,
                        const std::string &extension
);

template<>
error_t getFileNames<verbose::v>( std::vector<std::string> &fileNames,
                        const std::string &directory,
                        const std::string &prefix,
                        const std::string &substr,
                        const std::string &extension
);

template<>
error_t getFileNames<verbose::vv>( std::vector<std::string> &fileNames,
                        const std::string &directory,
                        const std::string &prefix,
                        const std::string &substr,
                        const std::string &extension
);

template<>
error_t getFileNames<verbose::vvv>( std::vector<std::string> &fileNames,
                        const std::string &directory,
                        const std::string &prefix,
                        const std::string &substr,
                        const std::string &extension
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
 * \test Verify creation of sequential file names \ref tests_ioutils_fileUtils_getSequentialFilename "[test doc]"
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

} // namespace ioutils
} // namespace mx

#endif // fileUtils_hpp
