/** \file mxlib.hpp
 * \author Jared R. Males
 * \brief Declarations of some libarary wide utilities
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

#ifndef mxlib_hpp
#define mxlib_hpp

#include <ios>

#include "mxlib_uncomp_version.h"

#include "error/error.hpp"
#include "error/exception.hpp"

namespace mx
{

/// Get the git repo name for the compiled portion of mxlib
/**
 * \returns the git repo name from when mxlib was compiled
 */
const char *mxlib_comp_git_repo();

/// Get the git repo url for the compiled portion of mxlib
/**
 * \returns the git repo url from when mxlib was compiled
 */
const char *mxlib_comp_git_url();

/// Get the git repo branch for the compiled portion of mxlib
/**
 * \returns the git repo branch from when mxlib was compiled
 */
const char *mxlib_comp_git_branch();

/// Get the source path for the compiled portion of mxlib
/**
 * \returns the path to the source from when mxlib was compiled
 */
const char *mxlib_comp_source_path();

/// Get the git repo sha1 hash for the compiled portion of mxlib
/**
 * \returns the git repo sha1 hash from when mxlib was compiled
 */
const char *mxlib_comp_git_sha1();

/// Get the git repo modified state the compiled portion of mxlib
/**
 * \returns the git repo modified state from when mxlib was compiled
 */
const bool mxlib_comp_git_modified();

/// Get the git repo untracked files flag for the compiled portion of mxlib
/**
 * \returns the git repo untracked files flag from when mxlib was compiled
 */
const bool mxlib_comp_git_untracked();

[[deprecated("use mxlib_comp_git_branch()")]]
const char *mxlib_comp_current_branch();

[[deprecated("use mxlib_comp_git_sha1()")]]
const char *mxlib_comp_current_sha1();

[[deprecated("use mxlib_comp_git_modified()")]]
const bool mxlib_comp_repo_modified();

/// Dump the git status of a repository to a stream
/** Prints the provided SHA1 hash and whether or not the library
 * has been modified since that commit.
 *
 * \tparam iosT a std::ostream-like type
 * \tparam comment a character to print at the beginning of each line.  Default is '#'.
 */
template <typename iosT, char comment = '#'>
iosT &dumpGitStatus(
    iosT &ios,                      ///< [in] a std::ostream-like stream.
    const std::string &repoName,    ///< [in] The name of the repository
    const std::string &branch,      ///< [in] The name of the branch
    const std::string &sha1,        ///< [in] The sha1 hash for the current commit
    const bool &modified,           ///< [in] Whether or not the repository is currently modified
    const std::string &section = "" ///< [in] [optional] Descriptive sub-section name (e.g. headers vs. compiled)
)
{
    // This causes the stream to not output a '\0' if comment = '\0'
    char c[] = { comment, '\0' };

    bool sect = !( section == "" );
    std::string space = "  ";
    if( sect )
        space += " ";

    ios << c << "--------------------------------------------\n";
    ios << c << " " << repoName << " git status:\n";
    if( sect )
        ios << c << "  " << section << ":\n";
    ios << c << space << "branch: " << branch << "\n";
    ios << c << space << "SHA1: " << sha1 << "\n";
    ios << c << space << "modified flag: " << std::boolalpha << modified << "\n";
    ios << c << "--------------------------------------------\n";

    return ios;
}

/// Dump the current git status of the mxlib library to a stream
/** Prints the current SHA1 hash and whether or not the library
 * has been modified since that commit.
 *
 * \tparam iosT a std::ostream-like type
 * \tparam comment a character to print at the beginning of each line.  Default is '#'.
 */
template <typename iosT, char comment = '#'>
iosT &dumpGitStatus( iosT &ios /**< [in] a std::ostream-like stream. */ )
{
    return dumpGitStatus( ios,
                          "mxlib",
                          mxlib_comp_git_branch(),
                          mxlib_comp_git_sha1(),
                          mxlib_comp_repo_modified() );
}

} // namespace mx

#endif // mxlib_hpp
