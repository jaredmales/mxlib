/** \file gitRepo.cpp
 * \author Jared R. Males
 * \brief Interrogate the current state of a git repository (definitions)
 * \ingroup utils_files
 */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#include "mxError.hpp"
#include "sys/gitRepo.hpp"

#include "ipc/processInterface.hpp"
#include "ioutils/stringUtils.hpp"
#include "ioutils/fileUtils.hpp"

namespace mx
{
namespace sys
{

gitRepo::gitRepo()
{
}

gitRepo::gitRepo( const std::string &d )
{
    dir( d );
}

void gitRepo::dir( const std::string &d )
{
    m_dir = d;

    getGitName();
    getGitHash();
    getGitModified();
    getGitFileState();
}

std::string gitRepo::dir()
{
    return m_dir;
}

std::string gitRepo::gitDir()
{
    return m_dir + "/.git";
}

int gitRepo::getGitName()
{
    int retVal;
    std::vector<std::string> stdOut, stdErr;

    if( ipc::runCommand( retVal, stdOut, stdErr, { "git", "--git-dir=" + gitDir(), "rev-parse", "--show-toplevel" } ) <
        0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error running git command" );
        return -1;
    }

    if( stdErr.size() > 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error returned by git" );

        for( size_t n = 0; n < stdErr.size(); ++n )
        {
            std::cerr << "err: " << stdErr[n] << "\n";
        }

        return -1;
    }

    if( stdOut.size() < 1 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "nothing returned by git" );
        return -1;
    }

    if( stdOut.size() > 1 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "too much returned by git" );
        return -1;
    }

    m_name = ioutils::pathFilename( stdOut[0] );

    return 0;
}

int gitRepo::getGitHash()
{
    int retVal;
    std::vector<std::string> stdOut, stdErr;

    if( ipc::runCommand( retVal,
                         stdOut,
                         stdErr,
                         { "git", "--git-dir=" + gitDir(), "--work-tree=" + dir(), "log", "-1", "--format=%H" } ) < 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error running git command" );
        return -1;
    }

    if( stdErr.size() > 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error returned by git" );

        for( size_t n = 0; n < stdErr.size(); ++n )
        {
            std::cerr << "err: " << stdErr[n] << "\n";
        }

        return -1;
    }

    if( stdOut.size() < 1 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "nothing returned by git" );
        return -1;
    }

    if( stdOut.size() > 1 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "too much returned by git" );
        return -1;
    }

    m_hash = stdOut[0];

    return 0;
}

int gitRepo::getGitModified()
{
    int retVal;
    std::vector<std::string> stdOut, stdErr;

    if( ipc::runCommand( retVal,
                         stdOut,
                         stdErr,
                         { "git", "--git-dir=" + gitDir(), "--work-tree=" + dir(), "diff-index", "HEAD", "--" } ) < 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error running git command" );
        return -1;
    }

    if( stdErr.size() > 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error returned by git" );
        for( size_t n = 0; n < stdErr.size(); ++n )
        {
            std::cerr << "err: " << stdErr[n] << "\n";
        }

        return -1;
    }

    m_modified = ( stdOut.size() > 0 ); // rv;

    return 0;
}

int gitRepo::getGitFileState()
{
    int retVal;
    std::vector<std::string> stdOut, stdErr;

    if( ipc::runCommand(
            retVal, stdOut, stdErr, { "git", "--git-dir=" + gitDir(), "--work-tree=" + dir(), "status" } ) < 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error running git command" );
        return -1;
    }

    if( stdErr.size() > 0 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "error returned by git" );

        for( size_t n = 0; n < stdErr.size(); ++n )
        {
            std::cerr << "err: " << stdErr[n] << "\n";
        }

        return -1;
    }

    if( stdOut.size() < 1 )
    {
        mxError( "mx::sys::gitRepo::getGitName", MXE_PROCERR, "nothing returned by git" );
        return -1;
    }

    m_modifiedFiles.clear();
    m_deletedFiles.clear();
    m_renamedFiles.clear();
    m_renamedFiles2.clear();
    m_untrackedFiles.clear();

    for( size_t n = 0; n < stdOut.size(); ++n )
    {
        if( stdOut[n].find( "On b" ) != std::string::npos )
        {
            m_branch = stdOut[n].substr( sizeof( "On branch" ) );
        }
        else if( stdOut[n].find( "Changes not" ) != std::string::npos )
        {
            // Changes not staged for commit

            ++n;
            while( stdOut[n].size() > 0 )
            {
                if( n >= stdOut.size() )
                    break;

                if( stdOut[n].find( "mod" ) != std::string::npos )
                {
                    m_modifiedFiles.insert( ioutils::removeWhiteSpace( stdOut[n].substr( sizeof( "modified:" ) ) ) );
                }
                else if( stdOut[n].find( "del" ) != std::string::npos )
                {
                    m_deletedFiles.insert( ioutils::removeWhiteSpace( stdOut[n].substr( sizeof( "deleted:" ) ) ) );
                }
                else if( stdOut[n].find( "ren" ) != std::string::npos )
                {
                    size_t a = stdOut[n].find( " -> " );

                    m_renamedFiles.insert( ioutils::removeWhiteSpace(
                        stdOut[n].substr( sizeof( "renamed:" ), a - sizeof( "renamed:" ) ) ) );
                    m_renamedFiles2.insert( ioutils::removeWhiteSpace( stdOut[n].substr( a + sizeof( " ->" ) ) ) );
                }
                ++n;
            }
        }
        else if( stdOut[n].find( "Changes to" ) != std::string::npos )
        {
            // Changes to be committed

            ++n;
            while( stdOut[n].size() > 0 )
            {
                if( n >= stdOut.size() )
                    break;

                if( stdOut[n].find( "mod" ) != std::string::npos )
                {
                    m_modifiedFiles.insert( ioutils::removeWhiteSpace( stdOut[n].substr( sizeof( "modified:" ) ) ) );
                }
                else if( stdOut[n].find( "del" ) != std::string::npos )
                {
                    m_deletedFiles.insert( ioutils::removeWhiteSpace( stdOut[n].substr( sizeof( "deleted:" ) ) ) );
                }
                else if( stdOut[n].find( "ren" ) != std::string::npos )
                {
                    size_t a = stdOut[n].find( " -> " );

                    m_renamedFiles.insert( ioutils::removeWhiteSpace(
                        stdOut[n].substr( sizeof( "renamed:    " ), a - sizeof( "renamed:" ) ) ) );
                    m_renamedFiles2.insert( ioutils::removeWhiteSpace( stdOut[n].substr( a + sizeof( " ->" ) ) ) );
                }
                ++n;
            }
        }
        else if( stdOut[n].find( "Untracked" ) != std::string::npos )
        {
            // Untracked files:

            ++n;
            ++n;
            while( stdOut[n].size() > 0 )
            {
                if( n >= stdOut.size() )
                    break;
                m_untrackedFiles.insert( ioutils::removeWhiteSpace( stdOut[n] ) );
                ++n;
            }
        }
    }

    return 0;
}

std::string gitRepo::name()
{
    return m_name;
}

std::string gitRepo::branch()
{
    return m_branch;
}

std::string gitRepo::hash()
{
    return m_hash;
}

bool gitRepo::modified()
{
    return m_modified;
}

bool gitRepo::isNotCommitted( const std::string &file )
{
    return ( m_modifiedFiles.count( file ) + m_deletedFiles.count( file ) + m_renamedFiles.count( file ) +
                 m_renamedFiles2.count( file ) + m_untrackedFiles.count( file ) >
             0 );
}

} // namespace sys
} // namespace mx
