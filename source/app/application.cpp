/** \file application.cpp
 * \author Jared R. Males
 * \brief Implementation of a class for managing applications
 *
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

#include "app/application.hpp"
#include "sys/environment.hpp"

namespace mx
{
namespace app
{

application::~application()
{
    return;
}

int application::main( int argc, char **argv )
{
    m_argc = argc;
    m_argv = argv;

    setup( argc, argv );

    if( doHelp )
    {
        help();
        return 1;
    }

    if( !m_preserveConfig )
    {
        config.clear();
    }

    return execute();
}

void application::setConfigPathGlobal( const std::string &s )
{
    m_configPathGlobal = s;
}

void application::setConfigPathUser( const std::string &s )
{
    m_configPathUser = s;
}

void application::setConfigPathLocal( const std::string &s )
{
    m_configPathLocal = s;
}

void application::setupConfig() // virtual
{
    return;
}

void application::loadConfig() // virtual
{
    return;
}

int application::execute() // virtual
{
    return 0;
}

void application::setup( int argc, char **argv )
{
    invokedName = argv[0];

    setupStandardConfig();
    setupStandardHelp();

    setupBasicConfig();
    setupConfig();

    setDefaults( argc, argv );

    config.readConfig( m_configPathGlobal, m_requireConfigPathGlobal );
    config.readConfig( m_configPathUser, m_requireConfigPathUser );
    config.readConfig( m_configPathLocal, m_requireConfigPathLocal );

    // Parse CL just to get the CL config.
    config.parseCommandLine( argc, argv, "config" );

    // And now get the value of it and parse it.
    loadStandardConfig();
    config.readConfig( m_configPathCL );

    // Now parse the command line for real.
    config.parseCommandLine( argc, argv );

    loadStandardHelp();

    loadBasicConfig();
    loadConfig();

    checkConfig();
}

int application::reReadConfig()
{
    config.readConfig( m_configPathGlobal );
    config.readConfig( m_configPathUser );
    config.readConfig( m_configPathLocal );

    config.readConfig( m_configPathCL );

    // Now parse the command line for real.
    config.parseCommandLine( m_argc, m_argv );

    return 0;
}

void application::setDefaults( int argc, char **argv )
{
    static_cast<void>( argc );
    static_cast<void>( argv );

    std::string tmp;

    if( m_configPathGlobal_env != "" )
    {
        m_configPathGlobal = sys::getEnv( m_configPathGlobal_env );
    }

    if( m_configPathUser_env != "" )
    {
        m_configPathUser = sys::getEnv( m_configPathUser_env );
    }

    if( m_configPathUser != "" )
    {
        // If it's a relative path, add it to the HOME directory
        if( m_configPathUser[0] != '/' && m_configPathUser[0] != '~' )
        {
            tmp = sys::getEnv( "HOME" );
            tmp += "/" + m_configPathUser;
            m_configPathUser = tmp;
        }
    }

    if( m_configPathLocal_env != "" )
    {
        m_configPathLocal = sys::getEnv( m_configPathLocal_env.c_str() );
    }

    if( m_configPathCLBase_env != "" )
    {
        m_configPathCLBase = sys::getEnv( m_configPathCLBase_env.c_str() );

        if( m_configPathCLBase.size() > 0 )
            if( m_configPathCLBase[m_configPathCLBase.size() - 1] != '/' )
                m_configPathCLBase += '/';
    }

    return;
}

void application::setupStandardConfig() // virtual
{
    config.add( "config", "c", "config", argType::Required, "", "config", false, "string", "A local config file" );
}

void application::setupStandardHelp() // virtual
{
    config.add( "help", "h", "help", argType::True, "", "", false, "none", "Print this message and exit" );
}

void application::loadStandardConfig() // virtual
{
    config( m_configPathCL, "config" );
    m_configPathCL = m_configPathCLBase + m_configPathCL;
}

void application::loadStandardHelp() // virtual
{
    config( doHelp, "help" );
}

void application::setupBasicConfig() // virtual
{
    return;
}

void application::loadBasicConfig() // virtual
{
    return;
}

void application::checkConfig() // virtual
{
    return;
}

void application::optionHelp( configTarget &tgt,
                              ioutils::textTable &tt ) // virtual
{
    std::string tmp;
    int row = tgt.orderAdded;

    if( tgt.shortOpt != "" )
    {
        tmp = "-" + tgt.shortOpt;
        tt.addCell( row, 0, tmp );
    }

    if( tgt.longOpt != "" )
    {
        tmp = "--" + tgt.longOpt;
        tt.addCell( row, 1, tmp );
    }

    tmp = "";
    if( tgt.section != "" )
    {
        tmp = tgt.section + ".";
    }

    if( tgt.keyword != "" )
    {
        tmp += tgt.keyword;
        tt.addCell( row, 2, tmp );
    }

    if( tgt.helpType != "" )
    {
        tmp = "<" + tgt.helpType + "> ";
        tt.addCell( row, 3, tmp );
    }

    tt.addCell( row, 4, tgt.helpExplanation );
}

void application::help() // virtual
{
    appConfigurator::targetIterator targIt;
    appConfigurator::clOnlyTargetIterator clOnlyTargIt;

    ioutils::textTable tt;

    int otherColWidth = m_helpSOColWidth + m_helpLOColWidth + m_helpCFColWidth + m_helpTypeColWidth;

    tt.m_colWidths = {
        m_helpSOColWidth, m_helpLOColWidth, m_helpCFColWidth, m_helpTypeColWidth, m_helpWidth - otherColWidth - 4 - 4 };
    tt.m_lineStart = "    ";
    tt.m_colSep = " ";
    tt.m_rowSep = "";

    std::cerr << "usage: " << invokedName << " [options] " << m_nonOptionHelp << "\n";
    std::cerr << "\n";
    std::cerr << "  Required arguments:\n";

    for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt != config.clOnlyTargets.end(); ++clOnlyTargIt )
    {
        if( clOnlyTargIt->isRequired == true )
        {
            optionHelp( *clOnlyTargIt, tt );
        }
    }

    for( targIt = config.m_targets.begin(); targIt != config.m_targets.end(); ++targIt )
    {
        if( targIt->second.isRequired == true )
        {
            optionHelp( targIt->second, tt );
        }
    }

    tt.outPut( std::cerr );
    tt.m_rows.clear();

    // row = 0;
    std::cerr << "\n  Optional arguments:\n";

    for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt != config.clOnlyTargets.end(); ++clOnlyTargIt )
    {
        if( clOnlyTargIt->isRequired == false )
        {
            optionHelp( *clOnlyTargIt, tt );
        }
    }

    for( targIt = config.m_targets.begin(); targIt != config.m_targets.end(); ++targIt )
    {
        if( targIt->second.isRequired == false )
        {
            optionHelp( targIt->second, tt );
        }
    }

    tt.outPut( std::cerr );
}

} // namespace app
} // namespace mx
