/** \file clOptions.cpp
 * \author Jared R. Males
 * \brief Implementatino of a command line parser
 *
 * \ingroup mxApp_files
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

#include "app/clOptions.hpp"

namespace mx
{
namespace app
{

static const char *falseStr = "false";
static const char *trueStr = "true";
static const char *blankStr = "";

static option::ArgStatus Arg_Required( const option::Option &option, bool msg )
{
    static_cast<void>( msg );

    if( option.arg != 0 )
        return option::ARG_OK;

    return option::ARG_ILLEGAL;
}

clOptions::~clOptions()
{
    if( options )
        delete[] options;
    if( buffer )
        delete[] buffer;
}

void clOptions::clear()
{
    if( options )
        delete[] options;
    options = 0;

    if( buffer )
        delete[] buffer;
    buffer = 0;

    map.clear();
    typeMap.clear();
    descriptions.clear();
}

void clOptions::add( const std::string &optName, const char *const shortOpt, const char *const longOpt, int argT )
{
    mapIterator it;
    it = map.find( optName );

    if( it == map.end() )
    {
        if( argT == argType::Optional )
        {
            descriptions.push_back( { nOpts, argT, shortOpt, longOpt, option::Arg::Optional, "" } );
        }
        else if( argT == argType::Required )
        {
            descriptions.push_back( { nOpts, argT, shortOpt, longOpt, Arg_Required, "" } );
        }
        else
        {
            descriptions.push_back( { nOpts, argT, shortOpt, longOpt, option::Arg::None, "" } );
        }

        map.insert( { optName, nOpts } );
        typeMap.insert( { optName, argT } );

        ++nOpts;
        return;
    }
}

void clOptions::parse( int argc, char **argv, std::vector<std::string> *nonOptions )
{
    argc -= ( argc > 0 );
    argv += ( argc > 0 ); // skip program name argv[0] if present

    // If not already done, we push the unknown catcher and the termination descriptor.
    if( descriptions.back().index != 0 )
    {
        descriptions.push_back(
            { nOpts, 0, "", "", option::Arg::None, "" } ); // This is inserted to catch unknown options
        descriptions.push_back( { 0, 0, 0, 0, 0, 0 } );
    }

    // Now allocate.
    option::Stats stats( descriptions.data(), argc, argv );
    options = new option::Option[stats.options_max];
    buffer = new option::Option[stats.buffer_max];

    option::Parser parse( false, descriptions.data(), argc, argv, options, buffer );

    if( nonOptions )
    {
        nonOptions->resize( parse.nonOptionsCount() );

        for( int i = 0; i < parse.nonOptionsCount(); ++i )
        {
            ( *nonOptions )[i] = parse.nonOption( i );
        }
    }
}

const char *clOptions::operator[]( const std::string &key )
{
    mapIterator it = map.find( key );

    if( it == map.end() )
    {
        std::cerr << "oh no\n";
        return 0;
    }

    int typeDesc = typeMap[key];

    // If this option is not pressent, either return the opposite T/F condition or blank
    if( options[it->second].type() == argType::None )
    {
        if( typeDesc == argType::False )
            return trueStr;
        if( typeDesc == argType::True )
            return falseStr;
        return blankStr;
    }

    if( typeDesc == argType::False || typeDesc == argType::True )
    {
        if( options[it->second].last()->type() == argType::False )
            return falseStr;
        else
            return trueStr;
    }

    if( options[it->second].arg == 0 )
        return blankStr;

    return options[it->second].last()->arg;
}

int clOptions::count( const std::string &key )
{
    mapIterator it = map.find( key );

    if( it == map.end() )
        return -1;

    return options[it->second].count();
}

void clOptions::getAll( std::vector<std::string> &args, const std::string &key )
{
    mapIterator it = map.find( key );

    if( it == map.end() )
    {
        std::cerr << "oh no\n";
        return;
    }

    int typeDesc = typeMap[key];

    // If this option is not present, either return the opposite T/F condition or blank
    if( options[it->second].type() == argType::None )
    {

        if( typeDesc == argType::False )
        {
            args.resize( 1 );
            args[0] = trueStr;
            return;
        }

        if( typeDesc == argType::True )
        {
            args.resize( 1 );
            args[0] = falseStr;
            return;
        }

        args.clear();

        return;
    }

    if( typeDesc == argType::False || typeDesc == argType::True )
    {
        args.resize( 1 );

        if( options[it->second].last()->type() == argType::False )
            args[0] = falseStr;
        else
            args[0] = trueStr;
        return;
    }

    if( options[it->second].arg == 0 )
    {
        args.clear();
        return;
    }

    int N = options[it->second].count();

    args.resize( N );

    int i = 0;

    for( option::Option *opt = options[it->second]; opt != NULL && i < N; opt = opt->next() )
    {
        args[i] = opt->arg;
        ++i;
    }
}

bool clOptions::optSet( const std::string &key )
{
    mapIterator it = map.find( key );

    if( it == map.end() )
        return false; // Not found --> e.g. if neither command line short nor long option set.

    if( options[it->second].type() != argType::None )
        return true;
    return false;
};

int clOptions::numUnknown()
{
    // The dummy description to catch unknown options is inserted at position nOpts.
    return options[nOpts].count();
}

int clOptions::unknown( std::vector<std::string> &unk )
{
    unk.clear();

    if( numUnknown() == 0 )
        return 0;

    for( option::Option *opt = options[nOpts]; opt != NULL; opt = opt->next() )
    {
        unk.push_back( opt->name );
    }

    return 0;
}

} // namespace app
} // namespace mx
