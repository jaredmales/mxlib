/** \file environment.cpp
 * \author Jared R. Males
 * \brief Utilities for working with the environment
 * \ingroup utils_files
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

#include "sys/environment.hpp"

namespace mx
{
namespace sys
{

std::string getEnv( const std::string &estr )
{
    char *e = getenv( estr.c_str() );

    if( e )
        return std::string( e );
    else
        return std::string( "" );
}

} // namespace sys
} // namespace mx
