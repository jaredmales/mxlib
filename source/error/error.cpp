/** \file error.cpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Definitions for the mxlib error reporting system.
 * \ingroup error_handling_files
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
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

#include "error/error.hpp"
#include <format>

namespace mx
{

namespace internal
{

template <>
std::string mxlib_error_message<verbose::o>( [[maybe_unused]] const error_t &code,
                                             [[maybe_unused]] const std::string &expl,
                                             [[maybe_unused]] const std::source_location &loc )
{
    return "";
}

template <>
std::string mxlib_error_message<verbose::v>( const error_t &code,
                                             const std::string &expl,
                                             [[maybe_unused]] const std::source_location &loc )
{
    return std::format( "{}: {}.", errorName( code ), expl );
}

template <>
std::string
mxlib_error_message<verbose::vv>( const error_t &code, const std::string &expl, const std::source_location &loc )
{
    return std::format( "{} ({}): {}. [{} {}]",
                        errorMessage( code ),
                        errorName( code ),
                        expl,
                        loc.file_name(),
                        loc.line() );
}

template <>
std::string
mxlib_error_message<verbose::vvv>( const error_t &code, const std::string &expl, const std::source_location &loc )
{
    return std::format( "An error has occurred in mxlib:\n"
                        "        error: {} ({})\n"
                        "  explanation: {}\n"
                        "      in file: {}\n"
                        "      at line: {}\n"
                        "       source: {}\n",
                        errorMessage( code ),
                        errorName( code ),
                        expl,
                        loc.file_name(),
                        loc.line(),
                        loc.function_name() );
}

template <>
std::string mxlib_error_message<verbose::o>( const error_t &code, [[maybe_unused]] const std::source_location &loc )
{
    return "";
}

template <>
std::string mxlib_error_message<verbose::v>( const error_t &code, [[maybe_unused]] const std::source_location &loc )
{
    return errorName( code );
}

template <>
std::string mxlib_error_message<verbose::vv>( const error_t &code, const std::source_location &loc )
{
    return std::format( "{} ({}). [{} {}]", errorMessage( code ), errorName( code ), loc.file_name(), loc.line() );
}

template <>
std::string mxlib_error_message<verbose::vvv>( const error_t &code, const std::source_location &loc )
{
    return std::format( "An error has occurred in mxlib:\n"
                        "        error: {} ({})\n"
                        "      in file: {}\n"
                        "      at line: {}\n"
                        "       source: {}\n",
                        errorMessage( code ),
                        errorName( code ),
                        loc.file_name(),
                        loc.line(),
                        loc.function_name() );
}

template <>
error_t mxlib_error_report<verbose::o>( const error_t &code,
                                        [[maybe_unused]] const std::string &expl,
                                        [[maybe_unused]] const std::source_location &loc )
{
    return code;
}

template <>
error_t mxlib_error_report<verbose::o>( const error_t &code, [[maybe_unused]] const std::source_location &loc )
{
    return code;
}


} // namespace internal

template <>
std::string error_message<verbose::o>( const error_t &code,
                                       [[maybe_unused]] const std::string &expl,
                                       [[maybe_unused]] const std::source_location &loc )
{
    return "";
}

template <>
std::string error_message<verbose::v>( const error_t &code, const std::string &expl, const std::source_location &loc )
{
    return internal::mxlib_error_message<verbose::v>( code, expl, loc );
}

template <>
std::string error_message<verbose::vv>( const error_t &code, const std::string &expl, const std::source_location &loc )
{
    return internal::mxlib_error_message<verbose::vv>( code, expl, loc );
}

template <>
std::string error_message<verbose::vvv>( const error_t &code, const std::string &expl, const std::source_location &loc )
{
    return std::format( "An error has occurred:\n"
                        "        error: {} ({})\n"
                        "  explanation: {}\n"
                        "      in file: {}\n"
                        "      at line: {}\n"
                        "       source: {}\n",
                        errorMessage( code ),
                        errorName( code ),
                        expl,
                        loc.file_name(),
                        loc.line(),
                        loc.function_name() );
}

template <>
std::string error_message<verbose::o>( [[maybe_unused]] const error_t &code,
                                       [[maybe_unused]] const std::source_location &loc )
{
    return "";
}

template <>
std::string error_message<verbose::v>( const error_t &code, const std::source_location &loc )
{
    return internal::mxlib_error_message<verbose::v>( code, loc );
}

template <>
std::string error_message<verbose::vv>( const error_t &code, const std::source_location &loc )
{
    return internal::mxlib_error_message<verbose::vv>( code, loc );
}

template <>
std::string error_message<verbose::vvv>( const error_t &code, const std::source_location &loc )
{
    return std::format( "An error has occurred:\n"
                        "    error: {} ({})\n"
                        "  in file: {}\n"
                        "  at line: {}\n"
                        "   source: {}\n",
                        errorMessage( code ),
                        errorName( code ),
                        loc.file_name(),
                        loc.line(),
                        loc.function_name() );
}

template <>
error_t error_report<verbose::o>( const error_t &code,
                                  [[maybe_unused]] const std::string &expl,
                                  [[maybe_unused]] const std::source_location &loc )
{
    return code;
}

template <>
error_t error_report<verbose::o>( const error_t &code, [[maybe_unused]] const std::source_location &loc )
{
    return code;
}

} // namespace mx
