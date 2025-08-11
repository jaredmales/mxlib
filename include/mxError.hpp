/** \file mxError.hpp
 * \brief The mxlib error reporting system.
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

#ifndef __mxError__
#define __mxError__

#include <cerrno>
#include <cstring>
#include <sstream>
#include <source_location>
#include <iostream>
#include <string>

#include "error/mxErrorOld.hpp"
#include "error/error_t.hpp"

namespace mx
{

/** \defgroup error_verbosity Error Report Verbosity
 * \ingroup error_handling
 *
 * Error reports can be controlled with a verbosity template parameter,
 * usually called `verboseT`.  The types in this namespace are the standard
 * ones supported by `mxlib`.
 *
 */

/// Namespace for the error reporting verbosity levels
/**
 * \ingroup error_verbosity
 */
namespace verbose
{

/// Verbosity level 0, no reports are generated or printed to stderr.
/** o stands for off
 *
 * \ingroup error_verbosity
 */
struct o
{
    static constexpr int level = 0;
};

/// Verbosity level 1.  Minimal reports with no source location information.
/** A typical example:
 * \verbatim
  dirnotfound: /tmp/fileUtils_test/dirnf was not found.
 * \endverbatim
 *
 * \ingroup error_verbosity
 */
struct v
{
    static constexpr int level = 1;
};

/// Verbosity level 2.  Additional information is provided, including source file and line
/** A typical example:
 * \verbatim
  The directory was not found (dirnotfound): /tmp/fileUtils_test/dirnf was not found. [ioutils/fileUtils.cpp 186]
 * \endverbatim
 *
 * \ingroup error_verbosity
 */
struct vv
{
    static constexpr int level = 2;
};

/// Verbosity level 3.  A full report is provided.
/** A typical example:
 * \verbatim
  An error has occurred in mxlib:
        error: The directory was not found (dirnotfound)
  explanation: /tmp/fileUtils_test/dirnf was not found
      in file: ioutils/fileUtils.cpp
      at line: 186
       source: mx::error_t mx::ioutils::impl::getFileNames(std::vector<std::__cxx11::basic_string<char> >&,
  const std::string&, const std::string&, const std::string&, const std::string&) [with verboseT =
  mx::verbose::vvv; std::string = std::__cxx11::basic_string<char>]
  \endverbatim
 *
 * \ingroup error_verbosity
 */
struct vvv
{
    static constexpr int level = 3;
};

} // namespace verbose

/** \defgroup error_internal Internal Error Reporting
 * \ingroup error_handling
 *
 * These functions are meant to be used by mxlib itself.  You should use the non-`internal` versions,
 * which don't start with `mxlib_`, in most
 * other cases.
 */

/// Namespace for the internal error reporting functions.
/**
 * \ingroup error_internal
 */
namespace internal
{

/// Format a report given an mxlib \ref error_t code and explanation.
/** What is included depends on the verbosity level set by the template parameter
 *  This is for internal `mxlib` use, it includes `mxlib` in the vvv report.
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the formatted message
 *
 * \ingroup error_internal
 */
template <class verboseT>
std::string mxlib_error_message( const error_t &code,     /**< [in] is an mx::error_t error code*/
                                 const std::string &expl, /**< [in] [optional] if more information can be provided,
                                                                               use this to inform the user.*/
                                 const std::source_location &loc /**< [in] [optional] source location */
                                 = std::source_location::current() );

/// Specialization of \ref mxlib_error_message for \ref verbose::o
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref mxlib_error_message for \ref verbose::v
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<verbose::v>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref mxlib_error_message for \ref verbose::vv
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<verbose::vv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref mxlib_error_message for \ref verbose::vvv
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<verbose::vvv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Format a report given an mxlib \ref error_t code.
/** What is included depends on the verbosity level set by the template parameter
 *  This is for internal mxlib use, it includes `mxlib` in the vvv report.
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the formatted message
 *
 * \ingroup error_internal
 */
template <class verboseT>
std::string mxlib_error_message( const error_t &code,            /**< [in] is an mx::error_t error code*/
                                 const std::source_location &loc /**< [in] [optional] source location */
                                 = std::source_location::current() );

/** \brief Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for
 * \ref verbose::o
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<verbose::o>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for \ref verbose::v
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<verbose::v>( const error_t &code, const std::source_location &loc );

/** \brief  Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for \ref verbose::vv
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<verbose::vv>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for \ref verbose::vvv
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<verbose::vvv>( const error_t &code, const std::source_location &loc );

/// Print a report to stderr given an mxlib \ref error_t \p code and explanation and return the \p code.
/** What is printed depends on the verbosity level set by the template parameter
 * This is for internal mxlib use, it includes `mxlib` in the vvv report.
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the provided \ref error_t \p code
 *
 * \ingroup error_internal
 */
template <class verboseT>
error_t mxlib_error_report( const error_t &code,            /**< [in] is an mx::error_t error code*/
                            const std::string &expl,        /**< [in] [optional] if more information can be provided,
                                                                                 use this to inform the user.*/
                            const std::source_location &loc /**< [in] [optional] source location */
                            = std::source_location::current() )
{
    if( verboseT::level > 0 )
    {
        std::cerr << mxlib_error_message<verboseT>( code, expl, loc ) << '\n';
    }

    return code;
}

/// Specialization of \ref mxlib_error_report for \ref verbose::o
/**
 * \ingroup error_internal
 */
template <>
error_t mxlib_error_report<verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Print a report to stderr given an mxlib \ref error_t \p code and return the \p code.
/** What is printed depends on the verbosity level set by the template parameter
 * This is for internal mxlib use, it includes `mxlib` in the vvv report.
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the provided \ref error_t \p code
 *
 * \ingroup error_internal
 */
template <class verboseT>
error_t mxlib_error_report( const error_t &code /**< [in] is an mx::error_t error code*/,
                            const std::source_location &loc /**< [in] [optional] source location */
                            = std::source_location::current() )
{
    if( verboseT::level > 0 )
    {
        std::cerr << mxlib_error_message<verboseT>( code, loc ) << '\n';
    }

    return code;
}

/** \brief  Specialization of \ref mxlib_error_report(const error_t &, const std::source_location&) "mxlib_error_report"
 * for \ref verbose::o
 *
 * \ingroup error_internal
 */
template <>
error_t mxlib_error_report<verbose::o>( const error_t &code, const std::source_location &loc );

} // namespace internal

/// Format a report given an mxlib \ref error_t \p code and explanation.
/** What is included depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the formatted message
 *
 * \ingroup error_handling
 */
template <class verboseT>
std::string error_message( const error_t &code,            /**< [in] is an mx::error_t error code*/
                           const std::string &expl,        /**< [in] [optional] if more information can be provided,
                                                                     use this to inform the user.*/
                           const std::source_location &loc /**< [in] [optional] source location */
                           = std::source_location::current() );

/// Specialization of \ref error_message for \ref verbose::o
/**
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref error_message for \ref verbose::v
/**
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::v>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref error_message for \ref verbose::vv
/**
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::vv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref error_message for \ref verbose::vvv
/**
 * \ingroup error_handling
 */
template <>
std::string
error_message<verbose::vvv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Format a report given an mxlib \ref error_t \p code.
/** What is included depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the formatted message
 *
 * \ingroup error_handling
 */
template <class verboseT>
std::string error_message( const error_t &code,            /**< [in] is an mx::error_t error code*/
                           const std::source_location &loc /**< [in] [optional] source location */
                           = std::source_location::current() );

/** \brief Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for  \ref
 * verbose::o
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::o>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for \ref
 * verbose::v
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::v>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for \ref
 * verbose::vv
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::vv>( const error_t &code, const std::source_location &loc );

/** \brief  Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for \ref
 * verbose::vvv
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<verbose::vvv>( const error_t &code, const std::source_location &loc );

/// Print a report to stderr given an mxlib \ref error_t \p code and explanation and return the \p code.
/** What is printed depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the provided \ref error_t \p code
 *
 * \ingroup error_handling
 */
template <class verboseT>
error_t error_report( const error_t &code,            /**< [in] is an mx::error_t error code*/
                      const std::string &expl,        /**< [in] [optional] if more information can be provided,
                                                                           use this to inform the user.*/
                      const std::source_location &loc /**< [in] [optional] source location */
                      = std::source_location::current() )
{
    if( verboseT::level > 0 )
    {
        std::cerr << error_message<verboseT>( code, expl, loc ) << '\n';
    }

    return code;
}

/** \brief  Specialization of \ref error_report for \ref verbose::o
 *
 * \ingroup error_handling
 */
template <>
error_t error_report<verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Print a report to stderr given an mxlib \ref error_t \p code and return the \p code.
/** What is printed depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the provided \ref error_t \p code
 *
 * \ingroup error_handling
 */
template <class verboseT>
error_t error_report( const error_t &code,            /**< [in] is an mx::error_t error code*/
                      const std::source_location &loc /**< [in] [optional] source location */
                      = std::source_location::current() )
{
    if( verboseT::level > 0 )
    {
        std::cerr << error_message<verboseT>( code, loc ) << '\n';
    }

    return code;
}

/** \brief  Specialization of \ref error_report(const error_t &, const std::source_location&) "error_report"
 * for \ref verbose::o
 *
 * \ingroup error_handling
 */
template <>
error_t error_report<verbose::o>( const error_t &code, const std::source_location &loc );

} // namespace mx

#endif //__mxError__
