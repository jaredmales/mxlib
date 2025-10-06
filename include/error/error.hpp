/** \file error.hpp
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

#ifndef error_error_hpp
#define error_error_hpp

#include <source_location>
#include <iostream>
#include <string>

#include "error_t.hpp"

namespace mx
{

/// Check if an error_t code is an error
/** Treats errors as `true/false`, where `false` means equal to `error_t::noerror` and all other values are `true`.
 *  So `not` returns `true` if \p errc is equal to `error_t::noerror`.  So one can
 *  use `!!` to check if an code is an error.  Example:
 *  \code
 *   error_t errc = error_t::noerror;
 *   if(!!errc) // will be false, no error
 *   {
 *       //handle errors
 *   }
 *
 *   errc = error_t::error;
 *   if(!!errc) // will be true, error
 *   {
 *       //handle errors
 *   }
 * \endcode
 *
 * \ingroup error_handling
 */
inline bool operator!( const error_t &errc /**< [in] the error_t code to check */ )
{
    return ( errc == error_t::noerror );
}

/// Check if an error_t code is an error
/** Treats errors as `true/false`, where `false` means equal to `error_t::noerror` and all other values are `true`.
 *  Example:
 *  \code
 *   error_t errc = error_t::noerror;
 *   if(errc == true) // will be false, no error
 *   {
 *       //handle errors
 *   }
 *
 *   errc = error_t::error;
 *   if(errc == true) // will be true, error
 *   {
 *       //handle errors
 *   }
 * \endcode
 *
 * \ingroup error_handling
 */
inline bool operator==( const error_t &errc, /**< [in] the error_t code to check */
                        bool tf              /**< [in] the value to compare to */
)
{
    return ( ( errc != error_t::noerror ) == tf );
}

/// Check if an error_t code is not an error
/** Treats errors as `true/false`, where `false` means equal to `error_t::noerror` and all other values are `true`.
 *  Example:
 *  \code
 *   error_t errc = error_t::noerror;
 *   if(errc != false) // will be false, no error
 *   {
 *       //handle errors
 *   }
 *
 *   errc = error_t::error;
 *   if(errc != false) // will be true, error
 *   {
 *       //handle errors
 *   }
 * \endcode
 *
 * \ingroup error_handling
 */
inline bool operator!=( const error_t &errc, /**< [in] the error_t code to check */
                        bool tf              /**< [in] the value to compare to */
)
{
    return ( ( errc != error_t::noerror ) != tf );
}

/// Check if an error_t code is an error
/** Treats errors as `true/false`, where `false` means equal to `error_t::noerror` and all other values are `true`.
 *
 * \ingroup error_handling
 */
inline bool isError( const error_t &errc /**< [in] the error_t code to check */ )
{
    return ( errc != error_t::noerror );
}

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
 * \endverbatim:
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

#ifndef MXLIB_DEFAULT_VERBOSITY
    #define MXLIB_DEFAULT_VERBOSITY vv
#endif

/// The default verbosity
typedef MXLIB_DEFAULT_VERBOSITY d;

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
template <class verboseT=verbose::d>
std::string mxlib_error_message( const error_t &code,     /**< [in] is an mx::error_t error code*/
                                 const std::string &expl, /**< [in] [opt] if more information can be provided,
                                                                               use this to inform the user.*/
                                 const std::source_location &loc /**< [in] [opt] source location */
                                 = std::source_location::current() );

/// Specialization of \ref mxlib_error_message for \ref verbose::o
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<mx::verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

// note: keep the mx:: in mx::verbose:: so doxygen matches things up in test cases

/// Specialization of \ref mxlib_error_message for \ref verbose::v
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<mx::verbose::v>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref mxlib_error_message for \ref verbose::vv
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<mx::verbose::vv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref mxlib_error_message for \ref verbose::vvv
/**
 * \ingroup error_internal
 */
template <>
std::string
mxlib_error_message<mx::verbose::vvv>( const error_t &code, const std::string &expl, const std::source_location &loc );

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
template <class verboseT=verbose::d>
std::string mxlib_error_message( const error_t &code,            /**< [in] is an mx::error_t error code*/
                                 const std::source_location &loc /**< [in] [opt] source location */
                                 = std::source_location::current() );

/** \brief Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for
 * \ref verbose::o
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<mx::verbose::o>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for \ref verbose::v
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<mx::verbose::v>( const error_t &code, const std::source_location &loc );

/** \brief  Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for \ref verbose::vv
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<mx::verbose::vv>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref mxlib_error_message(const error_t &, const std::source_location&)
 * "mxlib_error_message" for \ref verbose::vvv
 *
 * \ingroup error_internal
 */
template <>
std::string mxlib_error_message<mx::verbose::vvv>( const error_t &code, const std::source_location &loc );

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
template <class verboseT=verbose::d>
error_t mxlib_error_report( const error_t &code,            /**< [in] is an mx::error_t error code*/
                            const std::string &expl,        /**< [in] [opt] if more information can be provided,
                                                                                 use this to inform the user.*/
                            const std::source_location &loc /**< [in] [opt] source location */
                            = std::source_location::current() )
{
    if( code == error_t::noerror ) // avoid output at any level
    {
        return code;
    }

    if constexpr( verboseT::level > 0 )
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
error_t
mxlib_error_report<mx::verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

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
template <class verboseT=verbose::d>
error_t mxlib_error_report( const error_t &code /**< [in] is an mx::error_t error code*/,
                            const std::source_location &loc /**< [in] [opt] source location */
                            = std::source_location::current() )
{
    if( code == error_t::noerror )
    {
        return code;
    }

    if constexpr( verboseT::level > 0 )
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
error_t mxlib_error_report<mx::verbose::o>( const error_t &code, const std::source_location &loc );

/** \brief Perform an error check, if an error occurs report it and return the error.  Does not return on no error.
 *
 * Scope protected so the error_t value does not interfere with other values.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param fxn is the function to call and check the return value of
 *
 * \ingroup error_internal
 *
 */
#define mxlib_error_check( fxn )                                                                                       \
    {                                                                                                                  \
        mx::error_t __mxlib_error_check_errc = fxn;                                                                    \
        if( __mxlib_error_check_errc != mx::error_t::noerror )                                                         \
        {                                                                                                              \
            return internal::mxlib_error_report<verboseT>( __mxlib_error_check_errc );                                 \
        }                                                                                                              \
    }

/** \brief Perform an error check, if an error occurs report it, and return the error code even if no error.
 *
 * Scope protected so the error_t value does not interfere with other values.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param fxn is the function to call and check the return value of
 *
 * \ingroup error_internal
 */
#define mxlib_error_return( fxn )                                                                                      \
    {                                                                                                                  \
        mx::error_t __mxlib_error_return_errc = fxn;                                                                   \
        if( __mxlib_error_return_errc != mx::error_t::noerror )                                                        \
        {                                                                                                              \
            return internal::mxlib_error_report<verboseT>( __mxlib_error_return_errc );                                \
        }                                                                                                              \
        return mx::error_t::noerror;                                                                                   \
    }

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
template <class verboseT = verbose::d>
std::string error_message( const error_t &code,            /**< [in] is an mx::error_t error code*/
                           const std::string &expl,        /**< [in] [opt] if more information can be provided,
                                                                     use this to inform the user.*/
                           const std::source_location &loc /**< [in] [opt] source location */
                           = std::source_location::current() );

/// Specialization of \ref error_message for \ref verbose::o
/**
 * \ingroup error_handling
 */
template <>
std::string
error_message<mx::verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref error_message for \ref verbose::v
/**
 * \ingroup error_handling
 */
template <>
std::string
error_message<mx::verbose::v>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref error_message for \ref verbose::vv
/**
 * \ingroup error_handling
 */
template <>
std::string
error_message<mx::verbose::vv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Specialization of \ref error_message for \ref verbose::vvv
/**
 * \ingroup error_handling
 */
template <>
std::string
error_message<mx::verbose::vvv>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Format a report given an mxlib \ref error_t \p code.
/** What is included depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the formatted message
 *
 * \ingroup error_handling
 */
template <class verboseT = verbose::vvv>
std::string error_message( const error_t &code,            /**< [in] is an mx::error_t error code*/
                           const std::source_location &loc /**< [in] [opt] source location */
                           = std::source_location::current() );

/** \brief Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for  \ref
 * verbose::o
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<mx::verbose::o>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for \ref
 * verbose::v
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<mx::verbose::v>( const error_t &code, const std::source_location &loc );

/** \brief Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for \ref
 * verbose::vv
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<mx::verbose::vv>( const error_t &code, const std::source_location &loc );

/** \brief  Specialization of \ref error_message(const error_t &, const std::source_location&) "error_message" for \ref
 * verbose::vvv
 *
 * \ingroup error_handling
 */
template <>
std::string error_message<mx::verbose::vvv>( const error_t &code, const std::source_location &loc );

/// Print a report to stderr given an mxlib \ref error_t \p code and explanation and return the \p code.
/** What is printed depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the provided \ref error_t \p code
 *
 * \ingroup error_handling
 */
template <class verboseT = verbose::vvv>
error_t error_report( const error_t &code,            /**< [in] is an mx::error_t error code*/
                      const std::string &expl,        /**< [in] [opt] if more information can be provided,
                                                                           use this to inform the user.*/
                      const std::source_location &loc /**< [in] [opt] source location */
                      = std::source_location::current() )
{
    if( code == error_t::noerror ) // avoid output at any level
    {
        return code;
    }

    if constexpr( verboseT::level > 0 )
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
error_t error_report<mx::verbose::o>( const error_t &code, const std::string &expl, const std::source_location &loc );

/// Print a report to stderr given an mxlib \ref error_t \p code and return the \p code.
/** What is printed depends on the verbosity level set by the template parameter
 *
 * \tparam verboseT sets the verbosity level based on its `level` member.
 *
 * \returns the provided \ref error_t \p code
 *
 * \ingroup error_handling
 */
template <class verboseT = verbose::vvv>
error_t error_report( const error_t &code,            /**< [in] is an mx::error_t error code*/
                      const std::source_location &loc /**< [in] [opt] source location */
                      = std::source_location::current() )
{
    if( code == error_t::noerror ) // avoid output at any level
    {
        return code;
    }

    if constexpr( verboseT::level > 0 )
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
error_t error_report<mx::verbose::o>( const error_t &code, const std::source_location &loc );

/** \brief Perform an error check on the output of a function.
 *
 * If an error occurs report it and return the error.  Do not return on no error.
 *
 * Scope protected so the error_t value does not interfere with other values.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param fxn is the function to call and check the return value of
 *
 * \ingroup error_handling
 */
#define mx_error_check( fxn )                                                                                          \
    {                                                                                                                  \
        mx::error_t __mxlib_error_check_errc = fxn;                                                                    \
        if( __mxlib_error_check_errc != mx::error_t::noerror )                                                         \
        {                                                                                                              \
            return error_report<verboseT>( __mxlib_error_check_errc );                                                 \
        }                                                                                                              \
    }

/** \brief Perform an error check on an error_t code.
 *
 * If an error is indicated report it and return the error.  Do not return on no error.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param errc is the code to check for an error value
 *
 * \ingroup error_handling
 */
#define mx_error_check_code( errc )                                                                                    \
    if( errc != mx::error_t::noerror )                                                                                 \
    {                                                                                                                  \
        return error_report<verboseT>( errc );                                                                         \
    }

/** \brief Perform an error check on the output of a function.
 *
 * If an error occurs report it and return an arbitrary type. Does not return on no error. This is intended to be used
 * in `int main()`, or any other case where the return value is something other than \ref mx::error_t.
 *
 * Scope protected so the \ref mx::error_t value does not interfere with other values.
 *
 * \note this requires that the verbosity template parameter is \b verboseT.  You may need to typedef it.
 *
 * \param fxn is the function to call and check the return value of (can also just be a code).
 *
 * \ingroup error_handling
 *
 */
#define mx_error_check_rv( fxn, rv )                                                                                   \
    {                                                                                                                  \
        mx::error_t __mxlib_error_check_errc = fxn;                                                                    \
        if( __mxlib_error_check_errc != mx::error_t::noerror )                                                         \
        {                                                                                                              \
            error_report<verboseT>( __mxlib_error_check_errc );                                                        \
            return rv;                                                                                                 \
        }                                                                                                              \
    }

/** \brief Perform an error check on an error_t code.
 *
 * If an error occurs report it and return an arbitrary type. Does not return on no error. This is intended to be used
 * in `int main()`, or any other case where the return value is something other than \ref mx::error_t.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param errc is the code to check for an error value
 *
 * \ingroup error_handling
 */
#define mx_error_check_code_rv( errc, rv )                                                                             \
    if( errc != mx::error_t::noerror )                                                                                 \
    {                                                                                                                  \
        error_report<verboseT>( errc );                                                                                \
        return rv;                                                                                                     \
    }

/** \brief Perform an error check on a function and return the result
 *
 * If an error occurs it is reported. If no error occurs error_t::noerror is returned.
 *
 * Scope protected so the error_t value does not interfere with other values.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param fxn is the function to call and check the return value of
 *
 * \ingroup error_handling
 */
#define mx_error_return( fxn )                                                                                         \
    {                                                                                                                  \
        mx::error_t __mx_error_return_errc = fxn;                                                                      \
        if( __mx_error_return_errc != mx::error_t::noerror )                                                           \
        {                                                                                                              \
            return error_report<verboseT>( __mx_error_return_errc );                                                   \
        }                                                                                                              \
        return mx::error_t::noerror;                                                                                   \
    }

/** \brief Perform an error check on an error_t code and return the result
 *
 * If an error occurs it is reported and returned. If no error occurs error_t::noerror is returned.
 *
 * \note this requires that the verbosity template parameter is \b verboseT
 *
 * \param errc is the function to call and check the return value of
 *
 * \ingroup error_handling
 */
#define mx_error_return_code( errc )                                                                                   \
    if( errc != mx::error_t::noerror )                                                                                 \
    {                                                                                                                  \
        return error_report<verboseT>( errc );                                                                         \
    }                                                                                                                  \
    return mx::error_t::noerror;

} // namespace mx

#endif // error_error_hpp
