/** \file exception.hpp
 * \brief The mxlib exception class.
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

#ifndef error_exception_hpp
#define error_exception_hpp

#include <vector>

#include "error.hpp"

namespace mx
{

/// Augments an exception with the source file and line
/**
 * \tparam baseexcept is the base class exception which takes a string as constructor argument
 */
template <class verboseT = verbose::d>
class exception : public std::exception
{
  protected:
    std::string m_what;                   ///< The full what message (message + file information).

    error_t m_code{ error_t::exception }; ///< The error_t code

    std::string m_message;                ///< The explanatory message

    std::source_location m_location;

  public:
    /// Constructor with location only.
    /**
     */
    explicit exception( const std::source_location loc = /**< [in] [opt] the location where this was thrown*/
                        std::source_location::current() )
        : m_location{ loc }
    {
        m_what = error_message<verboseT>( m_code, m_location );
    }

    /// Constructor with message.
    /**
     */
    explicit exception( const std::string &msg,          /**< [in] the error description message) */
                        const std::source_location loc = /**< [in] [opt] the location where this was thrown*/
                        std::source_location::current() )
        : m_message{ msg }, m_location{ loc }
    {
        m_what = error_message<verboseT>( m_code, m_message, m_location );
    }

    /// Constructor with message and cod
    /**
     */
    exception( error_t code,                    /**< [in] a descriptive error code */
               const std::string &msg,          /**< [in] the error description (what message) */
               const std::source_location loc = /**< [in] [opt] the location where this was thrown*/
               std::source_location::current() )
        : m_code{ code }, m_message{ msg }, m_location{ loc }
    {
        m_what = error_message<verboseT>( m_code, m_message, m_location );
    }

    /// Constructor with code
    /** The message is filled in using \ref errorMessage.
     *
     * The what() message becomes "code_message (code) [file line]".
     */
    exception( error_t code,                    /**< [in] a descriptive error code */
               const std::source_location loc = /**< [in] [opt] the location where this was thrown*/
               std::source_location::current() )
        : m_code{ code }, m_location{ loc }
    {
        m_what = error_message<verboseT>( m_code, m_location );
    }

    /// Get the what string
    /** \returns the value of m_what.c_str()
     *
     */
    virtual const char *what() const noexcept
    {
        return m_what.c_str();
    }

    /// Get the message
    /** \returns the value of m_message
     *
     */
    const std::string &message() const
    {
        return m_message;
    }

    /// Get the source file
    /** \returns the value of m_location.file_name()
     *
     */
    const std::string file_name() const
    {
        return m_location.file_name();
    }

    /// Get the source line
    /** \returns the value of m_location.line()
     *
     */
    int line() const
    {
        return m_location.line();
    }

    /// Get the error code
    /** \returns the value of m_code
     *
     */
    error_t code() const
    {
        return m_code;
    }
};

/// Extract the explanatory string of nested exceptions, placing them in a vector.
/** If the exception is nested, this recurses to extract the explanatory string of the
 * next exception it holds.
 */
inline
void unwind_exceptions( std::vector<std::string> &whats, /**< [out] the vector of what messages*/
                        const std::exception &e          /**< [in] the exception */
)
{
    whats.push_back( e.what() );

    try
    {
        std::rethrow_if_nested( e );
    }
    catch( const std::exception &nestedException )
    {
        unwind_exceptions( whats, nestedException );
    }
    catch( ... )
    {
    }
}

/// Print nested exceptions to stderr, from the earliest to latest
/** Formats in a ladder
 *
*/
inline
void print_exceptions( std::vector<std::string> &whats, /**< [out] a vector of what messages*/
                       const std::string &message =     /**< [in] [optional] the top message to print */
                       "exception(s) thrown" )
{
    std::cerr << message << ":\n";
    std::cerr << "  " << whats.back() << '\n';
    for( size_t n = 1; n < whats.size(); ++n )
    {
        std::cerr << std::string( 2, ' ' ) << std::string( ( n - 1 ) * 4, ' ' ) << "|-->" << whats[whats.size() - 1 - n]
                  << '\n';
    }
}

} // namespace mx

#endif // error_exception_hpp
