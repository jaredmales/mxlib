/** \file exception_test.cpp
 * \brief Tests mxlib exception construction, inspection, and nested-exception reporting.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/error/exception.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/** \defgroup exception_unit_tests Exception Unit Tests
 * \ingroup error_unit_tests
 */

namespace unitTest::error_exception_test
{

/** \brief Verifies all mx::exception constructors and accessors used by hciReduce error paths.
 *
 * \ingroup exception_unit_tests
 */
TEST_CASE( "exception preserves error context", "[error::exception]" )
{
    const int locationLine = __LINE__ + 1;
    const mx::exception<mx::verbose::vv> locationOnly;

    REQUIRE( locationOnly.code() == mx::error_t::exception );
    REQUIRE( locationOnly.message().empty() );
    REQUIRE( locationOnly.line() == locationLine );
    REQUIRE( locationOnly.file_name().find( "exception_test.cpp" ) != std::string::npos );
    REQUIRE( std::string( locationOnly.what() ).find( "exception" ) != std::string::npos );

    const mx::exception<mx::verbose::vv> messageOnly( "custom explanation" );
    REQUIRE( messageOnly.code() == mx::error_t::exception );
    REQUIRE( messageOnly.message() == "custom explanation" );
    REQUIRE( std::string( messageOnly.what() ).find( "custom explanation" ) != std::string::npos );

    const mx::exception<mx::verbose::vv> codeAndMessage( mx::error_t::invalidarg, "invalid input" );
    REQUIRE( codeAndMessage.code() == mx::error_t::invalidarg );
    REQUIRE( codeAndMessage.message() == "invalid input" );
    REQUIRE( std::string( codeAndMessage.what() ).find( "invalid input" ) != std::string::npos );

    const mx::exception<mx::verbose::vv> codeOnly( mx::error_t::sizeerr );
    REQUIRE( codeOnly.code() == mx::error_t::sizeerr );
    REQUIRE( codeOnly.message().empty() );
    REQUIRE( std::string( codeOnly.what() ).find( "sizeerr" ) != std::string::npos );
}

/** \brief Verifies recursive extraction of standard and non-standard nested exceptions.
 *
 * \ingroup exception_unit_tests
 */
TEST_CASE( "unwind_exceptions extracts nested messages", "[error::exception]" )
{
    std::vector<std::string> messages;

    try
    {
        try
        {
            throw std::runtime_error( "inner failure" );
        }
        catch( ... )
        {
            std::throw_with_nested( mx::exception<mx::verbose::vv>( mx::error_t::exception, "outer failure" ) );
        }
    }
    catch( const std::exception &exception )
    {
        mx::unwind_exceptions( messages, exception );
    }

    REQUIRE( messages.size() == 2 );
    REQUIRE( messages[0].find( "outer failure" ) != std::string::npos );
    REQUIRE( messages[1] == "inner failure" );

    messages.clear();
    try
    {
        try
        {
            throw 17;
        }
        catch( ... )
        {
            std::throw_with_nested(
                mx::exception<mx::verbose::vv>( mx::error_t::exception, "non-standard nested failure" ) );
        }
    }
    catch( const std::exception &exception )
    {
        mx::unwind_exceptions( messages, exception );
    }

    REQUIRE( messages.size() == 1 );
    REQUIRE( messages[0].find( "non-standard nested failure" ) != std::string::npos );
}

/** \brief Verifies the ladder format produced by mx::print_exceptions.
 *
 * \ingroup exception_unit_tests
 */
TEST_CASE( "print_exceptions formats nested messages", "[error::exception]" )
{
    std::vector<std::string> messages{ "outer failure", "inner failure" };
    std::ostringstream captured;
    std::streambuf *original = std::cerr.rdbuf( captured.rdbuf() );

    mx::print_exceptions( messages, "processing failed" );

    std::cerr.rdbuf( original );
    REQUIRE( captured.str() == "processing failed:\n  inner failure\n  |-->outer failure\n" );
}

} // namespace unitTest::error_exception_test
