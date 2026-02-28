/** \file exception_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <sstream>

#include "../../../include/error/exception.hpp"

namespace unitTest
{
namespace errorTest
{
namespace exceptionTest
{

/// Test exception constructors and accessors
/**
 * \ingroup error_exception_unit_tests
 */
TEST_CASE( "Test exception constructors and accessors", "[error::exception]" )
{
    #if MXLIBTEST_DOXYGEN_REF // doxygen links
    mx::exception<> ex0;
    mx::exception<> ex1( "test message" );
    mx::exception<> ex2( mx::error_t::error, "test message" );
    mx::exception<> ex3( mx::error_t::error );
    #endif

    SECTION( "default constructor" )
    {
        mx::exception<> ex;

        REQUIRE( ex.code() == mx::error_t::exception );
        REQUIRE( ex.message() == "" );
        REQUIRE( std::string( ex.what() ) != "" );
        REQUIRE( ex.line() > 0 );
        REQUIRE( ex.file_name().find( "exception_test.cpp" ) != std::string::npos );
    }

    SECTION( "message constructor" )
    {
        mx::exception<> ex( "test message" );

        REQUIRE( ex.code() == mx::error_t::exception );
        REQUIRE( ex.message() == "test message" );
        REQUIRE( std::string( ex.what() ).find( "test message" ) != std::string::npos );
    }

    SECTION( "code+message constructor" )
    {
        mx::exception<> ex( mx::error_t::sizeerr, "size failed" );

        REQUIRE( ex.code() == mx::error_t::sizeerr );
        REQUIRE( ex.message() == "size failed" );
        REQUIRE( std::string( ex.what() ).find( "size failed" ) != std::string::npos );
        REQUIRE( std::string( ex.what() ).find( mx::errorName( mx::error_t::sizeerr ) ) != std::string::npos );
    }

    SECTION( "code-only constructor" )
    {
        mx::exception<> ex( mx::error_t::notimpl );

        REQUIRE( ex.code() == mx::error_t::notimpl );
        REQUIRE( ex.message() == "" );
        REQUIRE( std::string( ex.what() ).find( mx::errorName( mx::error_t::notimpl ) ) != std::string::npos );
    }
}

/// Test exception location forwarding
/**
 * \ingroup error_exception_unit_tests
 */
TEST_CASE( "Test exception location forwarding", "[error::exception]" )
{
    const std::source_location loc = std::source_location::current();
    mx::exception<> ex( mx::error_t::dirnotfound, "missing dir", loc );

    REQUIRE( ex.file_name() == loc.file_name() );
    REQUIRE( ex.line() == static_cast<int>( loc.line() ) );
}

/// Test unwind_exceptions
/**
 * \ingroup error_exception_unit_tests
 */
TEST_CASE( "Test unwind_exceptions", "[error::exception]" )
{
    #if MXLIBTEST_DOXYGEN_REF // doxygen link
    std::vector<std::string> whats;
    mx::exception<> ex;
    mx::unwind_exceptions( whats, ex );
    #endif

    SECTION( "std::exception-only nested chain recurses through all levels" )
    {
        try
        {
            try
            {
                try
                {
                    throw mx::exception<>( "inner exception" );
                }
                catch( ... )
                {
                    std::throw_with_nested( mx::exception<>( "middle exception" ) );
                }
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<>( "outer exception" ) );
            }
        }
        catch( const std::exception &e )
        {
            std::vector<std::string> whats;
            mx::unwind_exceptions( whats, e );

            REQUIRE( whats.size() == 3 );
            REQUIRE( whats[0].find( "outer exception" ) != std::string::npos );
            REQUIRE( whats[1].find( "middle exception" ) != std::string::npos );
            REQUIRE( whats[2].find( "inner exception" ) != std::string::npos );
        }
        catch( ... )
        {
            FAIL( "Expected std::exception chain" );
        }
    }

    SECTION( "mixed nested chain ends in non-std::exception" )
    {
        try
        {
            try
            {
                try
                {
                    throw 42;
                }
                catch( ... )
                {
                    std::throw_with_nested( mx::exception<>( "middle exception" ) );
                }
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<>( "outer exception" ) );
            }
        }
        catch( const std::exception &e )
        {
            const auto *nested = dynamic_cast<const std::nested_exception *>( &e );
            REQUIRE( nested != nullptr );
            try
            {
                nested->rethrow_nested();
                FAIL( "Expected nested exception to be rethrown" );
            }
            catch( const std::exception &middle )
            {
                const auto *middleNested = dynamic_cast<const std::nested_exception *>( &middle );
                REQUIRE( middleNested != nullptr );
                try
                {
                    middleNested->rethrow_nested();
                    FAIL( "Expected nested non-std::exception to be rethrown" );
                }
                catch( int value )
                {
                    REQUIRE( value == 42 );
                }
                catch( ... )
                {
                    FAIL( "Expected nested int exception" );
                }
            }
            catch( ... )
            {
                FAIL( "Expected nested middle std::exception" );
            }

            std::vector<std::string> whats;
            mx::unwind_exceptions( whats, e );

            REQUIRE( whats.size() == 2 );
            REQUIRE( whats[0].find( "outer exception" ) != std::string::npos );
            REQUIRE( whats[1].find( "middle exception" ) != std::string::npos );
        }
    }
}

/// Test print_exceptions
/**
 * \ingroup error_exception_unit_tests
 */
TEST_CASE( "Test print_exceptions", "[error::exception]" )
{
    #if MXLIBTEST_DOXYGEN_REF // doxygen link
    std::vector<std::string> whats{ "outer", "inner" };
    mx::print_exceptions( whats, "top message" );
    #endif

    std::vector<std::string> whats{ "outer", "inner" };

    std::stringstream ss;
    auto *old = std::cerr.rdbuf( ss.rdbuf() );
    mx::print_exceptions( whats, "top message" );
    std::cerr.rdbuf( old );

    std::string out = ss.str();

    REQUIRE( out.find( "top message" ) != std::string::npos );
    REQUIRE( out.find( "inner" ) != std::string::npos );
    REQUIRE( out.find( "outer" ) != std::string::npos );
    REQUIRE( out.find( "|-->" ) != std::string::npos );
}

} // namespace exceptionTest
} // namespace errorTest
} // namespace unitTest
