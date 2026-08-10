/** \file fitsHeaderCard_test.cpp
 * \brief Tests FITS header-card parsing and formatting.
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ioutils/fits/fitsHeaderCard.hpp"
using namespace mx::fits;

namespace mx
{
namespace unitTest
{
namespace fitsTest
{
namespace fitsHeaderCardTest
{

/// Verify conversion of types
/**
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "fitsHeaderCard setting types", "[ioutils::fits::fitsHeaderCard]" )
{
    SECTION( "a fitsHeaderCard constructed with char type" )
    {
        SECTION( "setting a char from a string" )
        {
            fitsHeaderCard fhc( "KEYWORD", "39", fitsType<char>(), "this comment" );

            REQUIRE( fhc.keyword() == "KEYWORD" );
            std::string s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.type() == fitsType<char>() );
            REQUIRE( fhc.valueGood() == false );
            REQUIRE( fhc.valueStrGood() == true );
            REQUIRE( fhc.comment() == "this comment" );

            char c = fhc.value<char>();
            REQUIRE( c == 39 );
            REQUIRE( fhc.valueGood() == true );

            int i = fhc.value<int>();
            REQUIRE( i == 39 );
            REQUIRE( fhc.type() == fitsType<char>() );

            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == true );

            fhc.type( fitsType<int>() );
            REQUIRE( fhc.type() == fitsType<int>() );
            REQUIRE( fhc.Int() == 39 );
            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == false );

            s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.valueStrGood() == true );
        }
        SECTION( "setting a char from a char" )
        {
            fitsHeaderCard fhc( "KEYWORD", static_cast<char>( 39 ), "this comment" );

            REQUIRE( fhc.keyword() == "KEYWORD" );
            REQUIRE( fhc.type() == fitsType<char>() );
            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == false );
            REQUIRE( fhc.comment() == "this comment" );

            char c = fhc.value<char>();
            REQUIRE( c == 39 );
            REQUIRE( fhc.valueGood() == true );

            // Test that valueStr stayed false and then read it
            REQUIRE( fhc.valueStrGood() == false );
            std::string s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.valueStrGood() == true );

            int i = fhc.value<int>();
            REQUIRE( i == 39 );
            REQUIRE( fhc.type() == fitsType<char>() );

            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == true );

            mx::error_t errc = fhc.type( fitsType<int>() );
            REQUIRE( !errc );

            REQUIRE( fhc.type() == fitsType<int>() );
            errc = mx::error_t::error;
            REQUIRE( fhc.Int( &errc ) == 39 );
            REQUIRE( errc == mx::error_t::noerror );

            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == false );

            s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.valueStrGood() == true );
        }
    }
}

/// Removing white space around string values
/**
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Removing white space around string values", "[ioutils::fits::fitsHeaderCard]" )
{
    // clang-format off
    #ifdef MXLIBTEST_DOXYGEN_REF
        fitsHeaderCard fhc;
        mx::error_t errc;
        fhc.value( mx::meta::tagT<std::string>(), errc );
    #endif
    // clang-format on

    SECTION( "typical case" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'simple         '", "comment" );
        REQUIRE( fhc.String() == "simple         " );
    }

    SECTION( "no ' and no space" )
    {
        fitsHeaderCard fhc( "KEYTEST", "simple", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "no spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'simple'", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "no '" )
    {
        fitsHeaderCard fhc( "KEYTEST", "simple         ", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "space at beginning" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'   simple         '", "comment" );
        REQUIRE( fhc.String() == "   simple         " );
    }

    SECTION( "spaces at beginning, no spaces at end" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'     simple'", "comment" );
        REQUIRE( fhc.String() == "     simple" );
    }

    SECTION( "spaces at beginning, no '" )
    {
        fitsHeaderCard fhc( "KEYTEST", "   simple         ", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "empty" )
    {
        fitsHeaderCard fhc( "KEYTEST", "", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "one '" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "two ''" )
    {
        fitsHeaderCard fhc( "KEYTEST", "''", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "two '', spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'  '", "comment" );
        REQUIRE( fhc.String() == "  " );
    }

    SECTION( "spaces only" )
    {
        fitsHeaderCard fhc( "KEYTEST", "    ", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "' part of value" )
    {
        fitsHeaderCard fhc( "KEYTEST", "z'", "comment" );
        REQUIRE( fhc.String() == "z'" );
    }

    SECTION( "' part of value at end with spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'z'   '", "comment" );
        REQUIRE( fhc.String() == "z'   " );
    }

    SECTION( "' part of value at beginning with spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "''z   '", "comment" );
        REQUIRE( fhc.String() == "'z   " );
    }

    SECTION( "spaces before and after '' with ' in value surrounded by spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "   '  'z   '  ", "comment" );
        REQUIRE( fhc.String() == "  'z   " );
    }
}

/// CONTINUE-ing a card
/**
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "CONTINUE-ing a card", "[ioutils::fits::fitsHeaderCard]" )
{
    SECTION("normal CONTINUE, trailing space")
    {
        fitsHeaderCard fhc( "KEYTEST", "'one,two,three,&'", "" );
        fitsHeaderCard fhcc( "CONTINUE", "", "'four,five,six' / the comment" );

        REQUIRE(fhc.appendContinue(fhcc) == error_t::noerror);

        REQUIRE(fhc.keyword() == "KEYTEST");
        REQUIRE(fhc.String() == "one,two,three,four,five,six");
        REQUIRE(fhc.comment() == "the comment");

    }

    SECTION("normal CONTINUE, trailing space")
    {
        fitsHeaderCard fhc( "KEYTEST", "'one two three &'", "" );
        fitsHeaderCard fhcc( "CONTINUE", "", "'four five six' / the comment" );

        REQUIRE(fhc.appendContinue(fhcc) == error_t::noerror);

        REQUIRE(fhc.keyword() == "KEYTEST");
        REQUIRE(fhc.String() == "one two three four five six");
        REQUIRE(fhc.comment() == "the comment");

    }

    SECTION("normal CONTINUE, leading space")
    {
        fitsHeaderCard fhc( "KEYTEST", "'one two three&'", "" );
        fitsHeaderCard fhcc( "CONTINUE", "", "' four five six' / the comment" );

        REQUIRE(fhc.appendContinue(fhcc) == error_t::noerror);

        REQUIRE(fhc.keyword() == "KEYTEST");
        REQUIRE(fhc.String() == "one two three four five six");
        REQUIRE(fhc.comment() == "the comment");

    }
}

} // namespace fitsHeaderCardTest
} // namespace fitsTest
} // namespace unitTest
} // namespace mx
