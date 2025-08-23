/** \file fitsHeaderCard_test.cpp
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ioutils/fits/fitsHeaderCard.hpp"
using namespace mx::fits;

namespace unitTest
{
namespace ioutilsTest
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

            fhc.type( fitsType<int>() );
            REQUIRE( fhc.type() == fitsType<int>() );
            mx::error_t errc = mx::error_t::error;
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

} // namespace fitsHeaderCardTest
} // namespace fitsTest
} // namespace ioutilsTest
} // namespace unitTest
