/** \file fitsFile_test.cpp
 */
#include "../../../catch2/catch.hpp"

using namespace Catch::Matchers;

#include "../../../../include/ioutils/fits/fitsHeader.hpp"
using namespace mx::fits;

#include "../../../../include/improc/eigenImage.hpp"

namespace mx
{
namespace unitTest
{
namespace fitsTest
{
namespace fitsHeaderTest
{
/**
 *
 *//// Test functionality of headersToValues

/// Test functionality of headersToValues
/**
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Test functionality of headersToValues", "[ioutils::fits::fitsHeader]" )
{
// to force doxygen linking
#ifdef MXLIB_DOXYGEN_PROTECTED_REF
#endif

    SECTION( "a vector of floating point values, no errors" )
    {
        std::vector<fitsHeader<mx::verbose::d>> fhs( 5 );

        fhs[0].append<float>( "KEYTEST", 1, "test comment" );
        fhs[1].append<float>( "KEYTEST", 2, "test comment" );
        fhs[2].append<float>( "KEYTEST", 3, "test comment" );
        fhs[3].append<float>( "KEYTEST", 4, "test comment" );
        fhs[4].append<float>( "KEYTEST", 5, "test comment" );

        std::vector<float> vals;
        std::vector<size_t> bad;
        mx::error_t errc = headersToValues( vals, bad, fhs, "KEYTEST" );

        REQUIRE( !errc );

        REQUIRE( vals[0] == 1 );
        REQUIRE( vals[1] == 2 );
        REQUIRE( vals[2] == 3 );
        REQUIRE( vals[3] == 4 );
        REQUIRE( vals[4] == 5 );
    }

    SECTION( "a vector of floating point values, no errors" )
    {
        std::vector<fitsHeader<mx::verbose::d>> fhs( 5 );

        fhs[0].append<float>( "KEYTEST", 1, "test comment" );
        fhs[1].append<float>( "KEYTEST", 2, "test comment" );
        // fhs[2] is not set
        fhs[3].append<float>( "KEYTEST", 4, "test comment" );
        fhs[4].append<float>( "KEYTEST", 5, "test comment" );

        std::vector<float> vals;
        std::vector<size_t> bad;
        mx::error_t errc = headersToValues( vals, bad, fhs, "KEYTEST" );

        REQUIRE( !!errc );

        REQUIRE( vals[0] == 1 );
        REQUIRE( vals[1] == 2 );
        REQUIRE( vals[2] == std::numeric_limits<float>::max() );
        REQUIRE( vals[3] == 4 );
        REQUIRE( vals[4] == 5 );
    }
}

} // namespace fitsHeaderTest
} // namespace fitsTest
} // namespace unitTest
} // namespace mx
