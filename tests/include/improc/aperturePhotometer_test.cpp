/** \file aperturePhotometer_test.cpp
 * \brief Tests for the aperturePhotometer production API.
 */

#include "../../catch2/catch.hpp"

#include "../../../include/improc/aperturePhotometer.hpp"

/** \defgroup aperturePhotometer_unit_tests aperturePhotometer Unit Tests
 * \ingroup improc_unit_tests
 */

namespace unitTest::improc_aperturePhotometer_test
{

/** \brief Verifies cumulative photometry through mx::improc::aperturePhotometer::cumPhot.
 *
 * \ingroup aperturePhotometer_unit_tests
 * \todo Add mask and background-annulus cases for mx::improc::aperturePhotometer::cumPhot.
 */
TEST_CASE( "aperturePhotometer accumulates a centered image", "[improc::aperturePhotometer]" )
{
    mx::improc::aperturePhotometer<double> photometer;
    mx::improc::eigenImage<double> image( 3, 3 );
    std::vector<double> cumulative;

    image.setOnes();

    REQUIRE( photometer.resize( image.rows(), image.cols() ) == 0 );
    REQUIRE( photometer.cumPhot( cumulative, image ) == 0 );
    REQUIRE( cumulative.size() == 3 );
    CHECK( cumulative[0] == 1.0 );
    CHECK( cumulative[1] == 5.0 );
    CHECK( cumulative[2] == 9.0 );
}

} // namespace unitTest::improc_aperturePhotometer_test
