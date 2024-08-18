/** \file geo_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/geo.hpp"

/** Verify compilation and calculations of math::angleMod.
 * Tests that various angle modulos are calculated correctly for both radians and degrees.
 *
 * \anchor tests_math_geo_angleMod
 */
SCENARIO( "Verify compilation and calculations of math::angleMod", "[math::angleMod]" )
{
    GIVEN( "angles in degrees" )
    {
        WHEN( "positive angle, no changes" )
        {
            double q = mx::math::angleMod<mx::math::degreesT<double>>( 43.2 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 43.2 ) );
        }

        WHEN( "positive angle, no changes" )
        {
            double q = mx::math::angleMod<mx::math::degreesT<double>>( 353.2 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 353.2 ) );
        }

        WHEN( "positive angle, exactly 360" )
        {
            double q = mx::math::angleMod<mx::math::degreesT<double>>( 360.0 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 0.0 ) );
        }

        WHEN( "positive angle, mod needed" )
        {
            double q = mx::math::angleMod<mx::math::degreesT<double>>( 362.0 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 2.0 ) );
        }
    }

    GIVEN( "angles in radians" )
    {
        WHEN( "positive angle, no changes" )
        {
            double q = mx::math::angleMod<mx::math::radiansT<double>>( mx::math::dtor( 43.2 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor( 43.2 ) ) );
        }

        WHEN( "positive angle, no changes" )
        {
            double q = mx::math::angleMod<mx::math::radiansT<double>>( mx::math::dtor( 353.2 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor( 353.2 ) ) );
        }

        WHEN( "positive angle, exactly 2pi" )
        {
            double q = mx::math::angleMod<mx::math::radiansT<double>>( mx::math::dtor( 360.0 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor( 0.0 ) ) );
        }

        WHEN( "positive angle, mod needed" )
        {
            double q = mx::math::angleMod<mx::math::radiansT<double>>( mx::math::dtor( 362.0 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor( 2.0 ) ) );
        }
    }
}

/** Verify compilation and calculations of math::angleDiff.
 * Tests that various angle differences are calculated correctly for both radians and degrees.
 *
 * \anchor tests_math_geo_angleDiff
 */
SCENARIO( "Verify compilation and calculations of math::angleDiff", "[math::angleDiff]" )
{
    GIVEN( "angles in degrees" )
    {
        WHEN( "positive, first angle is 0, not crossing 180/360" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 0.0, 43.2 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 43.2 ) );
        }

        WHEN( "negative, second angle is 0, not crossing 180/360" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 43.2, 0.0 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( -43.2 ) );
        }

        WHEN( "positive, first angle is 360, not crossing 180/360" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 360.0, 43.2 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 43.2 ) );
        }

        WHEN( "negative, second angle is 3600, not crossing 180/360" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 43.2, 360.0 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( -43.2 ) );
        }

        WHEN( "positive, crossing 360" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 340.0, 23.2 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 43.2 ) );
        }

        WHEN( "negative, crossing 180/360" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 23.2, 340.0 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( -43.2 ) );
        }

        WHEN( "positive, crossing 180" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 160.0, 206.2 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( 46.2 ) );
        }

        WHEN( "negative, crossing 180" )
        {
            double q = mx::math::angleDiff<mx::math::degreesT<double>>( 206.2, 160.0 );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( -46.2 ) );
        }
    }

    GIVEN( "angles in radians" )
    {
        WHEN( "positive, first angle is 0, not crossing pi/2pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 0.0 ),
                                                                        mx::math::dtor<double>( 43.2 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( 43.2 ) ) );
        }

        WHEN( "negative, second angle is 0, not crossing pi/2pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 43.2 ),
                                                                        mx::math::dtor<double>( 0.0 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( -43.2 ) ) );
        }

        WHEN( "positive, first angle is 360, not crossing pi/2pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 360.0 ),
                                                                        mx::math::dtor<double>( 43.2 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( 43.2 ) ) );
        }

        WHEN( "negative, second angle is 3600, not crossing pi/2pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 43.2 ),
                                                                        mx::math::dtor<double>( 360.0 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( -43.2 ) ) );
        }

        WHEN( "positive, crossing 2pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 340.0 ),
                                                                        mx::math::dtor<double>( 23.2 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( 43.2 ) ) );
        }

        WHEN( "negative, crossing 2pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 23.2 ),
                                                                        mx::math::dtor<double>( 340.0 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( -43.2 ) ) );
        }

        WHEN( "positive, crossing pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 160.0 ),
                                                                        mx::math::dtor<double>( 206.2 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( 46.2 ) ) );
        }

        WHEN( "negative, crossing pi" )
        {
            double q = mx::math::angleDiff<mx::math::radiansT<double>>( mx::math::dtor<double>( 206.2 ),
                                                                        mx::math::dtor<double>( 160.0 ) );

            REQUIRE_THAT( q, Catch::Matchers::WithinRel( mx::math::dtor<double>( -46.2 ) ) );
        }
    }
}
