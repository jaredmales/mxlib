/** \file moffat_test.cpp
 * \brief Tests Moffat profile functions.
 */
#include "../../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../../include/math/func/moffat.hpp"

/** compiling 1D Moffat function.
 *
 *
 */
/**
 * \ingroup moffat_unit_tests
 */
TEST_CASE( "compiling 1D Moffat function", "[math::func::moffat]" )
{
    GIVEN( "the 1D Moffat function" )
    {
        WHEN( "compiling" )
        {
            double mv = mx::math::func::moffat<double>( 0., 0., 1., 0., 1., 1. );
            REQUIRE_THAT( mv, Catch::Matchers::WithinAbs( 1.0, 1e-12 ) );
        }
    }
}

/** compiling 2D Moffat function.
 *
 *
 */
/**
 * \ingroup moffat_unit_tests
 */
TEST_CASE( "compiling 2D Moffat function", "[math::func::moffat]" )
{
    GIVEN( "the 2D Moffat function" )
    {
        WHEN( "compiling" )
        {
            double mv = mx::math::func::moffat2D<double>( 0., 0., 0., 1., 0., 0., 1., 1. );
            REQUIRE_THAT( mv, Catch::Matchers::WithinAbs( 1.0, 1e-12 ) );
        }
    }
}

/** compiling Moffat FWHM
 *
 *
 */
/**
 * \ingroup moffat_unit_tests
 */
TEST_CASE( "compiling Moffat FWHM", "[math::func::moffat]" )
{
    GIVEN( "the Moffat FWHM" )
    {
        WHEN( "compiling" )
        {
            double mv = mx::math::func::moffatFWHM<double>( 1., 1. );
            REQUIRE_THAT( mv, Catch::Matchers::WithinAbs( 2.0, 1e-12 ) );
        }
    }
}
