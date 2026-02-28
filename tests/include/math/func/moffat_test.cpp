/** \file moffat_test.cpp
 */
#include "../../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../../include/math/func/moffat.hpp"/// compiling 1D Moffat function


/// compiling 1D Moffat function
/**
 * \ingroup moffat_unit_tests
 */
TEST_CASE( "compiling 1D Moffat function", "[math::func::moffat]" )
{
    SECTION( "the 1D Moffat function" )
    {
        SECTION( "compiling" )
        {
            double mv = mx::math::func::moffat<double>( 0., 0., 1., 0., 1., 1. );
            REQUIRE_THAT( mv, Catch::Matchers::WithinAbs( 1.0, 1e-12 ) );
        }
    }
}/// compiling 2D Moffat function


/// compiling 2D Moffat function
/**
 * \ingroup moffat_unit_tests
 */
TEST_CASE( "compiling 2D Moffat function", "[math::func::moffat]" )
{
    SECTION( "the 2D Moffat function" )
    {
        SECTION( "compiling" )
        {
            double mv = mx::math::func::moffat2D<double>( 0., 0., 0., 1., 0., 0., 1., 1. );
            REQUIRE_THAT( mv, Catch::Matchers::WithinAbs( 1.0, 1e-12 ) );
        }
    }
}/// compiling Moffat FWHM


/// compiling Moffat FWHM
/**
 * \ingroup moffat_unit_tests
 */
TEST_CASE( "compiling Moffat FWHM", "[math::func::moffat]" )
{
    SECTION( "the Moffat FWHM" )
    {
        SECTION( "compiling" )
        {
            double mv = mx::math::func::moffatFWHM<double>( 1., 1. );
            REQUIRE_THAT( mv, Catch::Matchers::WithinAbs( 2.0, 1e-12 ) );
        }
    }
}
