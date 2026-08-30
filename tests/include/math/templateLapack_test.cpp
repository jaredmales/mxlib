/** \file templateLapack_test.cpp
 *  \brief Tests the typed LAPACK wrappers.
 *  \ingroup templateLapack_unit_tests
 */

#include "../../catch2/catch.hpp"

#include "../../../include/math/templateLapack.hpp"

#include <array>
#include <cmath>
#include <limits>

/// Verify typed lamch wrappers for supported real precisions.
/**
 * \ingroup templateLapack_unit_tests
 */
TEST_CASE( "getting lamch values", "[templateLapack]" )
{
    GIVEN( "a precision" )
    {

        WHEN( "precision is float" )
        {
            // This is just a compilation test.
            float ch = mx::math::lamch<float>( 'E' );
            REQUIRE( ch > 0 );
        }
        WHEN( "precision is double" )
        {
            // This is just a compilation test.
            double ch = mx::math::lamch<double>( 'E' );
            REQUIRE( ch > 0 );
        }
    }
}

/// solving a rank-one secular equation with laed9
/**
 * The mx::math::laed9<float>() and mx::math::laed9<double>() wrappers are applied to a two-dimensional
 * diagonal-plus-rank-one eigenproblem with an analytic spectrum.
 *
 * \ingroup templateLapack_unit_tests
 */
TEST_CASE( "solving a rank-one secular equation with laed9", "[templateLapack]" )
{
    SECTION( "float" )
    {
        std::array<float, 2> eigenvalues{};
        std::array<float, 4> workspace{};
        std::array<float, 4> eigenvectors{};
        std::array<float, 2> poles{ 1.0F, 3.0F };
        std::array<float, 2> update{ 0.6F, 0.8F };

        const MXLAPACK_INT info = mx::math::laed9<float>( eigenvalues.data(),
                                                          workspace.data(),
                                                          eigenvectors.data(),
                                                          2,
                                                          1,
                                                          2,
                                                          2,
                                                          2,
                                                          2.0F,
                                                          poles.data(),
                                                          update.data(),
                                                          2 );

        const float tolerance = 64.0F * std::numeric_limits<float>::epsilon();
        REQUIRE( info == 0 );
        REQUIRE( std::abs( eigenvalues[0] - 1.4F ) <= tolerance );
        REQUIRE( std::abs( eigenvalues[1] - 4.6F ) <= tolerance );
        for( std::size_t column = 0; column < 2; ++column )
        {
            const float first = eigenvectors[2 * column];
            const float second = eigenvectors[2 * column + 1];
            const float firstResidual = 1.72F * first + 0.96F * second - eigenvalues[column] * first;
            const float secondResidual = 0.96F * first + 4.28F * second - eigenvalues[column] * second;
            CAPTURE( column, firstResidual, secondResidual );
            REQUIRE( std::hypot( firstResidual, secondResidual ) <= tolerance );
            REQUIRE( std::abs( std::hypot( first, second ) - 1.0F ) <= tolerance );
        }
        REQUIRE( std::abs( eigenvectors[0] * eigenvectors[2] + eigenvectors[1] * eigenvectors[3] ) <= tolerance );
    }

    SECTION( "double" )
    {
        std::array<double, 2> eigenvalues{};
        std::array<double, 4> workspace{};
        std::array<double, 4> eigenvectors{};
        std::array<double, 2> poles{ 1.0, 3.0 };
        std::array<double, 2> update{ 0.6, 0.8 };

        const MXLAPACK_INT info = mx::math::laed9<double>( eigenvalues.data(),
                                                           workspace.data(),
                                                           eigenvectors.data(),
                                                           2,
                                                           1,
                                                           2,
                                                           2,
                                                           2,
                                                           2.0,
                                                           poles.data(),
                                                           update.data(),
                                                           2 );

        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
        REQUIRE( info == 0 );
        REQUIRE( std::abs( eigenvalues[0] - 1.4 ) <= tolerance );
        REQUIRE( std::abs( eigenvalues[1] - 4.6 ) <= tolerance );
        for( std::size_t column = 0; column < 2; ++column )
        {
            const double first = eigenvectors[2 * column];
            const double second = eigenvectors[2 * column + 1];
            const double firstResidual = 1.72 * first + 0.96 * second - eigenvalues[column] * first;
            const double secondResidual = 0.96 * first + 4.28 * second - eigenvalues[column] * second;
            CAPTURE( column, firstResidual, secondResidual );
            REQUIRE( std::hypot( firstResidual, secondResidual ) <= tolerance );
            REQUIRE( std::abs( std::hypot( first, second ) - 1.0 ) <= tolerance );
        }
        REQUIRE( std::abs( eigenvectors[0] * eigenvectors[2] + eigenvectors[1] * eigenvectors[3] ) <= tolerance );
    }
}
