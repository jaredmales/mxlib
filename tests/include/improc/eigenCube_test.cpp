/** \file eigenCube_test.cpp
 * \brief Tests of masked image combination in eigenCube.
 * \author OpenAI Codex
 */

#include "../../catch2/catch.hpp"

#include <limits>
#include <utility>
#include <vector>

#include "../../../include/improc/eigenCube.hpp"

/** \cond
 * Explicit instantiation compile-checks the representative double-precision cube and emits its header-defined
 * methods for coverage accounting in this test translation unit.
 */
template class mx::improc::eigenCube<double>;
/** \endcond */

namespace unitTest
{
namespace improcTest
{
namespace eigenCubeTest
{

/// Verify masked mean combinations accept an inclusive good-pixel fraction and use the cube scalar sentinel.
/** Exercises mx::improc::eigenCube::mean for masked and weighted-masked combinations. */
/**
 * \ingroup eigenCube_unit_tests
 */
TEST_CASE( "Masked eigenCube means use an inclusive good-pixel threshold", "[improc::eigenCube]" )
{
    mx::improc::eigenCube<double> cube( 1, 1, 4 );
    cube.image( 0 )( 0, 0 ) = 2.0;
    cube.image( 1 )( 0, 0 ) = 4.0;
    cube.image( 2 )( 0, 0 ) = 8.0;
    cube.image( 3 )( 0, 0 ) = 16.0;

    mx::improc::eigenCube<int> mask( 1, 1, 4 );
    mask.image( 0 )( 0, 0 ) = 1;
    mask.image( 1 )( 0, 0 ) = 1;
    mask.image( 2 )( 0, 0 ) = 0;
    mask.image( 3 )( 0, 0 ) = 0;

    mx::improc::eigenImage<double> result;

    cube.mean( result, mask, 0.49 );
    REQUIRE( result( 0, 0 ) == 3.0 );

    cube.mean( result, mask, 0.5 );
    REQUIRE( result( 0, 0 ) == 3.0 );

    cube.mean( result, mask, 0.51 );
    REQUIRE( result( 0, 0 ) == mx::improc::invalidNumber<double>() );
    REQUIRE( result( 0, 0 ) != static_cast<double>( mx::improc::invalidNumber<float>() ) );

    std::vector<double> weights{ 1.0, 3.0, 5.0, 7.0 };
    cube.mean( result, weights, mask, 0.5 );
    REQUIRE( result( 0, 0 ) == 3.5 );

    mask.setZero();
    cube.mean( result, mask, 0.0 );
    REQUIRE( result( 0, 0 ) == mx::improc::invalidNumber<double>() );

    mask.image( 0 )( 0, 0 ) = 1;
    cube.mean( result, mask, 0.0 );
    REQUIRE( result( 0, 0 ) == 2.0 );

    for( int plane = 0; plane < mask.planes(); ++plane )
    {
        mask.image( plane )( 0, 0 ) = 1;
    }

    cube.mean( result, mask, 1.0 );
    REQUIRE( result( 0, 0 ) == 7.5 );
}

/** \brief Verifies unmasked and weighted image combinations used when HCI masking is disabled.
 *
 * \ingroup eigenCube_unit_tests
 */
TEST_CASE( "Unmasked eigenCube combinations use every plane", "[improc::eigenCube]" )
{
    mx::improc::eigenCube<double> cube( 1, 1, 4 );
    cube.image( 0 )( 0, 0 ) = 1.0;
    cube.image( 1 )( 0, 0 ) = 2.0;
    cube.image( 2 )( 0, 0 ) = 3.0;
    cube.image( 3 )( 0, 0 ) = 100.0;

    mx::improc::eigenImage<double> result;
    std::vector<double> weights{ 1.0, 1.0, 1.0, 3.0 };

    cube.mean( result );
    REQUIRE( result( 0, 0 ) == 26.5 );

    cube.mean( result, weights );
    REQUIRE( result( 0, 0 ) == 51.0 );

    cube.median( result );
    REQUIRE( result( 0, 0 ) == 2.5 );

    cube.sigmaMean( result, 100.0 );
    REQUIRE( result( 0, 0 ) == 26.5 );

    cube.sigmaMean( result, weights, 100.0 );
    REQUIRE( result( 0, 0 ) == 51.0 );
}

/// Verify masked median combination honors the mask and inclusive good-pixel fraction.
/** Exercises mx::improc::eigenCube::median, including its even-sample median behavior. */
/**
 * \ingroup eigenCube_unit_tests
 */
TEST_CASE( "Masked eigenCube median honors its mask and threshold", "[improc::eigenCube]" )
{
    mx::improc::eigenCube<double> cube( 1, 1, 4 );
    cube.image( 0 )( 0, 0 ) = 1.0;
    cube.image( 1 )( 0, 0 ) = 100.0;
    cube.image( 2 )( 0, 0 ) = 5.0;
    cube.image( 3 )( 0, 0 ) = 9.0;

    mx::improc::eigenCube<int> mask( 1, 1, 4 );
    mask.image( 0 )( 0, 0 ) = 1;
    mask.image( 1 )( 0, 0 ) = 0;
    mask.image( 2 )( 0, 0 ) = 1;
    mask.image( 3 )( 0, 0 ) = 0;

    mx::improc::eigenImage<double> result;

    cube.median( result, mask, 0.5 );
    REQUIRE( result( 0, 0 ) == 3.0 );

    cube.median( result, mask, 0.5001 );
    REQUIRE( result( 0, 0 ) == mx::improc::invalidNumber<double>() );

    mask.setZero();
    cube.median( result, mask, 0.0 );
    REQUIRE( result( 0, 0 ) == mx::improc::invalidNumber<double>() );

    for( int plane = 0; plane < mask.planes(); ++plane )
    {
        mask.image( plane )( 0, 0 ) = 1;
    }

    cube.median( result, mask, 1.0 );
    REQUIRE( result( 0, 0 ) == 7.0 );
}

/// Verify masked sigma means apply the inclusive threshold with and without weights.
/** Exercises mx::improc::eigenCube::sigmaMean for masked and weighted-masked combinations. */
/**
 * \ingroup eigenCube_unit_tests
 */
TEST_CASE( "Masked eigenCube sigma means use an inclusive good-pixel threshold", "[improc::eigenCube]" )
{
    mx::improc::eigenCube<double> cube( 1, 1, 4 );
    cube.image( 0 )( 0, 0 ) = 2.0;
    cube.image( 1 )( 0, 0 ) = 4.0;
    cube.image( 2 )( 0, 0 ) = 8.0;
    cube.image( 3 )( 0, 0 ) = 16.0;

    mx::improc::eigenCube<int> mask( 1, 1, 4 );
    mask.image( 0 )( 0, 0 ) = 1;
    mask.image( 1 )( 0, 0 ) = 1;
    mask.image( 2 )( 0, 0 ) = 0;
    mask.image( 3 )( 0, 0 ) = 0;

    mx::improc::eigenImage<double> result;
    std::vector<double> weights{ 1.0, 3.0, 5.0, 7.0 };

    cube.sigmaMean( result, mask, 100.0, 0.5 );
    REQUIRE( result( 0, 0 ) == 3.0 );

    cube.sigmaMean( result, weights, mask, 100.0, 0.5 );
    REQUIRE( result( 0, 0 ) == 3.5 );

    cube.sigmaMean( result, mask, 100.0, 0.51 );
    REQUIRE( result( 0, 0 ) == mx::improc::invalidNumber<double>() );
}

/// Verify masked combinations reject invalid masks, weights, and good-pixel thresholds.
/** Exercises argument validation in mx::improc::eigenCube::mean, median, and sigmaMean. */
/**
 * \ingroup eigenCube_unit_tests
 */
TEST_CASE( "eigenCube combinations validate masks weights and thresholds", "[improc::eigenCube]" )
{
    mx::improc::eigenCube<double> cube( 1, 1, 4 );
    cube.setZero();

    mx::improc::eigenCube<int> mask( 1, 1, 4 );
    mask.setZero();

    mx::improc::eigenCube<int> wrongMask( 1, 2, 4 );
    wrongMask.setZero();

    mx::improc::eigenImage<double> result;
    std::vector<double> weights{ 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> wrongWeights{ 1.0, 1.0, 1.0 };

    try
    {
        cube.mean( result, wrongMask );
        FAIL( "mismatched mask dimensions were accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::sizeerr );
    }

    try
    {
        cube.mean( result, wrongWeights );
        FAIL( "mismatched weights were accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::sizeerr );
    }

    try
    {
        cube.sigmaMean( result, wrongWeights, 3.0 );
        FAIL( "mismatched sigma-mean weights were accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::sizeerr );
    }

    try
    {
        cube.mean( result, mask, -0.01 );
        FAIL( "a negative threshold was accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::invalidarg );
    }

    try
    {
        cube.median( result, mask, 1.01 );
        FAIL( "a threshold above one was accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::invalidarg );
    }

    try
    {
        cube.sigmaMean( result, mask, 3.0, std::numeric_limits<double>::quiet_NaN() );
        FAIL( "a NaN threshold was accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::invalidarg );
    }

    try
    {
        cube.sigmaMean( result, weights, mask, 3.0, std::numeric_limits<double>::infinity() );
        FAIL( "an infinite threshold was accepted" );
    }
    catch( const mx::exception<> &error )
    {
        REQUIRE( error.code() == mx::error_t::invalidarg );
    }
}

/// Verify eigenCube copies own independent storage and moves transfer ownership safely.
/** Exercises mx::improc::eigenCube copy construction, assignment, move construction, move assignment, and shallowCopy.
 */
/**
 * \ingroup eigenCube_unit_tests
 */
TEST_CASE( "eigenCube copy and move operations preserve ownership", "[improc::eigenCube]" )
{
    mx::improc::eigenCube<double> source( 1, 2, 2 );
    source.image( 0 ) << 1.0, 2.0;
    source.image( 1 ) << 3.0, 4.0;

    mx::improc::eigenCube<double> copied( source );
    REQUIRE( copied.data() != source.data() );
    REQUIRE( copied.image( 1 )( 0, 1 ) == 4.0 );
    copied.image( 0 )( 0, 0 ) = 10.0;
    REQUIRE( source.image( 0 )( 0, 0 ) == 1.0 );

    mx::improc::eigenCube<double> assigned;
    assigned = source;
    REQUIRE( assigned.data() != source.data() );
    REQUIRE( assigned.image( 1 )( 0, 0 ) == 3.0 );

    mx::improc::eigenCube<double> moved( std::move( copied ) );
    REQUIRE( copied.data() == nullptr );
    REQUIRE( copied.rows() == 0 );
    REQUIRE( moved.image( 0 )( 0, 0 ) == 10.0 );

    mx::improc::eigenCube<double> moveAssigned;
    moveAssigned = std::move( assigned );
    REQUIRE( assigned.data() == nullptr );
    REQUIRE( assigned.planes() == 0 );
    REQUIRE( moveAssigned.image( 1 )( 0, 1 ) == 4.0 );

    mx::improc::eigenCube<double> transferred;
    transferred.shallowCopy( source, true );
    REQUIRE( source.data() == nullptr );
    REQUIRE( source.cols() == 0 );
    REQUIRE( transferred.image( 0 )( 0, 1 ) == 2.0 );

    transferred.shallowCopy( transferred, true );
    REQUIRE( transferred.image( 1 )( 0, 0 ) == 3.0 );
}

} // namespace eigenCubeTest
} // namespace improcTest
} // namespace unitTest
