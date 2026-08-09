/** \file clGainOpt_test.cpp
 * \brief Tests of closed-loop gain optimization transfer functions.
 */

#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/clGainOpt.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace
{

using optimizerT = mx::AO::analysis::clGainOpt<double>;

/// Require two complex values to agree within floating-point precision.
void requireComplexEqual( const std::complex<double> &actual, /**< [in] value produced by the retimed optimizer */
                          const std::complex<double> &expected /**< [in] value produced by the fresh optimizer */ )
{
    REQUIRE( actual.real() == Approx( expected.real() ).epsilon( 1e-12 ).margin( 1e-14 ) );
    REQUIRE( actual.imag() == Approx( expected.imag() ).epsilon( 1e-12 ).margin( 1e-14 ) );
}

} // namespace

/// Verify that changing the sampling interval invalidates all sampling-interval-dependent cached values.
/** Exercises mx::AO::analysis::clGainOpt::Ti and the public transfer-function calculations after the trigonometric
 * cache has already been populated.
 */
TEST_CASE( "Gain optimizer recomputes transfer functions after changing Ti", "[ao::analysis::clGainOpt]" )
{
    constexpr double initialTi = 0.001;
    constexpr double updatedTi = 0.0017;
    constexpr double delay = 0.0013;
    constexpr double gain = 0.42;

    const std::vector<double> frequency{ 25.0, 50.0, 75.0, 100.0 };
    const std::vector<double> fir{ 0.75, 0.20, -0.05 };
    const std::vector<double> iir{ 0.45, -0.08, 0.025 };

    optimizerT retimed( initialTi, delay );
    retimed.f( frequency );
    retimed.b( fir );
    retimed.a( iir );
    retimed.remember( 0.97 );

    for( std::size_t index = 0; index < frequency.size(); ++index )
    {
        retimed.olXfer( index );
    }

    retimed.Ti( updatedTi );

    optimizerT fresh( updatedTi, delay );
    fresh.f( frequency );
    fresh.b( fir );
    fresh.a( iir );
    fresh.remember( 0.97 );

    for( std::size_t index = 0; index < frequency.size(); ++index )
    {
        CAPTURE( index, frequency[index] );

        optimizerT::complexT retimedDm;
        optimizerT::complexT retimedDelay;
        optimizerT::complexT retimedController;
        const optimizerT::complexT retimedOpenLoop =
            retimed.olXfer( index, retimedDm, retimedDelay, retimedController );

        optimizerT::complexT freshDm;
        optimizerT::complexT freshDelay;
        optimizerT::complexT freshController;
        const optimizerT::complexT freshOpenLoop = fresh.olXfer( index, freshDm, freshDelay, freshController );

        requireComplexEqual( retimedDm, freshDm );
        requireComplexEqual( retimedDelay, freshDelay );
        requireComplexEqual( retimedController, freshController );
        requireComplexEqual( retimedOpenLoop, freshOpenLoop );
        requireComplexEqual( retimed.clETF( index, gain ), fresh.clETF( index, gain ) );
        REQUIRE( retimed.clETFPhase( index, gain ) ==
                 Approx( fresh.clETFPhase( index, gain ) ).epsilon( 1e-12 ).margin( 1e-14 ) );
        REQUIRE( retimed.clETF2( index, gain ) ==
                 Approx( fresh.clETF2( index, gain ) ).epsilon( 1e-12 ).margin( 1e-14 ) );
        requireComplexEqual( retimed.clNTF( index, gain ), fresh.clNTF( index, gain ) );
        REQUIRE( retimed.clNTF2( index, gain ) ==
                 Approx( fresh.clNTF2( index, gain ) ).epsilon( 1e-12 ).margin( 1e-14 ) );

        double retimedEtf = 0;
        double retimedNtf = 0;
        retimed.clTF2( retimedEtf, retimedNtf, index, gain );

        double freshEtf = 0;
        double freshNtf = 0;
        fresh.clTF2( freshEtf, freshNtf, index, gain );

        REQUIRE( retimedEtf == Approx( freshEtf ).epsilon( 1e-12 ).margin( 1e-14 ) );
        REQUIRE( retimedNtf == Approx( freshNtf ).epsilon( 1e-12 ).margin( 1e-14 ) );
    }

    const std::vector<double> disturbancePsd{ 4.0, 2.0, 0.7, 0.2 };
    const std::vector<double> noisePsd{ 0.01, 0.01, 0.01, 0.01 };
    REQUIRE( retimed.clVariance( disturbancePsd, noisePsd, gain ) ==
             Approx( fresh.clVariance( disturbancePsd, noisePsd, gain ) ).epsilon( 1e-12 ).margin( 1e-14 ) );

    double retimedOptimalVariance = 0;
    double retimedOptimalGain = 0;
    REQUIRE(
        retimed.optGainOpenLoop( retimedOptimalGain, retimedOptimalVariance, disturbancePsd, noisePsd, 0.5, false ) ==
        mx::error_t::noerror );

    double freshOptimalVariance = 0;
    double freshOptimalGain = 0;
    REQUIRE( fresh.optGainOpenLoop( freshOptimalGain, freshOptimalVariance, disturbancePsd, noisePsd, 0.5, false ) ==
             mx::error_t::noerror );

    REQUIRE( retimedOptimalGain == Approx( freshOptimalGain ).epsilon( 1e-12 ).margin( 1e-14 ) );
    REQUIRE( retimedOptimalVariance == Approx( freshOptimalVariance ).epsilon( 1e-12 ).margin( 1e-14 ) );
}

/// Verify that maximum-stable-gain calculation interpolates a bracketed pure-integrator Nyquist crossing.
/** Exercises mx::AO::analysis::clGainOpt::maxStableGain and its crossing diagnostics. */
TEST_CASE( "Gain optimizer interpolates a pure-integrator stability crossing", "[ao::analysis::clGainOpt]" )
{
    optimizerT optimizer( 0.001, 0.0015 );
    const std::vector<double> frequency{ 100.0, 124.0, 126.0, 150.0 };
    optimizer.f( frequency );

    const double crossingPhase = std::numbers::pi_v<double> / 4.0;
    const double expectedGain = crossingPhase * crossingPhase / ( 2.0 * std::sin( crossingPhase / 2.0 ) );
    const double lowerSampleGain = -1.0 / optimizer.olXfer( 1 ).real();

    double maximumGain = 0;
    optimizerT::maxStableGainReport report;
    REQUIRE( optimizer.maxStableGain( maximumGain, &report ) == mx::error_t::noerror );

    REQUIRE( report.status == optimizerT::maxStableGainStatus::crossingFound );
    REQUIRE( report.lowerIndex == 1 );
    REQUIRE( report.upperIndex == 2 );
    REQUIRE( report.lowerFrequency == 124.0 );
    REQUIRE( report.upperFrequency == 126.0 );
    REQUIRE( report.crossingFrequency == Approx( 125.0 ).margin( 0.02 ) );
    REQUIRE( maximumGain == report.gain );
    REQUIRE( std::abs( maximumGain - expectedGain ) < std::abs( lowerSampleGain - expectedGain ) );
}

/// Verify that maximum-stable-gain calculation reports invalid grids and missing crossings.
/** Exercises mx::AO::analysis::clGainOpt::maxStableGain failure statuses without sentinel gain values. */
TEST_CASE( "Gain optimizer reports missing stability crossings", "[ao::analysis::clGainOpt]" )
{
    optimizerT optimizer( 0.001, 0.0015 );
    double maximumGain = 0;
    optimizerT::maxStableGainReport report;

    optimizer.f( std::vector<double>{ 1.0 } );
    REQUIRE( optimizer.maxStableGain( maximumGain, &report ) == mx::error_t::sizeerr );
    REQUIRE( report.status == optimizerT::maxStableGainStatus::invalidInput );
    REQUIRE( mx::math::isNan( maximumGain ) );

    optimizer.f( std::vector<double>{ 1.0, 2.0, 3.0 } );
    REQUIRE( optimizer.maxStableGain( maximumGain, &report ) == mx::error_t::notfound );
    REQUIRE( report.status == optimizerT::maxStableGainStatus::noCrossing );
    REQUIRE( mx::math::isNan( maximumGain ) );
}

/// Verify that optimum-gain calculation enforces its interval and reports termination state.
/** Exercises both mx::AO::analysis::clGainOpt::optGainOpenLoop overloads for small intervals, invalid intervals,
 * stability failure, and forced iteration exhaustion.
 */
TEST_CASE( "Gain optimizer reports minimizer termination", "[ao::analysis::clGainOpt]" )
{
    optimizerT optimizer( 0.001, 0.0015 );
    std::vector<double> frequency;
    std::vector<double> disturbance;
    std::vector<double> noise;
    for( int index = 1; index <= 500; ++index )
    {
        const double value = static_cast<double>( index );
        frequency.push_back( value );
        disturbance.push_back( 1.0 / ( 1.0 + value * value ) );
        noise.push_back( 1e-6 );
    }
    optimizer.f( frequency );

    double optimalGain = 0;
    double variance = 0;
    optimizerT::optGainReport report;
    constexpr double smallMaximumGain = 1e-3;
    REQUIRE( optimizer.optGainOpenLoop( optimalGain, variance, disturbance, noise, smallMaximumGain, true, &report ) ==
             mx::error_t::noerror );
    REQUIRE( ( report.status == optimizerT::optGainStatus::converged ||
               report.status == optimizerT::optGainStatus::boundaryLimited ) );
    REQUIRE( report.minimumEvaluatedGain >= optimizer.m_minFindMin );
    REQUIRE( report.maximumEvaluatedGain <= optimizer.m_minFindMaxFact * smallMaximumGain );
    REQUIRE( mx::math::isFinite( optimalGain ) );
    REQUIRE( mx::math::isFinite( variance ) );

    REQUIRE( optimizer.optGainOpenLoop( optimalGain, variance, disturbance, noise, 1e-10, false, &report ) ==
             mx::error_t::invalidconfig );
    REQUIRE( report.status == optimizerT::optGainStatus::invalidInput );
    REQUIRE( mx::math::isNan( optimalGain ) );
    REQUIRE( mx::math::isNan( variance ) );

    optimizer.f( std::vector<double>{ 1.0, 2.0, 3.0 } );
    disturbance.resize( 3 );
    noise.resize( 3 );
    REQUIRE( optimizer.optGainOpenLoop( optimalGain, variance, disturbance, noise, false, &report ) ==
             mx::error_t::notfound );
    REQUIRE( report.status == optimizerT::optGainStatus::stabilityFailure );
    REQUIRE( report.stability.status == optimizerT::maxStableGainStatus::noCrossing );

    optimizer.f( frequency );
    disturbance.resize( frequency.size(), 1e-3 );
    noise.resize( frequency.size(), 1e-6 );
    optimizer.m_minFindMaxIter = 1;
    REQUIRE( optimizer.optGainOpenLoop( optimalGain, variance, disturbance, noise, 0.5, false, &report ) ==
             mx::error_t::timeout );
    REQUIRE( report.status == optimizerT::optGainStatus::iterationLimit );
    REQUIRE( report.iterations == 1 );
}
