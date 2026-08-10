/** \file clAOLinearPredictor_test.cpp
 * \brief Tests of closed-loop linear-predictor regularization.
 */

#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/clAOLinearPredictor.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace
{

using predictorT = mx::AO::analysis::clAOLinearPredictor<double>;
using optimizerT = mx::AO::analysis::clGainOpt<double>;

/// Construct a compact, valid PSD fixture for regularization tests.
void makeFixture( optimizerT &optimizer,            /**< [out] configured gain optimizer */
                  std::vector<double> &disturbance, /**< [out] positive disturbance PSD */
                  std::vector<double> &noise /**< [out] nonnegative noise PSD */ )
{
    constexpr std::size_t sampleCount = 64;
    std::vector<double> frequency( sampleCount );
    disturbance.resize( sampleCount );
    noise.assign( sampleCount, 1e-3 );

    for( std::size_t index = 0; index < sampleCount; ++index )
    {
        frequency[index] = 0.5 * static_cast<double>( index + 1 ) / static_cast<double>( sampleCount );
        disturbance[index] = 1.0 / ( 1.0 + 100.0 * frequency[index] * frequency[index] );
    }

    optimizer.f( frequency );
}

/// Run one regularization search with telemetry enabled.
mx::error_t runSearch( predictorT &predictor,            /**< [in,out] predictor under test */
                       optimizerT &optimizer,            /**< [in,out] configured gain optimizer */
                       std::vector<double> &disturbance, /**< [in] disturbance PSD */
                       std::vector<double> &noise /**< [in] noise PSD */ )
{
    double maximumGain = 0;
    double optimalGain = 0;
    double variance = 0;
    double scale = 0;
    return predictor
        .regularizeCoefficients<true>( maximumGain, optimalGain, variance, scale, optimizer, disturbance, noise, 4 );
}

} // namespace

/// Verify that invalid linear-predictor regularization controls are rejected before evaluation.
/** Exercises validation of the regularization interval, spacing, refinement divisor, and iteration limit. */
/**
 * \ingroup clAOLinearPredictor_unit_tests
 */
TEST_CASE( "Linear-predictor regularization rejects invalid search controls", "[ao::analysis::clAOLinearPredictor]" )
{
    // clang-format off
#ifdef __DOXY_ONLY__
    mx::AO::analysis::clAOLinearPredictor<double>::regularizeCoefficients<true>();
#endif
    // clang-format on

    optimizerT optimizer( 1.0, 1.5 );
    std::vector<double> disturbance;
    std::vector<double> noise;
    makeFixture( optimizer, disturbance, noise );

    GIVEN( "a fresh predictor" )
    {
        predictorT predictor;

        WHEN( "the initial precision is below the minimum" )
        {
            predictor.m_precision0 = 0.5 * predictor.m_minPrecision;
            REQUIRE( runSearch( predictor, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );
        }

        WHEN( "the initial precision equals the minimum" )
        {
            predictor.m_precision0 = predictor.m_minPrecision;
            REQUIRE( runSearch( predictor, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );
        }

        WHEN( "the initial precision is zero, negative, NaN, or infinite" )
        {
            const std::vector<double> invalidValues{ 0,
                                                     -1,
                                                     std::numeric_limits<double>::quiet_NaN(),
                                                     std::numeric_limits<double>::infinity() };
            for( const double value : invalidValues )
            {
                predictorT invalidPredictor;
                invalidPredictor.m_precision0 = value;
                REQUIRE( runSearch( invalidPredictor, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );
                REQUIRE( invalidPredictor.m_regularizationReport.status ==
                         predictorT::regularizationStatus::invalidControls );
                REQUIRE( invalidPredictor.m_regularizationReport.evaluations == 0 );
                REQUIRE( invalidPredictor.m_regResults.empty() );
            }
        }

        WHEN( "the initial precision is wider than the initial interval" )
        {
            predictor.m_precision0 = predictor.m_max_sc0 - predictor.m_min_sc0 + 1;
            REQUIRE( runSearch( predictor, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );
        }

        WHEN( "the scale interval, refinement divisor, or iteration limit is invalid" )
        {
            predictorT reversedInterval;
            reversedInterval.m_max_sc0 = reversedInterval.m_min_sc0;
            REQUIRE( runSearch( reversedInterval, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );

            predictorT invalidDivisor;
            invalidDivisor.m_dPrecision = 1;
            REQUIRE( runSearch( invalidDivisor, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );

            predictorT invalidIterations;
            invalidIterations.m_maxIts = 0;
            REQUIRE( runSearch( invalidIterations, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );
        }

        WHEN( "another floating-point search control is nonfinite or nonpositive" )
        {
            predictorT nonfiniteMinimumScale;
            nonfiniteMinimumScale.m_min_sc0 = std::numeric_limits<double>::quiet_NaN();
            REQUIRE( runSearch( nonfiniteMinimumScale, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );

            predictorT nonfiniteMaximumScale;
            nonfiniteMaximumScale.m_max_sc0 = std::numeric_limits<double>::infinity();
            REQUIRE( runSearch( nonfiniteMaximumScale, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );

            predictorT zeroMinimumPrecision;
            zeroMinimumPrecision.m_minPrecision = 0;
            REQUIRE( runSearch( zeroMinimumPrecision, optimizer, disturbance, noise ) == mx::error_t::invalidconfig );

            predictorT nonfiniteRefinementDivisor;
            nonfiniteRefinementDivisor.m_dPrecision = std::numeric_limits<double>::infinity();
            REQUIRE( runSearch( nonfiniteRefinementDivisor, optimizer, disturbance, noise ) ==
                     mx::error_t::invalidconfig );
        }
    }
}

/// Verify that linear-predictor regularization reports how its search terminated.
/** Exercises completed, boundary-limited, and iteration-limited regularization searches. */
/**
 * \ingroup clAOLinearPredictor_unit_tests
 */
TEST_CASE( "Linear-predictor regularization reports search termination", "[ao::analysis::clAOLinearPredictor]" )
{
    // clang-format off
#ifdef __DOXY_ONLY__
    mx::AO::analysis::clAOLinearPredictor<double>::regularizeCoefficients<true>();
#endif
    // clang-format on

    optimizerT optimizer( 1.0, 1.5 );
    std::vector<double> disturbance;
    std::vector<double> noise;
    makeFixture( optimizer, disturbance, noise );

    GIVEN( "a precision immediately above the minimum and one allowed iteration" )
    {
        predictorT predictor;
        predictor.m_min_sc0 = 10;
        predictor.m_max_sc0 = 10.0025;
        predictor.m_precision0 = std::nextafter( predictor.m_minPrecision, std::numeric_limits<double>::infinity() );
        predictor.m_maxIts = 1;

        const mx::error_t result = runSearch( predictor, optimizer, disturbance, noise );

        REQUIRE( result != mx::error_t::invalidconfig );
        REQUIRE( predictor.m_regularizationReport.status != predictorT::regularizationStatus::invalidControls );
        REQUIRE( predictor.m_regularizationReport.iterations == 1 );
        REQUIRE( predictor.m_regularizationReport.evaluations > 0 );
        REQUIRE_FALSE( predictor.m_regResults.empty() );
    }

    GIVEN( "a valid search forced to exhaust its iteration limit" )
    {
        predictorT predictor;
        predictor.m_min_sc0 = 10;
        predictor.m_max_sc0 = 12;
        predictor.m_precision0 = 1;
        predictor.m_minPrecision = 1e-9;
        predictor.m_maxIts = 1;

        REQUIRE( runSearch( predictor, optimizer, disturbance, noise ) == mx::error_t::timeout );
        REQUIRE( predictor.m_regularizationReport.status == predictorT::regularizationStatus::iterationLimit );
        REQUIRE( predictor.m_regularizationReport.iterations == 1 );
        REQUIRE( predictor.m_regularizationReport.evaluations > 0 );
    }

    GIVEN( "a valid search whose initial precision equals the interval width" )
    {
        predictorT predictor;
        predictor.m_min_sc0 = 10;
        predictor.m_max_sc0 = 12;
        predictor.m_precision0 = 2;
        predictor.m_minPrecision = 0.1;

        REQUIRE( runSearch( predictor, optimizer, disturbance, noise ) == mx::error_t::noerror );
        REQUIRE( ( predictor.m_regularizationReport.status == predictorT::regularizationStatus::converged ||
                   predictor.m_regularizationReport.status == predictorT::regularizationStatus::boundaryLimited ) );
        REQUIRE( predictor.m_regularizationReport.iterations > 0 );
        REQUIRE( predictor.m_regularizationReport.evaluations > 0 );
    }
}
