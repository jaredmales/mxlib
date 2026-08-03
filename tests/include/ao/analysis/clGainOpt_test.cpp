/** \file clGainOpt_test.cpp
 * \brief Tests of closed-loop gain optimization transfer functions.
 */

#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/clGainOpt.hpp"

#include <complex>
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
    double retimedMaximumGain = 0.5;
    const double retimedOptimalGain =
        retimed.optGainOpenLoop( retimedOptimalVariance, disturbancePsd, noisePsd, retimedMaximumGain, false );

    double freshOptimalVariance = 0;
    double freshMaximumGain = 0.5;
    const double freshOptimalGain =
        fresh.optGainOpenLoop( freshOptimalVariance, disturbancePsd, noisePsd, freshMaximumGain, false );

    REQUIRE( retimedOptimalGain == Approx( freshOptimalGain ).epsilon( 1e-12 ).margin( 1e-14 ) );
    REQUIRE( retimedOptimalVariance == Approx( freshOptimalVariance ).epsilon( 1e-12 ).margin( 1e-14 ) );
}
