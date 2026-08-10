/** \file turbSubHarmonic_test.cpp
 * \brief Tests low-frequency turbulence subharmonic generation.
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include <cmath>

#include "../../../../include/ao/analysis/aoSystem.hpp"
#include "../../../../include/ao/sim/turbAtmosphere.hpp"
#include "../../../../include/ao/sim/turbSubHarmonic.hpp"

/** \defgroup turbSubHarmonic_unit_tests turbSubHarmonic Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::ao_sim_turbSubHarmonic_test
{

/** \brief Verifies real quadrature subharmonics generate isotropic nonzero tilt variance.
 *
 * Exercises mx::AO::sim::turbSubHarmonic::initGrid and mx::AO::sim::turbSubHarmonic::screen through both
 * pre-calculated and direct mode evaluation.
 *
 * \ingroup turbSubHarmonic_unit_tests
 */
TEST_CASE( "turbSubHarmonic generates cosine and sine quadratures", "[ao::sim::turbSubHarmonic]" )
{
    using realT = double;
    using aoSystemT = mx::AO::analysis::aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>>;
    using atmosphereT = mx::AO::sim::turbAtmosphere<aoSystemT>;

#ifdef MXLIB_BUILD_COVERAGE
    constexpr uint32_t wfSz = 32;
    constexpr uint32_t scrnSz = 64;
    constexpr size_t trials = 8;
    constexpr realT tiltRatioMargin = 0.5;
#else
    constexpr uint32_t wfSz = 64;
    constexpr uint32_t scrnSz = 192;
    constexpr size_t trials = 128;
    constexpr realT tiltRatioMargin = 0.35;
#endif

    aoSystemT aosys;
    aosys.D( 25.0 );
    aosys.lam_sci( 0.5e-6 );
    aosys.atm.setSingleLayer( 0.16, 0.5e-6, 0.0, 0.0, 0.0, 0.0, 0.0 );

    atmosphereT atmosphere;
    atmosphere.setup( wfSz, 0, &aosys, 1 );
    atmosphere.setLayers( scrnSz );

    mx::AO::sim::turbSubHarmonic<atmosphereT> preCalculated;
    preCalculated.turbAtmo( &atmosphere );
    preCalculated.level( 1 );
    preCalculated.outerSubHarmonics( true );
    preCalculated.preCalc( true );
    preCalculated.initGrid( 0 );

    mx::AO::sim::turbSubHarmonic<atmosphereT> direct;
    direct.turbAtmo( &atmosphere );
    direct.level( 1 );
    direct.outerSubHarmonics( true );
    direct.preCalc( false );
    direct.initGrid( 0 );

    mx::improc::eigenImage<realT> preScreen( scrnSz, scrnSz );
    mx::improc::eigenImage<realT> directScreen( scrnSz, scrnSz );

    atmosphere.normVar().seed( 12345 );
    preScreen.setZero();
    preCalculated.screen( preScreen );
    atmosphere.normVar().seed( 12345 );
    directScreen.setZero();
    direct.screen( directScreen );

    REQUIRE( ( preScreen - directScreen ).abs().maxCoeff() < 1e-12 );

    realT xTiltMeanSq{ 0 };
    realT yTiltMeanSq{ 0 };
    for( size_t trial = 0; trial < trials; ++trial )
    {
        preScreen.setZero();
        preCalculated.screen( preScreen );

        realT xTilt{ 0 };
        realT yTilt{ 0 };
        for( uint32_t jj = 0; jj < scrnSz; ++jj )
        {
            realT y = jj - 0.5 * ( scrnSz - 1 );
            for( uint32_t ii = 0; ii < scrnSz; ++ii )
            {
                realT x = ii - 0.5 * ( scrnSz - 1 );
                xTilt += x * preScreen( ii, jj );
                yTilt += y * preScreen( ii, jj );
            }
        }

        xTiltMeanSq += xTilt * xTilt;
        yTiltMeanSq += yTilt * yTilt;
    }

    xTiltMeanSq /= trials;
    yTiltMeanSq /= trials;

    REQUIRE( xTiltMeanSq > 1e-6 );
    REQUIRE( yTiltMeanSq > 1e-6 );
    REQUIRE( xTiltMeanSq / yTiltMeanSq == Approx( 1.0 ).margin( tiltRatioMargin ) );
}

} // namespace unitTest::ao_sim_turbSubHarmonic_test
