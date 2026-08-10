/** \file turbAtmosphere_test.cpp
 * \brief Tests simulated turbulent-atmosphere generation.
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/aoSystem.hpp"
#include "../../../../include/ao/sim/turbAtmosphere.hpp"

/** \defgroup turbAtmosphere_unit_tests turbAtmosphere Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::ao_sim_turbAtmosphere_test
{

/** \brief Verifies a configured one-layer atmosphere produces a finite phase screen.
 *
 * Exercises mx::AO::sim::turbAtmosphere::setup, mx::AO::sim::turbAtmosphere::setLayers, and
 * mx::AO::sim::turbAtmosphere::genLayers with subharmonics enabled.
 *
 * \ingroup turbAtmosphere_unit_tests
 */
TEST_CASE( "turbAtmosphere generates a one-layer subharmonic screen", "[ao::sim::turbAtmosphere]" )
{
    using realT = double;
    using aoSystemT = mx::AO::analysis::aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>>;
    using atmosphereT = mx::AO::sim::turbAtmosphere<aoSystemT>;

    aoSystemT aosys;
    aosys.D( 25.0 );
    aosys.lam_sci( 0.5e-6 );
    aosys.atm.setSingleLayer( 0.16, 0.5e-6, 0.0, 0.0, 0.0, 0.0, 0.0 );

    atmosphereT atmosphere;
    atmosphere.setup( 64, 0, &aosys, 1 );
    atmosphere.retain( true );
    atmosphere.setLayers( 192 );
    atmosphere.normVar().seed( 67890 );
    atmosphere.genLayers();

    REQUIRE( atmosphere.nLayers() == 1 );
    REQUIRE( atmosphere.layer( 0 ).m_phase.isFinite().all() );
    REQUIRE( atmosphere.layer( 0 ).m_phase.square().sum() > 0 );
}

} // namespace unitTest::ao_sim_turbAtmosphere_test
