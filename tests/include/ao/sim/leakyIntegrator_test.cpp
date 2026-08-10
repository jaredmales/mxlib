/** \file leakyIntegrator_test.cpp
 * \brief Tests for leaky adaptive-optics integration.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/leakyIntegrator.hpp"

/** \cond */
/// Compile the leaky integrator and its measurement type for a representative scalar type.
template struct mx::AO::sim::wfMeasurement<double>;
template class mx::AO::sim::leakyIntegrator<double>;
/** \endcond */

/** \defgroup leakyIntegrator_unit_tests leakyIntegrator Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::ao_sim_leakyIntegrator_test
{

/** \brief Verifies single-mode leak configuration through mx::AO::sim::leakyIntegrator::leak.
 *
 * \ingroup leakyIntegrator_unit_tests
 * \todo Add behavioral assertions for the remaining APIs declared in mx::AO::sim::leakyIntegrator.
 */
TEST_CASE( "leakyIntegrator sets a single-mode leak", "[ao::sim::leakyIntegrator]" )
{
    mx::AO::sim::leakyIntegrator<double> integrator;

    REQUIRE( integrator.initialize( 2 ) == 0 );
    CHECK( integrator.leak( 1, 0.25 ) == 0 );
}

} // namespace unitTest::ao_sim_leakyIntegrator_test
