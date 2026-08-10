/** \file generalIntegrator_test.cpp
 * \brief Placeholder tests for general adaptive-optics integration.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/generalIntegrator.hpp"

/** \cond */
/// Compile the general integrator and its measurement type for a representative scalar type.
template struct mx::AO::sim::wfMeasurement<double>;
template class mx::AO::sim::generalIntegrator<double>;
/** \endcond */

/** \defgroup generalIntegrator_unit_tests generalIntegrator Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::placeholder::ao_sim_generalIntegrator_test
{

/** \brief Verifies that the generalIntegrator API is represented in the unit-test target.
 *
 * \ingroup generalIntegrator_unit_tests
 * \todo Add behavioral assertions for the APIs declared in mx::AO::sim::generalIntegrator.
 */
TEST_CASE( "generalIntegrator API has a test placeholder", "[ao::sim::generalIntegrator][placeholder]" )
{
    SUCCEED( "generalIntegrator behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_sim_generalIntegrator_test
