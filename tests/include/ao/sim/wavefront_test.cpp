/** \file wavefront_test.cpp
 * \brief Placeholder tests for simulated wavefronts.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/wavefront.hpp"

/** \cond */
/// Compile the simulated wavefront representation for a representative scalar type.
template struct mx::AO::sim::wavefront<double>;
/** \endcond */

/** \defgroup wavefront_unit_tests AO Wavefront Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::placeholder::ao_sim_wavefront_test
{

/** \brief Verifies that the simulated wavefront API is represented in the unit-test target.
 *
 * \ingroup wavefront_unit_tests
 * \todo Add behavioral assertions for the APIs declared in mx::AO::sim::wavefront.
 */
TEST_CASE( "simulated wavefront API has a test placeholder", "[ao::sim::wavefront][placeholder]" )
{
    SUCCEED( "Simulated wavefront behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_sim_wavefront_test
