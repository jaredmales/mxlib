/** \file turbSequence_test.cpp
 * \brief Placeholder tests for simulated turbulence sequences.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/turbSequence.hpp"

/** \cond */
/// Compile the turbulence-sequence implementation for a representative scalar type.
template struct mx::AO::sim::turbSequence<double>;
/** \endcond */

/** \defgroup turbSequence_unit_tests turbSequence Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::placeholder::ao_sim_turbSequence_test
{

/** \brief Verifies that the turbSequence API is represented in the unit-test target.
 *
 * \ingroup turbSequence_unit_tests
 * \todo Add temporary FITS sequence fixtures to verify mx::AO::sim::turbSequence::turbFnames success and maximum-file
 * clamping behavior.
 */
TEST_CASE( "turbSequence API has a test placeholder", "[ao::sim::turbSequence][placeholder]" )
{
    SUCCEED( "turbSequence behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_sim_turbSequence_test
