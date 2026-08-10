/** \file ccdDetector_test.cpp
 * \brief Placeholder tests for simulated CCD detectors.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/ccdDetector.hpp"

/** \cond */
/// Compile the CCD detector implementation for a representative scalar type.
template class mx::AO::sim::ccdDetector<double>;
/** \endcond */

/** \defgroup ccdDetector_unit_tests ccdDetector Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::placeholder::ao_sim_ccdDetector_test
{

/** \brief Verifies that the ccdDetector API is represented in the unit-test target.
 *
 * \ingroup ccdDetector_unit_tests
 * \todo Add behavioral assertions for the APIs declared in mx::AO::sim::ccdDetector.
 */
TEST_CASE( "ccdDetector API has a test placeholder", "[ao::sim::ccdDetector][placeholder]" )
{
    SUCCEED( "ccdDetector behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_sim_ccdDetector_test
