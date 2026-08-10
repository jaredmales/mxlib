/** \file aoConstants_test.cpp
 * \brief Placeholder tests for adaptive-optics constants.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/analysis/aoConstants.hpp"

/** \defgroup aoConstants_unit_tests aoConstants Unit Tests
 * \ingroup ao_analysis_unit_tests
 */

namespace unitTest::placeholder::ao_analysis_aoConstants_test
{

/** \brief Verifies that the mx::AO::constants API is represented in the unit-test target.
 *
 * \ingroup aoConstants_unit_tests
 * \todo Add behavioral assertions for mx::AO::constants::calcConstants, mx::AO::constants::a_SF, and
 * mx::AO::constants::a_PSD.
 */
TEST_CASE( "aoConstants API has a test placeholder", "[ao::analysis::aoConstants][placeholder]" )
{
    SUCCEED( "Adaptive-optics constants behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_analysis_aoConstants_test
