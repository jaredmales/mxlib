/** \file aoWFS_test.cpp
 * \brief Placeholder tests for adaptive-optics wavefront-sensor analysis.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/analysis/aoWFS.hpp"

/** \cond */
/// Compile the analytic wavefront-sensor implementations for representative scalar and stream types.
template struct mx::AO::analysis::wfs<double>;
template struct mx::AO::analysis::pywfsUnmod<double>;
template struct mx::AO::analysis::pywfsModAsymptotic<double>;
template struct mx::AO::analysis::shwfs<double>;
template struct mx::AO::analysis::calculatedWFS<double>;
/** \endcond */

/** \defgroup aoWFS_unit_tests aoWFS Unit Tests
 * \ingroup ao_analysis_unit_tests
 */

namespace unitTest::placeholder::ao_analysis_aoWFS_test
{

/** \brief Verifies that the mx::AO::analysis wavefront-sensor APIs are represented in the unit-test target.
 *
 * \ingroup aoWFS_unit_tests
 * \todo Add behavioral assertions for mx::AO::analysis::wfs and its concrete sensor models.
 */
TEST_CASE( "aoWFS API has a test placeholder", "[ao::analysis::aoWFS][placeholder]" )
{
    SUCCEED( "aoWFS behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_analysis_aoWFS_test
