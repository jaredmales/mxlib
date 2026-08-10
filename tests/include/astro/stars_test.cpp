/** \file stars_test.cpp
 * \brief Placeholder tests for stellar properties.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/stars.hpp"

/// \cond
// Compile every non-template member of a representative stellar-sequence specialization.
template struct mx::astro::mainSequence<double>;
/// \endcond

/** \defgroup stars_unit_tests Stellar Properties Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_stars_test
{

/** \brief Verifies that the mx::astro stellar APIs are represented in the unit-test target.
 *
 * \ingroup stars_unit_tests
 * \todo Add behavioral assertions for mx::astro::mainSequence and related stellar APIs.
 */
TEST_CASE( "stellar APIs have a test placeholder", "[astro::stars][placeholder]" )
{
    SUCCEED( "Stellar-property behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_stars_test
