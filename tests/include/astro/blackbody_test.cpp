/** \file blackbody_test.cpp
 * \brief Placeholder tests for blackbody spectra.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/blackbody.hpp"

/** \cond
 * Explicit instantiation compile-checks the representative double-precision SI blackbody and emits its
 * header-defined methods for coverage accounting in this test translation unit.
 */
template struct mx::astro::blackbody<mx::astro::units::si<double>>;
/** \endcond */

/** \defgroup blackbody_unit_tests blackbody Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_blackbody_test
{

/** \brief Verifies that the mx::astro::blackbody API is represented in the unit-test target.
 *
 * \ingroup blackbody_unit_tests
 * \todo Add behavioral assertions for mx::astro::blackbody.
 */
TEST_CASE( "blackbody API has a test placeholder", "[astro::blackbody][placeholder]" )
{
    SUCCEED( "blackbody behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_blackbody_test
