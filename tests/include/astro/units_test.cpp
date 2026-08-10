/** \file units_test.cpp
 * \brief Placeholder tests for astronomical unit systems.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/units.hpp"

/** \cond
 * Explicit instantiations compile-check all five double-precision unit systems. These constant-only types provide
 * compile ownership even when they emit no executable coverage counters.
 */
template struct mx::astro::units::si<double>;
template struct mx::astro::units::cgs<double>;
template struct mx::astro::units::solar<double>;
template struct mx::astro::units::earth<double>;
template struct mx::astro::units::jupiter<double>;
/** \endcond */

/** \defgroup units_unit_tests Astronomy Units Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_units_test
{

/** \brief Verifies that the mx::astro::units API is represented in the unit-test target.
 *
 * \ingroup units_unit_tests
 * \todo Add behavioral assertions for the unit systems declared in mx::astro::units.
 */
TEST_CASE( "astronomy units API has a test placeholder", "[astro::units][placeholder]" )
{
    SUCCEED( "Astronomy units behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_units_test
