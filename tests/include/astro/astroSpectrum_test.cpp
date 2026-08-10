/** \file astroSpectrum_test.cpp
 * \brief Placeholder tests for astronomical spectrum manipulation.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/astroSpectrum.hpp"
#include "../../../include/astro/astroSpectra.hpp"

/** \cond
 * Explicit instantiations compile-check the double-precision spectrum base and a representative basic SI spectrum,
 * and emit their header-defined methods for coverage accounting in this test translation unit.
 */
using astroSpectrumSiDouble = mx::astro::units::si<double>;
using astroSpectrumBasicSiDouble = mx::astro::basicSpectrum<astroSpectrumSiDouble>;
template struct mx::astro::baseSpectrum<double>;
template struct mx::astro::astroSpectrum<astroSpectrumBasicSiDouble>;
/** \endcond */

/** \defgroup astroSpectrum_unit_tests astroSpectrum Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_astroSpectrum_test
{

/** \brief Verifies that the mx::astro::astroSpectrum API is represented in the unit-test target.
 *
 * \ingroup astroSpectrum_unit_tests
 * \todo Add behavioral assertions for mx::astro::astroSpectrum and mx::astro::baseSpectrum.
 */
TEST_CASE( "astroSpectrum API has a test placeholder", "[astro::astroSpectrum][placeholder]" )
{
    SUCCEED( "astroSpectrum behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_astroSpectrum_test
