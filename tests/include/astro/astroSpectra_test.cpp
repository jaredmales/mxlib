/** \file astroSpectra_test.cpp
 * \brief Placeholder tests for predefined astronomical spectra.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/astroSpectra.hpp"

/// \cond
// Compile every non-template member of each spectrum provider for representative SI units.
using astroSpectraSiDouble = mx::astro::units::si<double>;
template struct mx::astro::basicSpectrum<astroSpectraSiDouble>;
template struct mx::astro::astroFilter<astroSpectraSiDouble>;
template struct mx::astro::sqWaveFilter<astroSpectraSiDouble>;
template struct mx::astro::calspecSpectrum<astroSpectraSiDouble>;
template struct mx::astro::picklesSpectrum<astroSpectraSiDouble>;
template struct mx::astro::earthAlbedo<astroSpectraSiDouble>;
template struct mx::astro::venusAlbedo<astroSpectraSiDouble>;
/// \endcond

/** \defgroup astroSpectra_unit_tests astroSpectra Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_astroSpectra_test
{

/** \brief Verifies that the mx::astro spectrum-provider APIs are represented in the unit-test target.
 *
 * \ingroup astroSpectra_unit_tests
 * \todo Add behavioral assertions for the APIs declared in astroSpectra.hpp.
 */
TEST_CASE( "astroSpectra API has a test placeholder", "[astro::astroSpectra][placeholder]" )
{
    SUCCEED( "astroSpectra behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_astroSpectra_test
