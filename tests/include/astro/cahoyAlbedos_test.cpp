/** \file cahoyAlbedos_test.cpp
 * \brief Placeholder tests for Cahoy-model albedo spectra.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/cahoyAlbedos.hpp"

/// \cond
// Compile every non-template member of representative SI-unit specializations.
using cahoyAlbedosSiDouble = mx::astro::units::si<double>;
template struct mx::astro::cahoySpectrumRaw<cahoyAlbedosSiDouble>;
template struct mx::astro::cahoyGrid<cahoyAlbedosSiDouble>;
/// \endcond

/** \defgroup cahoyAlbedos_unit_tests Cahoy Albedo Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_cahoyAlbedos_test
{

/** \brief Verifies that the mx::astro Cahoy albedo APIs are represented in the unit-test target.
 *
 * \ingroup cahoyAlbedos_unit_tests
 * \todo Add behavioral assertions for mx::astro::cahoySpectrumRaw and mx::astro::cahoyGrid.
 */
TEST_CASE( "Cahoy albedo APIs have a test placeholder", "[astro::cahoyAlbedos][placeholder]" )
{
    SUCCEED( "Cahoy albedo behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_cahoyAlbedos_test
