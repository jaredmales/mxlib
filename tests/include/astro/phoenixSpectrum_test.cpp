/** \file phoenixSpectrum_test.cpp
 * \brief Placeholder tests for PHOENIX model spectra.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/astro/units.hpp"
#include "../../../include/astro/phoenixSpectrum.hpp"

/// \cond
// Compile the representative provider and both header-defined maintenance function templates.
template struct mx::astro::phoenixSpectrum<mx::astro::units::si<double>>;
template void mx::astro::rewritePhoenixSpectrum<double>( const std::string &, double, double, int, double );
template void mx::astro::rewritePhoenixSpectrumBatch<double>( const std::string &, double, double, double );
/// \endcond

/** \defgroup phoenixSpectrum_unit_tests phoenixSpectrum Unit Tests
 * \ingroup astro_unit_tests
 */

namespace unitTest::placeholder::astro_phoenixSpectrum_test
{

/** \brief Verifies that the mx::astro::phoenixSpectrum API is represented in the unit-test target.
 *
 * \ingroup phoenixSpectrum_unit_tests
 * \todo Add behavioral assertions for mx::astro::phoenixSpectrum.
 */
TEST_CASE( "phoenixSpectrum API has a test placeholder", "[astro::phoenixSpectrum][placeholder]" )
{
    SUCCEED( "phoenixSpectrum behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::astro_phoenixSpectrum_test
