/** \file fitGammaDistribution_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/fit/fitGammaDistribution.hpp.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/fit/fitGammaDistribution.hpp"

/** \cond */
/// Compile the four-parameter fitter's theta setter without instantiating incompatible setGuess overloads.
template void
mx::math::fit::fitGammaDistribution<mx::math::fit::gammaDistribution_4param_fitter<double>>::theta( double );
/** \endcond */

/** \defgroup fitGammaDistribution_unit_tests fitGammaDistribution Unit Tests
 * \ingroup math_fit_unit_tests
 */

namespace unitTest::placeholder::math_fit_fitGammaDistribution_test
{

/** \brief Verifies that APIs declared in include/math/fit/fitGammaDistribution.hpp remain available to the unit-test
 * target.
 *
 * \ingroup fitGammaDistribution_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/fit/fitGammaDistribution.hpp.
 */
TEST_CASE( "fitGammaDistribution API has a test placeholder", "[math::fit::fitGammaDistribution][placeholder]" )
{
    SUCCEED( "fitGammaDistribution behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_fit_fitGammaDistribution_test
