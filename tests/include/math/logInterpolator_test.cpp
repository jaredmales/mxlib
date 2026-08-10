/** \file logInterpolator_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/logInterpolator.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/math/logInterpolator.hpp"

/** \cond */
/// Compile logarithmic interpolation with the supported linear GSL policy.
template class mx::math::logInterpolator<mx::math::gsl_interp_linear<double>>;
/** \endcond */

/** \defgroup logInterpolator_unit_tests logInterpolator Unit Tests
 * \ingroup math_unit_tests
 */

namespace unitTest::placeholder::math_logInterpolator_test
{

/** \brief Verifies that APIs declared in include/math/logInterpolator.hpp remain available to the unit-test target.
 *
 * \ingroup logInterpolator_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/logInterpolator.hpp.
 */
TEST_CASE( "logInterpolator API has a test placeholder", "[math::logInterpolator][placeholder]" )
{
    SUCCEED( "logInterpolator behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_logInterpolator_test
