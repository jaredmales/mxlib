/** \file expModGaussian_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/func/expModGaussian.hpp.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/func/expModGaussian.hpp"

/** \cond */
/// Compile the exponentially modified Gaussian mode functor for a representative scalar type.
template struct mx::math::func::emgModeFunc<double>;
/** \endcond */

/** \defgroup expModGaussian_unit_tests expModGaussian Unit Tests
 * \ingroup math_func_unit_tests
 */

namespace unitTest::placeholder::math_func_expModGaussian_test
{

/** \brief Verifies that APIs declared in include/math/func/expModGaussian.hpp remain available to the unit-test target.
 *
 * \ingroup expModGaussian_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/func/expModGaussian.hpp.
 */
TEST_CASE( "expModGaussian API has a test placeholder", "[math::func::expModGaussian][placeholder]" )
{
    SUCCEED( "expModGaussian behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_func_expModGaussian_test
