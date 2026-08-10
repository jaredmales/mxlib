/** \file array2FitGaussian2D_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/fit/array2FitGaussian2D.hpp.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/fit/array2FitGaussian2D.hpp"

/** \cond */
/// Compile the two-dimensional Gaussian array adapter for a representative scalar type.
template struct mx::math::fit::array2FitGaussian2D<double>;
/** \endcond */

/** \defgroup array2FitGaussian2D_unit_tests array2FitGaussian2D Unit Tests
 * \ingroup math_fit_unit_tests
 */

namespace unitTest::placeholder::math_fit_array2FitGaussian2D_test
{

/** \brief Verifies that APIs declared in include/math/fit/array2FitGaussian2D.hpp remain available to the unit-test
 * target.
 *
 * \ingroup array2FitGaussian2D_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/fit/array2FitGaussian2D.hpp.
 */
TEST_CASE( "array2FitGaussian2D API has a test placeholder", "[math::fit::array2FitGaussian2D][placeholder]" )
{
    SUCCEED( "array2FitGaussian2D behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_fit_array2FitGaussian2D_test
