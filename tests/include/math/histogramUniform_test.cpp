/** \file histogramUniform_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/histogramUniform.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/math/histogramUniform.hpp"

/** \cond */
/// Compile the uniform histogram implementation for a representative scalar type.
template class mx::math::histogramUniform<double>;
/** \endcond */

/** \defgroup histogramUniform_unit_tests histogramUniform Unit Tests
 * \ingroup math_unit_tests
 */

namespace unitTest::placeholder::math_histogramUniform_test
{

/** \brief Verifies that APIs declared in include/math/histogramUniform.hpp remain available to the unit-test target.
 *
 * \ingroup histogramUniform_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/histogramUniform.hpp.
 */
TEST_CASE( "histogramUniform API has a test placeholder", "[math::histogramUniform][placeholder]" )
{
    SUCCEED( "histogramUniform behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_histogramUniform_test
