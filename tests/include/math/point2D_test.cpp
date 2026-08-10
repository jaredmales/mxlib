/** \file point2D_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/point2D.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/math/geo.hpp"
#include "../../../include/math/point2D.hpp"

/** \cond */
/// Compile the two-dimensional point implementation with its documented radians policy.
template class mx::math::point2D<mx::math::radiansT<double>>;
/** \endcond */

/** \defgroup point2D_unit_tests point2D Unit Tests
 * \ingroup math_unit_tests
 */

namespace unitTest::placeholder::math_point2D_test
{

/** \brief Verifies that APIs declared in include/math/point2D.hpp remain available to the unit-test target.
 *
 * \ingroup point2D_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/point2D.hpp.
 */
TEST_CASE( "point2D API has a test placeholder", "[math::point2D][placeholder]" )
{
    SUCCEED( "point2D behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_point2D_test
