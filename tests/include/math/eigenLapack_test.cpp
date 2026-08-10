/** \file eigenLapack_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/eigenLapack.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/math/eigenLapack.hpp"

/** \cond */
/// Compile the symmetric-eigensolver workspace implementation for a representative scalar type.
template struct mx::math::syevrMem<double>;
/** \endcond */

/** \defgroup eigenLapack_unit_tests eigenLapack Unit Tests
 * \ingroup math_unit_tests
 */

namespace unitTest::placeholder::math_eigenLapack_test
{

/** \brief Verifies that APIs declared in include/math/eigenLapack.hpp remain available to the unit-test target.
 *
 * \ingroup eigenLapack_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/eigenLapack.hpp.
 */
TEST_CASE( "eigenLapack API has a test placeholder", "[math::eigenLapack][placeholder]" )
{
    SUCCEED( "eigenLapack behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_eigenLapack_test
