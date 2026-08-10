/** \file mftT_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/ft/mftT.hpp.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/ft/mftT.hpp"

/** \cond */
/// Compile the CPU rank-two matrix Fourier transform for representative scalar types.
template class mx::math::ft::mftT<double, std::complex<double>, 2, 0>;
/** \endcond */

/** \defgroup mftT_unit_tests mftT Unit Tests
 * \ingroup math_ft_unit_tests
 */

namespace unitTest::placeholder::math_ft_mftT_test
{

/** \brief Verifies that APIs declared in include/math/ft/mftT.hpp remain available to the unit-test target.
 *
 * \ingroup mftT_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/ft/mftT.hpp.
 */
TEST_CASE( "mftT API has a test placeholder", "[math::ft::mftT][placeholder]" )
{
    SUCCEED( "mftT behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_ft_mftT_test
