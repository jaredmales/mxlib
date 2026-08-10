/** \file fitEmpirical_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/fit/fitEmpirical.hpp.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/fit/fitEmpirical.hpp"

/** \cond */
/// Compile the empirical fitter and its array adapter for a representative specialization.
template struct mx::math::fit::array2FitEmpirical<double>;
template struct mx::math::fit::empirical2D_fitter<double>;
template class mx::math::fit::fitEmpirical2DGen<mx::math::fit::empirical2D_fitter<double>>;
/** \endcond */

/** \defgroup fitEmpirical_unit_tests fitEmpirical Unit Tests
 * \ingroup math_fit_unit_tests
 */

namespace unitTest::placeholder::math_fit_fitEmpirical_test
{

/** \brief Verifies that APIs declared in include/math/fit/fitEmpirical.hpp remain available to the unit-test target.
 *
 * \ingroup fitEmpirical_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/fit/fitEmpirical.hpp.
 */
TEST_CASE( "fitEmpirical API has a test placeholder", "[math::fit::fitEmpirical][placeholder]" )
{
    SUCCEED( "fitEmpirical behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_fit_fitEmpirical_test
