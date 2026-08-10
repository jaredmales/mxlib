/** \file fitMoffat_test.cpp
 * \brief Placeholder tests for APIs declared in include/math/fit/fitMoffat.hpp.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/fit/fitMoffat.hpp"

/** \cond */
/// Compile the one- and two-dimensional Moffat fitters for representative scalar specializations.
template struct mx::math::fit::array2FitMoffat1D<double>;
template struct mx::math::fit::moffat1D_fitter<double>;
template class mx::math::fit::fitMoffat1D<double>;
template struct mx::math::fit::array2FitMoffat<double>;
template struct mx::math::fit::moffat2D_sym_fitter<double>;
template class mx::math::fit::fitMoffat2D<mx::math::fit::moffat2D_sym_fitter<double>>;
/** \endcond */

/** \defgroup fitMoffat_unit_tests fitMoffat Unit Tests
 * \ingroup math_fit_unit_tests
 */

namespace unitTest::placeholder::math_fit_fitMoffat_test
{

/** \brief Verifies that APIs declared in include/math/fit/fitMoffat.hpp remain available to the unit-test target.
 *
 * \ingroup fitMoffat_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/math/fit/fitMoffat.hpp.
 */
TEST_CASE( "fitMoffat API has a test placeholder", "[math::fit::fitMoffat][placeholder]" )
{
    SUCCEED( "fitMoffat behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::math_fit_fitMoffat_test
