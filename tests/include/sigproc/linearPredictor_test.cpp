/** \file linearPredictor_test.cpp
 * \brief Placeholder tests for APIs declared in include/sigproc/linearPredictor.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sigproc/linearPredictor.hpp"

/** \cond */
/// Compile the linear-predictor implementation for a representative scalar type.
template struct mx::sigproc::linearPredictor<double>;
/** \endcond */

/** \defgroup linearPredictor_unit_tests linearPredictor Unit Tests
 * \ingroup sigproc_unit_tests
 */

namespace unitTest::placeholder::sigproc_linearPredictor_test
{

/** \brief Verifies that APIs declared in include/sigproc/linearPredictor.hpp remain available to the unit-test target.
 *
 * \ingroup linearPredictor_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/sigproc/linearPredictor.hpp.
 */
TEST_CASE( "linearPredictor API has a test placeholder", "[sigproc::linearPredictor][placeholder]" )
{
    SUCCEED( "linearPredictor behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::sigproc_linearPredictor_test
