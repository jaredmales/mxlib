/** \file autocorrelation_test.cpp
 * \brief Placeholder tests for APIs declared in include/sigproc/autocorrelation.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sigproc/autocorrelation.hpp"

/** \cond */
/// Compile the PSD-to-autocorrelation functor for a representative scalar type.
template struct mx::sigproc::autocorrelationFromPSD<double>;
/** \endcond */

/** \defgroup autocorrelation_unit_tests autocorrelation Unit Tests
 * \ingroup sigproc_unit_tests
 */

namespace unitTest::placeholder::sigproc_autocorrelation_test
{

/** \brief Verifies that APIs declared in include/sigproc/autocorrelation.hpp remain available to the unit-test target.
 *
 * \ingroup autocorrelation_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/sigproc/autocorrelation.hpp.
 */
TEST_CASE( "autocorrelation API has a test placeholder", "[sigproc::autocorrelation][placeholder]" )
{
    SUCCEED( "autocorrelation behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::sigproc_autocorrelation_test
