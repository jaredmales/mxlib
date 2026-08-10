/** \file psdVarMean_test.cpp
 * \brief Placeholder tests for APIs declared in include/sigproc/psdVarMean.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sigproc/psdVarMean.hpp"

/** \cond */
/// Compile the PSD variance-of-the-mean evaluator with its standard parameter type.
template struct mx::sigproc::psdVarMeanParams<double>;
template struct mx::sigproc::psdVarMean<mx::sigproc::psdVarMeanParams<double>>;
/** \endcond */

/** \defgroup psdVarMean_unit_tests psdVarMean Unit Tests
 * \ingroup sigproc_unit_tests
 */

namespace unitTest::placeholder::sigproc_psdVarMean_test
{

/** \brief Verifies that APIs declared in include/sigproc/psdVarMean.hpp remain available to the unit-test target.
 *
 * \ingroup psdVarMean_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/sigproc/psdVarMean.hpp.
 */
TEST_CASE( "psdVarMean API has a test placeholder", "[sigproc::psdVarMean][placeholder]" )
{
    SUCCEED( "psdVarMean behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::sigproc_psdVarMean_test
