/** \file averagePeriodogram_test.cpp
 * \brief Placeholder tests for APIs declared in include/sigproc/averagePeriodogram.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sigproc/averagePeriodogram.hpp"

/** \cond */
/// Compile the average-periodogram implementation for a representative scalar type.
template class mx::sigproc::averagePeriodogram<double>;
/** \endcond */

/** \defgroup averagePeriodogram_unit_tests averagePeriodogram Unit Tests
 * \ingroup sigproc_unit_tests
 */

namespace unitTest::placeholder::sigproc_averagePeriodogram_test
{

/** \brief Verifies that APIs declared in include/sigproc/averagePeriodogram.hpp remain available to the unit-test
 * target.
 *
 * \ingroup averagePeriodogram_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/sigproc/averagePeriodogram.hpp.
 */
TEST_CASE( "averagePeriodogram API has a test placeholder", "[sigproc::averagePeriodogram][placeholder]" )
{
    SUCCEED( "averagePeriodogram behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::sigproc_averagePeriodogram_test
