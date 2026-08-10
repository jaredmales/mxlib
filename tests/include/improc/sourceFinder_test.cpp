/** \file sourceFinder_test.cpp
 * \brief Placeholder tests for the sourceFinder production API.
 */

#include "../../catch2/catch.hpp"

#include "../../../include/improc/sourceFinder.hpp"

/** \cond
 * Explicit instantiation compile-checks the representative double-precision source finder and emits its
 * header-defined methods for coverage accounting in this test translation unit.
 */
template class mx::improc::sourceFinder<double>;
/** \endcond */

/** \defgroup sourceFinder_unit_tests sourceFinder Unit Tests
 * \ingroup improc_unit_tests
 */

namespace unitTest::placeholder::improc_sourceFinder_test
{

/** \brief Records pending behavioral verification of the sourceFinder production API.
 *
 * \ingroup sourceFinder_unit_tests
 * \todo Add behavioral assertions for the production API declared in \c include/improc/sourceFinder.hpp.
 */
TEST_CASE( "sourceFinder production API placeholder", "[improc::sourceFinder][placeholder]" )
{
    SUCCEED( "Behavioral API assertions are pending." );
}

} // namespace unitTest::placeholder::improc_sourceFinder_test
