/** \file exception_test.cpp
 * \brief Placeholder tests for mxlib exceptions.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/error/exception.hpp"

/** \cond
 * Explicit instantiation compile-checks the default verbose exception specialization and emits its header-defined
 * methods for coverage accounting in this test translation unit.
 */
template class mx::exception<mx::verbose::d>;
/** \endcond */

/** \defgroup exception_unit_tests Exception Unit Tests
 * \ingroup error_unit_tests
 */

namespace unitTest::placeholder::error_exception_test
{

/** \brief Verifies that the mx::exception API is represented in the unit-test target.
 *
 * \ingroup exception_unit_tests
 * \todo Add behavioral assertions for mx::exception.
 */
TEST_CASE( "exception API has a test placeholder", "[error::exception][placeholder]" )
{
    SUCCEED( "Exception behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::error_exception_test
