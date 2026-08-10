/** \file changeable_test.cpp
 * \brief Placeholder tests for change tracking.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/base/changeable.hpp"

/// \cond
namespace unitTest::compile
{
struct changeableDerived
{
};
} // namespace unitTest::compile

// Compile every non-template member of a representative CRTP specialization.
template class mx::base::changeable<unitTest::compile::changeableDerived>;
/// \endcond

/** \defgroup changeable_unit_tests changeable Unit Tests
 * \ingroup base_unit_tests
 */

namespace unitTest::placeholder::base_changeable_test
{

/** \brief Verifies that the mx::base::changeable API is represented in the unit-test target.
 *
 * \ingroup changeable_unit_tests
 * \todo Add behavioral assertions for mx::base::changeable.
 */
TEST_CASE( "changeable API has a test placeholder", "[base::changeable][placeholder]" )
{
    SUCCEED( "changeable behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::base_changeable_test
