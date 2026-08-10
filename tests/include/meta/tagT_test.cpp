/** \file tagT_test.cpp
 * \brief Placeholder tests for APIs declared in include/meta/tagT.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/meta/tagT.hpp"

/** \cond */
/// Compile the tag wrapper for a representative payload type.
template struct mx::meta::tagT<int>;
/** \endcond */

/** \defgroup tagT_unit_tests tagT Unit Tests
 * \ingroup meta_unit_tests
 */

namespace unitTest::placeholder::meta_tagT_test
{

/** \brief Verifies that APIs declared in include/meta/tagT.hpp remain available to the unit-test target.
 *
 * \ingroup tagT_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/meta/tagT.hpp.
 */
TEST_CASE( "tagT API has a test placeholder", "[meta::tagT][placeholder]" )
{
    SUCCEED( "tagT behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::meta_tagT_test
