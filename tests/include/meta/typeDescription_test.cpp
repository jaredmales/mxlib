/** \file typeDescription_test.cpp
 * \brief Placeholder tests for APIs declared in include/meta/typeDescription.hpp.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/meta/typeDescription.hpp"

/** \cond */
namespace unitTest::compileAnchor
{
struct typeDescriptionProbe
{
};
} // namespace unitTest::compileAnchor

/// Compile the generic type-description implementation without selecting an explicit specialization.
template struct mx::meta::typeDescription<unitTest::compileAnchor::typeDescriptionProbe>;
/** \endcond */

/** \defgroup typeDescription_unit_tests typeDescription Unit Tests
 * \ingroup meta_unit_tests
 */

namespace unitTest::placeholder::meta_typeDescription_test
{

/** \brief Verifies that APIs declared in include/meta/typeDescription.hpp remain available to the unit-test target.
 *
 * \ingroup typeDescription_unit_tests
 * \todo Add behavioral assertions for the APIs declared in include/meta/typeDescription.hpp.
 */
TEST_CASE( "typeDescription API has a test placeholder", "[meta::typeDescription][placeholder]" )
{
    SUCCEED( "typeDescription behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::meta_typeDescription_test
