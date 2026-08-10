/** \file pupil_test.cpp
 * \brief Placeholder tests for adaptive-optics pupil utilities.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/ao/pupil.hpp"

/** \defgroup pupil_unit_tests AO Pupil Unit Tests
 * \ingroup ao_unit_tests
 */

namespace unitTest::placeholder::ao_pupil_test
{

/** \brief Verifies that the adaptive-optics pupil API is represented in the unit-test target.
 *
 * \ingroup pupil_unit_tests
 * \todo Add behavioral assertions for mx::AO::circularPupil and mx::AO::circularApodizedPupil.
 */
TEST_CASE( "adaptive-optics pupil API has a test placeholder", "[ao::pupil][placeholder]" )
{
    SUCCEED( "Adaptive-optics pupil behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_pupil_test
