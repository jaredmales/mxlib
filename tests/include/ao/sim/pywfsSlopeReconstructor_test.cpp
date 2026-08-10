/** \file pywfsSlopeReconstructor_test.cpp
 * \brief Placeholder tests for pyramid-sensor slope reconstruction.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/pywfsSlopeReconstructor.hpp"

/** \cond */
/// Compile the pyramid-sensor slope reconstructor for a representative scalar type.
template class mx::AO::sim::pywfsSlopeReconstructor<double>;
/** \endcond */

/** \defgroup pywfsSlopeReconstructor_unit_tests pywfsSlopeReconstructor Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::placeholder::ao_sim_pywfsSlopeReconstructor_test
{

/** \brief Verifies that the pywfsSlopeReconstructor API is represented in the unit-test target.
 *
 * \ingroup pywfsSlopeReconstructor_unit_tests
 * \todo Add behavioral assertions for the APIs declared in mx::AO::sim::pywfsSlopeReconstructor.
 */
TEST_CASE( "pywfsSlopeReconstructor API has a test placeholder", "[ao::sim::pywfsSlopeReconstructor][placeholder]" )
{
    SUCCEED( "pywfsSlopeReconstructor behavioral assertions are pending." );
}

} // namespace unitTest::placeholder::ao_sim_pywfsSlopeReconstructor_test
