/** \file directPhaseSensor_test.cpp
 * \brief Tests for simulated direct phase sensors.
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/ao/sim/ccdDetector.hpp"
#include "../../../../include/ao/sim/directPhaseSensor.hpp"

/// \cond
// Compile every defined non-template member of a representative CPU specialization.
template class mx::AO::sim::directPhaseSensor<double, mx::AO::sim::ccdDetector<double>>;
/// \endcond

/** \defgroup directPhaseSensor_unit_tests directPhaseSensor Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::ao_sim_directPhaseSensor_test
{

/** \brief Verifies that the directPhaseSensor detector-dimension setters preserve both dimensions.
 *
 * \ingroup directPhaseSensor_unit_tests
 */
TEST_CASE( "directPhaseSensor detector dimensions remain synchronized", "[ao::sim::directPhaseSensor]" )
{
    mx::AO::sim::directPhaseSensor<double, mx::AO::sim::ccdDetector<double>> sensor;

    sensor.detSize( 4, 5 );
    REQUIRE( sensor.detRows() == 4 );
    REQUIRE( sensor.detCols() == 5 );

    sensor.detRows( 6 );
    REQUIRE( sensor.detRows() == 6 );
    REQUIRE( sensor.detCols() == 5 );

    sensor.detCols( 7 );
    REQUIRE( sensor.detRows() == 6 );
    REQUIRE( sensor.detCols() == 7 );
}

} // namespace unitTest::ao_sim_directPhaseSensor_test
