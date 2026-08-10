/** \file pywfsSlopeReconstructor_test.cpp
 * \brief Tests pyramid-sensor slope reconstruction.
 */
#include "../../../catch2/catch.hpp"

#include <filesystem>

#include "../../../../include/ao/sim/pywfsSlopeReconstructor.hpp"

/** \cond */
/// Compile the pyramid-sensor slope reconstructor for a representative scalar type.
template class mx::AO::sim::pywfsSlopeReconstructor<double>;
/** \endcond */

/** \defgroup pywfsSlopeReconstructor_unit_tests pywfsSlopeReconstructor Unit Tests
 * \ingroup ao_sim_unit_tests
 */

namespace unitTest::ao_sim_pywfsSlopeReconstructor_test
{

/** \cond */
class pywfsSlopeReconstructor_test : public mx::AO::sim::pywfsSlopeReconstructor<double>
{
  public:
    /// Exposes the loaded quadrant mask for verification.
    const imageT &quadMask() const
    {
        return _quadMask;
    }
};
/** \endcond */

/** \brief Verifies that mx::AO::sim::pywfsSlopeReconstructor::calcMask preserves a configured FITS mask.
 *
 * \ingroup pywfsSlopeReconstructor_unit_tests
 */
TEST_CASE( "pywfsSlopeReconstructor loads its configured quadrant mask", "[ao::sim::pywfsSlopeReconstructor]" )
{
    const std::filesystem::path testDirectory =
        std::filesystem::temp_directory_path() / "mxlib-pywfsSlopeReconstructor-test";
    std::filesystem::remove_all( testDirectory );
    std::filesystem::create_directories( testDirectory );

    const std::filesystem::path previousDirectory = std::filesystem::current_path();
    std::filesystem::current_path( testDirectory );

    mx::improc::eigenImage<double> expectedMask( 2, 3 );
    expectedMask << 1, 0, 1, 0, 1, 0;

    const std::filesystem::path maskFile = testDirectory / "quadrant-mask.fits";
    mx::fits::fitsFile<double> fitsFile;
    REQUIRE( fitsFile.write( maskFile.string(), expectedMask ) == mx::error_t::noerror );

    pywfsSlopeReconstructor_test reconstructor;
    reconstructor.detRows( 8 );
    reconstructor.detCols( 8 );
    reconstructor.maskFile( maskFile.string() );

    REQUIRE( reconstructor.measurementSize() == 6 );
    REQUIRE( reconstructor.quadMask().isApprox( expectedMask ) );

    std::filesystem::current_path( previousDirectory );
    std::filesystem::remove_all( testDirectory );
}

} // namespace unitTest::ao_sim_pywfsSlopeReconstructor_test
