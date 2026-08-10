/** \file fourierTemporalPSD_test.cpp
 * \brief Tests of Fourier-mode temporal power spectral densities.
 */

#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/fourierTemporalPSD.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <unistd.h>

using aoSystemBaseT = mx::AO::analysis::aoSystem<double, mx::AO::analysis::vonKarmanSpectrum<double>>;

/// AO-system test type that forces local template instantiation under sanitizers.
struct aoSystemT : public aoSystemBaseT
{
};

using temporalPsdT = mx::AO::analysis::fourierTemporalPSD<double, aoSystemT>;

namespace
{
/// \cond fourierTemporalPSD_test_detail

/// Sentinel GSL error handler used to verify restoration after PSD calculations.
void testGslErrorHandler( const char *reason, /**< [in] GSL diagnostic text */
                          const char *file,   /**< [in] GSL source file */
                          int line,           /**< [in] GSL source line */
                          int errorNumber /**< [in] GSL status code */ )
{
    static_cast<void>( reason );
    static_cast<void>( file );
    static_cast<void>( line );
    static_cast<void>( errorNumber );
}

/// Return a null workspace to exercise allocation-failure propagation.
gsl_integration_workspace *failWorkspaceAllocation( size_t size /**< [in] requested workspace size */ )
{
    static_cast<void>( size );
    return nullptr;
}

/// Test evaluator exposing the protected allocator-injection constructor.
struct allocationFailureTemporalPsdT : public temporalPsdT
{
    /// Construct an evaluator whose workspace allocation always fails.
    allocationFailureTemporalPsdT() : temporalPsdT( &failWorkspaceAllocation )
    {
    }
};

/// Test evaluator exposing workspace allocation and ownership state.
struct workspaceOwnershipTemporalPsdT : public temporalPsdT
{
    /// Preserve move construction for the derived test evaluator.
    workspaceOwnershipTemporalPsdT( workspaceOwnershipTemporalPsdT && ) noexcept = default;

    /// Construct a test evaluator without allocating a workspace.
    workspaceOwnershipTemporalPsdT() = default;

    /// Allocate the workspace through the production helper.
    mx::error_t allocateTestWorkspace()
    {
        return allocateWorkspace();
    }

    /// Return whether this evaluator currently owns a workspace.
    bool ownsWorkspace() const
    {
        return m_workspace != nullptr;
    }
};

/// \endcond

} // namespace

/// Verify that Fourier temporal PSD evaluators uniquely own movable GSL workspaces.
/** Exercises the ownership contract of mx::AO::analysis::fourierTemporalPSD. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD workspace ownership is move-only", "[ao::analysis::fourierTemporalPSD]" )
{
    STATIC_REQUIRE_FALSE( std::is_copy_constructible_v<temporalPsdT> );
    STATIC_REQUIRE_FALSE( std::is_copy_assignable_v<temporalPsdT> );
    STATIC_REQUIRE( std::is_nothrow_move_constructible_v<temporalPsdT> );
    STATIC_REQUIRE( std::is_nothrow_move_assignable_v<temporalPsdT> );

    workspaceOwnershipTemporalPsdT source;
    REQUIRE( source.allocateTestWorkspace() == mx::error_t::noerror );
    REQUIRE( source.ownsWorkspace() );

    workspaceOwnershipTemporalPsdT destination( std::move( source ) );
    REQUIRE_FALSE( source.ownsWorkspace() );
    REQUIRE( destination.ownsWorkspace() );
    REQUIRE( destination.m_aosys == nullptr );
}

/// Verify that public Fourier temporal PSD calculations reject malformed inputs before modifying output.
/** Exercises mx::AO::analysis::fourierTemporalPSD::singleLayerPSD and
 * mx::AO::analysis::fourierTemporalPSD::multiLayerPSD precondition handling. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD validates public calculation inputs", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.lam_sci( 1.0e-6 );
    aoSystem.lam_wfs( 0.8e-6 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    temporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = mx::AO::analysis::basis::basic;

    std::vector<double> frequency{ 0.0, 1.0 };
    std::vector<double> psd( frequency.size(), 7.0 );
    const std::vector<double> unchanged = psd;
    temporalPsdT::reportT report;
    report.record( GSL_SUCCESS, 0, 0.0, 0.0, 0.0, 1.0, 1.0 );

    SECTION( "null AO system" )
    {
        temporalPsd.m_aosys = nullptr;
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, frequency.back(), &report ) ==
                 mx::error_t::invalidconfig );
    }

    SECTION( "empty vectors" )
    {
        std::vector<double> empty;
        REQUIRE( temporalPsd.singleLayerPSD( empty, empty, 1.0, 0.0, 0, 1, 0, &report ) == mx::error_t::sizeerr );
    }

    SECTION( "mismatched vector sizes" )
    {
        std::vector<double> undersizedPsd{ 7.0 };
        REQUIRE( temporalPsd.singleLayerPSD( undersizedPsd, frequency, 1.0, 0.0, 0, 1, frequency.back(), &report ) ==
                 mx::error_t::sizeerr );
        REQUIRE( undersizedPsd == std::vector<double>{ 7.0 } );
    }

    SECTION( "nonmonotone frequency grid" )
    {
        frequency[1] = frequency[0];
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) == mx::error_t::invalidarg );
    }

    SECTION( "negative frequency" )
    {
        frequency[0] = -1.0;
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) == mx::error_t::invalidarg );
    }

    SECTION( "nonfinite frequency grid" )
    {
        frequency[1] = std::numeric_limits<double>::quiet_NaN();
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) == mx::error_t::invalidarg );
    }

    SECTION( "nonfinite mode" )
    {
        REQUIRE(
            temporalPsd
                .singleLayerPSD( psd, frequency, std::numeric_limits<double>::infinity(), 0.0, 0, 1, 0, &report ) ==
            mx::error_t::invalidarg );
    }

    SECTION( "invalid parity" )
    {
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 0, 0, &report ) == mx::error_t::invalidarg );
    }

    SECTION( "negative cutoff" )
    {
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, -1.0, &report ) ==
                 mx::error_t::invalidarg );
    }

    SECTION( "invalid layer index" )
    {
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 1, 1, 0, &report ) == mx::error_t::invalidarg );
    }

    SECTION( "invalid tolerance" )
    {
        temporalPsd.absTol( 0 );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) ==
                 mx::error_t::invalidconfig );
    }

    SECTION( "invalid relative tolerance" )
    {
        temporalPsd.relTol( 1 );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) ==
                 mx::error_t::invalidconfig );
    }

    SECTION( "nonpositive aperture" )
    {
        aoSystem.D( 0 );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) ==
                 mx::error_t::invalidconfig );
    }

    SECTION( "zero wind" )
    {
        aoSystem.atm.layer_v_wind( std::vector<double>{ 0.0 } );
        REQUIRE( temporalPsd.multiLayerPSD<false>( psd, frequency, 1.0, 0.0, 1, 0, &report ) ==
                 mx::error_t::invalidconfig );
    }

    SECTION( "mismatched atmosphere vectors" )
    {
        aoSystem.atm.layer_dir( std::vector<double>{} );
        REQUIRE( temporalPsd.multiLayerPSD<false>( psd, frequency, 1.0, 0.0, 1, 0, &report ) == mx::error_t::sizeerr );
    }

    REQUIRE( psd == unchanged );
    REQUIRE( report.integrationsAttempted == 0 );
}

/// Verify that GSL workspace allocation failure is returned without evaluating or modifying the PSD.
/** Exercises mx::AO::analysis::fourierTemporalPSD::singleLayerPSD with an injected failing workspace allocator. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD reports workspace allocation failure", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    allocationFailureTemporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = mx::AO::analysis::basis::basic;

    std::vector<double> frequency{ 1.0 };
    std::vector<double> psd{ 7.0 };
    temporalPsdT::reportT report;
    REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, 0, &report ) == mx::error_t::allocerr );
    REQUIRE( psd == std::vector<double>{ 7.0 } );
    REQUIRE( report.integrationsAttempted == 0 );
}

/// Verify that PSD-grid generation validates its public inputs before creating output.
/** Exercises mx::AO::analysis::fourierTemporalPSD::makePSDGrid precondition handling. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD validates grid-generation inputs", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    temporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = mx::AO::analysis::basis::basic;

    const std::filesystem::path outputDirectory =
        std::filesystem::temp_directory_path() /
        ( "mxlib_fourierTemporalPSD_invalid_grid_" + std::to_string( static_cast<long long>( getpid() ) ) );
    std::error_code filesystemError;
    std::filesystem::remove_all( outputDirectory, filesystemError );
    REQUIRE_FALSE( filesystemError );

    SECTION( "empty output directory" )
    {
        REQUIRE( temporalPsd.makePSDGrid( "", 1, 1.0, 1.0 ) == mx::error_t::invalidarg );
    }

    SECTION( "nonpositive spatial extent" )
    {
        REQUIRE( temporalPsd.makePSDGrid( outputDirectory.string(), 0, 1.0, 1.0 ) == mx::error_t::invalidarg );
    }

    SECTION( "nonpositive frequency spacing" )
    {
        REQUIRE( temporalPsd.makePSDGrid( outputDirectory.string(), 1, 0.0, 1.0 ) == mx::error_t::invalidarg );
    }

    SECTION( "nonfinite frequency spacing" )
    {
        REQUIRE(
            temporalPsd.makePSDGrid( outputDirectory.string(), 1, std::numeric_limits<double>::quiet_NaN(), 1.0 ) ==
            mx::error_t::invalidarg );
    }

    SECTION( "nonpositive maximum frequency" )
    {
        REQUIRE( temporalPsd.makePSDGrid( outputDirectory.string(), 1, 1.0, 0.0 ) == mx::error_t::invalidarg );
    }

    SECTION( "negative exact-calculation cutoff" )
    {
        REQUIRE( temporalPsd.makePSDGrid( outputDirectory.string(), 1, 1.0, 1.0, -1.0 ) == mx::error_t::invalidarg );
    }

    SECTION( "unrepresentable sample count" )
    {
        REQUIRE( temporalPsd.makePSDGrid( outputDirectory.string(), 1, std::numeric_limits<double>::min(), 1.0 ) ==
                 mx::error_t::sizeerr );
    }

    SECTION( "invalid AO system" )
    {
        temporalPsd.m_aosys = nullptr;
        REQUIRE( temporalPsd.makePSDGrid( outputDirectory.string(), 1, 1.0, 1.0 ) == mx::error_t::invalidconfig );
    }

    REQUIRE_FALSE( std::filesystem::exists( outputDirectory ) );
}

/// Verify that PSD-grid generation propagates calculation and output failures.
/** Exercises mx::AO::analysis::fourierTemporalPSD::makePSDGrid failure handling after successful preflight. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD reports grid-generation failures", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    const std::filesystem::path outputPath =
        std::filesystem::temp_directory_path() /
        ( "mxlib_fourierTemporalPSD_grid_failure_" + std::to_string( static_cast<long long>( getpid() ) ) );
    std::error_code filesystemError;
    std::filesystem::remove_all( outputPath, filesystemError );
    REQUIRE_FALSE( filesystemError );

    SECTION( "mode calculation failure" )
    {
        allocationFailureTemporalPsdT temporalPsd;
        temporalPsd.m_aosys = &aoSystem;
        temporalPsd._useBasis = mx::AO::analysis::basis::basic;

        REQUIRE( temporalPsd.makePSDGrid( outputPath.string(), 1, 1.0, 1.0 ) == mx::error_t::allocerr );
        REQUIRE( std::filesystem::exists( outputPath / "params.txt" ) );
        REQUIRE( std::filesystem::exists( outputPath / "psds" / "freq.binv" ) );
    }

    SECTION( "output path is a regular file" )
    {
        std::ofstream outputFile( outputPath );
        REQUIRE( outputFile.is_open() );
        outputFile.close();

        temporalPsdT temporalPsd;
        temporalPsd.m_aosys = &aoSystem;
        temporalPsd._useBasis = mx::AO::analysis::basis::basic;
        REQUIRE( temporalPsd.makePSDGrid( outputPath.string(), 1, 1.0, 1.0 ) == mx::error_t::enotdir );
    }

    std::filesystem::remove_all( outputPath, filesystemError );
    REQUIRE_FALSE( filesystemError );
}

/// Verify that Fourier PSD calculation requires a valid atmosphere and accepts every complete preset.
/** Exercises mx::AO::analysis::aoAtmosphere::validate,
 * mx::AO::analysis::aoAtmosphere::loadGuyon2005, mx::AO::analysis::aoAtmosphere::loadLCO,
 * mx::AO::analysis::aoAtmosphere::setSingleLayer, and
 * mx::AO::analysis::fourierTemporalPSD::singleLayerPSD. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD enforces atmosphere validity at calculation", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );

    temporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = mx::AO::analysis::basis::basic;

    std::vector<double> frequency{ 1.0 };
    std::vector<double> psd{ 7.0 };

    SECTION( "default atmosphere" )
    {
        REQUIRE( aoSystem.atm.validate() == mx::error_t::sizeerr );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, frequency.back() ) ==
                 mx::error_t::sizeerr );
        REQUIRE( psd == std::vector<double>{ 7.0 } );
    }

    SECTION( "Guyon 2005 preset with infinite outer scale" )
    {
        aoSystem.atm.loadGuyon2005();
        REQUIRE( aoSystem.atm.L_0( 0 ) == 0.0 );
        REQUIRE( aoSystem.atm.validate() == mx::error_t::noerror );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, frequency.back() ) ==
                 mx::error_t::noerror );
    }

    SECTION( "LCO preset" )
    {
        aoSystem.atm.loadLCO();
        REQUIRE( aoSystem.atm.validate() == mx::error_t::noerror );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, frequency.back() ) ==
                 mx::error_t::noerror );
    }

    SECTION( "single-layer preset" )
    {
        aoSystem.atm.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );
        REQUIRE( aoSystem.atm.validate() == mx::error_t::noerror );
        REQUIRE( temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, frequency.back() ) ==
                 mx::error_t::noerror );
    }
}

/// Verify temporal-PSD tail initialization at and below its nominal averaging width.
/** Exercises zero, one, 49, and 50 exactly integrated frequency bins. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD handles short exact tails", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    temporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = mx::AO::analysis::basis::basic;

    constexpr std::size_t frequencyCount = 55;
    std::vector<double> frequency( frequencyCount );
    for( std::size_t index = 0; index < frequency.size(); ++index )
    {
        frequency[index] = static_cast<double>( index + 1 );
    }

    for( const std::size_t exactCount : { 0U, 1U, 49U, 50U } )
    {
        DYNAMIC_SECTION( exactCount << " exact bins" )
        {
            std::vector<double> psd( frequency.size(), 0.0 );
            const double maximumExactFrequency = exactCount == 0 ? 0.5 : frequency[exactCount - 1];
            temporalPsdT::reportT report;

            const mx::error_t result =
                temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, maximumExactFrequency, &report );

            if( exactCount == 0 )
            {
                REQUIRE( result == mx::error_t::invalidarg );
                REQUIRE( report.integrationsAttempted == 0 );
                continue;
            }

            REQUIRE( result == mx::error_t::noerror );
            REQUIRE( report.integrationsAttempted == exactCount );

            const std::size_t averageCount = std::min<std::size_t>( exactCount, 50 );
            const double exponent = aoSystem.atm.alpha( 0 ) + 2.0;
            double expected = 0.0;
            for( std::size_t offset = averageCount; offset > 0; --offset )
            {
                const std::size_t index = exactCount - offset;
                expected += psd[index] * std::pow( frequency[index] / frequency[exactCount], exponent );
            }
            expected /= static_cast<double>( averageCount );

            REQUIRE( psd[exactCount] == Approx( expected ).epsilon( 1e-12 ) );
            REQUIRE( std::isfinite( psd[exactCount] ) );
            REQUIRE( psd[exactCount] >= 0.0 );

            const double expectedLast =
                psd[exactCount] * std::pow( frequency[exactCount] / frequency.back(), exponent );
            REQUIRE( psd.back() == Approx( expectedLast ).epsilon( 1e-12 ) );
        }
    }
}

/// Verify aggregation and formatting of Fourier temporal-PSD quadrature diagnostics.
/** Exercises `fourierTemporalPSDReport::record`, `fourierTemporalPSDReport::merge`, and
 * `fourierTemporalPSDReport::write` with multiple GSL statuses and atmospheric layers. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD aggregates quadrature diagnostics", "[ao::analysis::fourierTemporalPSD]" )
{
    temporalPsdT::reportT report;
    report.record( GSL_SUCCESS, 0, 1.0, 2.0, 0.1, 0.1, 0.01 );
    report.record( GSL_EROUND, 0, 2.0, 2.0, 0.4, 0.1, 0.01 );
    report.record( GSL_EROUND, 1, 3.0, 4.0, 1.0, 0.1, 0.01 );

    temporalPsdT::reportT other;
    other.record( GSL_EDIVERGE, 2, 4.0, 1.0, 0.5, 0.1, 0.01 );
    report.merge( other );

    REQUIRE( report.integrationsAttempted == 4 );
    REQUIRE( report.integrationsConverged == 1 );
    REQUIRE( report.failureCount() == 3 );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).count == 2 );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).countByLayer.at( 0 ) == 1 );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).countByLayer.at( 1 ) == 1 );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).maximumAbsoluteError == Approx( 1.0 ) );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).maximumToleranceRatio == Approx( 10.0 ) );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).worstLayer == 1 );
    REQUIRE( report.gslStatus.at( GSL_EROUND ).worstFrequency == Approx( 3.0 ) );
    REQUIRE( report.gslStatus.at( GSL_EDIVERGE ).count == 1 );

    std::ostringstream summary;
    report.write( summary );
    REQUIRE( summary.str().find( "1/4 converged" ) != std::string::npos );
    REQUIRE( summary.str().find( gsl_strerror( GSL_EROUND ) ) != std::string::npos );
    REQUIRE( summary.str().find( gsl_strerror( GSL_EDIVERGE ) ) != std::string::npos );
    REQUIRE( summary.str().find( "layers {0: 1, 1: 1}" ) != std::string::npos );
}

/// Verify permissive and strict handling of GSL quadrature statuses.
/** Exercises the status policy used by `fourierTemporalPSD::singleLayerPSD`. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD applies quadrature policy", "[ao::analysis::fourierTemporalPSD]" )
{
    using mx::AO::analysis::fourierTemporalPSDPolicy;
    using mx::AO::analysis::fourierTemporalPSD_detail::applyPolicy;

    REQUIRE( applyPolicy( GSL_SUCCESS, fourierTemporalPSDPolicy::permissive ) == mx::error_t::noerror );
    for( const int status : { GSL_EMAXITER, GSL_EROUND, GSL_ESING, GSL_EDIVERGE } )
    {
        REQUIRE( applyPolicy( status, fourierTemporalPSDPolicy::permissive ) == mx::error_t::noerror );
        REQUIRE( applyPolicy( status, fourierTemporalPSDPolicy::strict ) == mx::error_t::liberr );
    }
    REQUIRE( applyPolicy( GSL_EDOM, fourierTemporalPSDPolicy::permissive ) == mx::error_t::invalidconfig );
    REQUIRE( applyPolicy( GSL_EINVAL, fourierTemporalPSDPolicy::permissive ) == mx::error_t::invalidconfig );
    REQUIRE( applyPolicy( GSL_ENOMEM, fourierTemporalPSDPolicy::permissive ) == mx::error_t::allocerr );
    REQUIRE( applyPolicy( GSL_EFAILED, fourierTemporalPSDPolicy::permissive ) == mx::error_t::liberr );
}

/// Verify multilayer error propagation and restoration of the caller's GSL handler.
/** Exercises `fourierTemporalPSD::multiLayerPSD` with an invalid basis. */
/**
 * \ingroup fourierTemporalPSD_unit_tests
 */
TEST_CASE( "Fourier temporal PSD propagates multilayer errors", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    temporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = -1;

    std::vector<double> frequency{ 1.0 };
    std::vector<double> psd( frequency.size(), 0.0 );
    temporalPsdT::reportT report;

    gsl_error_handler_t *previousHandler = gsl_set_error_handler( &testGslErrorHandler );
    const mx::error_t result =
        temporalPsd.multiLayerPSD<false>( psd,
                                          frequency,
                                          1.0,
                                          0.0,
                                          1,
                                          frequency.back(),
                                          &report,
                                          mx::AO::analysis::fourierTemporalPSDPolicy::permissive );
    gsl_error_handler_t *observedHandler = gsl_set_error_handler( previousHandler );

    REQUIRE( result == mx::error_t::invalidarg );
    REQUIRE( report.integrationsAttempted == 0 );
    REQUIRE( observedHandler == &testGslErrorHandler );
}
