/** \file fourierTemporalPSD_test.cpp
 * \brief Tests of Fourier-mode temporal power spectral densities.
 */

#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/fourierTemporalPSD.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using aoSystemBaseT = mx::AO::analysis::aoSystem<double, mx::AO::analysis::vonKarmanSpectrum<double>>;

/// AO-system test type that forces local template instantiation under sanitizers.
struct aoSystemT : public aoSystemBaseT
{
};

using temporalPsdT = mx::AO::analysis::fourierTemporalPSD<double, aoSystemT>;

namespace
{

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

} // namespace

/// Verify temporal-PSD tail initialization at and below its nominal averaging width.
/** Exercises zero, one, 49, and 50 exactly integrated frequency bins. */
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
