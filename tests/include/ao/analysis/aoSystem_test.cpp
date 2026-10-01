/** \file aoSystem_test.cpp
 * \brief Tests adaptive-optics system configuration.
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/aoSystem.hpp"

#include <sstream>

typedef double realT;

using namespace mx::app;
using namespace mx::AO::analysis;

/// Verify parsing of AO-system configuration settings.
/** Exercises mx::AO::analysis::aoSystem::setupConfig and mx::AO::analysis::aoSystem::loadConfig. */
/**
 * \ingroup aoSystem_unit_tests
 */
TEST_CASE( "Loading aoSystem config settings", "[ao::analysis::aoSystem]" )
{
    GIVEN( "no config file" )
    {
        aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aosys; // This will be cumulative
        WHEN( "verifying the defaults" )
        {

            REQUIRE( aosys.D() == 0.0 );
            REQUIRE( aosys.d_min( 0 ) == 0.0 );
            REQUIRE( aosys.d_opt() == 1e-50 );
            REQUIRE( aosys.optd() == false );
            REQUIRE( aosys.optd_delta() == 1.0 );
            REQUIRE( aosys.wfsBeta()->_id == "Ideal WFS" );
            REQUIRE( aosys.lam_wfs() == 0 );
            REQUIRE( aosys.npix_wfs().size() == 1 );
            REQUIRE( aosys.npix_wfs( (size_t)0 ) == 0 );
            REQUIRE( aosys.ron_wfs().size() == 1 );
            REQUIRE( aosys.ron_wfs( (size_t)0 ) == 0 );
            REQUIRE( aosys.Fbg().size() == 1 );
            REQUIRE( aosys.Fbg( (size_t)0 ) == 0 );
            REQUIRE( aosys.minTauWFS().size() == 1 );
            REQUIRE( aosys.minTauWFS( (size_t)0 ) == 0 );
            REQUIRE( aosys.bin_npix() == true );
            REQUIRE( aosys.tauWFS() == 0 );
            REQUIRE( aosys.deltaTau() == 0 );
            REQUIRE( aosys.optTau() == true );
            REQUIRE( aosys.lam_sci() == 0 );
            REQUIRE( aosys.zeta() == 0 );
            REQUIRE( aosys.secZeta() == 1.0 );
            REQUIRE( aosys.fit_mn_max() == 100 );
            REQUIRE( aosys.circularLimit() == false );
            REQUIRE( aosys.spatialFilter_ku() == std::numeric_limits<realT>::max() );
            REQUIRE( aosys.spatialFilter_kv() == std::numeric_limits<realT>::max() );
            REQUIRE( aosys.ncp_wfe() == 0 );
            REQUIRE( aosys.ncp_alpha() == 2.0 );
            REQUIRE( aosys.vibration_rms() == 0.0 );
            REQUIRE( aosys.vibration_units() == "mas" );
            REQUIRE( aosys.central_obs_k() == 1.03 );
            REQUIRE( aosys.vibrationStrehl() == 1.0 );
            REQUIRE( aosys.F0() == 0.0 );
            REQUIRE( aosys.starMag() == 0 );
            REQUIRE( aosys.Fg() == 0 );
            REQUIRE( aosys.Fg( 2.5 ) == 0 );
        }
    }

    GIVEN( "a valid config file" )
    {
        aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aosys; // This will be cumulative

        WHEN( "setting wfs to ideal" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf", { "aosys" }, { "wfs" }, { "ideal" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.wfsBeta()->_id == "Ideal WFS" );
        }
        WHEN( "setting wfs to unmodPyWFS" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf", { "aosys" }, { "wfs" }, { "unmodPyWFS" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );

            aosys.loadConfig( config );

            REQUIRE( aosys.wfsBeta()->_id == "Unmodulated Pyramid" );
        }
        WHEN( "setting wfs to asympModPyWFS" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf", { "aosys" }, { "wfs" }, { "asympModPyWFS" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.wfsBeta()->_id == "Asymptotic Modulated Pyramid" );
        }
        WHEN( "setting wfs to SHWFS" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf", { "aosys" }, { "wfs" }, { "SHWFS" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.wfsBeta()->_id == "Shack Hartmann" );
        }
        WHEN( "typical settings part 1" ) // These are broken into parts just to keep things <  1 line
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf",
                             { "aosys", "aosys", "aosys", "aosys", "aosys" },
                             { "D", "d_min", "optd", "optd_delta", "lam_wfs" },
                             { "7.6", "0.122", "true", "0.0037", "800e-9" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.D() == Approx( 7.6 ) );
            REQUIRE( aosys.d_min( 0 ) == Approx( 0.122 ) );
            REQUIRE( aosys.optd() == true );
            REQUIRE( aosys.optd_delta() == Approx( 0.0037 ) );
            REQUIRE( aosys.lam_wfs() == Approx( 800e-9 ) );
        }
        WHEN( "typical settings part 2" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf",
                             { "aosys", "aosys", "aosys", "aosys", "aosys" },
                             { "npix_wfs", "ron_wfs", "Fbg", "minTauWFS", "bin_npix" },
                             { "100,1001", "0.5,0.25", "0.01,0.03", "0.0001,0.0005", "false" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.npix_wfs( 0 ) == 100 );
            REQUIRE( aosys.npix_wfs( 1 ) == 1001 );
            REQUIRE( aosys.ron_wfs( 0 ) == 0.5 );
            REQUIRE( aosys.ron_wfs( 1 ) == 0.25 );
            REQUIRE( aosys.Fbg( 0 ) == 0.01 );
            REQUIRE( aosys.Fbg( 1 ) == 0.03 );
            REQUIRE( aosys.minTauWFS( 0 ) == 0.0001 );
            REQUIRE( aosys.minTauWFS( 1 ) == 0.0005 );
            REQUIRE( aosys.bin_npix() == false );
        }
        WHEN( "typical settings part 3" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf",
                             { "aosys", "aosys", "aosys", "aosys", "aosys", "aosys", "aosys", "aosys", "aosys" },
                             { "tauWFS",
                               "deltaTau",
                               "optTau",
                               "lam_sci",
                               "zeta",
                               "fit_mn_max",
                               "circularLimit",
                               "spatialFilter_ku",
                               "spatialFilter_kv" },
                             { "0.0007", "0.0001", "false", "802e-9", "0.15", "12", "true", "50", "53" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.tauWFS() == Approx( 0.0007 ) );
            REQUIRE( aosys.deltaTau() == Approx( 0.0001 ) );
            REQUIRE( aosys.optTau() == false );
            REQUIRE( aosys.lam_sci() == Approx( 802e-9 ) );
            REQUIRE( aosys.zeta() == Approx( 0.15 ) );
            REQUIRE( aosys.secZeta() == Approx( 1 / cos( 0.15 ) ) );
            REQUIRE( aosys.fit_mn_max() == 12 );
            REQUIRE( aosys.circularLimit() == true );
            REQUIRE( aosys.spatialFilter_ku() == 50 );
            REQUIRE( aosys.spatialFilter_kv() == 53 );
        }
        WHEN( "typical settings part 4" )
        {
            appConfigurator config;

            writeConfigFile(
                "aoSystem.conf",
                { "aosys", "aosys", "aosys", "aosys", "atm" }, // last one verifies that atm config gets processed
                { "ncp_wfe", "ncp_alpha", "starMag", "F0", "r_0" },
                { "30e-9", "2.34", "8.2", "5e9", "0.25" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.ncp_wfe() == Approx( 30e-9 ) );
            REQUIRE( aosys.ncp_alpha() == Approx( 2.34 ) );
            REQUIRE( aosys.F0() == Approx( 5e9 ) );
            REQUIRE( aosys.starMag() == Approx( 8.2 ) );
            REQUIRE( aosys.Fg() == Approx( 5e9 * pow( 10, -0.4 * 8.2 ) ) );
            REQUIRE( aosys.Fg( 9.0 ) == Approx( 5e9 * pow( 10, -0.4 * 9.0 ) ) );
            REQUIRE( aosys.atm.r_0() == Approx( 0.25 ) );
        }
        WHEN( "setting image quality in milliarcseconds" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf",
                             { "aosys", "aosys", "aosys" },
                             { "vibration_rms", "vibration_units", "central_obs_k" },
                             { "5.5", "mas", "1.12" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.vibration_rms() == Approx( 5.5 ) );
            REQUIRE( aosys.vibration_units() == "mas" );
            REQUIRE( aosys.central_obs_k() == Approx( 1.12 ) );
        }
        WHEN( "setting image quality in radians" )
        {
            appConfigurator config;

            writeConfigFile( "aoSystem.conf",
                             { "aosys", "aosys", "aosys" },
                             { "vibration_rms", "vibration_units", "central_obs_k" },
                             { "2.5e-8", "rad", "1.08" } );

            aosys.setupConfig( config );
            config.readConfig( "aoSystem.conf" );
            aosys.loadConfig( config );

            REQUIRE( aosys.vibration_rms() == Approx( 2.5e-8 ) );
            REQUIRE( aosys.vibration_units() == "rad" );
            REQUIRE( aosys.central_obs_k() == Approx( 1.08 ) );
        }
    }
}

/// Verify that AO-system configuration loading propagates atmosphere validity.
/** Exercises mx::AO::analysis::aoSystem::loadConfig and mx::AO::analysis::aoAtmosphere::setSingleLayer. */
/**
 * \ingroup aoSystem_unit_tests
 */
TEST_CASE( "AO-system configuration reports atmosphere validity", "[ao::analysis::aoSystem]" )
{
    aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aoSystem;
    appConfigurator config;
    aoSystem.setupConfig( config );

    SECTION( "default atmosphere" )
    {
        REQUIRE( aoSystem.loadConfig( config ) == mx::error_t::sizeerr );
    }

    SECTION( "complete atmosphere" )
    {
        aoSystem.atm.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );
        REQUIRE( aoSystem.loadConfig( config ) == mx::error_t::noerror );
    }
}

/** Verify AO-system vibration conversion, image-quality formulas, validation, and setup output.
 * Exercises mx::AO::analysis::aoSystem::vibration_rms,
 * mx::AO::analysis::aoSystem::vibrationStrehl, mx::AO::analysis::aoSystem::central_obs_k,
 * mx::AO::analysis::aoSystem::fwhm, and mx::AO::analysis::aoSystem::dumpAOSystem.
 *
 * \ingroup aoSystem_unit_tests
 */
TEST_CASE( "AO-system vibration and FWHM calculations", "[ao::analysis::aoSystem]" )
{
    aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aosys;
    constexpr realT diameter = 6.5;
    constexpr realT wavelength = 656e-9;
    constexpr realT vibrationMas = 5.0;
    const realT vibrationRad = vibrationMas * mx::math::pi<realT>() / 648000000.0;

    aosys.D( diameter );
    aosys.lam_sci( wavelength );

    SECTION( "central obscuration factor is present without vibration" )
    {
        REQUIRE( aosys.vibrationStrehl() == 1.0 );
        REQUIRE( aosys.fwhm() == Approx( aosys.central_obs_k() * wavelength / diameter ) );
    }

    SECTION( "milliarcseconds and radians are equivalent" )
    {
        const realT lambdaD = wavelength / diameter;
        const realT expectedStrehl = 1.0 / ( 1.0 + 8.0 * mx::math::ln_two<realT>() / pow( aosys.central_obs_k(), 2 ) *
                                                       pow( vibrationRad / lambdaD, 2 ) );
        const realT expectedFwhm = sqrt( pow( aosys.central_obs_k() * lambdaD, 2 ) +
                                         8.0 * mx::math::ln_two<realT>() * pow( vibrationRad, 2 ) );

        aosys.vibration_rms( vibrationMas, "mas" );
        REQUIRE( aosys.vibration_rms() == Approx( vibrationMas ) );
        REQUIRE( aosys.vibration_units() == "mas" );
        REQUIRE( aosys.vibrationStrehl() == Approx( expectedStrehl ) );
        REQUIRE( aosys.fwhm() == Approx( expectedFwhm ) );

        aosys.vibration_rms( vibrationRad, "rad" );
        REQUIRE( aosys.vibration_rms() == Approx( vibrationRad ) );
        REQUIRE( aosys.vibration_units() == "rad" );
        REQUIRE( aosys.vibrationStrehl() == Approx( expectedStrehl ) );
        REQUIRE( aosys.fwhm() == Approx( expectedFwhm ) );
    }

    SECTION( "central obscuration factor is configurable" )
    {
        constexpr realT centralObsK = 1.2;
        aosys.central_obs_k( centralObsK );
        aosys.vibration_rms( vibrationMas, "mas" );

        const realT lambdaD = wavelength / diameter;
        const realT expectedStrehl =
            1.0 / ( 1.0 + 8.0 * mx::math::ln_two<realT>() / pow( centralObsK, 2 ) * pow( vibrationRad / lambdaD, 2 ) );
        const realT expectedFwhm =
            sqrt( pow( centralObsK * lambdaD, 2 ) + 8.0 * mx::math::ln_two<realT>() * pow( vibrationRad, 2 ) );

        REQUIRE( aosys.central_obs_k() == Approx( centralObsK ) );
        REQUIRE( aosys.vibrationStrehl() == Approx( expectedStrehl ) );
        REQUIRE( aosys.fwhm() == Approx( expectedFwhm ) );
    }

    SECTION( "invalid image-quality inputs are rejected" )
    {
        REQUIRE_THROWS_AS( aosys.vibration_rms( -1, "mas" ), mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.vibration_rms( std::numeric_limits<realT>::infinity(), "mas" ),
                           mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.vibration_rms( std::numeric_limits<realT>::quiet_NaN(), "rad" ),
                           mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.vibration_rms( 1, "arcsec" ), mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.central_obs_k( 0 ), mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.central_obs_k( -1 ), mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.central_obs_k( std::numeric_limits<realT>::infinity() ),
                           mx::exception<mx::verbose::d> );
        REQUIRE_THROWS_AS( aosys.central_obs_k( std::numeric_limits<realT>::quiet_NaN() ),
                           mx::exception<mx::verbose::d> );
    }

    SECTION( "setup output is round-trippable" )
    {
        aosys.atm.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );
        aosys.vibration_rms( vibrationMas, "mas" );
        aosys.central_obs_k( 1.12 );
        std::ostringstream output;
        aosys.dumpAOSystem( output );

        REQUIRE( output.str().find( "#    vibration_rms = 5" ) != std::string::npos );
        REQUIRE( output.str().find( "#    vibration_units = mas" ) != std::string::npos );
        REQUIRE( output.str().find( "#    central_obs_k = 1.12" ) != std::string::npos );
    }
}

/** Verify vibration multiplies both AO-system Strehl calculations without changing WFE.
 * Exercises mx::AO::analysis::aoSystem::strehl, mx::AO::analysis::aoSystem::vibrationStrehl,
 * and mx::AO::analysis::aoSystem::wfeVar.
 *
 * \ingroup aoSystem_unit_tests
 */
TEST_CASE( "AO-system vibration multiplies final Strehl", "[ao::analysis::aoSystem]" )
{
    aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aosys;
    aosys.loadMagAOX();
    aosys.fit_mn_max( 2 );
    aosys.optd( false );
    aosys.starMag( 0 );

    const realT actuatorPitch = aosys.d_min( 0 );
    const realT directStrehl = aosys.strehl( actuatorPitch, 0 );
    const realT cachedStrehl = aosys.strehl();
    const realT wfeVariance = aosys.wfeVar();
    const realT measurementVariance = aosys.measurementErrorTotal();
    const realT timeDelayVariance = aosys.timeDelayErrorTotal();
    const realT fittingVariance = aosys.fittingErrorTotal();
    const realT chromScintOPDVariance = aosys.chromScintOPDErrorTotal();
    const realT chromIndexVariance = aosys.chromIndexErrorTotal();
    const realT dispAnisoOPDVariance = aosys.dispAnisoOPDErrorTotal();
    const realT ncpVariance = aosys.ncpError();

    aosys.vibration_rms( 5, "mas" );
    const realT vibrationStrehl = aosys.vibrationStrehl();

    REQUIRE( aosys.strehl( actuatorPitch, 0 ) == Approx( directStrehl * vibrationStrehl ) );
    REQUIRE( aosys.strehl() == Approx( cachedStrehl * vibrationStrehl ) );
    REQUIRE( aosys.wfeVar() == Approx( wfeVariance ) );
    REQUIRE( aosys.measurementErrorTotal() == Approx( measurementVariance ) );
    REQUIRE( aosys.timeDelayErrorTotal() == Approx( timeDelayVariance ) );
    REQUIRE( aosys.fittingErrorTotal() == Approx( fittingVariance ) );
    REQUIRE( aosys.chromScintOPDErrorTotal() == Approx( chromScintOPDVariance ) );
    REQUIRE( aosys.chromIndexErrorTotal() == Approx( chromIndexVariance ) );
    REQUIRE( aosys.dispAnisoOPDErrorTotal() == Approx( dispAnisoOPDVariance ) );
    REQUIRE( aosys.ncpError() == Approx( ncpVariance ) );
}
