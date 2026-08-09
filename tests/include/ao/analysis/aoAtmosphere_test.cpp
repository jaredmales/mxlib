/** \file aoAtmosphere_test.cpp
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/aoAtmosphere.hpp"

#include <limits>

typedef double realT;

using namespace mx::app;
using namespace mx::AO::analysis;

/// Verify parsing and validation of atmosphere configuration settings.
/** Exercises mx::AO::analysis::aoAtmosphere::setupConfig and
 * mx::AO::analysis::aoAtmosphere::loadConfig. */
SCENARIO( "Loading aoAtmosphere config settings", "[ao::analysis::aoAtmosphere]" )
{
    GIVEN( "a valid config file" )
    {
        aoAtmosphere<realT> atm;
        atm.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );

        WHEN( "all normal settings" )
        {
            appConfigurator config;

            writeConfigFile(
                "aoAtmosphere.conf",
                { "atm", "atm", "atm", "atm", "atm", "atm", "atm", "atm", "atm", "atm" },
                { "r_0", "lam_0", "L_0", "l_0", "layer_z", "h_obs", "H", "layer_Cn2", "layer_v_wind", "layer_dir" },
                { "0.25",
                  "0.4e-9",
                  "10,15,22.5",
                  "0.1,0.01,0.001",
                  "1000,2000,5000",
                  "3002",
                  "2.3",
                  "1,1,3.0",
                  "10.5,11.2,23.7",
                  "0.1,0.6,1.6" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.r_0() == 0.25 );
            REQUIRE( atm.lam_0() == 0.4e-9 );

            REQUIRE( atm.L_0( 0 ) == 10 );
            REQUIRE( atm.L_0( 1 ) == 15 );
            REQUIRE( atm.L_0( 2 ) == 22.5 );

            REQUIRE( atm.l_0( 0 ) == 0.1 );
            REQUIRE( atm.l_0( 1 ) == 0.01 );
            REQUIRE( atm.l_0( 2 ) == 0.001 );

            REQUIRE( atm.layer_z( 0 ) == 1000 );
            REQUIRE( atm.layer_z( 1 ) == 2000 );
            REQUIRE( atm.layer_z( 2 ) == 5000 );

            REQUIRE( atm.h_obs() == 3002 );

            REQUIRE( atm.H() == 2.3 );

            REQUIRE( atm.layer_Cn2( 0 ) == Approx( 0.2 ) );
            REQUIRE( atm.layer_Cn2( 1 ) == Approx( 0.2 ) );
            REQUIRE( atm.layer_Cn2( 2 ) == Approx( 0.6 ) );

            REQUIRE( atm.layer_v_wind( 0 ) == Approx( 10.5 ) );
            REQUIRE( atm.layer_v_wind( 1 ) == Approx( 11.2 ) );
            REQUIRE( atm.layer_v_wind( 2 ) == Approx( 23.7 ) );

            REQUIRE( atm.layer_dir( 0 ) == Approx( 0.1 ) );
            REQUIRE( atm.layer_dir( 1 ) == Approx( 0.6 ) );
            REQUIRE( atm.layer_dir( 2 ) == Approx( 1.6 ) );

            REQUIRE( atm.v_wind() == Approx( 19.278644205193533 ) );
            REQUIRE( atm.z_mean() == Approx( 3886.4496373530733 ) );

            REQUIRE( atm.nonKolmogorov() == false );
        }
        WHEN( "setting v_wind and z_mean" )
        {
            atm.L_0( { 25.0, 25.0, 25.0 } );
            atm.l_0( { 0.0, 0.0, 0.0 } );
            appConfigurator config;

            writeConfigFile( "aoAtmosphere.conf",
                             { "atm", "atm", "atm", "atm", "atm", "atm" },
                             { "layer_z", "layer_Cn2", "layer_v_wind", "layer_dir", "v_wind", "z_mean" },
                             { "1000,2000,5000", "1,1,3.0", "10.5,11.2,23.7", "0.1,0.6,1.6", "15.0", "5001.0" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.v_wind() == Approx( 15.0 ) );
            REQUIRE( atm.z_mean() == Approx( 5001.0 ) );

            REQUIRE( atm.nonKolmogorov() == false );
        }
        WHEN( "setting nonKolmogorov to false" )
        {
            appConfigurator config;

            writeConfigFile( "aoAtmosphere.conf", { "atm" }, { "nonKolmogorov" }, { "false" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.nonKolmogorov() == false );
        }
        WHEN( "setting nonKolmogorov to false" )
        {
            appConfigurator config;

            writeConfigFile( "aoAtmosphere.conf", { "atm" }, { "nonKolmogorov" }, { "true" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.nonKolmogorov() == true );
        }

        WHEN( "setting nonKolmogorov alpha" )
        {
            atm.nonKolmogorov( false );
            REQUIRE( atm.nonKolmogorov() == false );

            appConfigurator config;

            writeConfigFile( "aoAtmosphere.conf", { "atm" }, { "alpha" }, { "4.7" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.alpha( 0 ) == Approx( 4.7 ) );
            REQUIRE( atm.nonKolmogorov() == true );
        }

        WHEN( "setting nonKolmogorov beta" )
        {
            atm.nonKolmogorov( false );
            REQUIRE( atm.nonKolmogorov() == false );

            appConfigurator config;

            writeConfigFile( "aoAtmosphere.conf", { "atm" }, { "beta" }, { "0.026" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.beta( 0 ) == Approx( 0.026 ) );
            REQUIRE( atm.nonKolmogorov() == true );
        }

        WHEN( "setting nonKolmogorov beta_0" )
        {
            atm.nonKolmogorov( false );
            REQUIRE( atm.nonKolmogorov() == false );

            appConfigurator config;

            writeConfigFile( "aoAtmosphere.conf", { "atm" }, { "beta_0" }, { "1e-7" } );

            atm.setupConfig( config );
            config.readConfig( "aoAtmosphere.conf" );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );

            REQUIRE( atm.beta_0( 0 ) == Approx( 1e-7 ) );
            REQUIRE( atm.nonKolmogorov() == true );
        }
    }
    GIVEN( "command line options" )
    {
        aoAtmosphere<realT> atm;
        atm.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );

        WHEN( "setting nonKolmogorov to true" )
        {
            appConfigurator config;
            std::vector<std::string> argvs( 2 );

            argvs[0] = "test";
            argvs[1] = "--atm.nonKolmogorov=true";

            char *argv[2];
            argv[0] = (char *)argvs[0].data();
            argv[1] = (char *)argvs[1].data();

            REQUIRE( atm.nonKolmogorov() == false );
            atm.setupConfig( config );
            config.parseCommandLine( argvs.size(), argv );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );
            REQUIRE( atm.nonKolmogorov() == true );
        }

        WHEN( "setting nonKolmogorov to true" )
        {
            appConfigurator config;
            std::vector<std::string> argvs( 2 );

            argvs[0] = "test";
            argvs[1] = "--atm.nonKolmogorov=false";

            char *argv[2];
            argv[0] = (char *)argvs[0].data();
            argv[1] = (char *)argvs[1].data();

            atm.nonKolmogorov( true );
            REQUIRE( atm.nonKolmogorov() == true );
            atm.setupConfig( config );
            config.parseCommandLine( argvs.size(), argv );
            REQUIRE( atm.loadConfig( config ) == mx::error_t::noerror );
            REQUIRE( atm.nonKolmogorov() == false );
        }
    }
}

/// Verify that layer-strength normalization rejects invalid inputs without changing state.
/** Exercises mx::AO::analysis::aoAtmosphere::layer_Cn2. */
TEST_CASE( "Atmosphere layer-strength normalization is transactional", "[ao::analysis::aoAtmosphere]" )
{
    aoAtmosphere<realT> atmosphere;
    atmosphere.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );
    const std::vector<realT> originalStrength = atmosphere.layer_Cn2();

    SECTION( "empty strengths" )
    {
        REQUIRE( atmosphere.layer_Cn2( std::vector<realT>{} ) == mx::error_t::invalidarg );
    }

    SECTION( "negative strength" )
    {
        REQUIRE( atmosphere.layer_Cn2( { 1.0, -1.0 } ) == mx::error_t::invalidarg );
    }

    SECTION( "nonfinite strength" )
    {
        REQUIRE( atmosphere.layer_Cn2( { 1.0, std::numeric_limits<realT>::quiet_NaN() } ) == mx::error_t::invalidarg );
    }

    SECTION( "nonfinite sum" )
    {
        REQUIRE( atmosphere.layer_Cn2( { std::numeric_limits<realT>::max(), std::numeric_limits<realT>::max() } ) ==
                 mx::error_t::invalidarg );
    }

    SECTION( "zero sum" )
    {
        REQUIRE( atmosphere.layer_Cn2( { 0.0, 0.0 } ) == mx::error_t::invalidarg );
    }

    SECTION( "invalid reference wavelength" )
    {
        REQUIRE( atmosphere.layer_Cn2( { 1.0 }, -1.0 ) == mx::error_t::invalidarg );
    }

    REQUIRE( atmosphere.layer_Cn2() == originalStrength );
    REQUIRE( atmosphere.validate() == mx::error_t::noerror );
}

/// Verify that complete atmosphere validation rejects every layer-vector mismatch.
/** Exercises mx::AO::analysis::aoAtmosphere::validate with each core and non-Kolmogorov vector setter. */
TEST_CASE( "Atmosphere validation rejects mismatched layer vectors", "[ao::analysis::aoAtmosphere]" )
{
    aoAtmosphere<realT> atmosphere;
    atmosphere.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );

    SECTION( "layer strengths" )
    {
        REQUIRE( atmosphere.layer_Cn2( { 1.0, 1.0 } ) == mx::error_t::noerror );
    }

    SECTION( "outer scales" )
    {
        atmosphere.L_0( { 25.0, 25.0 } );
    }

    SECTION( "inner scales" )
    {
        atmosphere.l_0( { 0.0, 0.0 } );
    }

    SECTION( "layer heights" )
    {
        atmosphere.layer_z( { 0.0, 0.0 } );
    }

    SECTION( "wind speeds" )
    {
        atmosphere.layer_v_wind( { 10.0, 10.0 } );
    }

    SECTION( "wind directions" )
    {
        atmosphere.layer_dir( { 0.0, 0.0 } );
    }

    SECTION( "non-Kolmogorov exponents" )
    {
        atmosphere.alpha( { 11.0 / 3.0, 11.0 / 3.0 } );
    }

    SECTION( "non-Kolmogorov normalizations" )
    {
        atmosphere.beta( { 1.0, 1.0 } );
    }

    SECTION( "non-Kolmogorov constants" )
    {
        atmosphere.beta_0( { 0.0, 0.0 } );
    }

    REQUIRE( atmosphere.validate() == mx::error_t::sizeerr );
}

/// Verify that atmosphere validation rejects nonphysical scalar and layer values.
/** Exercises mx::AO::analysis::aoAtmosphere::validate after mutation through its public setters. */
TEST_CASE( "Atmosphere validation rejects invalid physical values", "[ao::analysis::aoAtmosphere]" )
{
    aoAtmosphere<realT> atmosphere;
    atmosphere.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );

    SECTION( "invalid Fried parameter" )
    {
        atmosphere.r_0( 0.0, 0.5e-6 );
    }

    SECTION( "invalid observatory height" )
    {
        atmosphere.h_obs( -1.0 );
    }

    SECTION( "invalid atmospheric scale height" )
    {
        atmosphere.H( 0.0 );
    }

    SECTION( "nonfinite outer scale" )
    {
        atmosphere.L_0( std::vector<realT>{ std::numeric_limits<realT>::infinity() } );
    }

    SECTION( "negative inner scale" )
    {
        atmosphere.l_0( std::vector<realT>{ -1.0 } );
    }

    SECTION( "negative layer height" )
    {
        atmosphere.layer_z( std::vector<realT>{ -1.0 } );
    }

    SECTION( "negative wind speed" )
    {
        atmosphere.layer_v_wind( std::vector<realT>{ -1.0 } );
    }

    SECTION( "nonfinite wind direction" )
    {
        atmosphere.layer_dir( std::vector<realT>{ std::numeric_limits<realT>::quiet_NaN() } );
    }

    SECTION( "invalid non-Kolmogorov normalization" )
    {
        atmosphere.nonKolmogorov( true );
        atmosphere.beta( std::vector<realT>{ 0.0 } );
    }

    REQUIRE( atmosphere.validate() == mx::error_t::invalidconfig );
}

/// Verify that atmosphere validation preserves model-specific physical conventions.
/** Exercises mx::AO::analysis::aoAtmosphere::validate for static layers, infinite outer scale, and explicit
 * non-Kolmogorov normalization. */
TEST_CASE( "Atmosphere validation accepts supported physical conventions", "[ao::analysis::aoAtmosphere]" )
{
    aoAtmosphere<realT> atmosphere;
    atmosphere.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );

    SECTION( "static layer" )
    {
        atmosphere.layer_v_wind( std::vector<realT>{ 0.0 } );
    }

    SECTION( "infinite outer scale" )
    {
        atmosphere.L_0( std::vector<realT>{ 0.0 } );
    }

    SECTION( "non-Kolmogorov model without a Fried parameter" )
    {
        atmosphere.r_0( 0.0, 0.5e-6 );
        atmosphere.nonKolmogorov( true );
    }

    REQUIRE( atmosphere.validate() == mx::error_t::noerror );
}

/// Verify that configuration loading returns layer-strength validation failures.
/** Exercises mx::AO::analysis::aoAtmosphere::loadConfig and
 * mx::AO::analysis::aoAtmosphere::layer_Cn2. */
TEST_CASE( "Atmosphere configuration propagates invalid layer strengths", "[ao::analysis::aoAtmosphere]" )
{
    aoAtmosphere<realT> atmosphere;
    atmosphere.setSingleLayer( 0.2, 0.5e-6, 25.0, 0.0, 0.0, 10.0, 0.0 );

    appConfigurator config;
    writeConfigFile( "aoAtmosphere.conf", { "atm" }, { "layer_Cn2" }, { "0,0" } );
    atmosphere.setupConfig( config );
    config.readConfig( "aoAtmosphere.conf" );

    REQUIRE( atmosphere.loadConfig( config ) == mx::error_t::invalidarg );
    REQUIRE( atmosphere.validate() == mx::error_t::noerror );
    REQUIRE( atmosphere.layer_Cn2() == std::vector<realT>{ 1.0 } );
}
