/** \file appConfigurator_test.cpp
 *  \brief Tests application configuration parsing.
 */

#include "../../catch2/catch.hpp"

#include <fstream>

#include "../../../include/app/appConfigurator.hpp"

using namespace mx::app;

/// \cond appConfiguratorTestHelpers
namespace
{

std::string logName;
int logCode;
std::string logValue;
std::string logSource;

void captureConfigLog( const std::string &name, const int &code, const std::string &value, const std::string &source )
{
    logName = name;
    logCode = code;
    logValue = value;
    logSource = source;
}

void resetConfigLog()
{
    logName.clear();
    logCode = 0;
    logValue.clear();
    logSource.clear();
}

/// Return a deterministic parser allocation failure without changing parser state.
int allocationFailureParser( iniFile &parser /**< [in,out] parser ignored by the injected failure. */,
                             const std::string &filename /**< [in] filename ignored by the injected failure. */ )
{
    static_cast<void>( parser );
    static_cast<void>( filename );
    return -2;
}

/// Restore the production parser after a readConfig failure-injection test.
struct parserReset
{
    /// Reset the parser hook before the test body runs.
    parserReset()
    {
        mx::app::detail::appConfiguratorParser = nullptr;
    }

    /// Restore the production parser when the test scope exits.
    ~parserReset()
    {
        mx::app::detail::appConfiguratorParser = nullptr;
    }
};

} // namespace
/// \endcond

/// Verify that appConfigurator parses global, sectioned, and unused configuration entries.
/**
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "config file parsing", "[appConfigurator]" )
{
    GIVEN( "a basic config file" )
    {
        WHEN( "no sections" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "", "", "", "" },
                             { "key0", "key1", "key2", "key3", "key4", "key5" },
                             { "val0", "val1", "val2", "val3", "val4", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "", "key1", false, "", "" );
            config.add( "key2", "", "", 0, "", "key2", false, "", "" );
            config.add( "key3", "", "", 0, "", "key3", false, "", "" );
            config.add( "key4", "", "", 0, "", "key4", false, "", "" );
            config.add( "key5", "", "", 0, "", "key5", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0" );

            config( val, "key1" );
            REQUIRE( val == "val1" );

            config( val, "key2" );
            REQUIRE( val == "val2" );

            config( val, "key3" );
            REQUIRE( val == "val3" );

            config( val, "key4" );
            REQUIRE( val == "val4" );

            config( val, "key5" );
            REQUIRE( val == "val5" );
        }

        WHEN( "sections, unique keys" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "sect1", "sect1", "sect2", "sect2" },
                             { "key0", "key1", "key2", "key3", "key4", "key5" },
                             { "val0", "val1", "val2", "val3", "val4", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "", "key1", false, "", "" );
            config.add( "sect1.key2", "", "", 0, "sect1", "key2", false, "", "" );
            config.add( "sect1.key3", "", "", 0, "sect1", "key3", false, "", "" );
            config.add( "sect2.key4", "", "", 0, "sect2", "key4", false, "", "" );
            config.add( "sect2.key5", "", "", 0, "sect2", "key5", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0" );

            config( val, "key1" );
            REQUIRE( val == "val1" );

            config( val, "sect1.key2" );
            REQUIRE( val == "val2" );

            config( val, "sect1.key3" );
            REQUIRE( val == "val3" );

            config( val, "sect2.key4" );
            REQUIRE( val == "val4" );

            config( val, "sect2.key5" );
            REQUIRE( val == "val5" );
        }

        WHEN( "sections, repeated keys still unique within sections" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "sect1", "sect1", "sect2", "sect2" },
                             { "key0", "key1", "key2", "key3", "key2", "key3" },
                             { "val0", "val1", "val2", "val3", "val4", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "", "key1", false, "", "" );
            config.add( "sect1.key2", "", "", 0, "sect1", "key2", false, "", "" );
            config.add( "sect1.key3", "", "", 0, "sect1", "key3", false, "", "" );
            config.add( "sect2.key2", "", "", 0, "sect2", "key2", false, "", "" );
            config.add( "sect2.key3", "", "", 0, "sect2", "key3", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0" );

            config( val, "key1" );
            REQUIRE( val == "val1" );

            config( val, "sect1.key2" );
            REQUIRE( val == "val2" );

            config( val, "sect1.key3" );
            REQUIRE( val == "val3" );

            config( val, "sect2.key2" );
            REQUIRE( val == "val4" );

            config( val, "sect2.key3" );
            REQUIRE( val == "val5" );
        }
    }

    GIVEN( "a config file with unused entries" )
    {

        WHEN( "no sections, unused entry in middle" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "", "", "", "" },
                             { "key0", "key1", "key2", "key3", "key4", "key5" },
                             { "val0", "val1", "val2", "val3", "val4", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "", "key1", false, "", "" );
            config.add( "key2", "", "", 0, "", "key2", false, "", "" );
            config.add( "key4", "", "", 0, "", "key4", false, "", "" );
            config.add( "key5", "", "", 0, "", "key5", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            // Check normal parsing
            config( val, "key0" );
            REQUIRE( val == "val0" );

            config( val, "key1" );
            REQUIRE( val == "val1" );

            config( val, "key2" );
            REQUIRE( val == "val2" );

            config( val, "key4" );
            REQUIRE( val == "val4" );

            config( val, "key5" );
            REQUIRE( val == "val5" );

            // Check taht the unused on is unused
            REQUIRE( config.m_unusedConfigs[iniFile::makeKey( "", "key3" )].used == false );

            // Check that the unused one is available.
            config.configUnused( val, "", "key3" );
            REQUIRE( val == "val3" );

            // Check taht the unused on is now used
            REQUIRE( config.m_unusedConfigs[iniFile::makeKey( "", "key3" )].used == true );
        }

        WHEN( "sections, repeated keys, unused sections" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "sect1", "sect1", "sect2", "sect2", "sect3" },
                             { "key0", "key1", "key2", "key3", "key2", "key3", "key4" },
                             { "val0", "val1", "val2", "val3", "val4", "val5", "val6" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "", "key1", false, "", "" );
            config.add( "sect2.key2", "", "", 0, "sect2", "key2", false, "", "" );
            config.add( "sect2.key3", "", "", 0, "sect2", "key3", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0" );

            config( val, "key1" );
            REQUIRE( val == "val1" );

            config( val, "sect2.key2" );
            REQUIRE( val == "val4" );

            config( val, "sect2.key3" );
            REQUIRE( val == "val5" );

            config.configUnused( val, "sect1", "key2" );
            REQUIRE( val == "val2" );

            config.configUnused( val, "sect1", "key3" );
            REQUIRE( val == "val3" );

            config.configUnused( val, "sect3", "key4" );
            REQUIRE( val == "val6" );

            std::vector<std::string> sections;
            config.unusedSections( sections );
            REQUIRE( sections.size() == 2 );
            REQUIRE( ( sections[0] == "sect1" || sections[1] == "sect1" ) );
            REQUIRE( ( sections[0] == "sect3" || sections[1] == "sect3" ) );
        }
    }

#if 1
    GIVEN( "a config file with repeated keys within the same section" )
    {
        WHEN( "no sections" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "", "", "", "" },
                             { "key0", "key0", "key2", "key2", "key4", "key4" },
                             { "val0", "val1", "val2", "val3", "val4", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key2", "", "", 0, "", "key2", false, "", "" );
            config.add( "key4", "", "", 0, "", "key4", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0val1" );

            config( val, "key2" );
            REQUIRE( val == "val2val3" );

            config( val, "key4" );
            REQUIRE( val == "val4val5" );
        }

        WHEN( "repeated sections and keys" )
        {
            writeConfigFile( "/tmp/test7.conf",
                             { "", "sect1", "sect2", "sect1", "sect2", "sect2", "sect3" },
                             { "key0", "key1", "key2", "key1", "key2", "key2", "key3" },
                             { "val0", "val1", "val2", "val3", "val4", "val4.1", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "sect1", "key1", false, "", "" );
            config.add( "key2", "", "", 0, "sect2", "key2", false, "", "" );
            config.add( "key3", "", "", 0, "sect3", "key3", false, "", "" );

            config.readConfig( "/tmp/test7.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0" );

            config( val, "key1" );
            REQUIRE( val == "val1val3" );

            config( val, "key2" );
            REQUIRE( val == "val2val4val4.1" );

            config( val, "key3" );
            REQUIRE( val == "val5" );
        }

        WHEN( "multi-line keys" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "sect1", "sect2", "sect1", "sect2", "sect3" },
                             { "key0", "key1", "key2", "key1", "key2", "key3" },
                             { "val0\n    val0.1", "val1", "val2", "val3\n val3.1", "val4", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "sect1", "key1", false, "", "" );
            config.add( "key2", "", "", 0, "sect2", "key2", false, "", "" );
            config.add( "key3", "", "", 0, "sect3", "key3", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::string val;

            config( val, "key0" );
            REQUIRE( val == "val0val0.1" );

            config( val, "key1" );
            REQUIRE( val == "val1val3val3.1" );

            config( val, "key2" );
            REQUIRE( val == "val2val4" );

            config( val, "key3" );
            REQUIRE( val == "val5" );
        }
    }

    GIVEN( "a config file with vectors" )
    {
        WHEN( "no sections" )
        {
            writeConfigFile( "/tmp/test.conf",
                             { "", "", "" },
                             { "key0", "key1", "key2" },
                             { "val0,val1,val2", " val3, val4,    ", "val5" } );

            mx::app::appConfigurator config;
            config.add( "key0", "", "", 0, "", "key0", false, "", "" );
            config.add( "key1", "", "", 0, "", "key1", false, "", "" );
            config.add( "key2", "", "", 0, "", "key2", false, "", "" );

            config.readConfig( "/tmp/test.conf" );

            std::vector<std::string> vals;

            config( vals, "key0" );
            REQUIRE( vals[0] == "val0" );
            REQUIRE( vals[1] == "val1" );
            REQUIRE( vals[2] == "val2" );

            config( vals, "key1" );
            REQUIRE( vals[0] == "val3" );
            REQUIRE( vals[1] == "val4" );
            REQUIRE( vals[2] == "" );

            config( vals, "key2" );
            REQUIRE( vals[0] == "val5" );
        }
    }
#endif
}

/** \brief Verifies mx::app::appConfigurator::readConfig handles parser errors and records config sources.
 *
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "readConfig handles failures and tracks recognized sources", "[appConfigurator]" )
{
    mx::app::appConfigurator config;
    config.m_sources = true;
    config.add( "known", "", "", 0, "section", "known", false, "", "" );

    REQUIRE( config.readConfig( "" ) == 0 );
    REQUIRE( config.readConfig( "/tmp/mxlib-missing-config.conf", false ) == -1 );
    REQUIRE( config.readConfig( "/tmp/mxlib-missing-config.conf" ) == -1 );

    {
        std::ofstream badConfig( "/tmp/mxlib-invalid-config.conf" );
        badConfig << "[unterminated\n";
    }
    REQUIRE( config.readConfig( "/tmp/mxlib-invalid-config.conf" ) == -1 );

    {
        std::ofstream goodConfig( "/tmp/mxlib-read-config.conf" );
        goodConfig << "[section]\nknown = value\n[extra]\nunused = retained\n";
    }
    REQUIRE( config.readConfig( "/tmp/mxlib-read-config.conf" ) == 0 );
    REQUIRE( config.m_targets.at( "known" ).values.at( 0 ) == "value" );
    REQUIRE( config.m_targets.at( "known" ).sources.at( 0 ) == "/tmp/mxlib-read-config.conf" );
    REQUIRE( config.m_unusedConfigs.at( "extra=unused" ).values.at( 0 ) == "retained" );
    REQUIRE( config.m_unusedConfigs.at( "extra=unused" ).sources.at( 0 ) == "/tmp/mxlib-read-config.conf" );

    parserReset reset;
    mx::app::detail::appConfiguratorParser = allocationFailureParser;
    REQUIRE( config.readConfig( "/tmp/mxlib-allocation-failure.conf" ) == -1 );
}

/** \brief Verifies command-line target metadata and typed retrieval used by hciReduce configuration loading.
 *
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "command-line configuration targets preserve metadata and typed values", "[appConfigurator]" )
{
    mx::app::appConfigurator config;
    config.add( "images", "i", "images", mx::app::argType::Required, "files", "images", true, "string", "" );
    config.add( "threshold", "t", "threshold", mx::app::argType::Required, "quality", "threshold", false, "float", "" );
    config.add( "enabled", "", "enabled", mx::app::argType::True, "preProcess", "enabled", false, "bool", "" );

    char invokedName[] = "configTest";
    char imagesOption[] = "--images=target.fits";
    char thresholdOption[] = "--threshold=0.75";
    char enabledOption[] = "--enabled";
    char *argv[] = { invokedName, imagesOption, thresholdOption, enabledOption };
    config.parseCommandLine( 4, argv );

    REQUIRE( config.nAdded == 3 );
    REQUIRE( config.m_targets.at( "images" ).section == "files" );
    REQUIRE( config.m_targets.at( "images" ).keyword == "images" );
    REQUIRE( config.m_targets.at( "images" ).isRequired );
    REQUIRE( config.isSet( "images" ) );
    REQUIRE( config.isSet( "threshold" ) );
    REQUIRE( config.isSet( "enabled" ) );

    std::string images;
    float threshold{ 0 };
    bool enabled{ false };
    REQUIRE( config( images, "images" ) == 0 );
    REQUIRE( config( threshold, "threshold" ) == 0 );
    REQUIRE( config( enabled, "enabled" ) == 0 );
    REQUIRE( images == "target.fits" );
    REQUIRE( threshold == 0.75f );
    REQUIRE( enabled );

    int unset{ 19 };
    REQUIRE( config( unset, "not-added" ) == 0 );
    REQUIRE( unset == 19 );
}

/** \brief Verifies command-line filtering, indexed retrieval, unknown options, and reusable configuration state.
 *
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "command-line configuration tracks selected and unknown options", "[appConfigurator]" )
{
    mx::app::appConfigurator config;
    config.m_sources = true;
    config.add( "first", "", "first", mx::app::argType::Required, "", "first", false, "int", "" );
    config.add( "second", "", "second", mx::app::argType::Required, "", "second", false, "int", "" );
    char invokedName[] = "configTest";
    char firstOption[] = "--first=5";
    char secondOption[] = "--second=7";
    char unknownOption[] = "--unknown=value";
    char *argv[] = { invokedName, firstOption, secondOption, unknownOption };

    config.parseCommandLine( 4, argv, "first" );
    REQUIRE( config.isSet( "first" ) );
    REQUIRE_FALSE( config.isSet( "second" ) );
    REQUIRE( config.m_targets.at( "first" ).sources.at( 0 ) == "command line" );
    REQUIRE( config.m_unusedConfigs.size() == 1 );
    REQUIRE( config.m_unusedConfigs.begin()->second.set );
    REQUIRE( config.m_unusedConfigs.begin()->second.sources.at( 0 ) == "command line" );

    int first{ 0 };
    REQUIRE( config.get( first, "first", 0 ) == 0 );
    REQUIRE( first == 5 );
    REQUIRE( config.get( first, "first", 1 ) == -1 );

    config.parseCommandLine( 4, argv );
    int second{ 0 };
    REQUIRE( config( second, "second" ) == 0 );
    REQUIRE( second == 7 );
    REQUIRE( config.verbosity( "first" ) == 1 );
    REQUIRE( config.count( "first" ) == 2 );

    config.clear();
    REQUIRE( config.m_targets.empty() );
    REQUIRE( config.clOnlyTargets.empty() );
    REQUIRE( config.nonOptions.empty() );
    REQUIRE( config.m_unusedConfigs.empty() );
}

/** \brief Verifies that repeated command-line-only names retain their separate target registrations.
 *
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "command-line-only targets retain duplicate names", "[appConfigurator]" )
{
    mx::app::appConfigurator config;
    config.add( "verbosity", "v", "verbose", mx::app::argType::True, "", "", false, "", "" );
    config.add( "verbosity", "V", "very-verbose", mx::app::argType::True, "", "", false, "", "" );

    REQUIRE( config.m_targets.size() == 1 );
    REQUIRE( config.clOnlyTargets.size() == 1 );
    REQUIRE( config.m_targets.at( "verbosity" ).orderAdded == 0 );
    REQUIRE( config.clOnlyTargets.begin()->orderAdded == 1 );
    REQUIRE( config.clOnlyTargets.begin()->longOpt == "very-verbose" );
}

/** \brief Verifies scalar configuration defaults, typed conversions, indexed access, and source logging.
 *
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "scalar configuration access covers HCI value types", "[appConfigurator]" )
{
    mx::app::appConfigurator config;
    config.configLog = captureConfigLog;

    resetConfigLog();
    int defaultValue{ 42 };
    REQUIRE( config( defaultValue, "unset" ) == 0 );
    REQUIRE( defaultValue == 42 );
    REQUIRE( logName == "unset" );
    REQUIRE( logCode == mx::meta::typeDescription<int>::code() );
    REQUIRE( logValue == "42" );
    REQUIRE( logSource == "default" );

    resetConfigLog();
    REQUIRE( config.get( defaultValue, "unset-indexed", 0 ) == 0 );
    REQUIRE( defaultValue == 42 );
    REQUIRE( logName == "unset-indexed" );
    REQUIRE( logValue == "42" );
    REQUIRE( logSource == "default" );

    config.add( "integer", "", "", 0, "", "integer" );
    config.m_targets.at( "integer" ).set = true;
    config.m_targets.at( "integer" ).values = { "11", "23" };
    config.m_targets.at( "integer" ).sources = { "first.conf", "second.conf" };

    int integerValue{ 0 };
    resetConfigLog();
    REQUIRE( config.get( integerValue, "integer", 0 ) == 0 );
    REQUIRE( integerValue == 11 );
    REQUIRE( logValue == "11" );
    REQUIRE( logSource.empty() );
    REQUIRE( config.get( integerValue, "integer", 2 ) == -1 );

    config.m_sources = true;
    resetConfigLog();
    REQUIRE( config.get( integerValue, "integer" ) == 0 );
    REQUIRE( integerValue == 23 );
    REQUIRE( logValue == "23" );
    REQUIRE( logSource == "second.conf" );

    config.m_sources = false;

    config.add( "empty", "", "", 0, "", "empty" );
    config.m_targets.at( "empty" ).set = true;
    REQUIRE( config.get( integerValue, "empty" ) == -1 );

    config.add( "string", "", "", 0, "", "string" );
    config.add( "boolean", "", "", 0, "", "boolean" );
    config.add( "float", "", "", 0, "", "float" );
    config.add( "double", "", "", 0, "", "double" );
    config.m_targets.at( "string" ).set = true;
    config.m_targets.at( "boolean" ).set = true;
    config.m_targets.at( "float" ).set = true;
    config.m_targets.at( "double" ).set = true;
    config.m_targets.at( "string" ).values = { "target.fits" };
    config.m_targets.at( "boolean" ).values = { "true" };
    config.m_targets.at( "float" ).values = { "1.25" };
    config.m_targets.at( "double" ).values = { "2.5" };

    std::string stringValue;
    bool boolValue{ false };
    float floatValue{ 0 };
    double doubleValue{ 0 };
    REQUIRE( config( stringValue, "string" ) == 0 );
    REQUIRE( config( boolValue, "boolean" ) == 0 );
    REQUIRE( config( floatValue, "float" ) == 0 );
    REQUIRE( config( doubleValue, "double" ) == 0 );
    REQUIRE( stringValue == "target.fits" );
    REQUIRE( boolValue );
    REQUIRE( floatValue == 1.25f );
    REQUIRE( doubleValue == 2.5 );
}

/** \brief Verifies vector configuration defaults, empty values, delimiters, typed conversions, and source logging.
 *
 * \ingroup appConfigurator_unit_tests
 */
TEST_CASE( "vector configuration access covers HCI value types", "[appConfigurator]" )
{
    mx::app::appConfigurator config;
    config.configLog = captureConfigLog;

    std::vector<int> defaultValues{ 7, 8 };
    resetConfigLog();
    REQUIRE( config( defaultValues, "unset-vector" ) == 0 );
    REQUIRE( defaultValues == ( std::vector<int>{ 7, 8 } ) );
    REQUIRE( logName == "unset-vector" );
    REQUIRE( logValue == "[need a vector to string]" );
    REQUIRE( logSource == "default" );

    resetConfigLog();
    REQUIRE( config.get( defaultValues, "unset-vector-indexed", 0 ) == 0 );
    REQUIRE( defaultValues == ( std::vector<int>{ 7, 8 } ) );
    REQUIRE( logName == "unset-vector-indexed" );
    REQUIRE( logValue == "[need a vector to string]" );
    REQUIRE( logSource == "default" );

    config.add( "empty-vector", "", "", 0, "", "empty-vector" );
    config.m_targets.at( "empty-vector" ).set = true;
    config.m_targets.at( "empty-vector" ).values = { "" };
    config.m_targets.at( "empty-vector" ).sources = { "empty.conf" };
    config.m_sources = true;
    resetConfigLog();
    REQUIRE( config.get( defaultValues, "empty-vector", 0 ) == 0 );
    REQUIRE( defaultValues.empty() );
    REQUIRE( logValue.empty() );
    REQUIRE( logSource == "empty.conf" );

    config.m_sources = false;
    defaultValues = { 9 };
    resetConfigLog();
    REQUIRE( config.get( defaultValues, "empty-vector", 0 ) == 0 );
    REQUIRE( defaultValues.empty() );
    REQUIRE( logSource.empty() );
    REQUIRE( config.get( defaultValues, "empty-vector", 1 ) == -1 );

    config.add( "integer-vector", "", "", 0, "", "integer-vector" );
    config.m_targets.at( "integer-vector" ).set = true;
    config.m_targets.at( "integer-vector" ).values = { " 1,  2,3", "4,5,6" };
    config.m_targets.at( "integer-vector" ).sources = { "first.conf", "second.conf" };

    std::vector<int> integerValues;
    resetConfigLog();
    REQUIRE( config.get( integerValues, "integer-vector", 0 ) == 0 );
    REQUIRE( integerValues == ( std::vector<int>{ 1, 2, 3 } ) );
    REQUIRE( logValue == " 1,  2,3" );
    REQUIRE( logSource.empty() );

    config.m_sources = true;
    resetConfigLog();
    REQUIRE( config.get( integerValues, "integer-vector" ) == 0 );
    REQUIRE( integerValues == ( std::vector<int>{ 4, 5, 6 } ) );
    REQUIRE( logValue == "4,5,6" );
    REQUIRE( logSource == "second.conf" );

    config.m_sources = false;

    config.add( "float-vector", "", "", 0, "", "float-vector" );
    config.m_targets.at( "float-vector" ).set = true;
    config.m_targets.at( "float-vector" ).values = { "1.5,2.5" };
    std::vector<float> floatValues;
    REQUIRE( config( floatValues, "float-vector" ) == 0 );
    REQUIRE( floatValues == ( std::vector<float>{ 1.5f, 2.5f } ) );

    config.add( "string-vector", "", "", 0, "", "string-vector" );
    config.m_targets.at( "string-vector" ).set = true;
    config.m_targets.at( "string-vector" ).values = { "alpha, beta," };
    std::vector<std::string> stringValues;
    REQUIRE( config( stringValues, "string-vector" ) == 0 );
    REQUIRE( stringValues == ( std::vector<std::string>{ "alpha", "beta", "" } ) );

    config.add( "empty-values", "", "", 0, "", "empty-values" );
    config.m_targets.at( "empty-values" ).set = true;
    REQUIRE( config.get( integerValues, "empty-values" ) == -1 );
}

/*
} // namespace appConfiguratorTest
} // namespace appTest
} // namespace testUnit
 */
