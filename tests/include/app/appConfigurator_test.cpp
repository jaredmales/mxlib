/** \file appConfigurator_test.cpp
 *  \brief Tests application configuration parsing.
 */

#include "../../catch2/catch.hpp"

#include <fstream>

#include "../../../include/app/appConfigurator.hpp"

using namespace mx::app;

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

/*
} // namespace appConfiguratorTest
} // namespace appTest
} // namespace testUnit
 */
