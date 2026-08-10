/** \file application_test.cpp
 * \brief Tests for the application framework lifecycle.
 */
#include "../../catch2/catch.hpp"

#include <filesystem>
#include <fstream>

#include "../../../include/app/application.hpp"

/** \defgroup application_unit_tests application Unit Tests
 * \ingroup app_unit_tests
 */

/// \cond
namespace unitTest::app_application_test
{

class temporaryConfigFile
{
  public:
    explicit temporaryConfigFile( const std::string &contents )
        : m_path( std::filesystem::temp_directory_path() / "mxlib_application_test.conf" )
    {
        std::ofstream file{ m_path };
        file << contents;
    }

    ~temporaryConfigFile()
    {
        std::filesystem::remove( m_path );
    }

    std::string string() const
    {
        return m_path.string();
    }

  private:
    std::filesystem::path m_path;
};

class lifecycleApplication : public mx::app::application
{
  public:
    lifecycleApplication()
    {
        m_requireConfigPathGlobal = false;
        m_requireConfigPathUser = false;
        m_requireConfigPathLocal = false;
    }

    int m_setupCalls{ 0 };
    int m_loadCalls{ 0 };
    int m_executeCalls{ 0 };
    int m_helpCalls{ 0 };
    int m_value{ 0 };
    bool m_configAvailableInExecute{ false };

    void preserveConfig( bool preserve )
    {
        m_preserveConfig = preserve;
    }

    void configOnly( bool configOnly )
    {
        m_configOnly = configOnly;
    }

    void requireGlobalConfig( bool require )
    {
        m_requireConfigPathGlobal = require;
    }

  protected:
    void setupConfig() override
    {
        ++m_setupCalls;
        config.add( "value", "v", "value", mx::app::argType::Required, "test", "value", false, "int", "" );
    }

    void loadConfig() override
    {
        ++m_loadCalls;
        config( m_value, "value" );
    }

    int execute() override
    {
        ++m_executeCalls;
        m_configAvailableInExecute = config.isSet( "value" );
        return m_value;
    }

    void help() override
    {
        ++m_helpCalls;
    }
};

} // namespace unitTest::app_application_test
/// \endcond

namespace unitTest::app_application_test
{

/** \brief Verifies that mx::app::application runs setup, loading, and execution in order.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application runs the configured lifecycle", "[app::application]" )
{
    lifecycleApplication app;
    char invokedName[] = "lifecycleApplication";
    char valueOption[] = "--value=17";
    char *argv[] = { invokedName, valueOption };

    REQUIRE( app.main( 2, argv ) == 17 );
    REQUIRE( app.m_setupCalls == 1 );
    REQUIRE( app.m_loadCalls == 1 );
    REQUIRE( app.m_executeCalls == 1 );
    REQUIRE( app.m_value == 17 );
    REQUIRE_FALSE( app.m_configAvailableInExecute );
}

/** \brief Verifies that mx::app::application preserves configuration only when requested.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application optionally preserves configuration for execution", "[app::application]" )
{
    lifecycleApplication app;
    app.preserveConfig( true );
    char invokedName[] = "lifecycleApplication";
    char valueOption[] = "--value=23";
    char *argv[] = { invokedName, valueOption };

    REQUIRE( app.main( 2, argv ) == 23 );
    REQUIRE( app.m_value == 23 );
    REQUIRE( app.m_configAvailableInExecute );
}

/** \brief Verifies that mx::app::application short-circuits execution for help and config-only requests.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application short-circuits help and config-only requests", "[app::application]" )
{
    SECTION( "help" )
    {
        lifecycleApplication app;
        char invokedName[] = "lifecycleApplication";
        char helpOption[] = "--help";
        char *argv[] = { invokedName, helpOption };

        REQUIRE( app.main( 2, argv ) == 1 );
        REQUIRE( app.m_setupCalls == 1 );
        REQUIRE( app.m_loadCalls == 0 );
        REQUIRE( app.m_helpCalls == 1 );
        REQUIRE( app.m_executeCalls == 0 );
    }

    SECTION( "config-only" )
    {
        lifecycleApplication app;
        app.configOnly( true );
        char invokedName[] = "lifecycleApplication";
        char valueOption[] = "--value=41";
        char *argv[] = { invokedName, valueOption };

        REQUIRE( app.main( 2, argv ) == 1 );
        REQUIRE( app.m_setupCalls == 1 );
        REQUIRE( app.m_loadCalls == 1 );
        REQUIRE( app.m_executeCalls == 0 );
    }
}

/** \brief Verifies configuration-file loading, command-line precedence, and required-file failure behavior.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application loads configuration files", "[app::application]" )
{
    SECTION( "command-line config and value override" )
    {
        temporaryConfigFile configFile{ "[test]\nvalue=31\n" };
        lifecycleApplication app;
        char invokedName[] = "lifecycleApplication";
        std::string configOption = "--config=" + configFile.string();
        char valueOption[] = "--value=37";
        char *argv[] = { invokedName, configOption.data(), valueOption };

        REQUIRE( app.main( 3, argv ) == 37 );
        REQUIRE( app.m_loadCalls == 1 );
        REQUIRE( app.m_executeCalls == 1 );
        REQUIRE( app.m_value == 37 );
    }

    SECTION( "required global config failure" )
    {
        lifecycleApplication app;
        app.requireGlobalConfig( true );
        app.setConfigPathGlobal( "/tmp/mxlib_application_test_missing.conf" );
        char invokedName[] = "lifecycleApplication";
        char *argv[] = { invokedName };

        REQUIRE( app.main( 1, argv ) == 1 );
        REQUIRE( app.m_loadCalls == 0 );
        REQUIRE( app.m_helpCalls == 1 );
        REQUIRE( app.m_executeCalls == 0 );
    }
}

} // namespace unitTest::app_application_test
