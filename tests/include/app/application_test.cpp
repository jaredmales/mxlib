/** \file application_test.cpp
 * \brief Tests for the application framework lifecycle.
 */
#include "../../catch2/catch.hpp"

#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <sstream>

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

class scopedEnvironment
{
  public:
    scopedEnvironment( const std::string &name, const std::string &value ) : m_name( name )
    {
        if( const char *oldValue = std::getenv( m_name.c_str() ) )
        {
            m_oldValue = oldValue;
            m_hadOldValue = true;
        }

        setenv( m_name.c_str(), value.c_str(), 1 );
    }

    ~scopedEnvironment()
    {
        if( !m_hadOldValue )
        {
            unsetenv( m_name.c_str() );
        }
        else
        {
            setenv( m_name.c_str(), m_oldValue.c_str(), 1 );
        }
    }

  private:
    std::string m_name;
    std::string m_oldValue;
    bool m_hadOldValue{ false };
};

class scopedCerr
{
  public:
    explicit scopedCerr( std::ostream &stream ) : m_oldBuffer( std::cerr.rdbuf( stream.rdbuf() ) )
    {
    }

    ~scopedCerr()
    {
        std::cerr.rdbuf( m_oldBuffer );
    }

  private:
    std::streambuf *m_oldBuffer;
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

    void useDefaultHelp( bool useDefault )
    {
        m_useDefaultHelp = useDefault;
    }

    void setPathEnvironmentNames( const std::string &global,
                                  const std::string &user,
                                  const std::string &local,
                                  const std::string &configBase )
    {
        m_configPathGlobal_env = global;
        m_configPathUser_env = user;
        m_configPathLocal_env = local;
        m_configPathCLBase_env = configBase;
    }

    void applyDefaults()
    {
        setDefaults( 0, nullptr );
    }

    int reRead()
    {
        return reReadConfig();
    }

    int configuredValueCount()
    {
        return config.count( "value" );
    }

    const std::string &globalPath() const
    {
        return m_configPathGlobal;
    }

    const std::string &userPath() const
    {
        return m_configPathUser;
    }

    const std::string &localPath() const
    {
        return m_configPathLocal;
    }

    const std::string &configBasePath() const
    {
        return m_configPathCLBase;
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
        if( m_useDefaultHelp )
        {
            application::help();
        }
    }

  private:
    bool m_useDefaultHelp{ false };
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

/** \brief Verifies that mx::app::application resolves environment-selected configuration search paths.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application resolves configuration paths from the environment", "[app::application]" )
{
    scopedEnvironment global{ "MXLIB_APPLICATION_TEST_GLOBAL", "/tmp/global.conf" };
    scopedEnvironment user{ "MXLIB_APPLICATION_TEST_USER", "/tmp/user.conf" };
    scopedEnvironment local{ "MXLIB_APPLICATION_TEST_LOCAL", "/tmp/local-config" };
    scopedEnvironment configBase{ "MXLIB_APPLICATION_TEST_CONFIG_BASE", "/tmp/config-base" };

    lifecycleApplication app;
    app.setConfigPathGlobal( "ignored-global" );
    app.setConfigPathUser( "ignored-user" );
    app.setConfigPathLocal( "ignored-local" );
    app.setPathEnvironmentNames( "MXLIB_APPLICATION_TEST_GLOBAL",
                                 "MXLIB_APPLICATION_TEST_USER",
                                 "MXLIB_APPLICATION_TEST_LOCAL",
                                 "MXLIB_APPLICATION_TEST_CONFIG_BASE" );
    app.applyDefaults();

    REQUIRE( app.globalPath() == "/tmp/global.conf" );
    REQUIRE( app.userPath() == "/tmp/user.conf" );
    REQUIRE( app.localPath() == "/tmp/local-config/" );
    REQUIRE( app.configBasePath() == "/tmp/config-base/" );
}

/** \brief Verifies that mx::app::application emits its default help for registered options.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application renders default help", "[app::application]" )
{
    lifecycleApplication app;
    app.useDefaultHelp( true );
    char invokedName[] = "lifecycleApplication";
    char helpOption[] = "--help";
    char *argv[] = { invokedName, helpOption };
    std::ostringstream help;
    scopedCerr capture{ help };

    REQUIRE( app.main( 2, argv ) == 1 );
    REQUIRE( help.str().find( "usage: lifecycleApplication" ) != std::string::npos );
    REQUIRE( help.str().find( "--value" ) != std::string::npos );
}

/** \brief Verifies that mx::app::application re-reads registered configuration sources without rebuilding targets.
 *
 * \ingroup application_unit_tests
 */
TEST_CASE( "application re-reads its configuration stack", "[app::application]" )
{
    temporaryConfigFile configFile{ "[test]\nvalue=29\n" };
    lifecycleApplication app;
    app.preserveConfig( true );
    char invokedName[] = "lifecycleApplication";
    std::string configOption = "--config=" + configFile.string();
    char *argv[] = { invokedName, configOption.data() };

    REQUIRE( app.main( 2, argv ) == 29 );
    REQUIRE( app.configuredValueCount() == 1 );
    REQUIRE( app.reRead() == 0 );
    REQUIRE( app.configuredValueCount() == 2 );
    REQUIRE( app.m_setupCalls == 1 );
    REQUIRE( app.m_loadCalls == 1 );
}

} // namespace unitTest::app_application_test
