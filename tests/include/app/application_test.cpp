/** \file application_test.cpp
 * \brief Tests for the application framework lifecycle.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/app/application.hpp"

/** \defgroup application_unit_tests application Unit Tests
 * \ingroup app_unit_tests
 */

/// \cond
namespace unitTest::app_application_test
{

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
    int m_value{ 0 };
    bool m_configAvailableInExecute{ false };

    void preserveConfig( bool preserve )
    {
        m_preserveConfig = preserve;
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

} // namespace unitTest::app_application_test
