/** \file ompLoopWatcher_test.cpp
 * \brief Tests OpenMP loop progress tracking and output formatting.
 */

#include "../../catch2/catch.hpp"

#include <sstream>
#include <string>

#include "../../../include/ipc/ompLoopWatcher.hpp"

/** \defgroup ompLoopWatcher_unit_tests OpenMP Loop Watcher Unit Tests
 * \ingroup ipc_unit_tests
 */

namespace unitTest::ipc_ompLoopWatcher_test
{

/** \brief Verifies the default ompLoopWatcher batches status output and exposes each public counter wrapper.
 *
 * \ingroup ompLoopWatcher_unit_tests
 */
TEST_CASE( "Default OpenMP loop watcher tracks and reports progress", "[ipc::ompLoopWatcher][default]" )
{
    std::ostringstream output;
    mx::ipc::ompLoopWatcher<> watcher( 20, output );

    watcher.increment();
    watcher.advance( 2 );
    for( int iteration = 0; iteration < 9; ++iteration )
    {
        watcher.outputStatus();
    }
    REQUIRE( output.str().empty() );

    watcher.outputStatus();
    REQUIRE( output.str().find( "3 / 20 (15%)" ) != std::string::npos );
    REQUIRE( output.str().find( "s/loop" ) != std::string::npos );
    REQUIRE( output.str().find( "s left" ) != std::string::npos );
    REQUIRE( output.str().back() == '\r' );

    output.str( "" );
    output.clear();
    for( int iteration = 0; iteration < 10; ++iteration )
    {
        watcher.incrementAndOutputStatus();
    }
    REQUIRE( output.str().find( "13 / 20 (65%)" ) != std::string::npos );

    output.str( "" );
    output.clear();
    for( int iteration = 0; iteration < 10; ++iteration )
    {
        watcher.advanceAndOutputStatus( 1 );
    }
    REQUIRE( output.str().find( "23 / 20 (115%)" ) != std::string::npos );

    output.str( "" );
    output.clear();
    watcher.outputFinalStatus();
    REQUIRE( output.str().find( "s elapsed" ) != std::string::npos );
    REQUIRE( output.str().find( "s/loop" ) != std::string::npos );
    REQUIRE( output.str().back() == '\r' );

    output.str( "" );
    output.clear();
    watcher.clearOutput();
    REQUIRE( output.str().find_first_not_of( ' ' ) == output.str().size() - 1 );
    REQUIRE( output.str().back() == '\r' );
}

/** \brief Verifies ompLoopWatcher newline and unformatted output policies.
 *
 * \ingroup ompLoopWatcher_unit_tests
 */
TEST_CASE( "OpenMP loop watcher honors formatting policies", "[ipc::ompLoopWatcher][format]" )
{
    std::ostringstream prettyOutput;
    mx::ipc::ompLoopWatcher<std::ostream, true, true, true, true, true> prettyWatcher( 4, prettyOutput );
    for( int iteration = 0; iteration < 10; ++iteration )
    {
        prettyWatcher.incrementAndOutputStatus();
    }
    REQUIRE( prettyOutput.str().find( "10 / 4 (250%)" ) != std::string::npos );
    REQUIRE( prettyOutput.str().back() == '\n' );

    prettyOutput.str( "" );
    prettyOutput.clear();
    prettyWatcher.outputFinalStatus();
    REQUIRE( prettyOutput.str().find( "s elapsed" ) != std::string::npos );
    REQUIRE( prettyOutput.str().back() == '\n' );

    std::ostringstream rawOutput;
    mx::ipc::ompLoopWatcher<std::ostream, false, true, true, true, true> rawWatcher( 8, rawOutput );
    for( int iteration = 0; iteration < 10; ++iteration )
    {
        rawWatcher.incrementAndOutputStatus();
    }
    REQUIRE( rawOutput.str().find( "108125" ) == 0 );
    REQUIRE( rawOutput.str().back() == '\n' );

    rawOutput.str( "" );
    rawOutput.clear();
    rawWatcher.outputFinalStatus();
    REQUIRE_FALSE( rawOutput.str().empty() );
    REQUIRE( rawOutput.str().front() == ' ' );
    REQUIRE( rawOutput.str().back() == '\n' );

    std::ostringstream terseOutput;
    mx::ipc::ompLoopWatcher<std::ostream, false, false, false, false, false> terseWatcher( 2, terseOutput );
    for( int iteration = 0; iteration < 10; ++iteration )
    {
        terseWatcher.incrementAndOutputStatus();
    }
    REQUIRE( terseOutput.str() == "10" );

    terseOutput.str( "" );
    terseOutput.clear();
    terseWatcher.outputFinalStatus();
    REQUIRE( terseOutput.str().empty() );
    terseWatcher.clearOutput();
    REQUIRE( terseOutput.str().empty() );
}

} // namespace unitTest::ipc_ompLoopWatcher_test
