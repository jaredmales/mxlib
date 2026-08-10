/** \file ompLoopWatcher_test.cpp
 * \brief Placeholder tests for the ompLoopWatcher production API.
 */

#include "../../catch2/catch.hpp"

#include "../../../include/ipc/ompLoopWatcher.hpp"

/** \cond
 * Explicit instantiation compile-checks the default OpenMP loop watcher specialization and emits its header-defined
 * methods for coverage accounting in this test translation unit.
 */
template class mx::ipc::ompLoopWatcher<>;
/** \endcond */

/** \defgroup ompLoopWatcher_unit_tests OpenMP Loop Watcher Unit Tests
 * \ingroup ipc_unit_tests
 */

namespace unitTest::placeholder::ipc_ompLoopWatcher_test
{

/** \brief Records pending behavioral verification of the ompLoopWatcher production API.
 *
 * \ingroup ompLoopWatcher_unit_tests
 * \todo Add behavioral assertions for the production API declared in \c include/ipc/ompLoopWatcher.hpp.
 */
TEST_CASE( "ompLoopWatcher production API placeholder", "[ipc::ompLoopWatcher][placeholder]" )
{
    SUCCEED( "Behavioral API assertions are pending." );
}

} // namespace unitTest::placeholder::ipc_ompLoopWatcher_test
