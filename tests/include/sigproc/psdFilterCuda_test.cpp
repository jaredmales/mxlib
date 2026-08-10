/** \file psdFilterCuda_test.cpp
 * \brief Documentation-only ownership placeholder for the legacy psdFilterCuda header.
 */
#include "../../catch2/catch.hpp"

#include "../../../include/sigproc/psdFilterCuda.hpp"

/** \defgroup psdFilterCuda_unit_tests psdFilterCuda Unit Tests
 * \ingroup sigproc_unit_tests
 */

namespace unitTest::placeholder::sigproc_psdFilterCuda_test
{

/** \brief Keeps ownership of the legacy psdFilterCuda header visible in the documented test inventory.
 *
 * The header predates the rank-aware psdFilter API and depends on removed CUDA helper APIs, so this TEST_CASE is
 * intentionally not compiled.
 *
 * \ingroup psdFilterCuda_unit_tests
 * \todo Decide whether to remove or modernize the header before enabling compilation and behavioral tests for the
 * CUDA PSD filter.
 */
TEST_CASE( "psdFilterCuda header has a documentation-only ownership placeholder",
           "[sigproc::psdFilterCuda][placeholder]" )
{
    SUCCEED( "psdFilterCuda remains excluded pending removal or modernization." );
}

} // namespace unitTest::placeholder::sigproc_psdFilterCuda_test
