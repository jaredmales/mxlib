/** \file psdFilter_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/circularBuffer.hpp"

/** Scenario: using a circular buffer with branch-wrapping
 *
 * Verify circular buffer operation with branch-wrapping
 *
 * \anchor tests_sigproc_circularBuffer_branch
 */
SCENARIO( "creating a circular buffer with branching", "[sigproc::circularBuffer::circularBufferBranch]" )
{
    GIVEN( "a circular buffer" )
    {
        WHEN( "adding exactly max entries worth" )
        {
            mx::sigproc::circularBufferBranch<int, int> cb;
            cb.maxEntries(5);
            cb.nextEntry(0);
            cb.nextEntry(1);
            cb.nextEntry(2);
            cb.nextEntry(3);
            cb.nextEntry(4);

            REQUIRE( cb[0] == 0 );
            REQUIRE( cb[1] == 1 );
            REQUIRE( cb[2] == 2 );
            REQUIRE( cb[3] == 3 );
            REQUIRE( cb[4] == 4 );

            //reverse
            REQUIRE( cb[-1] == 4 );
            REQUIRE( cb[-2] == 3 );
            REQUIRE( cb[-3] == 2 );
            REQUIRE( cb[-4] == 1 );

            REQUIRE( cb.at(cb.earliest(),0) == 0 );
            REQUIRE( cb.at(cb.latest(),0) == 4 );
        }

        WHEN( "adding new values past the end" )
        {
            mx::sigproc::circularBufferBranch<int, int> cb;
            cb.maxEntries(5);
            cb.nextEntry(0);
            cb.nextEntry(1);
            cb.nextEntry(2);
            cb.nextEntry(3);
            cb.nextEntry(4);
            cb.nextEntry(5);
            cb.nextEntry(6);

            REQUIRE( cb[0] == 2 );
            REQUIRE( cb[1] == 3 );
            REQUIRE( cb[2] == 4 );
            REQUIRE( cb[3] == 5 );
            REQUIRE( cb[4] == 6 );

            //reverse
            REQUIRE( cb[-1] == 6 );
            REQUIRE( cb[-2] == 5 );
            REQUIRE( cb[-3] == 4 );
            REQUIRE( cb[-4] == 3 );

            REQUIRE( cb.at(cb.earliest(),0) == 2 );
            REQUIRE( cb.at(cb.latest(),0) == 6 );
        }

        WHEN( "wrapping when not full" )
        {
            mx::sigproc::circularBufferBranch<int, int> cb;
            cb.maxEntries(5);
            cb.nextEntry(0);
            cb.nextEntry(1);
            cb.nextEntry(2);
            cb.nextEntry(3);

            REQUIRE( cb[0] == 0 );
            REQUIRE( cb[1] == 1 );
            REQUIRE( cb[2] == 2 );
            REQUIRE( cb[3] == 3 );
            REQUIRE( cb[4] == 0 );

            //reverse
            REQUIRE( cb[-1] == 3 );
            REQUIRE( cb[-2] == 2 );
            REQUIRE( cb[-3] == 1 );
            REQUIRE( cb[-4] == 0 );
        }
    }
}

/** Scenario: using a circular buffer with index-wrapping
 *
 * Verify circular buffer operation with index-wrapping
 *
 * \anchor tests_sigproc_circularBuffer_branch
 */
SCENARIO( "creating a circular buffer with indexing", "[sigproc::circularBuffer::circularBufferIndex]" )
{
    GIVEN( "a circular buffer" )
    {
        WHEN( "adding exactly max entries worth" )
        {
            mx::sigproc::circularBufferIndex<int, int> cb;
            cb.maxEntries(5);
            cb.nextEntry(0);
            cb.nextEntry(1);
            cb.nextEntry(2);
            cb.nextEntry(3);
            cb.nextEntry(4);

            REQUIRE( cb[0] == 0 );
            REQUIRE( cb[1] == 1 );
            REQUIRE( cb[2] == 2 );
            REQUIRE( cb[3] == 3 );
            REQUIRE( cb[4] == 4 );

            //reverse
            REQUIRE( cb[-1] == 4 );
            REQUIRE( cb[-2] == 3 );
            REQUIRE( cb[-3] == 2 );
            REQUIRE( cb[-4] == 1 );

            REQUIRE( cb.at(cb.earliest(),0) == 0 );
            REQUIRE( cb.at(cb.latest(),0) == 4 );
        }

        WHEN( "adding new values past the end" )
        {
            mx::sigproc::circularBufferIndex<int, int> cb;
            cb.maxEntries(5);
            cb.nextEntry(0);
            cb.nextEntry(1);
            cb.nextEntry(2);
            cb.nextEntry(3);
            cb.nextEntry(4);
            cb.nextEntry(5);
            cb.nextEntry(6);

            REQUIRE( cb[0] == 2 );
            REQUIRE( cb[1] == 3 );
            REQUIRE( cb[2] == 4 );
            REQUIRE( cb[3] == 5 );
            REQUIRE( cb[4] == 6 );

            //reverse
            REQUIRE( cb[-1] == 6 );
            REQUIRE( cb[-2] == 5 );
            REQUIRE( cb[-3] == 4 );
            REQUIRE( cb[-4] == 3 );

            REQUIRE( cb.at(cb.earliest(),0) == 2 );
            REQUIRE( cb.at(cb.latest(),0) == 6 );
        }
    }
}

