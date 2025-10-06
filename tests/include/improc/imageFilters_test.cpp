/** \file imageUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#include "../../../include/improc/imageFilters.hpp"

using namespace Catch::Matchers;

namespace unitTest
{
namespace improcTest
{
namespace imageFiltersTest
{

/// Verify precalcKernel filter
/**
 * \ingroup imageFilters_unit_tests
 */
TEST_CASE( "Verify precalcKernel filter", "[improc::imageFilters]" )
{
    SECTION("with Gaussian kernel")
    {
        typedef mx::improc::gaussKernel<mx::improc::eigenImage<float>> kernelT;

        kernelT kernel(2);

        mx::improc::precalcKernel<kernelT> pck(kernel, 64,64, 31.5, 31.5);

        REQUIRE(pck.m_kernels.size() == 64*64);

        mx::improc::eigenImage<float> karr, karrpc;

        bool alleq = true;
        for(int cc = 0; cc < 64; ++cc)
        {
            for(int rr = 0; rr < 64; ++rr)
            {
                kernel.setKernel(rr-31.5, cc-31.5, karr);

                pck.setKernel(rr-31.5, cc-31.5, karrpc);

                alleq = ((karr - karrpc).sum() == 0);

                if(!alleq) break;
            }

            if(!alleq) break;
        }

        REQUIRE(alleq);
    }

    SECTION("with azBoxKernel kernel")
    {
        typedef mx::improc::azBoxKernel<mx::improc::eigenImage<float>> kernelT;

        kernelT kernel(3,5);

        mx::improc::precalcKernel<kernelT> pck(kernel, 64,64, 31.5, 31.5);

        REQUIRE(pck.m_kernels.size() == 64*64);

        mx::improc::eigenImage<float> karr, karrpc;

        bool alleq = true;
        for(int cc = 0; cc < 64; ++cc)
        {
            for(int rr = 0; rr < 64; ++rr)
            {
                kernel.setKernel(rr-31.5, cc-31.5, karr);

                pck.setKernel(rr-31.5, cc-31.5, karrpc);

                alleq = ((karr - karrpc).sum() == 0);

                if(!alleq) break;
            }

            if(!alleq) break;
        }

        REQUIRE(alleq);
    }
}

} // namespace imageFiltersTest
} // namespace improcTest
} // namespace unitTest
