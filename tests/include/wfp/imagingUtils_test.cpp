/** \file imagingUtils_test.cpp
 * \brief Tests CPU and CUDA imaging utilities.
 */
#include "../../catch2/catch.hpp"
#ifdef MXLIB_CUDA
#include "../../cudaTestUtils.hpp"
#endif

#define MX_NO_ERROR_REPORTS

#include <iostream>

#include "../../../include/wfp/imagingUtils.hpp"

/// Make a complex wavefront and extract it as an intensity image on CPU
/**
 * \ingroup imagingUtils_unit_tests
 */
TEST_CASE( "Make a complex wavefront and extract it as an intensity image on CPU", "[wfp]" )
{
    typedef float realT;
    typedef std::complex<realT> complexT;

    typedef mx::improc::eigenImage<realT> realImageT;
    typedef mx::improc::eigenImage<complexT> complexFieldT;


    int wfSz = 128;

    realImageT im(wfSz, wfSz);
    im.setZero();
    complexFieldT wf(wfSz, wfSz);

    for(int cc = 0; cc < wfSz; ++cc)
    {
        for(int rr = 0; rr < wfSz; ++rr)
        {
            wf(rr,cc) = std::complex(rr,cc);
        }
    }

    mx::wfp::extractIntensityImageAccum<realImageT, complexFieldT,0>(im,0,wfSz,0,wfSz,wf,0,0);

    bool fail = false;
    for(int cc = 0; cc < wfSz; ++cc)
    {
        for(int rr = 0; rr < wfSz; ++rr)
        {
            realT val = rr*rr + cc*cc;
            if(im(rr,cc) != val)
            {
                fail=true;
                break;
            }

        }

        if(fail)
        {
            break;
        }
    }

    REQUIRE(fail == false);
}

#if defined( MXLIB_CUDA ) || defined( __DOXY_ONLY__ )
/// Make a complex wavefront and extract it as an intensity image on GPU
/**
 * \ingroup imagingUtils_unit_tests
 */
TEST_CASE( "Make a complex wavefront and extract it as an intensity image on GPU", "[wfp]" )
{
    if( !mxlibTest::cudaDeviceAvailable() )
    {
        WARN( "CUDA runtime is available but no CUDA device is present" );
        return;
    }

    typedef float realT;
    typedef std::complex<realT> complexT;

    typedef mx::improc::eigenImage<realT> realImageT;
    typedef mx::improc::eigenImage<complexT> complexFieldT;


    int wfSz = 128;

    realImageT im(wfSz, wfSz);
    im.setZero();
    complexFieldT wf(wfSz, wfSz);

    for(int cc = 0; cc < wfSz; ++cc)
    {
        for(int rr = 0; rr < wfSz; ++rr)
        {
            wf(rr,cc) = std::complex(rr,cc);
        }
    }

    mx::cuda::cudaPtr<complexT> devWf;
    devWf.upload(wf.data(), wfSz, wfSz);

    mx::cuda::cudaPtr<realT> devIm;
    devIm.resize(wfSz, wfSz);

    mx::wfp::extractIntensityImageAccum<mx::cuda::cudaPtr<realT>, mx::cuda::cudaPtr<complexT>,1>(devIm,0,wfSz,0,wfSz,devWf,0,0);

    devIm.download(im.data());

    bool fail = false;
    for(int cc = 0; cc < wfSz; ++cc)
    {
        for(int rr = 0; rr < wfSz; ++rr)
        {
            realT val = rr*rr + cc*cc;
            if(im(rr,cc) != val)
            {
                fail=true;
                break;
            }

        }

        if(fail)
        {
            break;
        }
    }

    REQUIRE(fail == false);
}
#endif // MXLIB_CUDA
