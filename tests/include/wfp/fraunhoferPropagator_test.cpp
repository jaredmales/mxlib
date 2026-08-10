/** \file fraunhoferPropagator_test.cpp
 * \brief Tests CPU and CUDA Fraunhofer propagation.
 */
#include "../../catch2/catch.hpp"
#ifdef MXLIB_CUDA
#include "../../cudaTestUtils.hpp"
#endif

#define MX_NO_ERROR_REPORTS

#include <iostream>

// #define DEBUG

#include "../../../include/wfp/fraunhoferPropagator.hpp"
#include "../../../include/improc/eigenImage.hpp"
#include "../../../include/improc/imageMasks.hpp"
#include "../../../include/math/constants.hpp"

/// Make an Airy pattern and go back to pupil on CPU
/**
 * \ingroup fraunhoferPropagator_unit_tests
 */
TEST_CASE( "Make an Airy pattern and go back to pupil on CPU", "[wfp]" )
{
    typedef float realT;
    typedef std::complex<realT> complexT;
    typedef mx::improc::eigenImage<std::complex<realT>> complexFieldT;

    int wfSz = 512;
    int pupSz = 128;
    mx::improc::eigenImage<realT> pupil;
    pupil.resize( pupSz, pupSz );
    pupil.setZero();
    mx::improc::maskCircle( pupil, 63.5, 63.5, 63.5, 1 );

    complexFieldT complexPupil, complexFocal;
    mx::improc::eigenImage<realT> realFocal, realPupil;

    complexPupil.resize( wfSz, wfSz );
    complexFocal.resize( wfSz, wfSz );
    realFocal.resize( wfSz, wfSz );
    realPupil.resize( wfSz, wfSz );

    mx::wfp::fraunhoferPropagator<complexFieldT> fi;
    mx::wfp::makeComplexPupil( complexPupil, pupil, wfSz );

    fi.propagatePupilToFocal( complexFocal, complexPupil );

    mx::wfp::extractIntensityImage( realFocal, 0, complexFocal.rows(), 0, complexFocal.cols(), complexFocal, 0, 0 );

    realT expwr = pupil.square().sum();
    realT pwr = realFocal.sum();

    REQUIRE_THAT( expwr, Catch::Matchers::WithinAbs( pwr, ( pwr ) * ( 1e-4 ) ) );

    complexPupil.setZero();
    fi.propagateFocalToPupil( complexPupil, complexFocal );

    mx::wfp::extractIntensityImage( realPupil, 0, complexPupil.rows(), 0, complexPupil.cols(), complexPupil, 0, 0 );

    REQUIRE_THAT( realPupil.sum(), Catch::Matchers::WithinAbs( pupil.sum(), pupil.sum() * 1e-4 ) );
}

#if defined( MXLIB_CUDA ) || defined( __DOXY_ONLY__ )
/// Make an Airy pattern and go back to pupil on GPU
/**
 * \ingroup fraunhoferPropagator_unit_tests
 */
TEST_CASE( "Make an Airy pattern and go back to pupil on GPU", "[wfp]" )
{
    if( !mxlibTest::cudaDeviceAvailable() )
    {
        WARN( "CUDA runtime is available but no CUDA device is present" );
        return;
    }

    typedef float realT;
    typedef std::complex<realT> complexT;
    typedef mx::improc::eigenImage<std::complex<realT>> complexFieldT;

    int wfSz = 512;
    int pupSz = 128;
    mx::improc::eigenImage<realT> pupil;
    pupil.resize( pupSz, pupSz );
    pupil.setZero();
    mx::improc::maskCircle( pupil, 63.5, 63.5, 63.5, 1 );

    complexFieldT complexPupil, complexFocal;
    mx::improc::eigenImage<realT> realFocal, realPupil;

    complexPupil.resize( wfSz, wfSz );
    complexFocal.resize( wfSz, wfSz );
    realFocal.resize( wfSz, wfSz );
    realPupil.resize( wfSz, wfSz );

    mx::wfp::fraunhoferPropagator<complexFieldT, 1> fi;
    mx::wfp::makeComplexPupil( complexPupil, pupil, wfSz );

    mx::cuda::cudaPtr<complexT> dev_complexPupil, dev_complexFocal;

    dev_complexPupil.upload( complexPupil.data(), complexPupil.rows(), complexPupil.cols() );
    dev_complexFocal.resize( complexPupil.rows(), complexPupil.cols() );

    BREAD_CRUMB;

    fi.propagatePupilToFocal( dev_complexFocal, dev_complexPupil );

    BREAD_CRUMB;

    dev_complexFocal.download( complexFocal.data() );

    BREAD_CRUMB;

    mx::wfp::extractIntensityImage( realFocal, 0, complexFocal.rows(), 0, complexFocal.cols(), complexFocal, 0, 0 );

    BREAD_CRUMB;

    realT expwr = pupil.square().sum();
    realT pwr = realFocal.sum();

    REQUIRE_THAT( expwr, Catch::Matchers::WithinAbs( pwr, ( pwr ) * ( 1e-4 ) ) );

    BREAD_CRUMB;

    fi.propagateFocalToPupil( dev_complexPupil, dev_complexFocal );

    BREAD_CRUMB;

    complexPupil.setZero(); // make sure
    dev_complexPupil.download( complexPupil.data() );

    mx::wfp::extractIntensityImage( realPupil, 0, complexPupil.rows(), 0, complexPupil.cols(), complexPupil, 0, 0 );

    REQUIRE_THAT( realPupil.sum(), Catch::Matchers::WithinAbs( pupil.sum(), pupil.sum() * 1e-4 ) );
}
#endif // MXLIB_CUDA
