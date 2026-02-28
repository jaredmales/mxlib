/** \file pyramidSensor_test.cpp
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

// #define DEBUG

#include "../../../../include/ao/sim/pyramidSensor.hpp"
#include "../../../../include/ao/sim/ccdDetector.hpp"
using namespace mx::AO::sim;

#include "../../../../include/improc/eigenImage.hpp"
#include "../../../../include/improc/imageMasks.hpp"
// #include <mx/improc/milkImage.hpp>
// #include <mx/improc/eigenCube.hpp>
using namespace mx::improc;
/// Simulate a pyramid sensor on CPU
/**
 * \ingroup pyramidSensor_unit_tests
 */
TEST_CASE( "Simulate a pyramid sensor on CPU", "[ao::sim]" )
{
    typedef float realT;

    uint32_t pupSz = 56.0;
    uint32_t wfSz = 256.0;

    pyramidSensor<realT, ccdDetector<realT>> pwfs;

    eigenImage<realT> pupil, fm;
    pupil.resize( pupSz, pupSz );
    pupil.setConstant( 0 );
    maskCircle( pupil, 0.5 * ( 1.0 * pupSz - 1 ), 1 );

    wavefront<realT> wf;
    realT nphperpix = 1e10 / pupil.sum();

    wf.setAmplitude( pupil * sqrt( nphperpix ) );

    fm.resize( pupSz, pupSz );

    pwfs.D( 6.5 );
    pwfs.pupilSz( pupSz );
    pwfs.wfSz( wfSz );
    pwfs.pupilSep( 1 + 4.1 / 56. );

    pwfs.perStep( 0.5 );
    pwfs.modRadius( 3 );

    pwfs.detSize( 120, 120 );
    pwfs.detector.noNoise( true );
    pwfs.detector.expTime( 1 );

    pwfs.lambda( 850e-9 );

    wf.setPhase( fm * 0 );

    pwfs.senseWavefrontCal( wf );

    size_t N = 10;
    double t0 = mx::sys::get_curr_time();
    for( size_t n = 0; n < N; ++n )
    {
        pwfs.senseWavefrontCal( wf );
    }
    double t1 = mx::sys::get_curr_time();

    std::cerr << "\nCPU: " << 1.0 * N / ( t1 - t0 ) << " fps\n\n";

    eigenImage<realT> ref = pwfs.detectorImage.image();

    realT s = ref.sum();
    std::cout << s << '\n';

    REQUIRE( s > 9.9e9 );
}
/// Simulate a pyramid sensor on GPU
/**
 * \ingroup pyramidSensor_unit_tests
 */
TEST_CASE( "Simulate a pyramid sensor on GPU", "[ao::sim]" )
{
    typedef float realT;

    uint32_t pupSz = 56.0;
    uint32_t wfSz = 256.0;

    omp_set_num_threads( 4 );

    pyramidSensor<realT, ccdDetector<realT>, 1> pwfs;

    eigenImage<realT> pupil, fm;
    pupil.resize( pupSz, pupSz );
    pupil.setConstant( 0 );
    maskCircle( pupil, 0.5 * ( 1.0 * pupSz - 1 ), 1 );

    wavefront<realT> wf;
    realT nphperpix = 1e10 / pupil.sum();

    wf.setAmplitude( pupil * sqrt( nphperpix ) );

    fm.resize( pupSz, pupSz );

    pwfs.D( 6.5 );
    pwfs.pupilSz( pupSz );
    pwfs.wfSz( wfSz );
    pwfs.pupilSep( 1 + 4.1 / 56. );

    pwfs.perStep( 0.5 );
    pwfs.modRadius( 3 );

    pwfs.detSize( 120, 120 );
    pwfs.detector.noNoise( true );
    pwfs.detector.expTime( 1 );

    pwfs.lambda( 850e-9 );

    wf.setPhase( fm * 0 );

    pwfs.senseWavefrontCal( wf );

    size_t N = 10;
    double t0 = mx::sys::get_curr_time();
    for( size_t n = 0; n < N; ++n )
    {
        pwfs.senseWavefrontCal( wf );
    }
    double t1 = mx::sys::get_curr_time();

    std::cerr << "\nGPU: " << 1.0 * N / ( t1 - t0 ) << " fps\n\n";

    eigenImage<realT> ref = pwfs.detectorImage.image();

    realT s = ref.sum();
    std::cout << s << '\n';

    REQUIRE( s > 9.9e9 );
}
