
#define EIGEN_DONT_PARALLELIZE ( 1 )

#include <omp.h>

#include "../../../../include/ao/sim/pyramidSensor.hpp"
#include "../../../../include/ao/sim/ccdDetector.hpp"

#include "../../../../include/improc/eigenImage.hpp"
#include "../../../../include/improc/imageMasks.hpp"

#include "../../../../include/ioutils/fits/fitsFile.hpp"

#include "../../../../include/sys/timeUtils.hpp"

using namespace mx::sys;
using namespace mx::AO::sim;
using namespace mx::improc;
using namespace mx::fits;

int main()
{
    typedef float realT;
    typedef ccdDetector<realT> detectorT;

    pyramidSensor<realT, detectorT> pwfs;

    int wfSz = 768;
    int pupSz = 336;
    realT D = 6.5;

    eigenImage<realT> pupilMask( wfSz, wfSz );
    pupilMask.setZero();
    maskCircle( pupilMask, 0.5 * ( wfSz - 1.0 ), 0.5 * ( wfSz - 1.0 ), 0.5 * pupSz, 1.0 );

    pyramidSensor<realT, detectorT>::wavefrontT wf;
    wf.phase.resize( wfSz, wfSz );
    wf.phase.setZero();
    wf.amplitude = pupilMask;

    pwfs.wfSz( wfSz );
    pwfs.detSize( 120., 120. );

    pwfs.quadSz( pupSz * 60. / 56. );
    pwfs.wfPS( D / pupSz );
    pwfs.lambda( 0.8e-6 );
    pwfs.D( D );
    pwfs.perStep( 1 );
    pwfs.modRadius( 3.0 );
    pwfs.wholePixelModulation( true );
    pwfs.preAllocate();

    timespec ts;
    double t0 = get_curr_time( ts );
    for( int i = 0; i < 100; ++i )
        pwfs.doSenseWavefront( wf );
    double t1 = get_curr_time( ts );

    std::cerr << ( t1 - t0 ) / 100. << "\n";
    fitsFile<realT> ff;
    ff.write( "wfsImage.fits", pwfs.m_wfsImage.image );
}
