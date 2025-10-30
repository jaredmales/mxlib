/** \file pyramidSensor_test.cpp
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/sim/pyramidSensor.hpp"
#include "../../../../include/ao/sim/ccdDetector.hpp"
using namespace mx::AO::sim;


#include "../../../../include/ioutils/fits/fitsFile.hpp"
using namespace mx::fits;

#include "../../../../include/improc/eigenImage.hpp"
#include "../../../../include/improc/imageMasks.hpp"
//#include <mx/improc/milkImage.hpp>
//#include <mx/improc/eigenCube.hpp>
using namespace mx::improc;


/// Simulate a pyramid sensor on CPU
/**
 * \ingroup ao_sim_pyramidSensor_tests
 */
TEST_CASE( "Simulate a pyramid sensor on CPU", "[ao::sim]" )
{
    typedef float realT;

    uint32_t pupSz = 56.0;
    uint32_t wfSz = 256.0;


    pyramidSensor<realT, ccdDetector<realT>> pwfs;

    eigenImage<realT> pupil, fm;
    pupil.resize(pupSz, pupSz);
    pupil.setConstant(0);
    maskCircle(pupil, 0.5 * (1.0 * pupSz - 1), 1);

    wavefront<realT> wf;
    wf.setAmplitude(pupil * 1e6 * pow(64.0 / pupSz, 1.));

    fm.resize(pupSz, pupSz);

    pwfs.D(6.5);
    pwfs.pupilSz(pupSz);
    pwfs.wfSz(wfSz);
    pwfs.pupilSep(1 + 4.1 / 56.);

    pwfs.perStep(0.5);
    pwfs.modRadius(3);

    pwfs.detSize(120, 120);
    pwfs.detector.noNoise(true);
    pwfs.detector.expTime(1);

    pwfs.lambda(800e-9);

    wf.setPhase(fm * 0);
    pwfs.senseWavefrontCal(wf);
    eigenImage<realT> ref = pwfs.detectorImage.image;

    fitsFile<realT> ff;
    ff.write("ref.fits", ref);



}

