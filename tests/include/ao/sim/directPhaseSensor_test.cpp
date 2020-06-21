
#include "../../../../include/math/vectorUtils.hpp"
using namespace mx::math;

#include "../../../../include/improc/fitsFile.hpp"
#include "../../../../include/improc/eigenImage.hpp"
using namespace mx::improc;

#include "../../../../include/sigproc/fourierModes.hpp"
using namespace mx::sigproc;

#include "../../../../include/ao/sim/directPhaseSensor.hpp"
#include "../../../../include/ao/sim/ccdDetector.hpp"

using namespace mx::AO::sim;

int main()
{
   typedef double realT;
   
   typedef wavefront<realT> wavefrontT;
    
   int Dpix = 128;
   
   directPhaseSensor<realT, ccdDetector<realT>> dps;
   
   
   dps.lambda(0.8e-6);
   dps.detSize(Dpix, Dpix);
   dps.m_detector.ron(0.5);
   dps.m_detector.gain(500.0);
   dps.m_detector.qe(0.25);

   dps.filterWidth(24.0);
   dps.applyFilter(false);

   dps.simStep(1);
   dps.iTime(1);
   dps.roTime(1);
   
   dps.beta_p(1.0);
   
   eigenImage<realT> pupil;
   eigenImage<realT> mode;
   
   fitsFile<realT> ff;
   
   ff.read(pupil, "/home/jrmales/Data/mxAO/pupil/circular_0percent_128_0os/pupil.fits");
   dps.pupil(pupil);
   
   mode.resize(Dpix, Dpix);
   Eigen::Map<eigenImage<realT>> mm(mode.data(), mode.rows(), mode.cols());
   makeModifiedFourierMode( mm, 10, 10, -1);

   mode *= pupil;
   
   wavefrontT wf;
   wf.setAmplitude(pupil);
   wf.setPhase(mode);
   
   realT psum = pupil.sum();
   
   realT amp0 = (mode*wf.phase).sum()/psum;
   
   wf.iterNo = 0;
   bool newCV = dps.senseWavefront(wf);
   wf.iterNo = 1;
   std::vector<realT> amps(1000);
   
   realT F0 = 10000000;
   
   wf.amplitude = pupil * sqrt(F0)/sqrt(psum);
   
   newCV = dps.senseWavefront(wf);
   wf.iterNo = 2;
   
   for(size_t n=0; n<amps.size();++n)
   {
      newCV = dps.senseWavefront(wf);
      ++wf.iterNo;
   
      amps[n] = (mode*dps.m_detectorImage.image).sum()/psum;
   
      //std::cout << amps[n] << "\n";
   }
   
   realT mnAmp = vectorMean(amps);
   std::cerr << amp0 << " " << mnAmp << " " << sqrt(vectorVariance(amps))/mnAmp << " " <<  1./sqrt(F0)<< "\n";
   
   return 0;
}
