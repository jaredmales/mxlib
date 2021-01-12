


#ifndef wavefront_hpp
#define wavefront_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include "../../wfp/imagingArray.hpp"
#include "../../wfp/imagingUtils.hpp"


namespace mx
{
namespace AO
{
namespace sim
{
   
///Structure containing the phase and amplitude of a wavefront   
/** The phase and amplitude are stored separately in real valued arrays.
  */
template<typename _realT>
struct wavefront
{
   ///The floating point type used for all calculations
   typedef _realT realT;
   
   ///The data type of the real image, used to hold the phase and amplitude
   typedef  Eigen::Array<_realT, Eigen::Dynamic, Eigen::Dynamic>  realImageT; 

   ///The wavefront data type
   typedef wfp::imagingArray<std::complex<realT>, wfp::fftwAllocator<std::complex<realT> >, 0> complexAmplitudeT;
   
   ///The wavefront amplitude
   realImageT amplitude;
   
   ///The wavefront phase
   realImageT phase;

   ///The wavelength at which the wavefront is specified
   realT lambda;
   
   ///The iteration number of this wavefront
   realT iterNo;
   
   ///Zero the wavefront
   void setZero()
   {
      amplitude.setZero();
      phase.setZero();
   }

   void setAmplitude(const realImageT & amp)
   {
      amplitude = amp;
   }
   
   void setPhase(const realImageT & ph)
   {
      phase = ph;
   }
   
   void getWavefront(complexAmplitudeT & wf, int wfSz)
   {
      wfp::makeComplexPupil( wf, amplitude, phase, wfSz); 
   }
   
   void getWavefront( complexAmplitudeT & wf, 
                      realT dlambda, 
                      int wfSz
                    )
   {
      realImageT dphase = phase*(lambda/dlambda);
      wfp::makeComplexPupil( wf, amplitude, dphase, wfSz); 
   }
   
   
};

} //namespace sim
} //namespace AO
} //namespace mx


#endif

