/** \file fraunhoferPropagator.hpp
  * \brief Declares and defines a class for Fraunhofer propagation of optical wavefronts
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */


#ifndef wfp_fraunhoferPropagator_hpp
#define wfp_fraunhoferPropagator_hpp

#include "../imagingArray.hpp"
#include "imagingUtils.hpp"

#include "../fft/fft.hpp"

namespace mx
{
  
namespace wfp 
{
   
template<typename _wavefrontT, int _cudaGPU = 0>
class fraunhoferPropagator; 


/// Class to perform Fraunhofer propagation between pupil and focal planes using the CPU
/** This class uses the FFT to propagate between planes, and normalizes so that flux
  * is conserved.  For propagation from pupil to focal plane, the pupil wavefront is tilted so that 
  * the focal-plane image is centered at the geometric center of the array.  After propagation from 
  * focal plane to pupil plane, the pupil plane wavefront is un-tilted to restore the 
  * pupil to its original position.  
  * 
  * \tparam _wavefrontT is an Eigen::Array-like type, with std::complex values.
  * 
  * \ingroup imaging
  */ 
template<typename _wavefrontT>
class fraunhoferPropagator<_wavefrontT, 0> 
{
   
public:
   
   ///The wavefront data type
   typedef _wavefrontT wavefrontT;
   
   ///The complex data type
   typedef typename wavefrontT::Scalar complexT;
   
   ///The real data type
   typedef typename wavefrontT::Scalar::value_type realT;
   
protected:

   ///The size of the wavefront in pixels
   int wavefrontSizePixels;

   realT xcen; ///<x-coordinate of focal plane center, in pixels
   realT ycen; ///<x-coordinate of focal plane center, in pixels
     
   ///Phase screen for tilting the pupil plane so that the focal plane image is centered.  
   wavefrontT centerFocal;  
   
   ///Phase screen for un-tilting the pupil plane after propagating from a centered focal plane.
   wavefrontT centerPupil;
   
   mx::fftT<complexT, complexT,2,0> fft_fwd;
   mx::fftT<complexT, complexT,2,0> fft_back;
   
   ///Initialize members
   void initialize();
      
   
public:
   
   ///Constructor
   fraunhoferPropagator();
   
   ///Destructor
   ~fraunhoferPropagator();
   
   ///Apply the shift to a pupil wavefront which will center the resultant focal plane image.
   /** You must have allocated the shift screens first, by calling propagatePupilToFocal, propagateFocalToPupil, or setWavefrontSizePixels.
     */
   void shiftPupil( wavefrontT & complexPupil /**< [in/out] the complex pupil plane wavefront to shift*/);
   
   ///Propagate the wavefront from the pupil plane to the focal plane
   /** The pupil plane wavefront (complexPupil) is multiplied by a tilt to place the
     * image in the geometric center of the focal plane. 
     * 
     */ 
   void propagatePupilToFocal( wavefrontT & complexFocal, ///< [out] the focal plane wavefront.  Must be pre-allocated to same size as complexPupil. 
                               wavefrontT & complexPupil  ///< [in] the pupil plane wavefront. Modified due to application of centering tilt.
                             );
   
   ///Propagate the wavefront from Focal plane to Pupil plane
   /** 
     * After the fourier transform, the output pupil plane wavefront is de-tilted, restoring it
     * to the state prior to calling \ref propagatePupilToFocal
     * 
     */ 
   void propagateFocalToPupil( wavefrontT & complexPupil, ///< [out] the pupil plane wavefront. Must be pre-allocated to same size as complexFocal.
                               wavefrontT & complexFocal  ///< [in] the focal plane wavefront.
                             );
   
   ///Set the size of the wavefront, in pixels
   /** Checks if the size changes, does nothing if no change.  Otherwise,  calls
     * \ref makeShiftPhase to pre-calculate the tilt arrays and plans the FFTs.
     *
     */ 
   void setWavefrontSizePixels( int wfsPix /**< [in] the desired new size of the wavefront */ );
   
protected:
   
   ///Calculate the complex tilt arrays for centering and normalizing the wavefronts
   /** 
     */
   void makeShiftPhase();
   
};//class fraunhoferPropagator

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::initialize()
{
   wavefrontSizePixels = 0;
      
   xcen = 0;
   ycen = 0;
}
   
template<typename wavefrontT>
fraunhoferPropagator<wavefrontT>::fraunhoferPropagator()
{
   initialize();
}
   
template<typename wavefrontT>
fraunhoferPropagator<wavefrontT>::~fraunhoferPropagator()
{
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::shiftPupil( wavefrontT & complexPupil )
{
   complexPupil *= centerFocal;
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::propagatePupilToFocal( wavefrontT & complexFocal, 
                                                              wavefrontT & complexPupil  
                                                            )
{
   //First setup the tilt screens (does nothing if there's no change in size)
   setWavefrontSizePixels(complexPupil.rows());
      
   //Apply the centering shift -- this adjusts by 0.5 pixels and normalizes
   complexPupil *= centerFocal;
        
   //fft_fwd.fft(complexPupil.data(), complexFocal.data() );
   fft_fwd( complexFocal.data(), complexPupil.data() );      
}
 
template<typename wavefrontT> 
void fraunhoferPropagator<wavefrontT>::propagateFocalToPupil(wavefrontT & complexPupil, wavefrontT & complexFocal)
{     
   //First setup the tilt screens (does nothing if there's no change in size)
   setWavefrontSizePixels(complexPupil.rows());
   
   //fft_back.fft( complexFocal.data(), complexPupil.data());
   fft_back( complexPupil.data(), complexFocal.data() );      
   
   //Unshift the wavefront and normalize
   complexPupil *= centerPupil;
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::setWavefrontSizePixels(int wfsPix)
{
   //If no change in size, do nothing
   if(wfsPix == centerFocal.rows()) return;
         
   wavefrontSizePixels = wfsPix;
   
   xcen = 0.5*(wfsPix - 1.0);
   ycen = 0.5*(wfsPix - 1.0);
   
   makeShiftPhase();
   
   fft_fwd.plan(wfsPix, wfsPix);
   
   fft_back.plan(wfsPix, wfsPix, MXFFT_BACKWARD);      
}   

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::makeShiftPhase()
{      
   realT pi = boost::math::constants::pi<realT>();
   
   //The normalization is included in the tilt.
   realT norm = 1./(wavefrontSizePixels*sqrt(2));
   complexT cnorm = complexT(norm, norm);
   
   //Resize the center phases
   centerFocal.resize(wavefrontSizePixels, wavefrontSizePixels);
   centerPupil.resize(wavefrontSizePixels, wavefrontSizePixels);

   //Shift by 0.5 pixels
   realT arg = -2.0*pi*0.5*wavefrontSizePixels/(wavefrontSizePixels-1);

   for(int ii=0; ii < wavefrontSizePixels; ++ii)
   {
      for(int jj=0; jj < wavefrontSizePixels; ++jj)
      {     
         centerFocal(ii,jj) = cnorm*exp(complexT(0.,arg*((ii-xcen)+(jj-ycen))));
         centerPupil(ii,jj) = cnorm*exp(complexT(0., 0.5*pi - arg*((ii-xcen)+(jj-ycen))));
      }
   }
}

} //namespace wfp
} //namespace mx

#endif //wfp_fraunhoferPropagator_hpp

