/** \file fraunhoferPropagatorCuda.hpp
  * \brief Declares and defines a class for Fraunhofer propagation of optical wavefronts with Cuda
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */


#ifndef wfp_fraunhoferPropagatorCuda_hpp
#define wfp_fraunhoferPropagatorCuda_hpp

#include "../math/constants.hpp"

#include "fraunhoferPropagator.hpp"
#include "../cuda/templateCuda.hpp"
#include "../cuda/templateCudaPtr.hpp"

namespace mx
{
namespace wfp 
{

/// Class to perform Fraunhofer propagation between pupil and focal planes using a GPU
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
class fraunhoferPropagator<_wavefrontT,1>
{
public:

   ///The wavefront data type
   typedef _wavefrontT wavefrontT;
   
   ///The complex data type
   typedef typename wavefrontT::Scalar complexT;
   
   ///The real data type
   typedef typename wavefrontT::Scalar::value_type realT;
   
   typedef complexT devicePtrT;
   
protected:

   ///The size of the wavefront in pixels
   int wavefrontSizePixels {0};

   realT xcen {0}; ///<x-coordinate of focal plane center, in pixels
   realT ycen {0}; ///<x-coordinate of focal plane center, in pixels
     
   ///Phase screen for tilting the pupil plane so that the focal plane image is centered [GPU memory].
   //complexT * m_centerFocal {nullptr};
   mx::cuda::cudaPtr<complexT> m_centerFocal;
   
   ///Phase screen for un-tilting the pupil plane after propagating from a centered focal plane [GPU memory].
   //complexT * m_centerPupil {nullptr};
   mx::cuda::cudaPtr<complexT> m_centerPupil;

   ///Cuda FFT plan.  We only need one since the forward/inverse is part of execution.
   cufftHandle m_fftPlan {0};
   
public:
   ///Constructor
   fraunhoferPropagator();
   
   ///Destructor
   ~fraunhoferPropagator();
   
   ///Set the size of the wavefront, in pixels
   /** Checks if the size changes, does nothing if no change.  Otherwise,  calls
     * \ref makeShiftPhase to pre-calculate the tilt arrays and plans the FFTs.
     *
     */ 
   void setWavefrontSizePixels( int wfsPix /**< [in] the desired new size of the wavefront */ );
   
   ///Propagate the wavefront from the pupil plane to the focal plane
   /** The pupil plane wavefront (complexPupil) is multiplied by a tilt to place the
     * image in the geometric center of the focal plane. 
     * 
     */ 
   void propagatePupilToFocal( devicePtrT * complexFocal, ///< [out] the focal plane wavefront.  Must be pre-allocated to same size as complexPupil. 
                               devicePtrT * complexPupil  ///< [in] the pupil plane wavefront. Modified due to application of centering tilt.
                             );
   
   ///Propagate the wavefront from Focal plane to Pupil plane
   /** 
     * After the fourier transform, the output pupil plane wavefront is de-tilted, restoring it
     * to the state prior to calling \ref propagatePupilToFocal
     * 
     */ 
   void propagateFocalToPupil( devicePtrT * complexPupil, ///< [out] the pupil plane wavefront. Must be pre-allocated to same size as complexFocal.
                               devicePtrT * complexFocal  ///< [in] the focal plane wavefront.
                             );
   
   
protected:
   
   ///Calculate the complex tilt arrays for centering and normalizing the wavefronts
   /** 
     */
   void makeShiftPhase();
   
};
   
template<typename wavefrontT>
fraunhoferPropagator<wavefrontT, 1>::fraunhoferPropagator()
{
}

template<typename wavefrontT>
fraunhoferPropagator<wavefrontT, 1>::~fraunhoferPropagator()
{
   if( m_fftPlan )
   {
      checkCudaErrors(cufftDestroy(m_fftPlan));
   }
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT,1>::setWavefrontSizePixels(int wfsPix)
{
   //If no change in size, do nothing
   if(wfsPix == wavefrontSizePixels) return;
         
   wavefrontSizePixels = wfsPix;
   
   xcen = 0.5*(wfsPix - 1.0);
   ycen = 0.5*(wfsPix - 1.0);
   
   makeShiftPhase();
   
   //Plan the FFT 
   checkCudaErrors(cufftPlan2d(&m_fftPlan, wavefrontSizePixels, wavefrontSizePixels, CUFFT_C2C));
   
}   

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT,1>::propagatePupilToFocal( devicePtrT * complexFocal, 
                                                                devicePtrT * complexPupil  
                                                              )
{
   //Apply the centering shift -- this adjusts by 0.5 pixels and normalizes
   mx::cuda::pointwiseMul<cuComplex><<<32, 256>>>( (cuComplex *) complexPupil, (cuComplex *) m_centerFocal.m_devicePtr, wavefrontSizePixels*wavefrontSizePixels);
      
   cufftExecC2C(m_fftPlan, (cufftComplex *) complexPupil, (cufftComplex *) complexFocal, CUFFT_FORWARD);     

}
 
template<typename wavefrontT> 
void fraunhoferPropagator<wavefrontT, 1>::propagateFocalToPupil( devicePtrT * complexPupil, 
                                                                 devicePtrT * complexFocal
                                                               )
{     
   cufftExecC2C(m_fftPlan, (cufftComplex *) complexFocal, (cufftComplex *) complexPupil, CUFFT_INVERSE);
   
   //Unshift the wavefront and normalize
   mx::cuda::pointwiseMul<cuComplex><<<32, 256>>>( (cuComplex*) complexPupil, (cuComplex*) m_centerPupil.m_devicePtr, wavefrontSizePixels*wavefrontSizePixels);   
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT, 1>::makeShiftPhase()
{      
   constexpr realT pi = math::pi<realT>();
   
   //The normalization is included in the tilt.
   realT norm = 1./(wavefrontSizePixels*sqrt(2));
   complexT cnorm = complexT(norm, norm);
   
   ///Host memory to build the shift screens
   complexT * centerFocal = new complexT[wavefrontSizePixels*wavefrontSizePixels];
   
   complexT * centerPupil = new complexT[wavefrontSizePixels*wavefrontSizePixels];
   
   //Shift by 0.5 pixels
   realT arg = -2.0*pi*0.5*wavefrontSizePixels/(wavefrontSizePixels-1);

   for(int ii=0; ii < wavefrontSizePixels; ++ii)
   {
      for(int jj=0; jj < wavefrontSizePixels; ++jj)
      {
         centerFocal[ii*wavefrontSizePixels + jj] = cnorm*exp(complexT(0.,arg*((ii-xcen)+(jj-ycen))));
         centerPupil[ii*wavefrontSizePixels + jj] = cnorm*exp(complexT(0., 0.5*pi - arg*((ii-xcen)+(jj-ycen))));
      }
   }

   m_centerFocal.upload(centerFocal, wavefrontSizePixels*wavefrontSizePixels);
   
   m_centerPupil.upload(centerPupil, wavefrontSizePixels*wavefrontSizePixels);
   
   delete[] centerFocal;
   delete[] centerPupil;
   
}
   

} //namespace wfp
} //namespace mx

#endif //wfp_fraunhoferPropagatorCuda_hpp

