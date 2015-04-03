#ifndef __fourierImager_hpp__
#define __foutierImager_hpp__

namespace mx
{

#include "imagingUtils.hpp"
   
/// Class to perform Fraunhofer propagation between pupil and focal planes
/** This class uses the FFT to propagate between planes, and normalizes so that flux
  * is conserved.  The wavefronts are tilted so that the focal-plane image is centered
  * at the geometric center of the array.
  * 
  * \tparam arithT is the arithmetic type
  * 
  */ 
template<typename arithT>
class fourierImager 
{
   
public:
   ///The complex data type
   typedef std::complex<arithT> complexT;
   
   ///The wavefront data type
   typedef Eigen::Array<std::complex<arithT>, Eigen::Dynamic, Eigen::Dynamic> wavefrontT;
   
protected:

   ///The size of the wavefront in pixels
   int wavefrontSizePixels;

   ///The physical size of the wavefront
   /** This is used for calculating the plate scale
     */
   arithT wavefrontSizePhysical;

   arithT xcen; ///<x-coordinate of focal plane center, in pixels
   arithT ycen; ///<x-coordinate of focal plane center, in pixels
     
   ///Phase screen for shifting the Focal plane by 0.5 pixels  
   wavefrontT centerFocal;  
   
   ///Phase screen for shifting the Pupil planet by -0.5 pixels
   wavefrontT centerPupil;

   ///Initialize members
   void initialize()
   {
      wavefrontSizePixels = 0;
      wavefrontSizePhysical = 0;
      
      xcen = 0;
      ycen = 0;
   }
   
public:
   
   ///c'tor
   fourierImager()
   {
      initialize();
   }

   void reset()
   {
      centerFocal.resize(0,0);
   }
   
   ///Propagate the wavefront from the pupil plane to the focal plane
   /** The pupil plane wavefront (complexPupil) is multiplied by the centerFocal
     * tilt to shift by 0.5 pixels.  The image is not shifted to the center of the frame.
     * 
     * Note that complexPupil is modified.
     * 
     */ 
   void propagatePupilToFocal(wavefrontT & complexFocal, wavefrontT & complexPupil)
   {
      //First setup the tilt screens (does nothing if there's not change in size)
      setWavefrontSizePixels(complexPupil.rows());
      
      //DFT normalization, sqrt(2) for complex number
      arithT norm = wavefrontSizePixels/sqrt(2.);
      
      //Apply the centering shift -- this adjusts by 0.5 pixels
      complexPupil *= centerFocal;
            
      fft(complexFocal, complexPupil);
      
      //Normalize
      complexFocal = complexFocal / complexT(norm, norm);
      
   }
   
   ///Propagate the wavefront from Focal plane to Pupil plane
   /** After the fourier transform, the output Pupil plane wavefront is multiplied by
     * the centerPupil tilt to unshift by 0.5 pixels.  This should be called
     * after first propagating to the Focal plane from the initial Pupil plane.
     * 
     */ 
   void propagateFocalToPupil(wavefrontT & complexPupil, wavefrontT & complexFocal)
   {
      //DFT normalization, sqrt(2) for complex number
      arithT norm = wavefrontSizePixels/sqrt(2.);
      
      fft(complexPupil, complexFocal, FFTW_BACKWARD);
      
      //Unshift the wavefront
      complexPupil *= centerPupil;
      
      complexPupil /= complexT(norm, norm);
   }
   
   /// Set the physical size of the wavefront
   /**
     * \param wfsPhys is the new size of the wavefront
     */ 
   void setWavefrontSizePhysical(int wfsPhys)
   {      
      wavefrontSizePhysical = wfsPhys;
   }

   void setWavefrontSizePixels(int wfsPix)
   {
      //If no change in size, do nothing
      if(wfsPix == centerFocal.rows()) return;
            
      wavefrontSizePixels = wfsPix;
      
      xcen = 0.5*(wfsPix - 1.0);
      ycen = 0.5*(wfsPix - 1.0);
      
      makeShiftPhase();
   }   
   
protected:
   
   void makeShiftPhase()
   {      
      //Resize the center phases
      centerFocal.resize(wavefrontSizePixels, wavefrontSizePixels);
      centerPupil.resize(wavefrontSizePixels, wavefrontSizePixels);

      //Shift by 0.5 pixels
      arithT arg = -1.0*D2PI*0.5*wavefrontSizePixels/(wavefrontSizePixels-1);
  
      for(int ii=0; ii < wavefrontSizePixels; ++ii)
      {
         for(int jj=0; jj < wavefrontSizePixels; ++jj)
         {     
            centerFocal(ii,jj) = exp(complexT(0.,arg*((ii-xcen)+(jj-ycen))));
            centerPupil(ii,jj) = exp(complexT(0., 0.5*DPI - arg*((ii-xcen)+(jj-ycen))));
         }
      }
   }

   
};//class fourierImager

} //namespace mx

#endif //__fourierImager_hpp__

