/** \file fraunhoferPropagator.hpp
  * \brief Declares and defines a class for Fraunhofer propagation of optical wavefronts
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef wfp_fraunhoferPropagator_hpp
#define wfp_fraunhoferPropagator_hpp

#include "../math/constants.hpp"
#include "imagingArray.hpp"
#include "imagingUtils.hpp"

#include "../math/fft/fft.hpp"

namespace mx
{

namespace wfp
{

/// Class to perform Fraunhofer propagation between pupil and focal planes
/** This class uses the FFT to propagate between planes, and normalizes so that flux
  * is conserved.  For propagation from pupil to focal plane, the pupil wavefront is tilted so that
  * the focal-plane image is centered at the geometric center of the array.  After propagation from
  * focal plane to pupil plane, the pupil plane wavefront is un-tilted to restore the
  * pupil to its original position.
  *
  * \tparam _wavefrontT is an Eigen::Array-like type, with std::complex values.
  * 
  * \todo check if we should reverse the FFT orders to be proper
  * 
  * \ingroup imaging
  */
template<typename _wavefrontT>
class fraunhoferPropagator
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
   int m_wavefrontSizePixels {0};

   realT m_xcen {0}; ///<x-coordinate of focal plane center, in pixels
   realT m_ycen {0}; ///<x-coordinate of focal plane center, in pixels

   /// Determines how the image is centered.  
   /** If 0 (default) it is at 0.5*(wfSz-1), if true it is shifted by 0.5*m_wholePixel in each axis.
     */
   realT m_wholePixel {0};
   
   ///Phase screen for tilting the pupil plane so that the focal plane image is centered.
   wavefrontT m_centerFocal;

   ///Phase screen for un-tilting the pupil plane after propagating from a centered focal plane.
   wavefrontT m_centerPupil;

   ///FFT object for forward FFTs
   math::fft::fftT<complexT, complexT,2,0> m_fft_fwd;
   
   ///FFT object for backward FFTs
   math::fft::fftT<complexT, complexT,2,0> m_fft_back;

   ///Initialize members
   void initialize();


public:

   ///Constructor
   fraunhoferPropagator();

   ///Destructor
   ~fraunhoferPropagator();

   /// Get the value of the wholePixel parameter
   /** The wholePixel parameter determines how the image is centered.  If 0 (default) it is at 0.5*(wfSz-1), if true it is shifted by 0.5*m_wholePixel in each axis.
     * 
     * \returns the value of m_wholePixel 
     */
   int wholePixel();
   
   /// Set the value of the wholePixel parameter
   /** The wholePixel parameter determines how the image is centered.  If 0 (default) it is at 0.5*(wfSz-1), if true it is shifted by 0.5*m_wholePixel in each axis.
     */
   void wholePixel(realT wp /**< [in] the new wholePixel value */);
   
   ///Apply the shift to a pupil wavefront which will center the resultant focal plane image, and apply the normalization.
   /** You must have allocated the shift screens first, by calling propagatePupilToFocal, propagateFocalToPupil, or setWavefrontSizePixels.
     */
   void shiftPupil( wavefrontT & complexPupil /**< [in/out] the complex pupil plane wavefront to shift*/);

   ///Apply the shift to a pupil wavefront which will restore it to a centered pupil image, with correct flux.
   /** You must have allocated the shift screens first, by calling propagatePupilToFocal, propagateFocalToPupil, or setWavefrontSizePixels.
     */
   void unshiftPupil( wavefrontT & complexPupil /**< [in/out] the complex pupil plane wavefront to shift*/);
   
   ///Propagate the wavefront from the pupil plane to the focal plane
   /** The pupil plane wavefront (complexPupil) is multiplied by a tilt to place the
     * image in the geometric center of the focal plane.  This can be prevented by 
     * setting doCenter to false.
     * 
     */
   void propagatePupilToFocal( wavefrontT & complexFocal, ///< [out] the focal plane wavefront.  Must be pre-allocated to same size as complexPupil.
                               wavefrontT & complexPupil, ///< [in] the pupil plane wavefront. Modified due to application of centering tilt.
                               bool doCenter = true       ///< [in] [optional] set to false to not apply the centering shift
                             );

   ///Propagate the wavefront from Focal plane to Pupil plane
   /** After the fourier transform, the output pupil plane wavefront is de-tilted, restoring it
     * to the state prior to calling \ref propagatePupilToFocal.  This can be prevented by 
     * setting doCenter to false.
     *
     */
   void propagateFocalToPupil( wavefrontT & complexPupil, ///< [out] the pupil plane wavefront. Must be pre-allocated to same size as complexFocal.
                               wavefrontT & complexFocal, ///< [in] the focal plane wavefront.
                               bool doCenter = true       ///< [in] [optional] set to false to not apply the centering shift
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
   m_wavefrontSizePixels = 0;

   m_xcen = 0;
   m_ycen = 0;
}

template<typename wavefrontT>
fraunhoferPropagator<wavefrontT>::fraunhoferPropagator()
{
}

template<typename wavefrontT>
fraunhoferPropagator<wavefrontT>::~fraunhoferPropagator()
{
}

template<typename wavefrontT>
int fraunhoferPropagator<wavefrontT>::wholePixel()
{
   return m_wholePixel;
}
   
template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::wholePixel(realT wp)
{
   m_wholePixel = wp;
   
   //Re-make the shift phase if size already set.
   if(m_wavefrontSizePixels > 0) makeShiftPhase();
}
   

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::shiftPupil( wavefrontT & complexPupil )
{
   complexPupil *= m_centerFocal;
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::unshiftPupil( wavefrontT & complexPupil )
{
   complexPupil *= m_centerPupil;
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::propagatePupilToFocal( wavefrontT & complexFocal,
                                                              wavefrontT & complexPupil,
                                                              bool doCenter /*default = true*/
                                                            )
{
   //First setup the tilt screens (does nothing if there's no change in size)
   setWavefrontSizePixels(complexPupil.rows());

   //Apply the centering shift -- this adjusts by 0.5 pixels and normalizes
   if(doCenter) shiftPupil(complexPupil);
   
   //fft_fwd.fft(complexPupil.data(), complexFocal.data() );
   m_fft_fwd( complexFocal.data(), complexPupil.data() );
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::propagateFocalToPupil( wavefrontT & complexPupil, 
                                                              wavefrontT & complexFocal,
                                                              bool doCenter /*default = true*/
                                                            )
{
   //First setup the tilt screens (does nothing if there's no change in size)
   setWavefrontSizePixels(complexPupil.rows());

   //fft_back.fft( complexFocal.data(), complexPupil.data());
   m_fft_back( complexPupil.data(), complexFocal.data() );

   //Unshift the wavefront and normalize
   if(doCenter) unshiftPupil(complexPupil);
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::setWavefrontSizePixels(int wfsPix)
{
   //If no change in size, do nothing
   if(wfsPix == m_centerFocal.rows()) return;

   m_wavefrontSizePixels = wfsPix;

   m_xcen = 0.5*(wfsPix - 1.0);
   m_ycen = 0.5*(wfsPix - 1.0);

   makeShiftPhase();

   m_fft_fwd.plan(wfsPix, wfsPix);

   m_fft_back.plan(wfsPix, wfsPix, MXFFT_BACKWARD);
}

template<typename wavefrontT>
void fraunhoferPropagator<wavefrontT>::makeShiftPhase()
{
   constexpr realT pi = math::pi<realT>();

   //The normalization is included in the tilt.
   realT norm = 1./(m_wavefrontSizePixels*sqrt(2));
   complexT cnorm = complexT(norm, norm);

   //Resize the center phases
   m_centerFocal.resize(m_wavefrontSizePixels, m_wavefrontSizePixels);
   m_centerPupil.resize(m_wavefrontSizePixels, m_wavefrontSizePixels);

   //Shift by 0.5 pixels
   realT arg = -2.0*pi*0.5*(m_wavefrontSizePixels-m_wholePixel)/(m_wavefrontSizePixels-1);

   for(int ii=0; ii < m_wavefrontSizePixels; ++ii)
   {
      for(int jj=0; jj < m_wavefrontSizePixels; ++jj)
      {
         m_centerFocal(ii,jj) = cnorm*exp(complexT(0.,arg*((ii-m_xcen)+(jj-m_ycen))));
         m_centerPupil(ii,jj) = cnorm*exp(complexT(0., 0.5*pi - arg*((ii-m_xcen)+(jj-m_ycen))));
      }
   }
}

} //namespace wfp
} //namespace mx

#endif //wfp_fraunhoferPropagator_hpp
