/** \file imagingUtils.hpp
  * \brief Utilities for modeling image formation
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __imagingUtils_hpp__
#define __imagingUtils_hpp__

#include "imagingArray.hpp"

#include <cmath>
#include <boost/math/constants/constants.hpp>

namespace mx
{

/// Fill in an imagingArray with a circular pupil mask.
/** \ingroup imaging
  * \param m is the allocated Array
  * \param eps [optional] is the central obscuration.  0-1, default is 0.
  * \param rad [optional] is the desired radius. if 0 the maximum radius is used.
  */  
template<class arrayT> 
void circularPupil( arrayT & m, 
                     typename arrayT::dataT eps=0, 
                     typename arrayT::dataT rad=0 
                   )
{
   size_t l0 = m.szY();
   size_t l1 = m.szX();
   
   typename arrayT::dataT r;
   typename arrayT::dataT xc = 0.5*(l0-1);
   typename arrayT::dataT yc = 0.5*(l1-1);
   
   if(rad == 0) rad = 0.5*std::min(l0-1, l1-1);
   
   for(size_t i=0; i < l0; i++)
   {
      for(size_t j=0; j < l1; j++)
      {
         r = std::sqrt( std::pow(i-xc, 2) + std::pow(j-yc, 2) );
         
         if(r <= rad+0.5 && r >= eps*rad) m(i,j) = 1;
         else m(i,j) = 0;
      }
   }
}


///Create a complex pupil plane wavefront from a real amplitude mask.
/** The input real amplitude mask is placed in the center of a 0-padded complex array.
  *
  * \param [out] complexPupil the complex pupil plane wavefront
  * \param [in] realPupil a real amplitude mask.
  * \param [in] wavefrontSizePixels the desired size of the ouput wavefront, should be at least as big as realPupil
  * 
  * \ingroup imaging
  */ 
template<typename arrayOutT, typename arrayInT>
void makeComplexPupil( arrayOutT & complexPupil, 
                       const arrayInT & realPupil, 
                       int wavefrontSizePixels)
{
   
   complexPupil.resize(wavefrontSizePixels, wavefrontSizePixels);
   complexPupil.set(0);
     
   //Lower-left corner of insertion region
   int bl = 0.5*(complexPupil.szY()-1) - 0.5*(realPupil.szY()-1.);
   
   for(int i=0; i< realPupil.szX(); ++i)
   {
      for(int j=0; j < realPupil.szY(); ++j)
      {
         complexPupil(bl+i, bl+j) = realPupil(i,j)*( typename arrayOutT::dataT(1,0)); 
      }
   }
   //complexPupil.block(bl, bl, realPupil.szY(), realPupil.szY()) = realPupil*std::complex<realT>(1,0);

}


///Create a complex wavefront from a real amplitude and a real phase.
/** The wavefront is placed in the center of a 0-padded complex array.
  *
  * \param [out] complexWavefront the complex pupil plane wavefront
  * \param [in] realAmplitude is the real-valued amplitude.
  * \param [in] realPhase is the real-valued phase in radians, same size as realAmplitude
  * \param [in] wavefrontSizePixels the desired size of the ouput wavefront, should be at least as big as the real arrays
  * 
  * \ingroup imaging
  */ 
template<typename arrayOutT, typename arrayInT>
void makeComplexPupil( arrayOutT & complexWavefront, 
                       const arrayInT & realAmplitude, 
                       const arrayInT & realPhase, 
                       int wavefrontSizePixels)
{
   
   complexWavefront.resize(wavefrontSizePixels, wavefrontSizePixels);
   complexWavefront.set(0);
     
   //Lower-left corner of insertion region
   int bl = 0.5*(complexWavefront.szY()-1) - 0.5*(realAmplitude.szY()-1.);
   
   for(int i=0; i< realAmplitude.szX(); ++i)
   {
      for(int j=0; j < realAmplitude.szY(); ++j)
      {
         complexWavefront(bl+i, bl+j) = realAmplitude(i,j)*exp(  (typename arrayOutT::dataT(0,1)) * realPhase(i,j)); 
      }
   }
}


///Apply a tilt to a wavefront 
/**
  * \param complexWavefront [in/out] the complex wavefront to tilt, will be modified on output
  * \param xTilt [input] the amount of tilt in the x direction, in pixels 
  * \param yTilt [input] the amount of tilt in the y direction, in pixels 
  * 
  * \ingroup imaging
  */
template<typename wavefrontT>
void tiltWavefront( wavefrontT & complexWavefront, 
                    typename wavefrontT::dataT::value_type xTilt, 
                    typename wavefrontT::dataT::value_type yTilt)
{
   typedef typename wavefrontT::dataT complexT;
   typedef typename wavefrontT::dataT::value_type realT;
   
   realT pi = boost::math::constants::pi<realT>();
   
   int wfsSizeX = complexWavefront.szX();
   int wfsSizeY = complexWavefront.szY();
   
   realT xCen = 0.5*wfsSizeX;
   realT yCen = 0.5*wfsSizeY;

   realT argX = 2*pi/(wfsSizeX-1);
   realT argY = 2*pi/(wfsSizeY-1);
  
   for(int ii=0; ii < wfsSizeX; ++ii)
   {
      for(int jj=0; jj < wfsSizeY; ++jj)
      {     
         complexWavefront(ii,jj) = complexWavefront(ii,jj)*exp( complexT( (realT)0., argX*xTilt*(ii-xCen)+argY*yTilt*(jj-yCen)));
      }
   }
}

template< typename imageT>
void extractBlock(imageT & im,
                  int imX0,
                  int imXsz,
                  int imY0,
                  int imYsz,
                  imageT & wf,
                  int wfX0,
                  int wfY0)
{
   int im_idx;
   int im_rows = im.szX();
   
   int wf_idx;
   int wf_rows = wf.szX();
   
   typedef typename imageT::dataT dataT;
   
   dataT * im_data;
   dataT * wf_data;
   
   

   for(int j =0; j< imYsz; ++j)
   {
      im_data = &im.data()[imX0 + (imY0+j)*im_rows];
      wf_data = &wf.data()[wfX0 + (wfY0+j)*wf_rows];
      
      memcpy( im_data, wf_data, sizeof(dataT)*imXsz);
   }
}

template< typename realImageT,
          typename complexImageT>
void extractIntensityImage(realImageT & im,
                           int imX0,
                           int imXsz,
                           int imY0,
                           int imYsz,
                           complexImageT & wf,
                           int wfX0,
                           int wfY0)
{
   int im_idx;
   int im_rows = im.szX();
   
   int wf_idx;
   int wf_rows = wf.szX();
   
   typename realImageT::dataT * im_data;
   typename complexImageT::dataT * wf_data;
   
   for(int j =0; j< imXsz; ++j)
   {
      im_data = &im.data()[imX0 + (imY0+j)*im_rows];
      wf_data = &wf.data()[wfX0 + (wfY0+j)*wf_rows];
      for(int i=0; i<imYsz; ++i)
      {
         im_data[i] += norm(wf_data[i]);
      }
   }
}

} //namespace mx

#endif //__imagingUtils_hpp__

