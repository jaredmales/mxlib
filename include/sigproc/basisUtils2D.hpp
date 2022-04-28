/** \file basisUtils2D.hpp
  * \brief Utilities for a working with a 2D basis set
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef basisUtils_hpp
#define basisUtils_hpp

#include "../improc/eigenImage.hpp"
#include "../improc/eigenCube.hpp"

namespace mx
{
namespace sigproc 
{

/// Mask a basis set.
/** Multiplies each mode in a basis set by a mask.
  *
  * \returns 0 on success
  * \returns -1 on error
  *
  * \tparam realT a real floating point type.
  * 
  * \ingroup signal_processing
  */ 
template<typename realT>
int basisMask( improc::eigenCube<realT> & modes, ///< [in/out] the basis to normalize.
               improc::eigenImage<realT> & mask  ///< [in] 1/0 mask defining the domain of the basis
             )
{
   for(int i=0; i < modes.planes(); ++i)
   {
      modes.image(i) *= mask;
   }  
   
   return 0;
}

/// Mean-subtract a basis set.
/** Subtracts the mean of each mode calculated over a domain defined by a mask.
  *
  * \returns 0 on success
  * \returns -1 on error
  *
  * \tparam realT a real floating point type.
  * 
  * \ingroup signal_processing
  */
template<typename realT>
int basisMeanSub( improc::eigenCube<realT> & modes, ///< [in/out] the basis to normalize.
                  improc::eigenImage<realT> & mask, ///< [in] 1/0 mask defining the domain of the basis
                  bool postMult = true      ///< [in] [optional] if true, then each image is multiplied by the mask after subtraction.
                )
{
   realT maskSum = mask.sum();
   for(int i=0; i < modes.planes(); ++i)
   {
      float mean = (modes.image(i)*mask).sum()/maskSum;
      
      modes.image(i) -= mean;
      
      if(postMult) modes.image(i) *= mask;
   }  
   
   return 0;
}

/// Normalize a basis set.
/** RMS normalize each mode over a domain defined by a mask.
  *
  * \returns 0 on success
  * \returns -1 on error
  *
  * \tparam realT a real floating point type.
  * 
  * \ingroup signal_processing
  */
template<typename realT>
int basisNormalize( improc::eigenCube<realT> & modes, ///< [in/out] the basis to normalize.
                    improc::eigenImage<realT> & mask  ///< [in] 1/0 mask defining the domain of the normalization
                  )
{
   realT psum = mask.sum();
   for(int i=0; i < modes.planes(); ++i)
   {
      float norm = (modes.image(i)*mask).square().sum()/psum;
      
      modes.image(i)/=sqrt(norm);
   }  
   
   return 0;
}

/// Measure the amplitudes of a set of basis modes fit to an image.  Optionally subtract them.
/** Mode subtraction occurs one by one, so subtraction will work with non-orthogonal basis sets.
  * 
  * \returns 0 on success
  * \returns -1 on error
  * 
  * \tparam realT the floating point type.
  * 
  * \ingroup signal_processing
  */ 
template<typename realT>
int basisAmplitudes( std::vector<realT> & amps,        ///< [out] the amplitudes of each mode fit to the image (will be resized).
                     improc::eigenImage<realT> & im,   ///< [in/out] the image to fit.  Is subtracted in place if desired.
                     improc::eigenCube<realT> & modes, ///< [in] the modes to fit.
                     improc::eigenImage<realT> & mask, ///< [in] the 1/0 mask which defines the domain of the fit.
                     bool subtract = false,            ///< [in] [optional] if true then the modes are subtracted as they are fit to the image
                     int meanIgnore= 0,                ///< [in] [optional] if 1 then the mean, or if 2 the median, value is subtracted before fitting. If subtract  is false, this value is added back after the subtraction.
                     int N = -1                        ///< [in] [optional] the number of modes to actually fit.  If N < 0 then all modes are fit.
                   )
{
   if( N < 0) N = modes.planes();
   amps.resize(N);

   realT apertureNPix = mask.sum();

   realT mean;
   if(meanIgnore)
   {
      if(meanIgnore == 2) mean = improc::imageMedian(im,&mask);
      else mean = (im*mask).sum()/apertureNPix;

      im -= mean;
   }

   for(int i=0; i<N; ++i)
   {
      amps[i] = (im*modes.image(i)*mask).sum()/ apertureNPix;
      
      if(subtract)
      {
         im -= amps[i]*modes.image(i);
      }
   }

   if(meanIgnore && !subtract)
   {
      im += mean;
   }

   return 0;
}   


} //namespace sigproc 
} //namespace mx

#endif //basisUtils_hpp


