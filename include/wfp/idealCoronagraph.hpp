/** \file idealCoronagraph.hpp
  * \brief Declares and defines a class to describe the Ideal Coronagraph.
  * \ingroup imaging_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __idealCoronagraph_hpp__
#define __idealCoronagraph_hpp__

#include <iostream>

#include "fraunhoferPropagator.hpp"
#include "imagingUtils.hpp"

#include "../wfp/imagingArray.hpp"

#include "../ioutils/fits/fitsFile.hpp"
#include "../ioutils/fits/fitsUtils.hpp"
#include "../improc/imagePads.hpp"

namespace mx
{

namespace wfp
{

/// The Ideal Coronagraph
/** A simple toy coronagraph which operates only the pupil plane, subtracting the energy-minimizing wavefront \cite cavarroc_2006.  For an on-axis source with
  * a perfectly flat wavefront this results in perfect extinction.  This coronagraph is not physically realizable, but it is often useful
  * for modeling and analysis, and it has been shown that real-world optimized coronagraphs often approach ideal performance \cite sauvage_2010 [and Males in prep].
  *
  * \ingroup coronagraphs
  */
template<typename _realT>
struct idealCoronagraph
{
   ///The real floating point type
   typedef _realT realT;

   ///The complex floating point type
   typedef std::complex<realT> complexT;

   ///The wavefront complex field type
   typedef imagingArray<std::complex<realT>, fftwAllocator<std::complex<realT> >, 0> complexFieldT;

   ///The image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

   ///The directory where coronagraph files are stored
   std::string _fileDir;

   ///The linear size of the wavefront in pixels
   int _wfSz;

   imageT _realPupil;

   complexFieldT m_focalPlane;

   ///Fraunhofer propagator
   fraunhoferPropagator<complexFieldT> m_fi;

   /// Determines how the image is centered.  
   /** If 0 (default) it is at 0.5*(wfSz-1), if true it is shifted by 0.5*m_wholePixel in each axis.  This is passed to 
     * fraunhoferPropagator when it is resized.
     */ 
   realT m_wholePixel {0};
   
   idealCoronagraph();

   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */
   int wfSz();

   ///Set the wavefront size in pixels.
   /**
     */
   void wfSz(int sz /**< [in] is the new size */);

   /// Set the real pupil mask.
   /** The input mask does not have to be the same size as wfSz, as the stored mask will be padded.
     *
     * \returns 0 on success.
     * \returns -1 on error.
     */
   int setPupil( imageT & pupil /**< [in] the pupil mask. */);

   /// Load the real pupil mask from a FITS file.
   /** The input mask does not have to be the same size as wfSz, as the stored mask will be padded.
     *
     * \returns 0 on success.
     * \returns -1 on error.
     */
   int loadPupil( const std::string & pupilFile /**< [in] the FITS file containing the pupil mask. */);

   ///Load the components of the coronagraph (just a pupil) based in its base name
   /** Looks in _filDir for the files.
     *
     * \returns 0 on success.
     * \returns -1 on error.
     */
   int loadCoronagraph( const std::string & cName /**< The name of the coronagraph, without directory or file extensions */);

   /// Propagate the given pupil-plane wavefront through the coronagraph to the exit pupil plane
   int propagate( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);


   /// Propagate the given pupil-plane wavefront through the coronagraph to the exit pupil plane, and then to the final focal plane.
   int propagate( imageT & fpIntensity,      ///< [out] The intensity image in the focal plane.  This should be pre-allocated.
                  complexFieldT & pupilPlane ///< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph.
                );

   /// Propagate the given pupil-plane wavefront without the coronagraph.
   /** For the ideal coronagraph nothing is done.  This method is included for compliance with
     * with the coronagraph interface.
     */
   int propagateNC( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is un-modified. */);

   /// Propagate the given pupil-plane wavefront without the coronagraph to the exit pupil plane, and then to the final focal plane.
   /** For the ideal coronagraph nothing is done to the input wavefront.
     */
   int propagateNC( imageT & fpIntensity,      ///< [out] The intensity image in the focal plane.  This should be pre-allocated.
                    complexFieldT & pupilPlane ///< [in/out] The wavefront at the input pupil plane.  It is un-modified.
                  );
   
   bool apodize {false};
   imageT apodizer;

};

template<typename realT>
idealCoronagraph<realT>::idealCoronagraph()
{
   _wfSz = 0;
}

template<typename realT>
int idealCoronagraph<realT>::wfSz()
{
   return _wfSz;
}

template<typename realT>
void idealCoronagraph<realT>::wfSz(int sz)
{
   _wfSz = sz;

   m_fi.setWavefrontSizePixels(sz);
   m_fi.wholePixel(m_wholePixel);
   m_focalPlane.resize(sz, sz);

}

template<typename realT>
int idealCoronagraph<realT>::setPupil( imageT & pupil)
{
   if(_wfSz <= 0)
   {
      mxError("idealCoronagraph", MXE_PARAMNOTSET, "Must set wavefront size (wfSz) before setting up coronagraph.");
      return -1;
   }

   //Create Coronagraph pupil.
   improc::padImage(_realPupil, pupil, 0.5*(_wfSz-pupil.rows()),0);

   return 0;
}

template<typename realT>
int idealCoronagraph<realT>::loadPupil( const std::string & pupilFile)
{

   imageT pupil;

   fits::fitsFile<realT> ff;

   ff.read(pupil, pupilFile);

   return setPupil(pupil);

   return 0;
}

template<typename realT>
int idealCoronagraph<realT>::loadCoronagraph( const std::string & cName)
{

   if( _fileDir == "")
   {
      mxError("idealCoronagraph", MXE_PARAMNOTSET, "file directory (fileDir) not set.");
      return -1;
   }

   std::string pupilFile = _fileDir + "/" + cName + ".fits";

   return loadPupil(pupilFile);

}

template<typename realT>
int idealCoronagraph<realT>::propagate( complexFieldT & pupilPlane )
{
   if( pupilPlane.rows() != _realPupil.rows() || pupilPlane.cols() != _realPupil.cols())
   {
      mxError("idealCoronagraph", MXE_SIZEERR, "pupilPlane wavefront size does not match realPupil");
      return -1;
   }

   Eigen::Map<Eigen::Array<complexT,-1,-1> > eigWf(pupilPlane.data(), pupilPlane.cols(), pupilPlane.rows());
   Eigen::Map<Eigen::Array<realT,-1,-1> > eigPup(_realPupil.data(), _realPupil.cols(), _realPupil.rows());

   eigWf = (eigWf - ((eigWf*eigPup).sum()/(eigPup*eigPup).sum()))*eigPup;
   
   if(apodize)
   {
      eigWf *= apodizer;
   }

   return 0;
}

template<typename realT>
int idealCoronagraph<realT>::propagate( imageT & fpIntensity,
                                        complexFieldT & pupilPlane
                                      )
{
   propagate(pupilPlane);

   m_fi.propagatePupilToFocal(m_focalPlane, pupilPlane);

   int x0 = 0.5*(_wfSz-1) - 0.5*(fpIntensity.rows()-1);
   int y0 = 0.5*(_wfSz-1) - 0.5*(fpIntensity.cols()-1);

   extractIntensityImage(fpIntensity,0, fpIntensity.rows(),0,fpIntensity.cols(), m_focalPlane, x0,y0);

   return 0;
}

template<typename realT>
int idealCoronagraph<realT>::propagateNC( complexFieldT & pupilPlane )
{
   static_cast<void>(pupilPlane);
   return 0;
}

template<typename realT>
int idealCoronagraph<realT>::propagateNC( imageT & fpIntensity,
                                          complexFieldT & pupilPlane )
{
   propagateNC(pupilPlane);

   m_fi.propagatePupilToFocal(m_focalPlane, pupilPlane);

   int x0 = 0.5*(_wfSz-1) - 0.5*(fpIntensity.rows()-1);
   int y0 = 0.5*(_wfSz-1) - 0.5*(fpIntensity.cols()-1);

   extractIntensityImage(fpIntensity,0, fpIntensity.rows(),0,fpIntensity.cols(), m_focalPlane, x0,y0);

   return 0;
}


} //namespace wfp
} //namespace mx

#endif //__idealCoronagraph_hpp__
