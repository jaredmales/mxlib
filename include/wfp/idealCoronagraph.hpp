/** \file idealCoronagraph.hpp
  * \brief Declares and defines a class to describe the Ideal Coronagraph.
  * \ingroup imaging_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __idealCoronagraph_hpp__
#define __idealCoronagraph_hpp__

#include <iostream>

#include "../imagingArray.hpp"
#include "imagingUtils.hpp"

#include "../improc/fitsFile.hpp"
#include "../improc/fitsUtils.hpp"
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
   
   /// Propagate the given pupil-plane wavefront through the coronagraph
   int propagate( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);

   /// Propagate the given pupil-plane wavefront without the coronagraph. 
   /** For the ideal coronagraph nothing is done.  This method is included for compliance with
     * with the coronagraph interface.
     */ 
   int propagateNC( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is un-modified. */);
   
   
};

template<typename _realT>
idealCoronagraph<_realT>::idealCoronagraph()
{
   _wfSz = 0;
}

template<typename _realT>
int idealCoronagraph<_realT>::wfSz()
{
   return _wfSz;
}

template<typename _realT>
void idealCoronagraph<_realT>::wfSz(int sz)
{
   _wfSz = sz;   
}

template<typename _realT>
int idealCoronagraph<_realT>::setPupil( imageT & pupil)
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

template<typename _realT>
int idealCoronagraph<_realT>::loadPupil( const std::string & pupilFile)
{

   imageT pupil;
   
   improc::fitsFile<realT> ff;
   
   ff.read(pupil, pupilFile);
   
   return setPupil(pupil);
   
   return 0;
}

template<typename _realT>
int idealCoronagraph<_realT>::loadCoronagraph( const std::string & cName)
{

   if( _fileDir == "")
   {
      mxError("idealCoronagraph", MXE_PARAMNOTSET, "file directory (fileDir) not set.");
      return -1;
   }
   
   std::string pupilFile = _fileDir + "/" + cName + ".fits";
   
   return loadPupil(pupilFile);
   
}

template<typename _realT>
int idealCoronagraph<_realT>::propagate( complexFieldT & pupilPlane )
{
   if( pupilPlane.rows() != _realPupil.rows() || pupilPlane.cols() != _realPupil.cols())
   {
      mxError("idealCoronagraph", MXE_SIZEERR, "pupilPlane wavefront size does not match realPupil");
      return -1;
   }
   
   Eigen::Map<Eigen::Array<complexT,-1,-1> > eigWf(pupilPlane.data(), pupilPlane.cols(), pupilPlane.rows());
   Eigen::Map<Eigen::Array<realT,-1,-1> > eigPup(_realPupil.data(), _realPupil.cols(), _realPupil.rows());
      
   eigWf = eigWf - ((eigWf*eigPup).sum()/(eigPup*eigPup).sum())*eigPup;
   
   return 0;
}

template<typename _realT>
int idealCoronagraph<_realT>::propagateNC( complexFieldT & pupilPlane )
{
   return 0;
}




} //namespace wfp 
} //namespace mx

#endif //__idealCoronagraph_hpp__

