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
#include "../imagingUtils.hpp"

#include "../fitsFile.hpp"
#include "../fitsUtils.hpp"

//#include "../eigenImage.hpp"

namespace mx
{
   
namespace imaging 
{
   
/// The Ideal Coronagraph
/**
  * \ingroup coronagraphs
  */   
template<typename _realT>
struct idealCoronagraph
{
   typedef _realT realT;
   typedef std::complex<realT> complexT;
   
   ///The wavefront complex field type
   typedef imagingArray<std::complex<realT>, fftwAllocator<std::complex<realT> >, 0> complexFieldT;
   
   ///The image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

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
   
   /// Propagate the given pupil-plane wavefront through the coronagraph
   void propagate( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);

   /// Propagate the given pupil-plane wavefront without the coronagraph. 
   /** For the ideal coronagraph nothing is done.  This method is included for compliance with
     * with the coronagraph interface.
     */ 
   void propagateNC( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is un-modified. */);
   
   
};

template<typename _realT>
idealCoronagraph<_realT>::idealCoronagraph()
{
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
   //Create Coronagraph pupil.
   padImage(_realPupil, pupil, 0.5*(_wfSz-pupil.rows()),0);
   
   return 0;
}


template<typename _realT>
void idealCoronagraph<_realT>::propagate( complexFieldT & pupilPlane )
{
   Eigen::Map<Eigen::Array<complexT,-1,-1> > eigWf(pupilPlane.data(), pupilPlane.cols(), pupilPlane.rows());
   Eigen::Map<Eigen::Array<realT,-1,-1> > eigPup(_realPupil.data(), _realPupil.cols(), _realPupil.rows());
      
   eigWf = eigWf - ((eigWf*eigPup).sum()/(eigPup*eigPup).sum())*eigPup;
   
   return;
}

template<typename _realT>
void idealCoronagraph<_realT>::propagateNC( complexFieldT & pupilPlane )
{
   return;
}




} //namespace imaging    
} //namespace mx

#endif //__idealCoronagraph_hpp__

