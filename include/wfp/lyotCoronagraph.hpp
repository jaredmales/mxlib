/** \file lyotCoronagraph.hpp
  * \brief Declares and defines a class to describe and optimize a Lyot Coronagraph.
  * \ingroup imaging_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __lyotCoronagraph_hpp__
#define __lyotCoronagraph_hpp__

#include <iostream>

#include "imagingArray.hpp"
#include "imagingUtils.hpp"
#include "fraunhoferPropagator.hpp"

#include "../ioutils/fits/fitsFile.hpp"
#include "../ioutils/fits/fitsUtils.hpp"

#include "../improc/eigenImage.hpp"
#include "../improc/eigenCube.hpp"

namespace mx
{
   
namespace wfp
{
   
/// The Lyot Coronagraph
/** A generalized Lyot coronagraph, which can include a pupil apodization, a focal plane mask with complex transmission, and a 
  * Lyot stop.  Light is propagated through the coronagraph plane-by-plane, and the complex wavefront can be accessed at any plane.
  * Also provides functions for optimizing the pupil apodizer and focal plane mask.
  * 
  * \ingroup coronagraphs
  */   
template<typename _realT, typename _fpmaskFloatT>
struct lyotCoronagraph
{
public:
   ///The real floating point type
   typedef _realT realT;
   
   ///The real floating point type for mask calculations.
   typedef _fpmaskFloatT fpmaskFloatT;
   
   ///The complex floating point type
   typedef std::complex<realT> complexT;

   ///The wavefront complex field type
   typedef mx::wfp::imagingArray<std::complex<realT>, fftwAllocator<std::complex<realT> >, 0> complexFieldT;
   
   ///The focal plane mask type
   typedef Eigen::Array< std::complex<fpmaskFloatT>, Eigen::Dynamic, Eigen::Dynamic> fpMaskT;
   
   ///The image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

   
public:   
   
   ///The directory where coronagraph files are stored
   std::string m_fileDir {"coron"};
   
   ///The linear size of the wavefront in pixels
   int m_wfSz {0};
   
   ///Image containing the pupil apodization.
   imageT m_pupilApodizer;
   
   
   int m_maskSource; ///< 0= read from file, 1 = constructed by makeFocalMask, 2 = trans optimized.
   std::string m_maskFile; ///< Name of file from which mask was loaded.
   
   realT m_maskRad; ///<Radius of mask if it was constructed.
   realT m_maskTrans; ///<Transmission of mask if it was constructed.
   
   ///The focal plane  mask.
   fpMaskT m_focalMask;
   
   ///Image containing the lyot stop.
   imageT m_lyotStop;

   
   complexFieldT m_focalPlane;

   bool m_savePreMaskFocalPlane {false};
   complexFieldT m_preMaskFocalPlane;
   
   bool m_savePreLyotPupilPlane {false};
   complexFieldT m_preLyotPupilPlane;
   
   ///Fraunhofer propagator
   fraunhoferPropagator<complexFieldT> m_fp;

public:
   
   /// Default c'tor
   lyotCoronagraph();
   
   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */ 
   int wfSz();
   
   ///Set the wavefront size in pixels.
   /**
     */ 
   void wfSz(int sz /**< [in] is the new size */);

      
   ///Make the focal plane mask
   void makeFocalMask( realT rad, 
                       fpmaskFloatT trans = 0.0,
                       int sz = 0.0 );
   
   ///Load the apodizer from a FITS file 
   /**
     * \returns 0 on success.
     * \returns -1 on error. 
     */ 
   int loadApodizer( const std::string & apodName /**< [in] is the name of the FITS file containing the apodizer. */);
   
   ///Load the focal plane mask from a FITS file 
   /**
     * \returns 0 on success.
     * \returns -1 on error. 
     */ 
   int loadFocalMask( const std::string & fpmName /**< [in] is the name of the FITS file containing the focal plane mask. */);
   
   ///Load the Lyot stop from a FITS file
   /**
     * \returns 0 on success.
     * \returns -1 on error. 
     */ 
   int loadLyotStop( const std::string & lyotName /**< [in] is the name of the FITS file containing the lyot stop. */);

   ///Load the components of the coronagraph from FITS files 
   /**
     * \returns 0 on success.
     * \returns -1 on error.
     */ 
   int loadCoronagraph( const std::string & apodName, ///< [in] is the name of the FITS file containing the apodizer.
                        const std::string & fpmName, ///< [in] is the name of the FITS file containing the focal plane mask.
                        const std::string & lyotName ///< [in] is the name of the FITS file containing the lyot stop.
                      );
   
   ///Load the components of the coronagraph based in its base name
   /** Looks in _filDir for the files.
     * 
     * \returns 0 on success.
     * \returns -1 on error. 
     */ 
   int loadCoronagraph( const std::string & cName /**< [in] is the base name of the coronagraph without directory or file extension. */); 
   
   void applyApodizer( complexFieldT &pupilPlane );
   void applyFocalMask( complexFieldT & focalPlane );
   void applyLyotStop( complexFieldT &lyotPlane );
   
   ///Propagate the given pupil-plane wavefront through the coronagraph.
   int propagate( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);

   int propagate( imageT & fpIntensity,       ///< [out] The intensity image in the focal plane.  This should be pre-allocated.
                   complexFieldT & pupilPlane ///< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph.
                );
   
   /// Propagate the given pupil-plane wavefront without the coronagraph. 
   /** For a Lyot coronagraph, this applies the pupil apodization and Lyot stop, but not the FPM, to the given pupil-plane wavefront such 
     * that the result will produce the non-coronagraphic (off-axis) PSF.
     */
   int propagateNC( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);
   
   int propagateNC( imageT & fpIntensity,      ///< [out] The intensity image in the focal plane.  This should be pre-allocated.
                    complexFieldT & pupilPlane ///< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. 
                  );
   
   ///Optimize the pupil amplitude apodization and focal-plane mask complex transmission.
   /** Uses the algorithm in Guyon (2014) \cite guyon_2014 for optimizing an Apodized-Pupil Lyot Complex Mask Coronagraph (APLCMC).
     *
     */
   void optimizeAPLCMC( imageT & geomPupil,  ///< The geometric pupil mask, binary 1/0 transmission.
                        realT fpmRadPix, ///< The radius in pixels of the FPM.
                        realT relTol, ///< Relative tolerance for convergence.
                        realT absTol, ///< Absolute tolerance for convergence.
                        int maxIter, ///< Maximum number of iterations to allow.
                        const std::string & cname ///< Name of coronagraph, used as base-name for output files in m_fileDir.
                      );  
   
   void optimizeAPLCMC( imageT & geomPupil,  ///< The geometric pupil mask, binary 1/0 transmission.
                        realT fpmRadPix, ///< The radius in pixels of the FPM.
                        imageT & fpmPhase,
                        realT relTol, ///< Relative tolerance for convergence.
                        realT absTol, ///< Absolute tolerance for convergence.
                        int maxIter, ///< Maximum number of iterations to allow.
                        const std::string & cname ///< Name of coronagraph, used as base-name for output files in m_fileDir.
                      );  
   
};

template<typename _realT, typename _fpmaskFloatT>
lyotCoronagraph<_realT, _fpmaskFloatT>::lyotCoronagraph()
{
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::wfSz()
{
   return m_wfSz;
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::wfSz(int sz)
{
   m_wfSz = sz;
   
   m_fp.setWavefrontSizePixels(sz);
   m_focalPlane.resize(sz, sz);
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::makeFocalMask(_realT rad, 
                                                            _fpmaskFloatT trans, 
                                                            int sz )
{
   if(sz == 0)
   {
      sz = 2*rad+1;
      if(sz % 2 == 1) ++sz;
   }
   
   //Make a real circular mask
   improc::eigenImage<fpmaskFloatT> fpm(sz,sz);
   mx::wfp::circularPupil( fpm, 0., rad);

   //Convert to complex amplitude
   m_focalMask.resize(sz,sz);
   for(int i=0; i<sz; ++i)
   {
      for(int j=0; j<sz; ++j)
      {
         if(fpm(i,j) == 1) m_focalMask(i,j) = std::complex<fpmaskFloatT>(trans,0);
         else m_focalMask(i,j) = std::complex<fpmaskFloatT>(1,0);
      }
   }
   
   m_maskSource = 1;
   m_maskFile = "";
   m_maskRad = rad;
   m_maskTrans = 0.0;

}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::loadApodizer( const std::string & apodName)
{
   fits::fitsFile<_realT> ff;
   
   return ff.read(m_pupilApodizer, apodName);
   
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::loadFocalMask( const std::string & fpmName)
{
   fits::fitsFile<_realT> ff;
   improc::eigenCube<fpmaskFloatT> fpm;
   
   if(ff.read(fpm, fpmName) < 0) return -1;
   
   if(fpm.planes()==1)
   {
      m_focalMask = fpm.image(0);
   }
   else if(fpm.planes()==2)
   {
      m_focalMask.resize( fpm.rows(), fpm.cols());
      
      for( int r = 0; r < fpm.rows(); ++r)
      {
         for( int c = 0; c < fpm.cols(); ++c)
         {
            m_focalMask(r,c) = fpm.image(0)(r,c)*exp( std::complex<fpmaskFloatT>(0, fpm.image(1)(r,c)));
         }
      }
   }
   else
   {
      std::cerr << "too many planes in focal mask file\n";
      return -1;
   }
   
   m_maskSource = 0;
   m_maskFile = fpmName;
   m_maskRad = 0.0;
   m_maskTrans = 0.0;  
   
   return 0;
}
   
template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::loadLyotStop( const std::string & lyotName)
{
   fits::fitsFile<_realT> ff;
    
   return ff.read(m_lyotStop, lyotName);
   
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::loadCoronagraph( const std::string & apodName, 
                                                             const std::string & fpmName,
                                                             const std::string & lyotName)
{
   if(loadApodizer(apodName) < 0) return -1;
   if(loadFocalMask(fpmName) < 0) return -1;
   if(loadLyotStop(lyotName) < 0) return -1;
   
   return 0;
}
 
template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::loadCoronagraph( const std::string & cName)
{
   if( m_fileDir == "")
   {
      mxError("lyotCoronagraph", MXE_PARAMNOTSET, "file directory (fileDir) not set.");
      return -1;
   }

   std::string apodName= m_fileDir + "/" + cName + "_apod.fits";
   std::string fpmName = m_fileDir + "/" + cName + "_fpm.fits";
   std::string lyotName = m_fileDir + "/" + cName + "_lyot.fits";
      
   return loadCoronagraph(apodName, fpmName, lyotName);
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::applyApodizer( complexFieldT &pupilPlane )
{
   int sz = m_pupilApodizer.rows();
   
   realT w = 0.5*(sz - 1.0);
   
   realT xc = 0.5*(pupilPlane.rows() - 1);
   realT yc = 0.5*(pupilPlane.cols() - 1);
   
   for(int i=0; i<pupilPlane.rows(); ++i)
   {
      for(int j=0; j<pupilPlane.cols(); ++j)
      {
         if(i >= xc-w && i <= xc+w && j >= yc-w && j <= yc+w) pupilPlane( i, j) *= m_pupilApodizer((int)(i-(xc-w)),(int)(j-(yc-w)));
         else pupilPlane(i,j) *= 0;
      }
   }
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::applyFocalMask( complexFieldT &focalPlane )
{
   int sz = m_focalMask.rows();
   
   realT w = 0.5*(sz - 1.0);
   
   realT xc = 0.5*(focalPlane.rows() - 1);
   realT yc = 0.5*(focalPlane.cols() - 1);
   
   for(int i=0; i<sz; ++i)
   {
      for(int j=0; j<sz; ++j)
      {
         focalPlane( xc - w + i, yc - w + j ) *= m_focalMask(i,j);
      }
   }
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::applyLyotStop( complexFieldT & lyotPlane )
{
   int sz = m_lyotStop.rows();
   
   realT w = 0.5*(sz - 1.0);
   
   realT xc = 0.5*(lyotPlane.rows() - 1);
   realT yc = 0.5*(lyotPlane.cols() - 1);
   
   for(int i=0; i< lyotPlane.rows(); ++i)
   {
      for(int j=0; j< lyotPlane.cols(); ++j)
      {
         if(i >= xc-w && i <= xc+w && j >= yc-w && j <= yc+w) lyotPlane( i, j) *= m_lyotStop((int)(i - (xc-w)),(int)(j - (yc-w)));
         else lyotPlane(i,j) = 0;
      }
   }
   
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::propagate( complexFieldT & pupilPlane )
{
   applyApodizer( pupilPlane);
   
   m_fp.propagatePupilToFocal(m_focalPlane, pupilPlane); 
   
   if(m_savePreMaskFocalPlane)
   {
      m_preMaskFocalPlane = m_focalPlane;
   }
   
   applyFocalMask( m_focalPlane);
   
   m_fp.propagateFocalToPupil(pupilPlane, m_focalPlane);
   
   if(m_savePreLyotPupilPlane)
   {
      m_preLyotPupilPlane = pupilPlane;
   }
   
   applyLyotStop( pupilPlane);
   
   return 0;
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::propagate( imageT & fpIntensity,
                                                       complexFieldT & pupilPlane 
                                                     )
{
   propagate(pupilPlane);
   
   m_fp.propagatePupilToFocal(m_focalPlane, pupilPlane);
      
   int x0 = 0.5*(m_wfSz-1) - 0.5*(fpIntensity.rows()-1);
   int y0 = 0.5*(m_wfSz-1) - 0.5*(fpIntensity.cols()-1);
   
   extractIntensityImage(fpIntensity,0, fpIntensity.rows(),0,fpIntensity.cols(), m_focalPlane, x0,y0);
   
   return 0;
}
template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::propagateNC( complexFieldT & pupilPlane )
{
   applyApodizer( pupilPlane);
   
   m_fp.propagatePupilToFocal(m_focalPlane, pupilPlane); 
      
   m_fp.propagateFocalToPupil(pupilPlane, m_focalPlane);
   
   applyLyotStop( pupilPlane);
   
   return 0;
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::propagateNC( imageT & fpIntensity,
                                                         complexFieldT & pupilPlane 
                                                       )
{
   propagateNC(pupilPlane);
   
   m_fp.propagatePupilToFocal(m_focalPlane, pupilPlane);
      
   int x0 = 0.5*(m_wfSz-1) - 0.5*(fpIntensity.rows()-1);
   int y0 = 0.5*(m_wfSz-1) - 0.5*(fpIntensity.cols()-1);
   
   extractIntensityImage(fpIntensity,0, fpIntensity.rows(),0,fpIntensity.cols(), m_focalPlane, x0,y0);
   
   return 0;
}


template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::optimizeAPLCMC( imageT & geomPupil, 
                                                             realT fpmRadPix, 
                                                             realT relTol, 
                                                             realT absTol, 
                                                             int maxIter,
                                                             const std::string & cname )
{
   complexFieldT focalPlane, pupilPlane;

   pupilPlane.resize(m_wfSz, m_wfSz);
   focalPlane.resize(m_wfSz, m_wfSz);
   
   mx::wfp::imagingArray<realT,fftwAllocator<realT>, 0> mask(m_wfSz, m_wfSz);
   mx::wfp::circularPupil( mask, 0., fpmRadPix);   
   /*for(int c=0; c< mask.cols()*0.5; ++c)
   {
      for(int r=0; r < mask.rows(); ++r) mask(r,c) = 0;
   }*/

   //Initialize pupilImage
   mx::wfp::imagingArray<realT,fftwAllocator<realT>, 0> pupilImage(m_wfSz, m_wfSz);
   pupilImage.setZero();
   
   int gpLLi = 0.5*(m_wfSz-1) - 0.5*(geomPupil.rows()-1);
   int gpLLj = 0.5*(m_wfSz-1) - 0.5*(geomPupil.cols()-1);
   
   int gpURi = gpLLi + geomPupil.rows();
   int gpURj = gpLLj + geomPupil.cols();
   
   for(int i=0;i<geomPupil.rows();++i) 
   {
      for(int j=0;j< geomPupil.cols();++j) 
      {
         pupilImage(gpLLi + i, gpLLj + j) = geomPupil(i,j);
      }
   }
   
   realT lastLambdaA, LambdaA;
   
   lastLambdaA = 1;
   int n;
   std::string reason;
   for(n=0; n< maxIter; ++n)
   {
      mx::wfp::makeComplexPupil(pupilPlane, pupilImage, m_wfSz);
      
      m_fp.propagatePupilToFocal(focalPlane, pupilPlane);            

      for(int i=0; i < m_wfSz; ++i)
      {
         for(int j=0;j < m_wfSz; ++j)
         {
            focalPlane(i,j) *= mask(i,j);
         }
      }  
   
      m_fp.propagateFocalToPupil(pupilPlane, focalPlane);
   
      for(int i=0; i< m_wfSz; ++i) for(int j=0;j< m_wfSz; ++j) pupilImage(i,j) = abs(pupilPlane(i,j));
      
      
      LambdaA = 0;
      for(int i=0; i< m_wfSz; ++i)
      {
         for(int j=0;j< m_wfSz; ++j)
         {
            if( i >= gpLLi && i < gpURi && j >= gpLLj && j < gpURj)
            {
               pupilImage(i,j) *= geomPupil(i-gpLLi,j-gpLLj); 
               if(pupilImage(i,j) > LambdaA) LambdaA = pupilImage(i,j);
            }
            else
            {
               pupilImage(i,j) = 0;
            }
         }
      }
      //LambdaA = 1.0/LambdaA;
      for(int i=0; i< m_wfSz; ++i) for(int j=0;j<m_wfSz; ++j) pupilImage(i,j) /= LambdaA;
   
      std::cout <<  n << " " << LambdaA << "\n";
      if( fabs(LambdaA - lastLambdaA) < absTol)
      {
         std::cout << "Converged on absTol.\n";
         reason = "absTol";
         break;
      }
      if( fabs((LambdaA - lastLambdaA)/lastLambdaA) < relTol)
      {
         std::cout << "Converged on relTol.\n";
         reason = "relTol";
         break;
      }
      
      if(n == maxIter - 1)
      {
         std::cout << "maxIter reached.\n";
         reason = "maxIter";
      }
      
      lastLambdaA = LambdaA;
   }
   if(LambdaA > 1) LambdaA = 1;
      
   std::cout << "LambdaA: = " << LambdaA << "\n";
   
   realT trans = 1.0 - 1.0/LambdaA;

   int pupSize = geomPupil.rows();
   
   m_pupilApodizer.resize(pupSize, pupSize);
   
   extractBlock(m_pupilApodizer, 0, pupSize, 0, pupSize, pupilImage, 0.5*( (pupilImage.rows()-1) - (pupSize-1)), 0.5*( (pupilImage.rows()-1) - (pupSize-1))); 

   makeFocalMask(fpmRadPix, trans, pupSize);

   /*for(int c=0; c< m_focalMask.cols()*0.5; ++c)
   {
      for(int r=0; r < m_focalMask.rows(); ++r) 
      {
         if(m_focalMask(r,c).real()==1) m_focalMask(r,c) = 0;
      }
   }*/

   m_maskSource = 2;
   
   m_lyotStop = geomPupil;
   
   
   fits::fitsHeader head;
   
   head.append("", fits::fitsCommentType(), "----------------------------------------");
   head.append("", fits::fitsCommentType(), "lyotCoronagraph optimization Parameters:");
   head.append("", fits::fitsCommentType(), "----------------------------------------");
   head.append<int>("WFSZ", m_wfSz, "Size of wavefront used for FFTs (pixels)");
   head.append<realT>("FPMRADPX", fpmRadPix, "input radius of focal plane mask (pixels)");
   head.append<realT>("ABSTOL", absTol , "input absolute tolerance");
   head.append<realT>("RELTOL", relTol , "input relative tolerance");
   head.append<int>("MAXITER", maxIter , "input maximum iterations");
   head.append<int>("NITER", n , "actual number of iterations");
   head.append<std::string>("XREASON", reason , "reason for convergence");
   head.append<realT>("FPMTRANS", trans, "transmission of FPM");
   
   
   fits::fitsFile<double> ff;

   std::string fname = "coron/" + cname + "_apod.fits";
   
   ff.write(fname, m_pupilApodizer, head);
   
   fname = "coron/" + cname + "_fpm.fits";
   improc::eigenImage<fpmaskFloatT> fpm(m_focalMask.rows(), m_focalMask.cols());
   for(int r=0;r<m_focalMask.rows(); ++r)
   {
      for(int c=0; c< m_focalMask.cols(); ++c)
      {
         fpm(r,c) = m_focalMask(r,c).real();
      }
   }
   
   ff.write(fname, fpm,head);
   
   fname = "coron/" + cname + "_lyot.fits";
   ff.write(fname, m_lyotStop, head);
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::optimizeAPLCMC( imageT & geomPupil, 
                                                             realT fpmRadPix, 
                                                             imageT & fpmPhase,
                                                             realT relTol, 
                                                             realT absTol, 
                                                             int maxIter,
                                                             const std::string & cname )
{
   complexFieldT focalPlane, pupilPlane;

   pupilPlane.resize(m_wfSz, m_wfSz);
   focalPlane.resize(m_wfSz, m_wfSz);
   
   mx::wfp::imagingArray<realT,fftwAllocator<realT>, 0> mask(m_wfSz, m_wfSz);
   mx::wfp::circularPupil( mask, 0., fpmRadPix);   
  
   //Initialize pupilImage
   mx::wfp::imagingArray<realT,fftwAllocator<realT>, 0> pupilImage(m_wfSz, m_wfSz);
   pupilImage.setZero();
   
   int gpLLi = 0.5*(m_wfSz-1) - 0.5*(geomPupil.rows()-1);
   int gpLLj = 0.5*(m_wfSz-1) - 0.5*(geomPupil.cols()-1);
   
   int gpURi = gpLLi + geomPupil.rows();
   int gpURj = gpLLj + geomPupil.cols();
   
   for(int i=0;i<geomPupil.rows();++i) 
   {
      for(int j=0;j< geomPupil.cols();++j) 
      {
         pupilImage(gpLLi + i, gpLLj + j) = geomPupil(i,j);
      }
   }
   
   realT lastLambdaA, LambdaA;
   
   lastLambdaA = 1;
   int n;
   std::string reason;
   for(n=0; n< maxIter; ++n)
   {
      mx::wfp::makeComplexPupil(pupilPlane, pupilImage, m_wfSz);
      
      m_fp.propagatePupilToFocal(focalPlane, pupilPlane);            

      for(int i=0; i < m_wfSz; ++i)
      {
         for(int j=0;j < m_wfSz; ++j)
         {
            focalPlane(i,j) *= mask(i,j)*exp( std::complex<realT>(0, fpmPhase(i,j)));;
         }
      }  
   
      m_fp.propagateFocalToPupil(pupilPlane, focalPlane);
   
      for(int i=0; i< m_wfSz; ++i) for(int j=0;j< m_wfSz; ++j) pupilImage(i,j) = abs(pupilPlane(i,j));
      
      
      LambdaA = 0;
      for(int i=0; i< m_wfSz; ++i)
      {
         for(int j=0;j< m_wfSz; ++j)
         {
            if( i >= gpLLi && i < gpURi && j >= gpLLj && j < gpURj)
            {
               pupilImage(i,j) *= geomPupil(i-gpLLi,j-gpLLj); 
               if(pupilImage(i,j) > LambdaA) LambdaA = pupilImage(i,j);
            }
            else
            {
               pupilImage(i,j) = 0;
            }
         }
      }
      //LambdaA = 1.0/LambdaA;
      for(int i=0; i< m_wfSz; ++i) for(int j=0;j<m_wfSz; ++j) pupilImage(i,j) /= LambdaA;
   
      std::cout <<  n << " " << LambdaA << "\n";
      if( fabs(LambdaA - lastLambdaA) < absTol)
      {
         std::cout << "Converged on absTol.\n";
         reason = "absTol";
         break;
      }
      if( fabs((LambdaA - lastLambdaA)/lastLambdaA) < relTol)
      {
         std::cout << "Converged on relTol.\n";
         reason = "relTol";
         break;
      }
      
      if(n == maxIter - 1)
      {
         std::cout << "maxIter reached.\n";
         reason = "maxIter";
      }
      
      lastLambdaA = LambdaA;
   }
   if(LambdaA > 1) LambdaA = 1;
      
   std::cout << "LambdaA: = " << LambdaA << "\n";
   
   realT trans = 1.0 - 1.0/LambdaA;

   int pupSize = geomPupil.rows();
   
   m_pupilApodizer.resize(pupSize, pupSize);
   
   extractBlock(m_pupilApodizer, 0, pupSize, 0, pupSize, pupilImage, 0.5*( (pupilImage.rows()-1) - (pupSize-1)), 0.5*( (pupilImage.rows()-1) - (pupSize-1))); 

   makeFocalMask(fpmRadPix, trans, pupSize);
   m_maskSource = 2;
   
   m_lyotStop = geomPupil;
   
   
   fits::fitsHeader head;
   
   head.append("", fits::fitsCommentType(), "----------------------------------------");
   head.append("", fits::fitsCommentType(), "lyotCoronagraph optimization Parameters:");
   head.append("", fits::fitsCommentType(), "----------------------------------------");
   head.append<int>("WFSZ", m_wfSz, "Size of wavefront used for FFTs (pixels)");
   head.append<realT>("FPMRADPX", fpmRadPix, "input radius of focal plane mask (pixels)");
   head.append<realT>("ABSTOL", absTol , "input absolute tolerance");
   head.append<realT>("RELTOL", relTol , "input relative tolerance");
   head.append<int>("MAXITER", maxIter , "input maximum iterations");
   head.append<int>("NITER", n , "actual number of iterations");
   head.append<std::string>("XREASON", reason , "reason for convergence");
   head.append<realT>("FPMTRANS", trans, "transmission of FPM");
   
   
   fits::fitsFile<double> ff;

   std::string fname = "coron/" + cname + "_apod.fits";
   
   ff.write(fname, m_pupilApodizer, head);
   
   fname = "coron/" + cname + "_fpm.fits";
   improc::eigenImage<fpmaskFloatT> fpm(m_focalMask.rows(), m_focalMask.cols());
   for(int r=0;r<m_focalMask.rows(); ++r)
   {
      for(int c=0; c< m_focalMask.cols(); ++c)
      {
         fpm(r,c) = m_focalMask(r,c).real();
      }
   }
   
   ff.write(fname, fpm,head);
   
   fname = "coron/" + cname + "_lyot.fits";
   ff.write(fname, m_lyotStop, head);
}

} //namespace wfp   
} //namespace mx

#endif //__lyotCoronagraph_hpp__

