/** \file lyotCoronagraph.hpp
  * \brief Declares and defines a class to describe and optimize a Lyot Coronagraph.
  * \ingroup imaging_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __lyotCoronagraph_hpp__
#define __lyotCoronagraph_hpp__

#include <iostream>

#include "imagingArray.hpp"
#include "imagingUtils.hpp"
#include "fraunhoferImager.hpp"

#include "../fitsFile.hpp"
#include "../fitsUtils.hpp"

#include "../eigenImage.hpp"

namespace mx
{
   
namespace imaging 
{
   
/// The Lyot Coronagraph
/**
  * \ingroup coronagraphs
  */   
template<typename _realT, typename _fpmaskFloatT>
struct lyotCoronagraph
{
   typedef _realT realT;
   typedef _fpmaskFloatT fpmaskFloatT;
   
   typedef std::complex<realT> complexT;
   
   ///The wavefront complex field type
   typedef mx::imagingArray<std::complex<realT>, fftwAllocator<std::complex<realT> >, 0> complexFieldT;
   
   ///The focal plane mask type
   typedef Eigen::Array< fpmaskFloatT, Eigen::Dynamic, Eigen::Dynamic> fpMaskT;
   
   ///The image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

   std::string _fileDir;
   
   int _wfSz;
   
   imageT pupilApodizer;
   
   
   int maskSource; ///< 0= read from file, 1 = constructed by makeFocalMask, 2 = trans optimized
   std::string maskFile; ///< Name of file from which mask was loaded
   realT maskRad; ///<Radius of mask if it was constructed
   realT maskTrans; ///<Transmission of mask if it was constructed
   
   fpMaskT focalMask;
   
   
   
   imageT lyotStop;

   complexFieldT focalPlane;

   bool savePreMaskFocalPlane;
   complexFieldT preMaskFocalPlane;
   
   bool savePreLyotPupilPlane;
   complexFieldT preLyotPupilPlane;
   
   
   
   ///Fraunhofer propagator
   mx::fraunhoferImager<complexFieldT> fi;

   lyotCoronagraph();
   
   ///Get the wavefront size in pixels
   /**
     * \returns the wavefront size in pixels
     */ 
   int wfSz();
   
   ///Set the wavefront size in pixels.
   /**
     * \param sz is the new size
     */ 
   void wfSz(int sz);

      
   ///Make the focal plane mask
   void makeFocalMask( realT rad, 
                       fpmaskFloatT trans = 0.0,
                       int sz = 0.0 );
   
   ///Load the apodizer from a FITS file 
   /**
     * \param apodName is the name of the FITS file
     */ 
   void loadApodizer( const std::string & apodName);
   
   ///Load the focal plane mask from a FITS file 
   /**
     * \param fpmName is the name of the FITS file
     */ 
   void loadFocalMask( const std::string & fpmName);
   
   ///Load the Lyot stop from a FITS file
   /**
     * \param lyotName is the name of the FITS file
     */ 
   void loadLyotStop( const std::string & lyotName);
   
   ///Load the components of the coronagraph based in its base name
   /** Looks in "./coron/"
     * 
     * \param cName is the base name of the coronagraph
     */ 
   void loadCoronagraph( const std::string & cName); 
   
   ///Load the components of the coronagraph from FITS files 
   /**
     * \param apodName is the name of the FITS file
     * \param fpmName is the name of the FITS file
     * \param lyotName is the name of the FITS file
     */ 
   void loadCoronagraph( const std::string & apodName, 
                         const std::string & fpmName,
                         const std::string & lyotName);
   
   
   
   void applyApodizer( complexFieldT &pupilPlane );
   void applyFocalMask( complexFieldT & focalPlane );
   void applyLyotStop( complexFieldT &lyotPlane );
   
   ///Propagate the given pupil-plane wavefront through the coronagraph.
   void propagate( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);

   /// Propagate the given pupil-plane wavefront without the coronagraph. 
   /** For a Lyot coronagraph, this applies the pupil apodization and Lyot stop, but not the FPM, to the given pupil-plane wavefront such 
     * that the result will produce the non-coronagraphic (off-axis) PSF.
     */
   void propagateNC( complexFieldT & pupilPlane /**< [in/out] The wavefront at the input pupil plane.  It is modified by the coronagraph. */);
   
   
   void optimizeApodizer( imageT & geomPupil, 
                          realT fpmRadPix, 
                          realT relTol, 
                          realT absTol, 
                          int maxIter,
                          const std::string & cname );  
   
};

template<typename _realT, typename _fpmaskFloatT>
lyotCoronagraph<_realT, _fpmaskFloatT>::lyotCoronagraph()
{
   _fileDir = "coron/";
   
   savePreLyotPupilPlane = false;
   savePreMaskFocalPlane = false;
}

template<typename _realT, typename _fpmaskFloatT>
int lyotCoronagraph<_realT, _fpmaskFloatT>::wfSz()
{
   return _wfSz;
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::wfSz(int sz)
{
   _wfSz = sz;
   
   fi.setWavefrontSizePixels(sz);
   focalPlane.resize(sz, sz);
   
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
   
   focalMask.resize(sz, sz);
   mx::circularPupil( focalMask, 0., rad);
      
   for(int i=0; i<sz; ++i)
   {
      for(int j=0; j<sz; ++j)
      {
         if(focalMask(i,j) == 1) focalMask(i,j) *= trans;
         else focalMask(i,j) = 1;
      }
   }
   
   maskSource = 1;
   maskFile = "";
   maskRad = rad;
   maskTrans = 0.0;

}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::loadApodizer( const std::string & apodName)
{
   fitsFile<_realT> ff;
   
   ff.read(apodName, pupilApodizer);
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::loadFocalMask( const std::string & fpmName)
{
   fitsFile<_realT> ff;
   
   ff.read(fpmName, focalMask);
   
   maskSource = 0;
   maskFile = fpmName;
   maskRad = 0.0;
   maskTrans = 0.0;   
}
   
template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::loadLyotStop( const std::string & lyotName)
{
   fitsFile<_realT> ff;
   
   ff.read(lyotName, lyotStop);
}

 
template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::loadCoronagraph( const std::string & cName)
{
   std::string apodName= _fileDir + cName + "_apod.fits";
   std::string fpmName = _fileDir + cName + "_fpm.fits";
   std::string lyotName = _fileDir + cName + "_lyot.fits";
   
   loadApodizer(apodName);
   loadFocalMask(fpmName);
   loadLyotStop(lyotName);
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::loadCoronagraph( const std::string & apodName, 
                                                               const std::string & fpmName,
                                                               const std::string & lyotName)
{
   loadApodizer(apodName);
   loadFocalMask(fpmName);
   loadLyotStop(lyotName);
}


template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::applyApodizer( complexFieldT &pupilPlane )
{
   int sz = pupilApodizer.rows();
   
   realT w = 0.5*(sz - 1.0);
   
   realT xc = 0.5*(pupilPlane.rows() - 1);
   realT yc = 0.5*(pupilPlane.cols() - 1);
   
   for(int i=0; i<pupilPlane.rows(); ++i)
   {
      for(int j=0; j<pupilPlane.cols(); ++j)
      {
         if(i >= xc-w && i <= xc+w && j >= yc-w && j <= yc+w) pupilPlane( i, j) *= pupilApodizer(i-(xc-w),j-(yc-w));
         else pupilPlane(i,j) *= 0;
      }
   }
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::applyFocalMask( complexFieldT &focalPlane )
{
   int sz = focalMask.rows();
   
   realT w = 0.5*(sz - 1.0);
   
   realT xc = 0.5*(focalPlane.rows() - 1);
   realT yc = 0.5*(focalPlane.cols() - 1);
   
   for(int i=0; i<sz; ++i)
   {
      for(int j=0; j<sz; ++j)
      {
         focalPlane( xc - w + i, yc - w + j ) *= focalMask(i,j);
      }
   }
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::applyLyotStop( complexFieldT & lyotPlane )
{
   int sz = lyotStop.rows();
   
   realT w = 0.5*(sz - 1.0);
   
   realT xc = 0.5*(lyotPlane.rows() - 1);
   realT yc = 0.5*(lyotPlane.cols() - 1);
   
   for(int i=0; i< lyotPlane.rows(); ++i)
   {
      for(int j=0; j< lyotPlane.cols(); ++j)
      {
         if(i >= xc-w && i <= xc+w && j >= yc-w && j <= yc+w) lyotPlane( i, j) *= lyotStop(i - (xc-w),j - (yc-w));
         else lyotPlane(i,j) = 0;
      }
   }
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::propagate( complexFieldT & pupilPlane )
{
   applyApodizer( pupilPlane);
   
   fi.propagatePupilToFocal(focalPlane, pupilPlane); 
   
   if(savePreMaskFocalPlane)
   {
      preMaskFocalPlane = focalPlane;
   }
   
   applyFocalMask( focalPlane);
   
   fi.propagateFocalToPupil(pupilPlane, focalPlane);
   
   if(savePreLyotPupilPlane)
   {
      preLyotPupilPlane = pupilPlane;
   }
   
   applyLyotStop( pupilPlane);
   
}

template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::propagateNC( complexFieldT & pupilPlane )
{
   applyApodizer( pupilPlane);
   
   fi.propagatePupilToFocal(focalPlane, pupilPlane); 
      
   fi.propagateFocalToPupil(pupilPlane, focalPlane);
   
   applyLyotStop( pupilPlane);
   
}


template<typename _realT, typename _fpmaskFloatT>
void lyotCoronagraph<_realT, _fpmaskFloatT>::optimizeApodizer( imageT & geomPupil, 
                                                                realT fpmRadPix, 
                                                                realT relTol, 
                                                                realT absTol, 
                                                                int maxIter,
                                                                const std::string & cname )
{
   complexFieldT focalPlane, pupilPlane;

   pupilPlane.resize(_wfSz, _wfSz);
   focalPlane.resize(_wfSz, _wfSz);
   
   mx::imagingArray<realT,mx::fftwAllocator<realT>, 0> mask(_wfSz, _wfSz);
   mx::circularPupil( mask, 0., fpmRadPix);   
  
   //Initialize pupilImage
   mx::imagingArray<realT,mx::fftwAllocator<realT>, 0> pupilImage(_wfSz, _wfSz);
   pupilImage.setZero();
   
   int gpLLi = 0.5*(_wfSz-1) - 0.5*(geomPupil.rows()-1);
   int gpLLj = 0.5*(_wfSz-1) - 0.5*(geomPupil.cols()-1);
   
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
      mx::makeComplexPupil(pupilPlane, pupilImage, _wfSz);
      
      fi.propagatePupilToFocal(focalPlane, pupilPlane);            

      for(int i=0; i < _wfSz; ++i)
      {
         for(int j=0;j < _wfSz; ++j)
         {
            focalPlane(i,j) *= mask(i,j);
         }
      }  
   
      fi.propagateFocalToPupil(pupilPlane, focalPlane);
   
      for(int i=0; i< _wfSz; ++i) for(int j=0;j< _wfSz; ++j) pupilImage(i,j) = abs(pupilPlane(i,j));
      
      
      LambdaA = 0;
      for(int i=0; i< _wfSz; ++i)
      {
         for(int j=0;j< _wfSz; ++j)
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
      for(int i=0; i< _wfSz; ++i) for(int j=0;j<_wfSz; ++j) pupilImage(i,j) /= LambdaA;
   
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
   
   pupilApodizer.resize(pupSize, pupSize);
   
   extractBlock(pupilApodizer, 0, pupSize, 0, pupSize, pupilImage, 0.5*( (pupilImage.rows()-1) - (pupSize-1)), 0.5*( (pupilImage.rows()-1) - (pupSize-1))); 

   makeFocalMask(fpmRadPix, trans, pupSize);
   maskSource = 2;
   
   lyotStop = geomPupil;
   
   
   fitsHeader head;
   
   head.append("", fitsCommentType(), "----------------------------------------");
   head.append("", fitsCommentType(), "lyotCoronagraph optimization Parameters:");
   head.append("", fitsCommentType(), "----------------------------------------");
   head.append<int>("WFSZ", _wfSz, "Size of wavefront used for FFTs (pixels)");
   head.append<realT>("FPMRADPX", fpmRadPix, "input radius of focal plane mask (pixels)");
   head.append<realT>("ABSTOL", absTol , "input absolute tolerance");
   head.append<realT>("RELTOL", relTol , "input relative tolerance");
   head.append<int>("MAXITER", maxIter , "input maximum iterations");
   head.append<int>("NITER", n , "actual number of iterations");
   head.append<std::string>("XREASON", reason , "reason for convergence");
   head.append<realT>("FPMTRANS", trans, "transmission of FPM");
   
   
   mx::fitsFile<double> ff;

   std::string fname = "coron/" + cname + "_apod.fits";
   
   ff.write(fname, pupilApodizer, head);
   
   fname = "coron/" + cname + "_fpm.fits";
   ff.write(fname, focalMask,head);
   
   fname = "coron/" + cname + "_lyot.fits";
   ff.write(fname, lyotStop, head);
   
   

   
}

} //namespace imaging    
} //namespace mx

#endif //__lyotCoronagraph_hpp__

