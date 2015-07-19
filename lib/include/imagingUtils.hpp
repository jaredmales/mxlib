/** \file imagingUtils.hpp
  * \brief Utilities for modeling image formation
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __imagingUtils_hpp__
#define __imagingUtils_hpp__

#include <cmath>

namespace mx
{

/// Fill in an Eigen Array with a circular pupil mask.
/** \ingroup imaging
  * \param m is the allocated Array
  * \param eps [optional] is the central obscuration.  0-1, default is 0.
  * \param rad [optional] is the desired radius. if 0 the maximum radius is used.
  */  
template<class arrayT> 
void circularPupil( arrayT & m, 
                     typename arrayT::Scalar eps=0, 
                     typename arrayT::Scalar rad=0 
                   )
{
   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   typename arrayT::Scalar r;
   typename arrayT::Scalar xc = 0.5*(l0-1);
   typename arrayT::Scalar yc = 0.5*(l1-1);
   
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
template<typename arithT>
void makeComplexPupil(Eigen::Array<std::complex<arithT>, Eigen::Dynamic, Eigen::Dynamic> & complexPupil, 
                       Eigen::Array<arithT, Eigen::Dynamic, Eigen::Dynamic> & realPupil, int wavefrontSizePixels)
{
   
   complexPupil.resize(wavefrontSizePixels, wavefrontSizePixels);
   complexPupil.setZero();
   
   //complexPupil.bottomRightCorner(realPupil.rows(), realPupil.cols()) = realPupil*std::complex<arithT>(1,0);
   
   int bl = 0.5*(complexPupil.rows()-1) - 0.5*(realPupil.rows()-1.);
   //int ur = 0.5*(complexPupil.rows()-1) + 0.5*(realPupil.rows()-1.);
   
   complexPupil.block(bl, bl, realPupil.rows(), realPupil.rows()) = realPupil*std::complex<arithT>(1,0);

}

///Apply a tilt to a wavefront 
/**
  * \param complexWavefront [in/out] the complex wavefront to tilt, will be modified on output
  * \param xTilt [input] the amount of tilt in the x direction, in pixels 
  * \param yTilt [input] the amount of tilt in the y direction, in pixels 
  * 
  * \ingroup imaging
  */
template<typename arithT>
void tiltWavefront(Eigen::Array<std::complex<arithT>, Eigen::Dynamic, Eigen::Dynamic> & complexWavefront, 
                       arithT xTilt, arithT yTilt)
{
   typedef std::complex<arithT> complexT;
   
   int wfsSizeX = complexWavefront.cols();
   int wfsSizeY = complexWavefront.rows();
   
   arithT xCen = 0.5*wfsSizeX;
   arithT yCen = 0.5*wfsSizeY;

   arithT argX = D2PI/(wfsSizeX-1);
   arithT argY = D2PI/(wfsSizeY-1);
  
   for(int ii=0; ii < wfsSizeX; ++ii)
   {
      for(int jj=0; jj < wfsSizeY; ++jj)
      {     
         complexWavefront(ii,jj) = complexWavefront(ii,jj)*exp(complexT(0.,argX*xTilt*(ii-xCen)+argY*yTilt*(jj-yCen)));
      }
   }
}

} //namespace mx

#endif //__imagingUtils_hpp__

