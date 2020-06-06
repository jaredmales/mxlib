/** \file pupil.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Utilities for specifying pupils.
  * \ingroup mxAO_files
  * 
  */

#ifndef pupil_hpp
#define pupil_hpp

#include <mx/improc/fitsFile.hpp>
#include <mx/improc/eigenImage.hpp>

#include <mx/wfp/imagingUtils.hpp>

#include <mx/sigproc/signalWindows.hpp>

#include "aoPaths.hpp"

namespace mx
{
   
namespace AO
{

///Generate a circular pupil and saves it to disk
template<typename realT>
void circularPupil( const std::string & pupilName,
                    realT pupilDiamPixels,
                    realT pupilDiamMeters,
                    realT centralObs = 0,
                    realT overscan = 0
                  )
{
   using namespace mx::improc;
   using namespace mx::wfp;
   using namespace mx::sigproc;   

   /*Create pupil*/
   eigenImage<realT> pup;
   
   pup.resize(pupilDiamPixels, pupilDiamPixels);
   
   wfp::circularPupil( pup, centralObs,0,overscan);
   
   fitsHeader phead;
   phead.append("SCALE", pupilDiamMeters/pupilDiamPixels, "Scale in m/pix");
   phead.append("PUPILD", pupilDiamMeters, "Physical diameter of pupil image [m]");
   phead.append("CENTOBS", centralObs, "Central obscuration ratio");
   phead.append("OVERSCAN", overscan, "Fractional pixel overscan");
   
   fitsFile<realT> ff;
   
   std::string fName = mx::AO::path::pupil::pupilFile(pupilName, true);

   ff.write(fName, pup, phead);
   
}


///Generates a circular apodized pupil and saves it to disk
/** Apodization is with a Tukey window.
  */
template<typename realT>
void circularApodizedPupil( const std::string & pupilName,
                            int pupilDiamPixels,
                            realT pupilDiamMeters,
                            realT tukeyAlpha,
                            realT centralObs = 0,
                            realT overScan = 0)
{

   using namespace mx::improc;
   using namespace mx::wfp;
   using namespace mx::sigproc;   
   
   /*Create pupil*/
   Eigen::Array<realT, -1, -1> pup;
   
   pup.resize(pupilDiamPixels, pupilDiamPixels);
   
   
   realT cen = 0.5*(pupilDiamPixels - 1.0);
   
   if(centralObs == 0)
   {
      window::tukey2d<realT>(pup.data(), pupilDiamPixels, pupilDiamPixels + overScan, tukeyAlpha, cen,cen);
   }
   else
   {
      window::tukey2dAnnulus<realT>(pup.data(), pupilDiamPixels, pupilDiamPixels + overScan, centralObs, tukeyAlpha, cen,cen);
   }
   
   
   
   fitsHeader phead;
   phead.append("SCALE", pupilDiamMeters/pupilDiamPixels, "Scale in m/pix");
   phead.append("PUPILD", pupilDiamMeters, "Physical diameter of pupil image [m]");
   phead.append("CENTOBS", centralObs, "Central obscuration ratio");
   phead.append("TUKALPHA", tukeyAlpha, "Tukey window alpha parameter");
   phead.append("OVERSCAN", overScan, "Apodization overscan");
    
   fitsFile<realT> ff;
   
   std::string fName = mx::AO::path::pupil::pupilFile(pupilName, true);

   ff.write(fName, pup, phead);
   
}


} //namespace mx
} //namespace AO

   
#endif //__basis_hpp__
