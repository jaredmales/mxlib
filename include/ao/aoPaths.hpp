/** \file aoPaths.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Standardized paths for the mx::AO system.
  * \ingroup mxAO_files
  * 
  */

#ifndef __aoPaths_hpp__
#define __aoPaths_hpp__

#include "../sys/environment.hpp"

#include "../ioutils/fileUtils.hpp"

namespace mx
{

namespace AO
{




///Namespace for paths
/** \ingroup mxAO_paths
 */
namespace path
{

   
///Namespace for basis paths
/** \ingroup mxAO_paths
 */
namespace basis
{

   ///The root path for basis files
   /** 
     * \param[in] basisName the name of the basis  
     * \param[in] create [optional] create the directory if it noes not exist.
     * 
     * \returns the root path for basis files.
     * 
     * \ingroup mxAO_paths
     */
   std::string root(const std::string & basisName, bool create = false)
   {
      
      std::string path = mx::sys::getEnv("MX_AO_DATADIR");

      path += "/basis/" + basisName;

      if(create)
      {
         ioutils::createDirectories(path);
      }
      return path;
   }

   std::string modes(const std::string & basisName, bool create = false)
   {
      std::string path = root(basisName, create) + "/modes.fits";

      return path;
   }

   std::string spectrum(const std::string & basisName, bool create = false)
   {
      std::string path = root(basisName, create) + "/spectrum.fits";

      return path;
   }
} //namespace basis

//--------------------------------------------------------------------------------


namespace dm
{

   ///The root path for deformable mirror (DM) files
   /** 
     * \param[in] dmName the name of the DM 
     * \param[in] create [optional] create the directory if it noes not exist.
     * 
     * \returns the root path for DM files.
     */
   std::string root(const std::string & dmName, bool create = false)
   {
      std::string path = mx::sys::getEnv("MX_AO_DATADIR");

      path += "/dm/" + dmName;

      if(create)
      {
         ioutils::createDirectories(path);
      }

      return path;
   }

   ///The path for the deformable mirror (DM) influence functions.
   /** 
     * \param[in] dmName the name of the DM 
     * \param[in] create [optional] create the root directory if it noes not exist.
     * 
     * \returns the path to the FITS file containing the DM influence functions.
     */
   std::string influenceFunctions(const std::string & dmName, bool create = false)
   { 
      std::string path = root(dmName, create) + "/inf.fits";

      return path;
   }

   ///The path for the deformable mirror (DM) actuator positions
   /** 
     * \param[in] dmName the name of the DM 
     * \param[in] create [optional] create the root directory if it noes not exist.
     * 
     * \returns the path to the FITS file containing the DM actuator positions
     */
   std::string actuatorPositions(const std::string & dmName, bool create = false)
   { 
      std::string path = root(dmName, create) + "/actPos.fits";

      return path;
   }

   ///The path for the deformable mirror (DM) influence function pseudo-inverse.
   /** 
     * \param[in] dmName the name of the DM 
     * \param[in] create [optional] create the root directory if it noes not exist.
     * 
     * \returns the path to the FITS file containing the DM influence function pseudo-inverse matrix
     */
   std::string pseudoInverse(const std::string & dmName, bool create = false)
   {
      std::string path = root(dmName, create) + "/pinv.fits";

      return path;
   }

///The path for the deformable mirror (DM) influence function based mirror modes.
/** 
  * \param[in] dmName the name of the DM 
  * \param[in] create [optional] create the root directory if it noes not exist.
  * 
  * \returns the path to the FITS file containing the DM influence function based mirror modes
  */
std::string mirrorModes(const std::string & dmName, bool create = false)
{
   std::string path = basis::modes(dmName, create);

   return path;
}

///The path for the deformable mirror (DM) influence function pseudo-inverse singular values.
/** 
  * \param[in] dmName the name of the DM 
  * \param[in] create [optional] create the root directory if it noes not exist.
  * 
  * \returns the path to the FITS file containing the DM influence function pseudo-inverse singular values
  */
std::string singularValues(const std::string & dmName, bool create = false)
{
   std::string path = root(dmName, create) + "/singularValues.dat";

   return path;
}

///The root path for the deformable mirror (DM) basis related files.
/** 
  * \param[in] dmName the name of the DM 
  * \param[in] basisName the name of the basis set
  * \param[in] create [optional] create the root directory if it noes not exist.
  * 
  * \returns the root path for the DM basis related files
  */
std::string basisRoot(const std::string & dmName, const std::string & basisName, bool create = false)
{
   std::string path = root(dmName, create) + "/basis/" + basisName;
   
   if(create)
   {
      ioutils::createDirectories(path);
   }
   
   return path;
}

///The path for the modes-to-commands (M2c) matrix for a deformable mirror (DM) and a basis set.
/** 
  * \param[in] dmName the name of the DM 
  * \param[in] basisName the name of the basis set
  * \param[in] create [optional] create the root directory if it noes not exist.
  * 
  * \returns the path to the FITS file containting the M2c matrix for a DM and a basis set.
  */
std::string M2c(const std::string & dmName, const std::string & basisName, bool create = false)
{
   std::string path = basisRoot(dmName, basisName, create);
      
   if(create)
   {
      ioutils::createDirectories(path);
   }
   
   path += "/M2c.fits";

   return path;
}

// ///The path for the modes-to-commands (M2c) matrix for a deformable mirror (DM) and an orthogonalized basis set.
// /** 
//   * \param[in] dmName the name of the DM 
//   * \param[in] basisName the name of the basis set
//   * \param[in] create [optional] create the root directory if it noes not exist.
//   * 
//   * \returns the path to the FITS file containting the Mortho2c matrix for a DM and a basis set.
//   */
// std::string Mortho2c(const std::string & dmName, const std::string & basisName, const std::string & pupilName, bool create = false)
// {
// 
//    std::string path = basisRoot(dmName, basisName, create);
//    
//    path += "/" + pupilName;
//    
//    if(create)
//    {
//       ioutils::createDirectories(path);
//    }
//    
//    path += "/Mortho2c.fits";
// 
//    return path;
// }

///The path for the projected modes for a deformable mirror (DM) and a basis set.
/** These are the modes as reproduced by the DM actuators.
  * 
  * \param[in] dmName the name of the DM 
  * \param[in] basisName the name of the basis set
  * \param[in] pupilName the name of the pupil
  * \param[in] create [optional] create the root directory if it noes not exist.
  * 
  * \returns the path to the FITS file containting the M2c matrix for a DM and a basis set.
  */
std::string projectedModes(const std::string & dmName, const std::string & basisName, bool create = false)
{

   std::string path = basisRoot(dmName, basisName, create);
   
   if(create)
   {
      ioutils::createDirectories(path);
   }
   
   path += "/projectedModes.fits";

   return path;
}

// ///The path for the projected modes for a deformable mirror (DM) and an orthogonalized basis set.
// /** These are the modes as reproduced by the DM actuators.
//   * 
//   * \param[in] dmName the name of the DM 
//   * \param[in] basisName the name of the basis set
//   * \param[in] pupilName the name of the pupil
//   * \param[in] create [optional] create the root directory if it noes not exist.
//   * 
//   * \returns the path to the FITS file containting the M2c matrix for a DM and a basis set.
//   */
// std::string projectedOrthoModes(const std::string & dmName, const std::string & basisName, const std::string & pupilName, bool create = false)
// {
// 
//    std::string path = basisRoot(dmName, basisName, create);
//    
//    path += "/" + pupilName;
//    
//    if(create)
//    {
//       ioutils::createDirectories(path);
//    }
//    
//    path += "/projectedOrthoModes.fits";
// 
//    return path;
// }

} //namespace dm

//--------------------------------------------------------------------------------

namespace pupil
{
///The root path for pupil files
/** 
  * \param[in] pupilName the name of the basis  
  * \param[in] create [optional] create the directory if it noes not exist.
  * 
  * \returns the root path for pupil files.
  */
std::string root(const std::string & pupilName, bool create = false)
{
   
   std::string path = mx::sys::getEnv("MX_AO_DATADIR");

   path += "/pupil/" + pupilName;

   if(create)
   {
      ioutils::createDirectories(path);
   }

   
   return path;
}
   
   
///The path for the pupil FITS file.
/** 
  * \param[in] pupilName the name of the pupil file 
  * \param[in] create [optional] create the root directory if it noes not exist.
  * 
  * \returns the path to the FITS file containing the pupil as a 1/0 mask
  */
std::string pupilFile(const std::string & pupilName, bool create = false)
{
   std::string path = root(pupilName, create) + "/pupil.fits";

   return path;
}   
   
}// namespace pupil

//--------------------------------------------------------------------------------

namespace sys 
{
   
   ///The root path for system files
   /** 
     * \param[in] sysName the name of the system 
     * \param[in] create [optional] create the directory if it noes not exist.
     * 
     * \returns the root path for calibration files.
     */
   std::string root(const std::string & sysName, bool create = false)
   {
      std::string path = mx::sys::getEnv("MX_AO_DATADIR");

      path += "/system/" + sysName;

      if(create)
      {
         ioutils::createDirectories(path);
      }

      return path;
   }

   namespace cal
   {
   
      ///The root path for system calibration files
      /** 
        * \param[in] sysName the name of the system 
        * \param[in] create [optional] create the directory if it noes not exist.
        * 
        * \returns the root path for calibration files.
        */
      std::string root(const std::string & sysName, bool create = false)
      {

         std::string path = mx::AO::path::sys::root(sysName, create);

         path += "/cal";

         if(create)
         {
            ioutils::createDirectories(path);
         }
   
         return path;
      }

      ///Path for the system response calibration
      std::string sysResp( const std::string & sysName,
                           const std::string & dmName,
                           const std::string & wfsName, 
                           const std::string & pupilName, 
                           const std::string & basisName,
                           bool create = false)
      {
         std::string path = root(sysName, create);
   
         path += "/" + dmName + "/" + wfsName + "/" + pupilName + "/" + basisName;
   
         if(create)
         {
            ioutils::createDirectories(path);
         }

         return path;
      }


      ///Path for the system response calibration
      std::string rMat( const std::string & sysName, 
                        const std::string & dmName,
                        const std::string & wfsName, 
                        const std::string & pupilName,
                        const std::string & basisName,
                        const std::string & id,
                        bool create = false)
      {

         std::string path = sysResp(sysName, dmName, wfsName, pupilName, basisName, create);
   
         path += "/";
         path += "rMat_";
   
         path += id;
         path += ".fits";
   
         return path;
      }

      ///Path for the system response calibration
      std::string rImages( const std::string & sysName, 
                           const std::string & dmName,
                           const std::string & wfsName, 
                           const std::string & pupilName,
                           const std::string & basisName,
                           const std::string & id,
                           bool create = false)
      {

         std::string path = sysResp(sysName, dmName, wfsName, pupilName, basisName, create);
   
         path += "/";
         path += "rImages_";
   
         path += id;
         path += ".fits";
   
         return path;
      }

      ///Path for the system response interaction matrix
      std::string iMat( const std::string & sysName, 
                        const std::string & dmName,
                        const std::string & wfsName, 
                        const std::string & pupilName,
                        const std::string & basisName,
                        const std::string & id,
                        bool create = false)
      {

         std::string path = sysResp(sysName, dmName, wfsName, pupilName, basisName, create);
   
         path += "/";

         path += "iMat_";
   
         path += id;
         path += ".fits";
   
         return path;
      }  

      ///Path for the system response interaction matrix
      std::string U( const std::string & sysName, 
                     const std::string & dmName,
                     const std::string & wfsName, 
                     const std::string & pupilName,
                     const std::string & basisName,
                     const std::string & id,
                     bool create = false)
      {

         std::string path = sysResp(sysName, dmName, wfsName, pupilName, basisName, create);
   
         path += "/";

         path += "U_";
   
         path += id;
         path += ".fits";
   
         return path;
      }  
      
      ///Path for the system response interaction matrix
      std::string S( const std::string & sysName, 
                     const std::string & dmName,
                     const std::string & wfsName, 
                     const std::string & pupilName,
                     const std::string & basisName,
                     const std::string & id,
                     bool create = false)
      {

         std::string path = sysResp(sysName, dmName, wfsName, pupilName, basisName, create);
   
         path += "/";

         path += "S_";
   
         path += id;
         path += ".fits";
   
         return path;
      }
      
      ///Path for the system response interaction matrix
      std::string VT( const std::string & sysName, 
                      const std::string & dmName,
                      const std::string & wfsName, 
                      const std::string & pupilName,
                      const std::string & basisName,
                      const std::string & id,
                      bool create = false)
      {

         std::string path = sysResp(sysName, dmName, wfsName, pupilName, basisName, create);
   
         path += "/";

         path += "VT_";
   
         path += id;
         path += ".fits";
   
         return path;
      }
      
      std::string siGains( const std::string & sysName,
                           const std::string & dmName,
                           const std::string & wfsName, 
                           const std::string & pupilName,
                           const std::string & basisName,
                           double mag,
                           bool create = false
                         )
      {
         std::string path = root(sysName, create);
         path += "/" + dmName + "/" + wfsName + "/" + pupilName + "/" + basisName + "/gains/si";

         if(create)
         {
            ioutils::createDirectories(path);
         }

         path += "/optg_" + std::to_string(mag) + "mag.fits";
         
         return path;
      }
      
      std::string lpGains( const std::string & sysName,
                           const std::string & dmName,
                           const std::string & wfsName, 
                           const std::string & pupilName,
                           const std::string & basisName,
                           double mag,
                           int lpNc,
                           bool create = false
                         )
      {
         std::string path = root(sysName, create);
         path += "/" + dmName + "/" + wfsName + "/" + pupilName + "/" + basisName + "/gains/lp";

         if(create)
         {
            ioutils::createDirectories(path);
         }

         path += "/optg_" + std::to_string(mag) + "mag_lpNc_" + std::to_string(lpNc) + ".fits";
         
         return path;
      }
       
      std::string lpCoeff( const std::string & sysName,
                           const std::string & dmName,
                           const std::string & wfsName, 
                           const std::string & pupilName,
                           const std::string & basisName,
                           double mag,
                           int lpNc,
                           bool create = false
                         )
      {
         std::string path = root(sysName, create);
         path += "/" + dmName + "/" + wfsName + "/" + pupilName + "/" + basisName + "/gains/lp";

         if(create)
         {
            ioutils::createDirectories(path);
         }

         path += "/lpc_" + std::to_string(mag) + "mag_lpNc_" + std::to_string(lpNc) + ".fits";
         
         return path;
      }
   }// namespace cal

} //namespace sys 

//--------------------------------------------------------------------------------

} //namespace path



} //namespace AO

} //namespace mx


#endif //__mxao_hpp__
