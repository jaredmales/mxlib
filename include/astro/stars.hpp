/** \file stars.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Various utilities related to stars.
  * \ingroup astrofiles
  * 
  */

#include "../gslInterpolation.hpp"

#include "units.hpp"

#ifndef __mx_astro_stars_hpp__
#define __mx_astro_stars_hpp__

namespace mx
{
namespace astro
{



/** \ingroup stars
  * @{
  */
///Parse a main sequence spectral type string into a numeric code.
/** Expects standard spectral type strings such as "O5V" or "F2.5" or "G8.5V".  
  * The numeric code starts at 0 for "O0V", 10 for "A0V", through 90 for "Y0".  The subtypes are added to this value.  For instance,
  * "G8.5V" becomes 48.5. 
  * 
  * Only works for main sequence types.  Any other spectral type, such as MK class I-IV, will return -1.
  * 
  */ 
template<typename realT=float>
realT numSpType( std::string spType /**< [in] The spectral type string to parse.*/)
{
   
   spType = removeWhiteSpace(spType);

   spType = toUpper(spType);
   
   
   realT num = 0;
   
   switch( spType[0] )
   {
      case 'O':
         num = 0;
         break;
      case 'B':
         num = 10;
         break;
      case 'A':
         num = 20;
         break;
      case 'F':
         num = 30;
         break;
      case 'G':
         num = 40;
         break;
      case 'K':
         num = 50;
         break;
      case 'M':
         num = 60;
         break;
      case 'L':
         num = 70;
         break;
      case 'T':
         num = 80;
         break;
      case 'Y':
         num = 90;
         break;
      default:
         return -1;
   }
   
   if(spType.size() < 2) return -1;
   if(!isdigit(spType[1])) return -1;
   
   int i=1;
   while( i < spType.size() && ( isdigit(spType[i]) || spType[i] == '.')) ++i;
   
   std::string subType = spType.substr(1, i-1);
   
   num += convertFromString<realT>(subType);
   
   if(spType.size() == i) return num;
   if(spType[i] == 'V') return num;
   
   return -1;
}

/// Provide various characteristics of main sequence stars according to their spectral type.
/** Interpolates on "A Modern Mean Dwarf Stellar Color and Effective Temperature Sequence"
  * by Eric Mamajek, available at  http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt.
  * 
  * Loads values of the sequence at construction (hard coded), and then provides interpolation in between the points.
  * 
  * Version of the sequence used: 2017.03.14
  */
template<typename realT>
struct mainSequence
{
   std::vector<double> numTypes; ///< The numeric spectral type codes from the sequence
   std::vector<double> Teffs; ///< The effective temperatures from the sequence
   std::vector<double> logLs; ///< The log10 luminosities from the sequence
   std::vector<double> Mvs; ///< The absolute V magnitudes from the sequence
   std::vector<double> V_Rcs; ///< The V-Rc colors from the sequence
   std::vector<double> V_Ics; ///< The V-Ic colors from the sequence
   std::vector<double> V_Kss; ///< The V-Ks colors from the sequence
   std::vector<double> J_Hs; ///< The J-H colors from the sequence
   std::vector<double> H_Kss; ///< The H-Ks colors from the sequence
   
   gslInterpolator interpT; ///< The interpolator for effective Temperature
   double minT; ///< The minimum numeric spectral type for effective Temperature
   double maxT; ///< The maximum numeric spectral type for effective Temperature
   
   gslInterpolator interpL; ///< The interpolator for log luminosity
   double minL; ///< The minimum numeric spectral type for log luminosity
   double maxL; ///< The maximum numeric spectral type for log luminosity
   
   gslInterpolator interpMv; ///< The interpolator for absolute V magnitude
   double minMv; ///< The minimum numeric spectral type for absolute V magnitude
   double maxMv; ///< The maximum numeric spectral type for absolute V magnitude
   
   gslInterpolator interpVRc; ///< The interpolator for V-Rc color
   double minVRc; ///< The minimum numeric spectral type for V-Rc color
   double maxVRc; ///< The maximum numeric spectral type for V-Rc color
   
   gslInterpolator interpVIc; ///< The interpolator for V-Ic color
   double minVIc; ///< The minimum numeric spectral type for V-Ic color
   double maxVIc; ///< The maximum numeric spectral type for V-Ic color
   
   gslInterpolator interpVKs; ///< The interpolator for V-Ks color
   double minVKs; ///< The minimum numeric spectral type for V-Ks color
   double maxVKs; ///< The maximum numeric spectral type for V-Ks color
   
   gslInterpolator interpJH; ///< The interpolator for J-H color
   double minJH; ///< The minimum numeric spectral type for J-H color
   double maxJH; ///< The maximum numeric spectral type for J-H color
   
   gslInterpolator interpHKs; ///< The interpolator for H-Ks color
   double minHKs; ///< The minimum numeric spectral type for H-Ks color
   double maxHKs; ///< The maximum numeric spectral type for H-Ks color
   
   mainSequence()
   {
      numTypes = {3, 4, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 14, 15, 16, 17, 18, 19, 19.5, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 61.5, 62, 62.5, 63, 63.5, 64, 64.5, 65, 65.5, 66, 66.5, 67, 67.5, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 84.5, 85, 85.5, 86, 87, 87.5, 88, 88.5, 89, 90, 92};
      
      Teffs = {46000, 43000, 41500, 40000, 39000, 37300, 36500, 35000, 34500, 33000, 32500, 32000, 31500, 29000, 26000, 24500, 20600, 18500, 17000, 16700, 15700, 14500, 14000, 12500, 10700, 10400, 9700, 9200, 8840, 8550, 8270, 8080, 8000, 7800, 7500, 7440, 7200, 7030, 6810, 6720, 6640, 6510, 6340, 6240, 6170, 6040, 5920, 5880, 5770, 5720, 5680, 5660, 5590, 5530, 5490, 5340, 5280, 5170, 5040, 4840, 4620, 4450, 4200, 4050, 3970, 3880, 3850, 3700, 3650, 3550, 3500, 3410, 3250, 3200, 3100, 3030, 3000, 2850, 2710, 2650, 2600, 2500, 2450, 2250, 2100, 1960, 1830, 1700, 1590, 1490, 1410, 1350, 1300, 1260, 1230, 1200, 1160, 1120, 1090, 1050, 1010, 960, 840, 770, 700, 610, 530, 420, 250};
      
      interpT.setup( gsl_interp_linear, numTypes, Teffs );
      minT = 3.0;
      maxT = 92.0;
      
      logLs = {5.8, 5.67, 5.58, 5.46, 5.37, 5.24, 5.17, 5.05, 4.99, 4.87, 4.82, 4.76, 4.65, 4.47, 4.13, 3.89, 3.38, 3.1, 2.96, 2.91, 2.8, 2.54, 2.49, 2.27, 1.79, 1.73, 1.54, 1.41, 1.33, 1.29, 1.2, 1.16, 1.14, 1.06, 0.97, 0.97, 0.89, 0.77, 0.7, 0.67, 0.61, 0.54, 0.43, 0.36, 0.31, 0.26, 0.14, 0.13, 0.01, -0.01, -0.04, -0.05, -0.11, -0.12, -0.17, -0.25, -0.33, -0.37, -0.47, -0.58, -0.71, -0.75, -0.92, -1.04, -1.08, -1.22, -1.26, -1.42, -1.47, -1.57, -1.68, -1.78, -2.07, -2.2, -2.31, -2.52, -2.79, -3.02, -3.09, -3.21, -3.29, -3.34, -3.53, -3.57, -3.7, -3.84, -3.98, -4.11, -4.24, -4.36, -4.46, -4.55, -4.61, -4.66, -4.69, -4.73, -4.77, -4.84, -4.9, -4.95, -5.04, -5.12, -5.37, -5.54, -5.71, -5.93, -6.15, -6.52, -8.6};

      interpL.setup( gsl_interp_linear, numTypes, logLs );
      minL = 3.0;
      maxL = 92.0;
      
      Mvs = {-5.7, -5.5, -5.4, -5.2, -5.1, -4.9, -4.8, -4.6, -4.5, -4.3, -4.2, -4.1, -4, -3.6, -3.1, -2.8, -1.7, -1.4, -1.1, -1, -0.9, -0.5, -0.4, -0.2, 0.7, 0.8, 1.11, 1.34, 1.48, 1.55, 1.76, 1.84, 1.89, 2.07, 2.29, 2.3, 2.51, 2.79, 2.99, 3.08, 3.23, 3.4, 3.7, 3.87, 4.01, 4.15, 4.45, 4.5, 4.79, 4.86, 4.94, 4.98, 5.13, 5.18, 5.32, 5.55, 5.76, 5.91, 6.19, 6.57, 7.04, 7.25, 7.88, 8.3, 8.58, 9, 9.16, 9.8, 9.97, 10.3, 10.7, 11.14, 12.19, 12.8, 13.57, 14.3, 15.51, 16.62, 17.07, 17.81, 18.42, 18.88, 19.26, 19.8, 20.3, 20.8, 21.3, 22.2, 23.1, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};

      interpMv.setup( gsl_interp_linear, numTypes, Mvs );
      minMv = 3.0;
      maxMv = 75.0;
      
      V_Rcs = {-99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -0.115, -0.114, -0.094, -0.087, -0.08, -0.074, -0.07, -0.062, -0.058, -0.048, -0.028, -0.017, 0.001, 0.019, 0.042, 0.05, 0.078, 0.089, 0.094, 0.117, 0.14, 0.143, 0.166, 0.19, 0.213, 0.222, 0.236, 0.252, 0.276, 0.29, 0.3, 0.312, 0.336, 0.34, 0.363, 0.368, 0.374, 0.377, 0.388, 0.393, 0.404, 0.423, 0.443, 0.46, 0.487, 0.544, 0.64, 0.671, 0.771, 0.826, 0.857, 0.902, 0.913, 0.968, 0.978, 1.001, 1.041, 1.079, 1.178, 1.241, 1.345, 1.446, 1.656, 1.95, 2.003, 2.18, 2.16, 2.15, 1.89, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};

      interpVRc.setup( gsl_interp_linear, numTypes, V_Rcs );
      minVRc = 11.0;
      maxVRc = 69.0;
      
      V_Ics = {-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -0.369, -0.361, -0.355, -0.338, -0.325, -0.281, -0.23, -0.21, -0.192, -0.176, -0.165, -0.145, -0.133, -0.108, -0.061, -0.035, 0.004, 0.044, 0.091, 0.108, 0.164, 0.186, 0.197, 0.242, 0.288, 0.294, 0.339, 0.385, 0.432, 0.449, 0.476, 0.506, 0.553, 0.579, 0.599, 0.62, 0.664, 0.672, 0.713, 0.722, 0.733, 0.738, 0.758, 0.766, 0.786, 0.82, 0.853, 0.883, 0.929, 1.025, 1.19, 1.246, 1.448, 1.58, 1.664, 1.808, 1.848, 2.051, 2.089, 2.173, 2.306, 2.42, 2.68, 2.831, 3.073, 3.277, 3.664, 4.1, 4.284, 4.52, 4.56, 4.6, 4.37, 4.6, 4.8, 4.9, 5, 5.5, 6, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};

      interpVIc.setup( gsl_interp_linear, numTypes, V_Ics );
      minVIc = 9.0;
      maxVIc = 79.0;
      
      V_Kss = {-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -1, -0.977, -0.958, -0.913, -0.874, -0.752, -0.602, -0.544, -0.492, -0.447, -0.417, -0.358, -0.325, -0.254, -0.121, -0.048, 0.041, 0.101, 0.188, 0.228, 0.353, 0.403, 0.428, 0.528, 0.626, 0.638, 0.732, 0.828, 0.925, 0.961, 1.017, 1.079, 1.185, 1.244, 1.29, 1.34, 1.44, 1.458, 1.564, 1.59, 1.621, 1.635, 1.691, 1.712, 1.768, 1.861, 1.953, 2.034, 2.155, 2.41, 2.733, 2.835, 3.19, 3.418, 3.544, 3.737, 3.79, 4.065, 4.12, 4.24, 4.43, 4.6, 5, 5.25, 5.64, 5.94, 6.5, 7.3, 7.6, 8.05, 8.45, 8.73, 8.85, 9.3, 9.6, 9.9, 10.2, 10.8, 11.4, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};

      interpVKs.setup( gsl_interp_linear, numTypes, V_Kss );
      minVKs = 19.0;
      maxVKs = 79.0;
      
      J_Hs = {-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -0.164, -0.161, -0.159, -0.153, -0.148, -0.132, -0.113, -0.105, -0.098, -0.092, -0.089, -0.081, -0.077, -0.067, -0.05, -0.044, -0.032, -0.024, -0.01, -0.002, 0.022, 0.031, 0.036, 0.055, 0.075, 0.078, 0.098, 0.119, 0.14, 0.147, 0.159, 0.173, 0.199, 0.213, 0.225, 0.237, 0.262, 0.267, 0.293, 0.299, 0.307, 0.31, 0.324, 0.329, 0.342, 0.365, 0.387, 0.405, 0.432, 0.49, 0.544, 0.56, 0.605, 0.624, 0.631, 0.624, 0.622, 0.61, 0.607, 0.6, 0.589, 0.579, 0.558, 0.557, 0.564, 0.58, 0.588, 0.605, 0.609, 0.613, 0.65, 0.677, 0.749, 0.79, 0.8, 0.87, 1, 1.14, 1.13, 1.08, 1.1, 1.14, 1.1, 1.02, 1.02, 0.86, 0.68, 0.35, 0.2, 0.2, 0.2, 0.1, 0, 0.2, 0.2, 0.2, 0.1, -99, -99};

      interpJH.setup( gsl_interp_linear, numTypes, J_Hs );
      minJH = 19.0;
      maxJH = 89.0;
      
      H_Kss = {-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -0.071, -0.069, -0.067, -0.063, -0.059, -0.047, -0.032, -0.026, -0.021, -0.016, -0.013, -0.007, -0.004, 0.003, 0.016, 0.021, 0.028, 0.031, 0.034, 0.034, 0.037, 0.038, 0.038, 0.04, 0.042, 0.043, 0.045, 0.047, 0.05, 0.051, 0.052, 0.054, 0.057, 0.06, 0.061, 0.063, 0.067, 0.068, 0.073, 0.074, 0.075, 0.076, 0.079, 0.08, 0.082, 0.087, 0.091, 0.094, 0.099, 0.11, 0.125, 0.13, 0.151, 0.169, 0.18, 0.198, 0.203, 0.225, 0.228, 0.234, 0.244, 0.252, 0.269, 0.282, 0.301, 0.311, 0.329, 0.352, 0.364, 0.386, 0.422, 0.447, 0.481, 0.5, 0.54, 0.57, 0.63, 0.63, 0.65, 0.64, 0.62, 0.63, 0.63, 0.54, 0.45, 0.27, 0.08, -0.19, -0.06, -0.08, -0.1, -0.03, 0, -0.05, -0.05, -99, -0.2, -0.5, -99};

      interpHKs.setup( gsl_interp_linear, numTypes, H_Kss );
      minHKs = 19.0;
      maxHKs = 88.0;
      
   }
   
   /// Get the interpolated effective temperature
   /**
     * \returns the effective temperature. 
     * \returns -999 if the input type code is out of range.
     */ 
   realT Teff( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minT || numType > maxT) 
      {
         return -999;
      }
      
      return interpT.interpolate(numType);
   }
   
   /// Get the interpolated effective temperature
   /**
     * \returns the effective temperature. 
     * \returns -999 if the input type is out of range.
     */ 
   realT Teff( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return Teff(  numSpType(spType));
   }
   
   /// Get the interpolated log of luminosity
   /**
     * \returns the log of luminosity. 
     * \returns -999 if the input type code is out of range.
     */ 
   realT logL( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minL || numType > maxL) 
      {
         return -999;
      }
      
      return interpL.interpolate(numType);
   }
   
   /// Get the interpolated log of luminosity
   /**
     * \returns the log of luminosity. 
     * \returns -999 if the input type is out of range.
     */
   realT logL( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return logL(  numSpType(spType));
   }
   
   /// Get the interpolated absolute V magnitude
   /**
     * \returns the aboslute V magnitude
     * \returns -999 if the input type code is out of range.
     */ 
   realT Mv( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minMv || numType > maxMv) 
      {
         return -999;
      }
      
      return interpMv.interpolate(numType);
   }
   
   /// Get the interpolated absolute V magnitude
   /**
     * \returns the aboslute V magnitude
     * \returns -999 if the input type is out of range.
     */ 
   realT Mv( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return Mv(  numSpType(spType));
   }
   
   /// Get the interpolated absolute V-Rc color
   /**
     * \returns the aboslute V-Rc color
     * \returns -999 if the input type code is out of range.
     */ 
   realT V_Rc( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minVRc || numType > maxVRc) 
      {
         return -999;
      }
      
      return interpVRc.interpolate(numType);
   }
   
   /// Get the interpolated absolute V-Rc color
   /**
     * \returns the aboslute V-Rc color
     * \returns -999 if the input type is out of range.
     */ 
   realT V_Rc( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return V_Rc(  numSpType(spType));
   }
   
   /// Get the interpolated absolute V-Ic color
   /**
     * \returns the aboslute V-Ic color
     * \returns -999 if the input type code is out of range.
     */ 
   realT V_Ic( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minVIc || numType > maxVIc) 
      {
         return -999;
      }
      
      return interpVIc.interpolate(numType);
   }
   
   /// Get the interpolated absolute V-Ic color
   /**
     * \returns the aboslute V-Ic color
     * \returns -999 if the input type is out of range.
     */
   realT V_Ic( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return V_Ic(  numSpType(spType));
   }
   
   /// Get the interpolated absolute V-Ks color
   /**
     * \returns the aboslute V-Ks color
     * \returns -999 if the input type code is out of range.
     */ 
   realT V_Ks( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minVKs || numType > maxVKs) 
      {
         return -999;
      }
      
      return interpVKs.interpolate(numType);
   }
   
   /// Get the interpolated absolute V-Ks color
   /**
     * \returns the aboslute V-Ks color
     * \returns -999 if the input type is out of range.
     */ 
   realT V_Ks( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return V_Ks(  numSpType(spType));
   }
   
   /// Get the interpolated absolute J-H color
   /**
     * \returns the aboslute J-H color
     * \returns -999 if the input type code is out of range.
     */ 
   realT J_H( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minJH || numType > maxJH) 
      {
         return -999;
      }
      
      return interpJH.interpolate(numType);
   }
   
   /// Get the interpolated absolute J-H color
   /**
     * \returns the aboslute J-H color
     * \returns -999 if the input type is out of range.
     */
   realT J_H( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return J_H(  numSpType(spType));
   }
   
   /// Get the interpolated absolute H-Ks color
   /**
     * \returns the aboslute H-Ks color
     * \returns -999 if the input type code is out of range.
     */ 
   realT H_Ks( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < minHKs || numType > maxHKs) 
      {
         return -999;
      }
      
      return interpHKs.interpolate(numType);
   }
   
   /// Get the interpolated absolute H-Ks color
   /**
     * \returns the aboslute H-Ks color
     * \returns -999 if the input type is out of range.
     */
   realT H_Ks( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return H_Ks(  numSpType(spType));/// Get the interpolated absolute H-Ks color
   /**
     * \returns the aboslute H-Ks color
     * \returns -999 if the input type code is out of range.
     */
   }
};


///@}
} //namespace astro
} //namespace mx
#endif
