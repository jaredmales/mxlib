/** \file stars.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Various utilities related to stars.
  * \ingroup astrofiles
  * 
  */

#include "../math/gslInterpolator.hpp"

#include "units.hpp"

#ifndef mx_astro_stars_hpp
#define mx_astro_stars_hpp

#include "../ioutils/stringUtils.hpp"
#include "../ioutils/readColumns.hpp"

namespace mx
{
namespace astro
{

///Parse a main sequence spectral type string into a numeric code.
/** Expects standard spectral type strings such as "O5V" or "F2.5" or "G8.5V".  
  * The numeric code starts at 0 for "O0V", 10 for "A0V", through 90 for "Y0".  The subtypes are added to this value.  For instance,
  * "G8.5V" becomes 48.5. 
  * 
  * Only works for main sequence types.  Any other spectral type, such as MK class I-IV, will return -1.
  * 
  * \ingroup stars
  */ 
template<typename realT=float>
realT numSpType( std::string spType /**< [in] The spectral type string to parse.*/)
{
   
   spType = ioutils::removeWhiteSpace(spType);

   spType = ioutils::toUpper(spType);
   
   
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
         return -1000;
   }
   
   if(spType.size() < 2) return -1000;
   if(!isdigit(spType[1])) return -1000;
   
   int i=1;
   while( i < spType.size() && ( isdigit(spType[i]) || spType[i] == '.')) ++i;
   
   std::string subType = spType.substr(1, i-1);
   
   num += ioutils::convertFromString<realT>(subType);
   
   if(spType.size() == i) return num;
   if(spType[i] == 'V') return num;
   
   return -1*num;
}

/// Provide various characteristics of main sequence stars according to their spectral type.
/** Interpolates on "A Modern Mean Dwarf Stellar Color and Effective Temperature Sequence"
  * by Eric Mamajek, available at  http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt.
  * 
  * Loads values of the sequence at construction (hard coded), and then provides interpolation in between the points.
  * 
  * Version of the sequence used: 2017.03.14
  * 
  * \ingroup stars
  */
template<typename realT>
struct mainSequence
{
   std::vector<double> m_numTypes; ///< The numeric spectral type codes from the sequence
   std::vector<double> m_numTypesR; ///< The numeric spectral type codes from the sequence, in reversed order.
   
   std::vector<double> m_Teffs; ///< The effective temperatures from the sequence
   std::vector<double> m_TeffsR; ///< The effective temperatures from the sequence, reveresed order.
   
   std::vector<double> m_logLs; ///< The log10 luminosities from the sequence
   std::vector<double> m_rads; ///< The radii from the sequence
   std::vector<double> m_masses; ///< The masses from the sequence
   std::vector<double> m_Mvs; ///< The absolute V magnitudes from the sequence
   std::vector<double> m_V_Rcs; ///< The V-Rc colors from the sequence
   std::vector<double> m_V_Ics; ///< The V-Ic colors from the sequence
   std::vector<double> m_V_Kss; ///< The V-Ks colors from the sequence
   std::vector<double> m_J_Hs; ///< The J-H colors from the sequence
   std::vector<double> m_H_Kss; ///< The H-Ks colors from the sequence
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpT; ///< The interpolator for effective Temperature
   double m_minT; ///< The minimum numeric spectral type for effective Temperature
   double m_maxT; ///< The maximum numeric spectral type for effective Temperature
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpRad; ///< The interpolator for effective Temperature
   double m_minRad; ///< The minimum numeric spectral type for effective Temperature
   double m_maxRad; ///< The maximum numeric spectral type for effective Temperature
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpL; ///< The interpolator for log luminosity
   double m_minL; ///< The minimum numeric spectral type for log luminosity
   double m_maxL; ///< The maximum numeric spectral type for log luminosity
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpMv; ///< The interpolator for absolute V magnitude
   double m_minMv; ///< The minimum numeric spectral type for absolute V magnitude
   double m_maxMv; ///< The maximum numeric spectral type for absolute V magnitude
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpVRc; ///< The interpolator for V-Rc color
   double m_minVRc; ///< The minimum numeric spectral type for V-Rc color
   double m_maxVRc; ///< The maximum numeric spectral type for V-Rc color
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpVIc; ///< The interpolator for V-Ic color
   double m_minVIc; ///< The minimum numeric spectral type for V-Ic color
   double m_maxVIc; ///< The maximum numeric spectral type for V-Ic color
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpVKs; ///< The interpolator for V-Ks color
   double m_minVKs; ///< The minimum numeric spectral type for V-Ks color
   double m_maxVKs; ///< The maximum numeric spectral type for V-Ks color
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpJH; ///< The interpolator for J-H color
   double m_minJH; ///< The minimum numeric spectral type for J-H color
   double m_maxJH; ///< The maximum numeric spectral type for J-H color
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpHKs; ///< The interpolator for H-Ks color
   double m_minHKs; ///< The minimum numeric spectral type for H-Ks color
   double m_maxHKs; ///< The maximum numeric spectral type for H-Ks color
   
   
   math::gslInterpolator<math::gsl_interp_linear<double>> interpSpTfmT; ///< The interpolator for effective Temperature
   double m_minSpTfmT; ///< The minimum numeric spectral type for effective Temperature
   double m_maxSpTfmT; ///< The maximum numeric spectral type for effective Temperature
   
   void findMinMax( std::vector<double> & seq,
                    double & min,
                    double & max,
                    std::vector<double> & ref
                  )
   {
      if(seq.size()!=ref.size())
      {
         std::cerr << "size mismatch\n";
         min = 0;
         max = 0;
      }
      
      min = ref[0];
      size_t n;
      for(n=0; n < seq.size();++n)
      {
         if(seq[n] != -99) break;
      }

      if(n >= seq.size()-1)
      {
         min = 0;
         max = 0;
         return;
      }

      min = ref[n];
      
      
      for(;n<seq.size()-1;++n)
      {
         if(seq[n+1] == -99) break;
      }
   
      max = ref[n];
      ++n;
      for(;n<seq.size();++n)
      {
         seq[n] = -99;
      }
      
   }
   
   mainSequence()
   {
      #include "mamajek.inc"
      
      m_numTypesR.assign(m_numTypes.rbegin(), m_numTypes.rend());
      
      findMinMax(m_Teffs, m_minT, m_maxT, m_numTypes);
      m_TeffsR.assign(m_Teffs.rbegin(), m_Teffs.rend()); //get reversed vector after normalizing -99s.
      
      findMinMax(m_rads, m_minRad, m_maxRad, m_numTypes);
      
      findMinMax(m_logLs, m_minL, m_maxL, m_numTypes);
      findMinMax(m_Mvs, m_minMv, m_maxMv, m_numTypes);
      findMinMax(m_V_Rcs, m_minVRc, m_maxVRc, m_numTypes);
      findMinMax(m_V_Ics, m_minVIc, m_maxVIc, m_numTypes);
      findMinMax(m_V_Kss, m_minVKs, m_maxVKs, m_numTypes);
      findMinMax(m_J_Hs, m_minJH, m_maxJH, m_numTypes);
      findMinMax(m_H_Kss, m_minHKs, m_maxHKs, m_numTypes);
         
      findMinMax(m_numTypesR, m_minSpTfmT, m_maxSpTfmT, m_TeffsR);
      
      interpT.setup(  m_numTypes, m_Teffs );
      interpRad.setup( m_numTypes, m_rads);
      interpL.setup( m_numTypes, m_logLs );
      interpMv.setup( m_numTypes, m_Mvs );
      interpVRc.setup( m_numTypes, m_V_Rcs );
      interpVIc.setup( m_numTypes, m_V_Ics );
      interpVKs.setup( m_numTypes, m_V_Kss );
      interpJH.setup( m_numTypes, m_J_Hs );
      interpHKs.setup( m_numTypes, m_H_Kss );
      
      interpSpTfmT.setup( m_TeffsR, m_numTypesR);
   }
   
   /// Get the interpolated effective temperature
   /**
     * \returns the effective temperature. 
     * \returns -999 if the input type code is out of range.
     */ 
   realT Teff( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < m_minT || numType > m_maxT) 
      {
         return -999;
      }
      
      return interpT(numType);
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
   
   /// Get the interpolated radius
   /**
     * \returns the radius in Solar units. 
     * \returns -999 if the input type code is out of range.
     */ 
   realT radius( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < m_minRad || numType > m_maxRad) 
      {
         return -999;
      }
      
      return interpRad(numType);
   }
   
   /// Get the interpolated radius
   /**
     * \returns the radius in Solar units. 
     * \returns -999 if the input type is out of range.
     */ 
   realT radius( const std::string & spType /**< [in] the spectral type in standard format */)
   {
      return radius(  numSpType(spType));
   }
   
   /// Get the interpolated log of luminosity
   /**
     * \returns the log of luminosity. 
     * \returns -999 if the input type code is out of range.
     */ 
   realT logL( realT numType /**< [in] the numeric spectral type code, see \ref mx::astro::numSpType() */)
   {
      if( numType < m_minL || numType > m_maxL) 
      {
         return -999;
      }
      
      return interpL(numType);
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
      if( numType < m_minMv || numType > m_maxMv) 
      {
         return -999;
      }
      
      return interpMv(numType);
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
      if( numType < m_minVRc || numType > m_maxVRc) 
      {
         return -999;
      }
      
      return interpVRc(numType);
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
      if( numType < m_minVIc || numType > m_maxVIc) 
      {
         return -999;
      }
      
      return interpVIc(numType);
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
      if( numType < m_minVKs || numType > m_maxVKs) 
      {
         return -999;
      }
      
      return interpVKs(numType);
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
      if( numType < m_minJH || numType > m_maxJH) 
      {
         return -999;
      }
      
      return interpJH(numType);
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
      if( numType < m_minHKs || numType > m_maxHKs) 
      {
         return -999;
      }
      
      return interpHKs(numType);
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
   
   realT spTFromTeff( realT Teff )
   {
      if( Teff < m_minSpTfmT || Teff > m_maxSpTfmT)
      {
         return -999;
      }
      
      return floor(2*interpSpTfmT(Teff)+0.5)/2; //round to nearest 0.5 types.
   }
   
};


namespace maintenance
{
   

/// Read in the main sequence table of Mamajek, and construct the vectors for input into the mainSequence class.
/** This is used to create mamajek.inc.
  * The table should be copied to a text file, and all `...` replaced with -99, then any remaining . replaced with space.
  */
void makeMSTable( const std::string & fn )
{
   std::vector<std::string> SpT;
   std::vector<double> Teff;  
   std::vector<double> logT;   
   std::vector<double> logL;  
   std::vector<double> Mbol;  
   std::vector<double> BCv;   
   std::vector<double> Mv;    
   std::vector<double> B__V;    
   std::vector<double> Bt__Vt;  
   std::vector<double> G__V;    
   std::vector<double> Bp__Rp;  
   std::vector<double> G__Rp;   
   std::vector<double> M_G;    
   std::vector<double> b__y;    
   std::vector<double> U__B;   
   std::vector<double> V__Rc;   
   std::vector<double> V__Ic;   
   std::vector<double> V__Ks;   
   std::vector<double> J__H;    
   std::vector<double> H__Ks;    
   std::vector<double> Ks__W1;   
   std::vector<double> W1__W2;  
   std::vector<double> W1__W3;  
   std::vector<double> W1__W4;   
   std::vector<double> M_J;   
   std::vector<double> M_Ks;  
   std::vector<double> i__z;  
   std::vector<double> z__Y; 
   std::vector<double> R_Rsun; 
   std::vector<double> Msun;

   ioutils::readColumns(fn,SpT,Teff,logT,logL,Mbol,BCv,Mv,B__V,Bt__Vt,G__V,Bp__Rp,G__Rp,M_G,b__y,U__B,V__Rc,V__Ic,V__Ks,J__H,H__Ks,Ks__W1,W1__W2,W1__W3,
                  W1__W4,M_J,M_Ks,i__z,z__Y,R_Rsun, Msun);

   std::cout << "//" << fn << "\n";
   std::cout << "         m_numTypes = {" << numSpType(SpT[0]);
   for(int n=1;n<SpT.size();++n) std::cout << ',' << numSpType(SpT[n]);
   std::cout << "};\n";
   
   std::cout << "         m_Teffs = {" << Teff[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << Teff[n];
   std::cout << "};\n";
   
   std::cout << "         m_logLs = {" << logL[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << logL[n];
   std::cout << "};\n";
   
   std::cout << "         m_rads = {" << R_Rsun[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << R_Rsun[n];
   std::cout << "};\n";
   
   std::cout << "         m_masses = {" << Msun[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << Msun[n];
   std::cout << "};\n";
   
   std::cout << "         m_Mvs = {" << Mv[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << Mv[n];
   std::cout << "};\n";
   
   std::cout << "         m_V_Rcs = {" << V__Rc[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << V__Rc[n];
   std::cout << "};\n";
   
   std::cout << "         m_V_Ics = {" << V__Ic[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << V__Ic[n];
   std::cout << "};\n";
   
   std::cout << "         m_V_Kss = {" << V__Ks[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << V__Ks[n];
   std::cout << "};\n";
   
   std::cout << "         m_J_Hs = {" << J__H[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << J__H[n];
   std::cout << "};\n";
   
   std::cout << "         m_H_Kss = {" << H__Ks[0];
   for(int n=1;n<SpT.size();++n) std::cout << ',' << H__Ks[n];
   std::cout << "};\n";
}
} //namespace maintenance
} //namespace astro
} //namespace mx



#endif
