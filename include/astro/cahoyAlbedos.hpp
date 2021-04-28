/** \file cahoyAlbedos.hpp
  * \author Jared R. Males
  * \brief Type definitions for the Cahoy et al 2010 albedo spectra
  * \ingroup astrophot
  *
  */

#ifndef mx_astro_cahoyAlbedos_hpp
#define mx_astro_cahoyAlbedos_hpp

#include "units.hpp"
#include "../ioutils/readColumns.hpp"

namespace mx
{
namespace astro 
{

/// An albedo spectrum directly from the Cahoy et al (2010) \cite cahoy_2010 grid.
/** 
  * \ingroup astrophot_spectra
  */
template<typename _units>
struct cahoySpectrumRaw
{
   typedef _units units;
   typedef typename units::realT realT;
   
   typedef struct
   {
      realT sep; ///< The separation, in AU.  Can be 0.8, 2.0, 5.0, or 10.0.
      int metal; ///< The metallicty.  Can be 1, 3, 10, or 30.
      realT phase; ///< The phase. Must be one of the values in the grid: 0.0, 0.1745, 0.3491, 0.5236, 0.6981, 0.8727, 1.0472, 1.2217, 1.3963,1.5708,1.7453,1.9199, 2.0944, 2.2689, 2.4435, 2.6180, 2.7925,2.9671, 3.139.
   }
   paramsT;
   
   static const bool freq = false;
   
   ///Convert from um to SI m
   static constexpr realT wavelengthUnits = static_cast<realT>(1e6);
   
   ///Geometric albedos are dimensionless.
   static constexpr realT fluxUnits = static_cast<realT>(1); 
   
   static constexpr const char * dataDirEnvVar = "CAHOYALBEDO_DATADIR";
   
   ///The file name is constructed from the paramerters sep, metal, and phase.
   static std::string fileName( const paramsT & params )
   {
      std::string fname; 
   
      if(params.sep == 0.8)
      {
         fname = "0.8";
      }
      else if(params.sep == 2.0)
      {
         fname = "2.0";
      }
      else if(params.sep == 5.0)
      {
         fname = "5";
      }
      else if(params.sep == 10.0)
      {
         fname = "10";
      }
      else
      {
         mxError("cahoyAlbedos::fileName", MXE_INVALIDARG, "invalid separation (params.sep)");
         return "";
      }

      fname += "au_";
      fname += ioutils::convertToString<int>(params.metal);
      fname += "x_albedos/00_Spectral_Albedos   ";  
   
      char pstr[9];
      snprintf(pstr, 9, "%0.5f", params.phase);
   
      fname+=pstr;
   
      fname += ".txt";
      
      return fname;
   
   }
   
   ///Read the spectrum from the file specified by path.  Extra columns are discarded.
   static int readSpectrum( std::vector<realT> & rawLambda,
                            std::vector<realT> & rawSpectrum,
                            const std::string & path,
                            const paramsT & params
                          )
   {
      ioutils::skipCol sk;
      
      return ioutils::readColumns(path, sk, rawLambda, sk, rawSpectrum);
   }
   
};

/// A class to manage interpolation in separation and phase on the albedo spectrum grid from Cahoy et al (2010) \cite cahoy_2010.
/** Usage:
  \code
  std::vector<double> lambda;
  math::vectorScale(lambda, 100000, 4e-12, 600e-9); //Construct the wavelength scale
  
  cahoyGrid<units::si<double>> grid;
  grid.loadGrid(3, lambda); //Open the grid for 3x metallicity, pre-interpolating on the wavelength scale.
  
  grid.setParameters({1.0, 0.1}); //Set the separation to 1.0 AU, and the phase angle to 0.1 radians.
  grid.setSpectrum(lambda); //Perform the interpolation.  Note that lambda is only used for a size check here.
    
  \endcode
  * After the above steps you now have the Ag spectrum at 1.0 AU and 0.1 radians phase on the lambda scale.
  *
  * \ingroup astrophot
  */
template<typename _units>
struct cahoyGrid : public baseSpectrum<typename _units::realT>
{
   typedef _units units;
   typedef typename units::realT realT;
   
   typedef struct
   {
      realT sep;
      realT phase;
   } paramsT;
   
   paramsT _params;
   
   std::vector<realT> _sep;
   std::vector<realT> _phase;
   
   std::vector<realT> _lambda;
   std::vector<std::vector<std::vector<realT>>> _ag;
      
   ///Load the grid into memory for one of the metallicities and pre-interpolate onto the wavelength scale. 
   int loadGrid( int metal, ///< [in] the metalicity, can be 1, 3, 10 or 30. 
                 std::vector<realT> & lambda ///< [in] the wavelength scale.
               )
   {
      _sep = {0.8, 2.0, 2.0, 5.0, 10.0}; //2 is in there twice to prepare for transition distance.
      _phase = {0.0, 0.1745, 0.3491, 0.5236, 0.6981, 0.8727, 1.0472, 1.2217, 1.3963,1.5708,1.7453,1.9199, 2.0944, 2.2689, 2.4435, 2.6180, 2.7925,2.9671, 3.139};
      
      _ag.resize( _sep.size());
      for(int i=0; i< _ag.size(); ++i) _ag[i].resize(_phase.size());
      
      astroSpectrum<cahoySpectrumRaw<units>> rawSpect;
      
      for(int i=0; i< _sep.size(); ++i)
      {
         for(int j=0; j < _phase.size(); ++j)
         {
            rawSpect.setParameters({_sep[i], metal, _phase[j]});            
            if( rawSpect.setSpectrum(lambda) < 0)
            {
               mxError("cahoGrid::loadGrid", MXE_FILERERR, "Error reading albedo spectrum.");
               return -1;
            }
            
            _ag[i][j] = rawSpect._spectrum;
         }
      }

      _sep[1] = 1.0; //The transition from cloudy to not-cloudy occurs here.
      
      this->_spectrum.resize(lambda.size());
      return 0;
   }
   
   
   ///Set the separatio and phase of the spectrum.
   int setParameters( const paramsT & params /**< [in] Struct containting the separation and phase of the spectrum */)
   {
      _params = params;
      
      return 0;
   }
   
   ///Get an interpolated spectrum at arbitrary non-grid point using bilinear interpolation.
   /** If outside the endpoints of the grid, the grid is simply extended (i.e. constant albedo).
     */
   int setSpectrum( std::vector<realT> & lambda /**< [in] the wavelength grid.  Must be the same as used for construction and/or openGrid*/)
   {
      int i_l, i_u, j_l, j_u;
      realT phase, sep;
      
      
      if(lambda.size() != this->_spectrum.size())
      {
         mxError("cahoyGrid::setSpectrum", MXE_INVALIDARG, "wavelength grid (lambda) is not the same size as used for openGrid, or loadGrid not yet called.");
         return -1;
      }
      
      //Normalize phase
      phase = fmod(_params.phase, math::pi<realT>());
      if(phase < 0) phase += math::pi<realT>();
      
      sep = _params.sep;
      
      i_l = _sep.size()-1;
      while( _sep[i_l] > sep && i_l > 0) --i_l;
      
      i_u = 0;
      while( i_u < _sep.size()-1) 
      {
         if(_sep[i_u] >= sep) break;
         ++i_u;
      }
      
      j_l = _phase.size()-1;
      while( j_l > 0) 
      {
         if(_phase[j_l] <= phase) break;
         --j_l;
      }
      
      j_u = 0;
      while( j_u < _phase.size()-1) 
      {
         if(_phase[j_u] >= phase) break;
         ++j_u;
      }
      
      realT x = sep - _sep[i_l];
      if(x < 0) x = 0;
      
      realT y = (phase - _phase[j_l]);
      if(y < 0) y = 0;
      
      realT x0, x1;
            
      for(int i=0; i< this->_spectrum.size(); ++i)
      {
         x0 = _ag[i_l][j_l][i];
         x1 = _ag[i_u][j_l][i];

         if(y != 0 && j_u != j_l)
         {
            x0 += (_ag[i_l][j_u][i] - _ag[i_l][j_l][i])*y/(_phase[j_u] - _phase[j_l]);
            x1 += (_ag[i_u][j_u][i] - _ag[i_u][j_l][i])*y/(_phase[j_u] - _phase[j_l]);
         }
         if( x == 0 || i_u==i_l ) this->_spectrum[i] = x0;
         else this->_spectrum[i] = x0 + (x1-x0)*x/( _sep[i_u] - _sep[i_l]);
      }

      return 0;
   }

};

} //namespace mx

} //namespace astro

#endif //mx_astro_cahoyAlbedos_hpp
