/** \file aoAtmosphere.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to specify atmosphere parameters.
  * \ingroup mxAO_files
  * 
  */

#ifndef aoAtmosphere_hpp
#define aoAtmosphere_hpp



#include <iostream>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>


#include "../../mxlib.hpp"
#include "../../mxError.hpp"

#include "aoConstants.hpp"

#include "../../math/constants.hpp"

#include "../../app/appConfigurator.hpp"

namespace mx
{
namespace AO
{
namespace analysis 
{
   
///A class to specify atmosphere parameters and perform related calculations.
/** 
  * \todo Add layer outer scales.
  * 
  * \tparam realT is the real floating type in which all calculations are performed.
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename _realT>
class aoAtmosphere
{
public:
   
   typedef _realT realT; ///< The real floating type in which all calculations are performed.
   
   ///Constructor
   aoAtmosphere();

protected:
   realT m_r_0; ///< Fried's parameter, in m
   
   realT m_lam_0; ///<Wavelength of Fried's parameter, in m

   std::vector<realT> m_L_0; ///< The outer scale, in m
   
   std::vector<realT> m_l_0; ///< The inner scale of each layer, in m
   
   std::vector<realT> m_layer_z; ///< Vector of layer heights, in m, above the observatory.
   
   realT m_h_obs; ///< Height of the observatory above sea level, in m.  
   
   realT m_H; ///< The atmospheric scale height, in m.
   
   std::vector<realT> m_layer_Cn2; ///< Vector of layer strengths. 
   
   std::vector<realT> m_layer_v_wind; ///< Vector of layer wind speeds, in m/s.
   
   std::vector<realT> m_layer_dir; ///<Vector of layer wind directions, in radians.

   bool m_v_wind_updated; ///< whether or not m_v_wind has been updated after changes
   
   realT m_v_wind; ///< \f$ C_n^2 \f$ averaged windspeed
   
   realT m_dir_wind; ///< \f$ C_n^2 \f$ averaged direction
   
   bool m_z_mean_updated; ///< whether or not m_z_mean has been update after changes
   
   realT m_z_mean; ///< \f$ C_n^2 \f$ averaged layer height
   
   bool m_nonKolmogorov {false};
   
   realT m_alpha{2};
   
   realT m_beta{1};
   
   /// Checks if layer vectors have consistent length.
   /**
     * \returns 0 if all layer vectors are the same length
     * \returns -1 if not, and prints an error.
     */ 
   int checkLayers();
   
public:
   
   ///Get the value of Fried's parameter r_0 at the reference wavelength lam_0.
   /** 
     * \returns the curret value of m_r_0, in m.
     */ 
   realT r_0();
   
   ///Get the value of Fried's parameter r_0 at the specified wavelength.
   /**
     * \note This does not set the value of r_0.
     * 
     * \returns the current value of m_r_0 * pow(lam, 6/5), in m.
     */
   realT r_0( const realT & lam /**< [in] the wavelength, in m, at which to calculate r_0*/ );
   
   ///Set the value of Fried's parameter and the reference wavelength.
   /** If the provided reference wavelength is \<=0, then 0.5 microns is used.
     * 
     */ 
   void r_0( const realT & r0, ///< [in] is the new value of r_0, m  
             const realT & l0 ///< [in] is the new value of lam_0 (in m), if 0 then 0.5e-6 is the default.
           ); 
   
   ///Get the current value of the reference wavelength.
   /** This is the wavelength at which r_0 is specified.
     *
     * \returns the current value of m_lam_0, in m. 
     */
   realT lam_0();
   
   ///Get the value of the outer scale for a single layer.
   /**
     * \returns the current value of m_L_0[n], in m.
     */ 
   realT L_0( const size_t & n );
   
   ///Set the vector of layer outer scales.
   /**
     */ 
   void L_0( const std::vector<realT> & L0 /**< [in] is the new vector of outer scales, in m  */);

   /// Get the vector of outer scales.
   /**
     * \returns a copy of the vector of layer outer scales, in m
     */
   std::vector<realT> L_0();
   
   ///Get the value of the inner scale for a single layer.
   /**
     * \returns the current value of m_l_0[n], in m.
     */ 
   realT l_0( const size_t & n );
   
   ///Set the vector of layer inner scales.
   /**
     */ 
   void l_0( const std::vector<realT> & l0 /**< [in] is the new vector of inner scales, in m*/);

   /// Get the vector of inner scales.
   /**
     * \returns a copy of the vector of layer inner scales, in m
     */
   std::vector<realT> l_0();
   
   ///Get the height of a single layer.
   /**
     * \returns the height of layer n, in m.
     */
   realT layer_z( const size_t n /**< [in] specifies the layer. */);
   
   ///Get the vector layer heights.
   /** 
     * \returns a copy of the vector of layer heights:m_layer_z, in m.
     */ 
   std::vector<realT> layer_z();
   
   ///Set the vector of layer heights.
   /**
     */
   void layer_z(const std::vector<realT> & layz /**< [in] new vector of layer heights, in m. */ );
   
   /// Get the height of the observatory.
   /** 
     * \return the current value of m_h_obs, in m.
     */ 
   realT h_obs();
   
   /// Set the height of the observatory.
   /**
     */
   void h_obs( realT nh /**< [in] the new height of the observatory [m] */);
   
   /// Get the atmospheric scale height.
   /** 
     * \returns the current value of m_H, in m.
     */ 
   realT H();
   
   /// Set the atmospheric scale height.
   /**
     */
   void H( realT nH /**< [in] the new value of m_H [m] */);
   
   ///Get the strength of a single layer.
   /**
     * \returns the value of m_layer_Cn2[n].
     */
   realT layer_Cn2(const int n /**< [in] specifies the layer. */);
   
   ///Get the vector of layer strengths.
   /** 
     * \returns a copy of the vector of layer strengths: m_layer_Cn2.
     */ 
   std::vector<realT> layer_Cn2();
   
   ///Set the vector layer strengths, possibly calculating r_0.
   /**   
     * If a reference wavelength is specified (l0 > 0), then r_0 is set from the layer strengths
     * according to
     * 
     * Regardless of what units the strengths are specified in, they are stored normalized so that
     * \f$ \sum_n C_n^2 = 1 \f$.
     * 
     */  
   void layer_Cn2( const std::vector<realT> & cn2, ///<  [in] is a vector containing the layer strengths 
                   const realT l0 = 0  ///< [in] [optional] if l0 > 0, then r_0 is set from the layer strengths.
                 );
   
   ///Get the wind speed of a single layer.
   /**
     * \returns the value of m_layer_v_wind[n].
     */
   realT layer_v_wind(const int n /**< [in] specifies the layer. */);
   
   ///Get the vector of layer windspeeds
   /** 
     * \returns a copy of the vector of layer windspeeds: m_layer_v_wind.
     */ 
   std::vector<realT>  layer_v_wind();
   
   ///Set the vector of layer windspeeds.
   /**
     * \param 
     */
   void layer_v_wind(const std::vector<realT> & spd /**< [in] the new vector, which is copied to m_layer_v_wind. */);
   
   ///Get the wind direction of a single layer.
   /**
     * \returns the value of m_layer_dir[n], which is the wind direction in that layer in radians.
     */
   realT layer_dir(const int n /**< [in] specifies the layer. */);
   
   ///Get the vector of layer wind directions
   /** 
     * \returns a copy of the vector of layer wind directions: m_layer_dir.
     */
   std::vector<realT> layer_dir();
   
   ///Set the vector of layer wind directions.
   /**
     * \param 
     */
   void layer_dir(const std::vector<realT> & d /**< [in] the new vector of wind directions in radians, which is copied to m_layer_dir. */);
   
   ///Get the number of layers
   /**
     * \returns the size of the m_layer_Cn2 vector.
     */ 
   size_t n_layers();
   
   ///Get the 5/3 moment weighted mean wind speed
   /** Returns the weighted mean wind speed according to the 5/3's turbulence moment.  This is defined as
     * 
     \f[
      \bar{v} = \left[\sum_i C_N^2(z_i) v_i^{5/3} \right]^{3/5}
     \f]
     * See Hardy (1998) Section 3.3.6. \cite hardy_1998.
     * 
     * This is only re-calculated if either m_layer_Cn2 or m_layer_v_wind is changed, otherwise this just
     * returns the value of m_v_wind.
     * 
     * \returns the current value of m_v_wind.
     */ 
   realT v_wind();
   
   /// Get the mean wind speed
   /** Returns the weighted mean wind speed according to the 5/3's turbulence moment.  This is defined as
     * 
     \f[
      \bar{v} = \sum_i C_N^2(z_i) v_i 
     \f]
     * See Hardy (1998) Section 3.3.6. \cite hardy_1998.
     * 
     * This is only re-calculated if either m_layer_Cn2 or m_layer_v_wind is changed, otherwise this just
     * returns the value of m_v_wind.
     * 
     * \returns the current value of m_v_wind.
     */ 
   realT v_wind_mean();
   
   /// Get the mean-squared wind speed
   realT v_wind_mean2();
   
   realT v_max()
   {
      auto max = std::max_element(std::begin(m_layer_v_wind), std::end(m_layer_v_wind));
      
      return *max;
   }
   
   ///Get the weighted mean wind direction
   /** Returns the weighted mean wind direction according to the 5/3's turbulence moment.  This is defined as
     * 
     \f[
      \bar{\theta} = \left[\sum_i C_N^2(z_i) \theta_i^{5/3} \right]^{3/5}
     \f]
     * See Hardy (1998) Section 3.3.6. \cite hardy_1998.
     * 
     * This is only re-calculated if either m_layer_Cn2 or m_layer_v_wind is changed, otherwise this just
     * returns the value of m_dir_wind.
     * 
     * \returns the current value of m_dir_wind.
     */ 
   realT dir_wind();
   
protected:
   ///Recalculate m_v_wind
   /** Called by v_wind() whenever m_v_wind_updated is true.
     *
     * \todo handle dir_wind averaging across 0/2pi  
     */
   void update_v_wind();
  
public:   
   ///Set the weighted mean m_v_wind and renormalize the layer wind speeds
   /** Calling this function changes the values of m_layer_v_wind so that the 
     * Layer averaged 5/3 \f$C_n^2\f$ moment of wind speed is the new value specified by vw.
     * 
     */
   void v_wind(const realT & vw /**< [in] the new value of m_v_wind. */);
   
   
   ///Get the weighted mean layer height
   /** Returns the weighted layer height according to the 5/3's turbulence moment.  This is defined as
     \f[
      \bar{z} = \left[\sum_i C_N^2(z_i) z_i^{5/3} \right]^{3/5}
     \f]
     * See Hardy (1998) Section 3.3.6 and 3.7.2. \cite hardy_1998.
     * 
     * This is only re-calculated if either m_layer_Cn2 orm_layer_z is changed, otherwise this just
     * returns the value of m_z_mean.
     * 
     * \returns the current value of m_z_mean.
     */ 
   realT z_mean();
   
protected:
   ///Recalculate m_z_mean
   /** Called by z_mean() whenever m_z_mean_updated is true.
     */
   void update_z_mean();
  
public:   
   ///Set the weighted mean m_z_mean and renormalize the layer heights
   /** Calling this function changes the values ofm_layer_z so that the 
     * Layer averaged 5/3 \f$C_n^2\f$ moment of the height is the new value specified.
     * 
     */
   void z_mean(const realT & zm /**< [in] is the new value of m_v_wind. */);
   
   /// Set the value of m_nonKolmogorov 
   /** This flag indicates if non-Kolmogorov turbulence is being modeled.
     */
   void nonKolmogorov( const bool & nk /**< [in] the value of m_nonKolmogorov*/);

   /// Return the value of m_nonKolmogorov 
   /** This flag indicates if non-Kolmogorov turbulence is being modeled.
     * 
     * \returns the current value of m_nonKolmogorov 
     */
   bool nonKolmogorov();

   ///Return the PSD index for Kolmogorov turbulence.
   /** Satifies the requirements of psdParamsT.
     *
     * If m_nonKolmogorov is false, this returns
     * \[
       \alpha = \frac{11}{3}
       \]
     * Otherwise it returns the current value of m_alpha.
     *
     * \returns the PSD index.
     */
   realT alpha();
   
   ///Set the PSD index alpha.
   /** Setting alpha with this function also sets m_nonKolmogorov to true.
     */
   void alpha( realT a /**< [in] the new value of alpha */);
   
   ///Return the PSD normalization constant for Kolmogorov turbulence.
   /** Satifies the requirements of psdParamsT.
     *
     * If m_nonKolmogorov is false, this returns
     * \[
       \beta = \frac{0.0218}{r_0^{5/3}}
       \]
     *
     * Otherwise it returns the current value of m_beta.
     *
     * \returns the PSD normalization constant.
     */
   realT beta();
   
   ///Set the PSD normalization alpha.
   /** Setting beta with this function also sets m_nonKolmogorov to true.
     */
   void beta( realT b /**< [in] the new beta */);
   
   ///The fraction of the turbulence PSD in phase after Fresnel propagation.
   /** See Equation (14) of Guyon (2005) \cite guyon_2005. 
     *
     * \returns the value of the X function. 
     */
   realT X( realT k,  ///< [in] the spatial frequency, in inverse meters.
            realT lam_sci, ///< [in] is the science observation wavelength.
            realT secZ ///< [in] is the secant of the zenith distance.
          );
   
   ///The differential fraction of the turbulence PSD in phase after Fresnel propagation.
   /** See Equation (25) of Guyon (2005) \cite guyon_2005. 
     *
     * \returns the value of the dX function. 
     */
   realT dX( realT k,  ///< [in] the spatial frequency, in inverse meters.
             realT lam_sci, ///<  [in] is the science observation wavelength.
             realT lam_wfs ///< [in] is the wavefront sensor wavelength.
           );

   ///The fraction of the turbulence PSD in amplitude after Fresnel propagation.
   /** See Equation (15) of Guyon (2005) \cite guyon_2005. 
     *
     * \returns the value of the Y function. 
     */
   realT Y( realT k, ///< [in] the spatial frequency, in inverse meters.
            realT lam_sci, ///< [in] is the science observation wavelength.
            realT secZ ///< [in] is the secant of the zenith distance.
         );
   
   ///The differential fraction of the turbulence PSD in amplitude after Fresnel propagation.
   /** See Equation (27) of Guyon (2005) \cite guyon_2005. 
     *
     * \returns the value of the dY function. 
     */
   realT dY( realT k,       ///< [in] the spatial frequency, in inverse meters.
             realT lam_sci, ///< [in] is the science observation wavelength.
             realT lam_wfs  ///< [in] is the wavefront sensor wavelength.
           );
    
   
   realT n_air( realT lam /**< [in]The wavelength*/);
   
   realT X_Z( realT k, ///< [in] the spatial frequency, in inverse meters 
              realT lambda_sci, ///< [in] is the science observation wavelength.
              realT lambda_wfs, ///< [in] is the wavefront sensor wavelength.
              realT secZ ///< [in] is the secant of the zenith distance.
            );
   
   ///Calculate the full-width at half-maximum of a seeing limited image for this atmosphere for a small telescope (ignoring L_0)
   /** Calculate the FWHM of a seeing limited image with the current parameters according to Floyd et al. (2010) \cite floyd_2010
     \f[
      \epsilon_0 = 0.98\frac{\lambda_{sci}}{r_0(\lambda_sci)}.
      \f]
     *
     * \returns the value of the FWHM for the current atmosphere parameters.
     */ 
   realT fwhm0(realT lam_sci /**< [in] the wavelength of the science observation. */ );
   
   ///Calculate the full-width at half-maximum of a seeing limited image for this atmosphere for a large telescope (including L_0)
   /** Calculate the FWHM of a seeing limited image with the current parameters according to Floyd et al. (2010) \cite floyd_2010
     \f[
      \epsilon_0 = 0.98\frac{\lambda_{sci}}{r_0(\lambda_sci)}.
      \f]
     * If there is an outer scale (_L_0 > 0), then a correction is applied according to Tokovinin (2002) \cite tokovinin_2002
     \f[
      \left( \frac{\epsilon_{vK}}{\epsilon_0}\right)^2 = 1 - 2.183\left( \frac{r_0(\lambda_{sci}}{L_0}\right)^{0.356}
     \f]
     *
     *
     * \returns the value of the FWHM (\f$ \epsilon_{0/vK} \f$) for the current atmosphere parameters.
     */ 
   realT fwhm(realT lam_sci /**< [in] the wavelength of the science observation. */ );
   
   ///Get the greenwood frequency at the reference wavelength
   /**
     * \todo derive full value of the constant 
     */
   realT f_g();

   ///Get the greenwood frequency at a specified wavelength
   /**
     * 
     */
   realT f_g(realT lam_sci /**< [in] the wavelength of the science observation. */);
   
   ///Get tau_0 at the reference wavelength
   /**
     * \todo derive full value of the constant
     */ 
   realT tau_0();

   ///Get tau_0 at a specified wavelength.
   /**
     * 
     */ 
   realT tau_0(realT lam_sci /**< [in] the wavelength of the science observation. */);
   
   /// Scale v_wind so that tau_0 has the specified value at the specified wavelength.
   /** Does not modify r_0.
     */
   void tau_0( realT tau_0,  ///< [in] the desired tau_0
               realT lam_sci ///< [in] the wavelength of the science observation.
             );

   ///Load the default atmosphere model from Guyon (2005).
   /** Sets the parameters from Table 4 of Guyon (2005) \cite guyon_2005.
     */
   void loadGuyon2005();
   
   ///Load parameters corresponding to the median atmosphere of the GMT site survey at LCO.
   /**
    */
   void loadLCO();
   

   /// Set a single layer model.
   /** Sets all layer vectors to size=1 and populates their fields based on these arguments.
     * 
     */
   void setSingleLayer( realT r0,   ///< [in] is the new value of r_0
                        realT lam0, ///< [in] is the new value of lam_0, if 0 then 0.5 microns is the default.
                        realT L0,   ///< [in] the new outer scale
                        realT l0,   ///< [in] the new inner scale
                        realT lz,   ///< [in] the layer height
                        realT vw,   ///< [in] the layer wind-speed
                        realT dir   ///< [in] the layer wind direction.
                      );
      
   ///Output current parameters to a stream
   /** Prints a formatted list of all current parameters.
     *
     * \tparam iosT is a std::ostream-like type.
     * 
     * \todo update for new vector components (L_0, etc.)
     */ 
   template<typename iosT>
   iosT & dumpAtmosphere( iosT & ios /**< [in] a std::ostream-like stream. */);

   /// \name mx::application support
   /** @{
     */
   
   /// Setup the configurator to configure this class
   /**
     * \test Loading aoAtmosphere config settings \ref tests_ao_analysis_aoAtmosphere_config "[test doc]" 
     */
   void setupConfig( app::appConfigurator & config /**< [in] the app::configurator object*/);

   /// Load the configuration of this class from a configurator
   /**
     * \test Loading aoAtmosphere config settings \ref tests_ao_analysis_aoAtmosphere_config "[test doc]"
     */
   void loadConfig( app::appConfigurator & config /**< [in] the app::configurator object*/);


   /// @}
};

template<typename realT>
aoAtmosphere<realT>::aoAtmosphere()
{
   m_lam_0 = 0.5e-6;
   
   
   m_h_obs = 0;
   m_H = 8000;
   
   m_v_wind_updated = false;
   m_z_mean_updated = false;
}

template<typename realT>
int aoAtmosphere<realT>::checkLayers()
{
   size_t n = m_L_0.size();
   
   if( m_l_0.size() != n)
   {
      mxError("aoAtmosphere", MXE_SIZEERR, "mismatched layer numbers (inner scale vs. outer scale)");
      return -1;
   }
   
   if( m_layer_z.size() != n)
   {
      mxError("aoAtmosphere", MXE_SIZEERR, "mismatched layer numbers (layer_z  vs. outer scale)");
      return -1;
   }
   
   if( m_layer_Cn2.size() != n)
   {
      mxError("aoAtmosphere", MXE_SIZEERR, "mismatched layer numbers (layer_Cn2  vs. outer scale)");
      return -1;
   }
   
   if( m_layer_dir.size() != n)
   {
      mxError("aoAtmosphere", MXE_SIZEERR, "mismatched layer numbers (layer_dir  vs. outer scale)");
      return -1;
   }
   
   if( m_layer_v_wind.size() != n)
   {
      mxError("aoAtmosphere", MXE_SIZEERR, "mismatched layer numbers (layer_v_wind vs. outer scale)");
      return -1;
   }
   

   return 0;
}

template<typename realT>
realT aoAtmosphere<realT>::r_0()
{
   return m_r_0;
}

template<typename realT>
realT aoAtmosphere<realT>::r_0(const realT & lam)
{
   return m_r_0*pow(lam/m_lam_0, math::six_fifths<realT>());
}

template<typename realT>
void aoAtmosphere<realT>::r_0(const realT & r0, const realT  & l0)
{
   m_r_0 = r0;
   
   if( l0 > 0)
   {
      m_lam_0 = l0;
   }
   else
   {
      m_lam_0 = 0.5e-6;
   }
}

template<typename realT>
realT aoAtmosphere<realT>::lam_0()
{
   return m_lam_0;
}

template<typename realT>
realT aoAtmosphere<realT>::L_0(const size_t & n)
{
   return m_L_0[n];
}

template<typename realT>
void aoAtmosphere<realT>::L_0(const std::vector<realT> & L_0)
{
   m_L_0 = L_0;
}

template<typename realT>
std::vector<realT> aoAtmosphere<realT>::L_0()
{
   return m_L_0;
}

template<typename realT>
realT aoAtmosphere<realT>::l_0(const size_t & n)
{
   return m_l_0[n];
}

template<typename realT>
void aoAtmosphere<realT>::l_0(const std::vector<realT> & l_0)
{
   m_l_0 = l_0;
}

template<typename realT>
std::vector<realT> aoAtmosphere<realT>::l_0()
{
   return m_l_0;
}

template<typename realT>
realT aoAtmosphere<realT>::layer_z(const size_t n)
{
   return m_layer_z[n];
}

template<typename realT>
void aoAtmosphere<realT>::layer_z(const std::vector<realT> & z)
{
  m_layer_z = z;
  m_z_mean_updated = false;
}

template<typename realT>
std::vector<realT> aoAtmosphere<realT>::layer_z()
{
   return m_layer_z;
}

template<typename realT>
realT aoAtmosphere<realT>::h_obs()
{
   return m_h_obs;
}
   
template<typename realT>
void aoAtmosphere<realT>::h_obs( realT nh )
{
   m_h_obs = nh;
}

template<typename realT>
realT aoAtmosphere<realT>::H()
{
   return m_H;
}
   
template<typename realT>
void aoAtmosphere<realT>::H( realT nH )
{
   m_H = nH;
}
   

template<typename realT>
realT aoAtmosphere<realT>::layer_Cn2(const int n)
{
   return m_layer_Cn2[n];
}

template<typename realT>
std::vector<realT>  aoAtmosphere<realT>::layer_Cn2()
{
   return m_layer_Cn2;
}
   
template<typename realT>
void aoAtmosphere<realT>::layer_Cn2(const std::vector<realT> & cn2, const realT l0)
{
   m_layer_Cn2 = cn2; 
   
   realT layer_norm = 0;
   
   for(size_t i=0;i < m_layer_Cn2.size();++i) 
   {
      layer_norm += cn2[i];
   }
      
   for(size_t i=0;i< m_layer_Cn2.size();++i) m_layer_Cn2[i] = m_layer_Cn2[i]/layer_norm;
   
   if(l0 > 0)
   {   
      m_r_0 = 1.0 / pow(layer_norm * 5.520e13, math::three_fifths<realT>() );  
      m_lam_0 = l0;
   }

   m_v_wind_updated = false;
   m_z_mean_updated = false;
}

template<typename realT>
realT aoAtmosphere<realT>::layer_v_wind(const int n)
{
   return m_layer_v_wind[n];
}


template<typename realT>
std::vector<realT> aoAtmosphere<realT>::layer_v_wind()
{
   return m_layer_v_wind;
}

template<typename realT>
void aoAtmosphere<realT>::layer_v_wind(const std::vector<realT> & v)
{
   m_layer_v_wind = v;
   m_v_wind_updated = false;
}

template<typename realT>
realT  aoAtmosphere<realT>::layer_dir(const int n)
{
   return m_layer_dir[n];
}

template<typename realT>
std::vector<realT>  aoAtmosphere<realT>::layer_dir()
{
   return m_layer_dir;
}

template<typename realT>
void aoAtmosphere<realT>::layer_dir(const std::vector<realT> & dir)
{
   m_layer_dir = dir;
   m_v_wind_updated = false;
}

template<typename realT>
size_t aoAtmosphere<realT>::n_layers()
{
   return m_layer_Cn2.size();
}

template<typename realT>
realT aoAtmosphere<realT>::v_wind()
{
   if(m_v_wind_updated == false) update_v_wind();
   return m_v_wind;
}

template<typename realT>
realT aoAtmosphere<realT>::dir_wind()
{
   if(m_v_wind_updated == false) update_v_wind();
   return m_dir_wind;
}

template<typename realT>
void aoAtmosphere<realT>::update_v_wind()
{

   if( m_layer_v_wind.size() == 0 || m_layer_Cn2.size() == 0)
   {
      m_v_wind = 0;
      m_dir_wind = 0;
      m_v_wind_updated = true;
      return;
   }
   
   m_v_wind = 0;
   
   realT s = 0;
   realT c = 0;
   
   for(size_t i=0;i<m_layer_Cn2.size(); ++i)
   {
      m_v_wind += m_layer_Cn2[i] * pow(m_layer_v_wind[i], math::five_thirds<realT>() );
      s += pow(m_layer_v_wind[i], math::five_thirds<realT>())*sin(m_layer_dir[i]) ;//pow( sin(m_layer_dir[i]), math::five_thirds<realT>() );
      c += pow(m_layer_v_wind[i], math::five_thirds<realT>())*cos(m_layer_dir[i]) ;//pow( cos(m_layer_dir[i]), math::five_thirds<realT>() );
   }
      
   m_v_wind = pow(m_v_wind, math::three_fifths<realT>());
   
   //m_dir_wind = atan(pow(s, math::three_fifths<realT>()) / pow(c, math::three_fifths<realT>()));
   m_dir_wind = atan( s / c);
   if(m_dir_wind < 0) m_dir_wind += math::pi<realT>();
   
   m_v_wind_updated = true;
   
}

template<typename realT>
void aoAtmosphere<realT>::v_wind(const realT & vw)
{
   if(m_v_wind_updated == false) update_v_wind();
   
   realT vw_old = m_v_wind;
   
   m_v_wind = vw;
   
   //Now update the layers if needed
   if( m_layer_v_wind.size() > 0)
   {
      for(size_t i=0; i< m_layer_v_wind.size(); ++i)
      {
         m_layer_v_wind[i] = m_layer_v_wind[i]*(vw/vw_old);
      }
   }
}

template<typename realT>
realT aoAtmosphere<realT>::z_mean()
{
   if(m_z_mean_updated == false) update_z_mean();
   return m_z_mean;
}

template<typename realT>
void aoAtmosphere<realT>::update_z_mean()
{
   if(m_layer_z.size() == 0 || m_layer_Cn2.size() == 0)
   {
      m_z_mean = 0;
      m_z_mean_updated = true;
      return;
   }
   
   m_z_mean = 0;
   
   for(size_t i=0;i<m_layer_Cn2.size(); ++i)
   {
      m_z_mean += m_layer_Cn2[i] * pow(m_layer_z[i], math::five_thirds<realT>() );
   }
      
   m_z_mean = pow(m_z_mean, math::three_fifths<realT>());
   
   m_z_mean_updated = true;
   
}

template<typename realT>
void aoAtmosphere<realT>::z_mean(const realT & zm)
{
   if(m_z_mean_updated == false) update_z_mean();
   
   realT zh_old = m_z_mean;
   
   m_z_mean = zm;
   
   //Now update the layers if needed
   if(m_layer_z.size() > 0)
   {
      for(size_t i=0; i<m_layer_z.size(); ++i)
      {
        m_layer_z[i] =m_layer_z[i]*(zm/zh_old);
      }
   }
}

template<typename realT>
void aoAtmosphere<realT>::nonKolmogorov( const bool & nk )
{
   m_nonKolmogorov = nk;
}

template<typename realT>
bool aoAtmosphere<realT>::nonKolmogorov()
{
   return m_nonKolmogorov;
}

template<typename realT>
realT aoAtmosphere<realT>::alpha()
{
   if(!m_nonKolmogorov)
   {
      return math::eleven_thirds<realT>();
   }
   else
   {
      return m_alpha;
   }
}

template<typename realT>
void aoAtmosphere<realT>::alpha(realT a)
{
   m_nonKolmogorov = true;
   
   m_alpha = a;
}

template<typename realT>
realT aoAtmosphere<realT>::beta()
{
   if(!m_nonKolmogorov)
   {
      return constants::a_PSD<realT>()*pow(m_r_0, -math::five_thirds<realT>());
   }
   else
   {
      return m_beta;
   }
}
   
template<typename realT>
void aoAtmosphere<realT>::beta(realT b)
{
   m_nonKolmogorov = true;
   
   m_beta = b;
}

template<typename realT>
realT aoAtmosphere<realT>::X( realT k, 
                              realT lam_sci, 
                              realT secZ)
{
   realT c = 0;
    
   for(size_t i=0;i<m_layer_Cn2.size(); ++i)
   {
      c += m_layer_Cn2[i] * pow( cos(math::pi<realT>()*k*k*lam_sci *m_layer_z[i] * secZ), 2);
   }
   
   return c;
}
 
template<typename realT>
realT aoAtmosphere<realT>::dX(realT f, realT lam_sci, realT lam_wfs)
{
   realT c = 0;
 
   for(size_t i=0;i<m_layer_Cn2.size(); ++i)
   {   
      c += m_layer_Cn2[i]* pow((cos( math::pi<realT>()*f*f*lam_sci *m_layer_z[i]) - cos( math::pi<realT>()*f*f*lam_wfs *m_layer_z[i])), 2);
   }
   
   return c;
}

template<typename realT>
realT aoAtmosphere<realT>::Y( realT k, 
                              realT lam_sci,
                              realT secZ
                            )
{
   realT c = 0;
 
   for(size_t i=0;i<m_layer_Cn2.size(); ++i)
   {   
      c += m_layer_Cn2[i]*pow(sin( math::pi<realT>()*k*k*lam_sci *m_layer_z[i] * secZ), 2);
   }
   return c;
}

template<typename realT>
realT aoAtmosphere<realT>::dY(realT f, realT lam_sci, realT lam_wfs)
{
   realT c = 0;

   for(size_t i=0;i<m_layer_Cn2.size(); ++i)
   {   
      c += m_layer_Cn2[i]*pow( (sin(math::pi<realT>()*f*f*lam_sci *m_layer_z[i]) - sin( math::pi<realT>()*f*f*lam_wfs *m_layer_z[i])), 2);
   }
   
   return c;
}

template<typename realT>
realT aoAtmosphere<realT>::n_air( realT lambda)
{
   realT ll2 = static_cast<realT>(1)/pow(lambda/1e-6, 2);
   
   return 1.0 + 8.34213e-5 + 0.0240603/(130.0 - ll2) + 0.00015997/(38.9 - ll2);
}

template<typename realT>
realT aoAtmosphere<realT>::X_Z(realT k, realT lambda_i, realT lambda_wfs, realT secZ)
{
   realT c = 0;
   realT sinZ = sqrt(1.0 - pow(1.0/secZ,2));
   realT tanZ = sinZ*secZ;
   realT x0 = (n_air(lambda_wfs) - n_air(lambda_i)) * m_H*tanZ*secZ; 
   realT x;
   
   for(size_t i = 0; i < m_layer_Cn2.size(); ++i)
   {
      x = x0*(1-exp((m_layer_z[i]+m_h_obs)/m_H));
      c += m_layer_Cn2[i] * pow( cos(math::pi<realT>()*k*k*lambda_i *m_layer_z[i]*secZ), 2) * pow( sin(math::pi<realT>()*x*k*cos(0.*3.14/180.)), 2);
   }
   
   return 4*c;
}

template<typename realT>
realT aoAtmosphere<realT>::fwhm0(realT lam_sci)
{
   realT r0lam = r_0(lam_sci);
   
   realT fwhm = 0.98*(lam_sci/r0lam);

   return fwhm;
}


template<typename realT>
realT aoAtmosphere<realT>::fwhm(realT lam_sci)
{
   realT r0lam = r_0(lam_sci);
   
   realT fwhm = 0.98*(lam_sci/r0lam);

   ///\todo this needs to handle layers with different L_0
   if( L_0(0) > 0) fwhm *= sqrt( 1 - 2.183*pow(r0lam/L_0(0), 0.356));
   
   return fwhm;
}

template<typename realT>
realT aoAtmosphere<realT>::f_g()
{
   return 0.428*v_wind()/m_r_0;
}

template<typename realT>
realT aoAtmosphere<realT>::f_g(realT lam_sci)
{
   return 0.428*pow(m_lam_0/lam_sci, math::six_fifths<realT>())*v_wind()/m_r_0;
}

template<typename realT>
realT aoAtmosphere<realT>::tau_0()
{
   return 0.134/f_g();
}

template<typename realT>
realT aoAtmosphere<realT>::tau_0(realT lam_sci)
{
   return 0.134/f_g(lam_sci);
}

template<typename realT>
void aoAtmosphere<realT>::tau_0( realT tau_0,  
                                 realT lam_sci 
                               )
{
   realT vw = (0.134/tau_0) / 0.428 * m_r_0* pow(lam_sci/m_lam_0, math::six_fifths<realT>());
   v_wind(vw);
}

template<typename realT>
void aoAtmosphere<realT>::loadGuyon2005()
{
   layer_Cn2({0.2283, 0.0883, 0.0666, 0.1458, 0.3350, 0.1350});
   layer_z({500, 1000, 2000, 4000, 8000, 16000});
   layer_v_wind({10., 10.,  10.,  10.,  10.,  10.});
   layer_dir({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
   
   r_0(0.2, 0.5e-6);
   
   h_obs(4200);
}

template<typename realT>
void aoAtmosphere<realT>::loadLCO()
{
   layer_Cn2(      {0.42,      0.029,      0.062,    0.16,     0.11,     0.10,     0.12});
   layer_z(        {250.,       500.,       1000.,   2000.,    4000.,    8000.,    16000. });
   layer_v_wind(   {10.0,      10.0,         20.0,    20.0,    25.0,     30.0,      25.0});
   layer_dir(      {1.05,      1.05,        1.31,     1.31,    1.75,     1.92,     1.75});
   
   r_0(0.17, 0.5e-6);
   
   L_0({25.0,25.0,25.0,25.0, 25.0, 25.0, 25.0});
   l_0({0.0,0,0,0,0,0,0});
   
   h_obs(2400);
}

template<typename realT>
void aoAtmosphere<realT>::setSingleLayer( realT r0,
                                          realT lam0,
                                          realT L0,
                                          realT l0,
                                          realT lz,
                                          realT vw,
                                          realT dir
                                        )
{      
   r_0(r0, lam0);
   L_0(std::vector<realT>({L0}));
   l_0(std::vector<realT>({l0}));
   layer_Cn2(std::vector<realT>({1}));
   layer_z(std::vector<realT>({lz}));
   layer_v_wind(std::vector<realT>({vw}));
   layer_dir(std::vector<realT>({dir}));
}
   
template<typename realT>
template<typename iosT>
iosT & aoAtmosphere<realT>::dumpAtmosphere( iosT & ios)
{
   ios << "# Atmosphere Parameters:\n";
   ios << "#    r_0 = " << r_0() << '\n';
   ios << "#    lam_0 = " << lam_0() << '\n';
   ios << "#    tau_0 = " << tau_0(lam_0()) << '\n';
   ios << "#    L_0 = ";
   for(size_t i=0;i < n_layers()-1;++i) ios << L_0()[i] << ", ";
   ios <<  L_0()[n_layers()-1] << '\n';
   ios << "#    FWHM = " << fwhm(lam_0()) << '\n';
   ios << "#    n_layers = " << n_layers() << '\n';
   ios << "#    layer_z = ";
   for(size_t i=0;i < n_layers()-1;++i) ios << layer_z()[i] << ", ";
   ios <<  layer_z()[ n_layers()-1] << '\n';
   ios << "#    h_obs = " << h_obs() << '\n';
   ios << "#    H = " << H() << '\n';
   ios << "#    layer_Cn2 = ";
   for(size_t i=0;i< n_layers()-1;++i) ios << layer_Cn2()[i] << ", ";
   ios << layer_Cn2()[n_layers()-1] << '\n';
   ios << "#    layer_v_wind = ";
   for(size_t i=0;i< n_layers()-1;++i) ios << layer_v_wind()[i] << ", ";
   ios << layer_v_wind()[ n_layers()-1] << '\n';
   ios << "#    layer_dir = ";
   for(size_t i=0;i< n_layers()-1;++i) ios << layer_dir()[i] << ", ";
   ios << layer_dir()[ n_layers()-1] << '\n';
   ios << "#    mean v_wind = " << v_wind() << '\n';
   ios << "#    mean dir_wind = " << dir_wind() << '\n';
   ios << "#    mean z = " << z_mean() << '\n';
   
   
   return ios;
}
 
template<typename realT>
void aoAtmosphere<realT>::setupConfig( app::appConfigurator & config )
{
   using namespace mx::app;

   config.add("atm.r_0"          ,"", "atm.r_0"          , argType::Required, "atm", "r_0",           false, "real"        , "Fried's parameter [m]");
   config.add("atm.lam_0"        ,"", "atm.lam_0"        , argType::Required, "atm", "lam_0",         false, "real"        , "The reference wavlength for r_0 [m]");
   config.add("atm.L_0"          ,"", "atm.L_0"          , argType::Required, "atm", "L_0",           false, "vector<real>", "Layer outer scales [m]");
   config.add("atm.l_0"          ,"", "atm.l_0"          , argType::Required, "atm", "l_0",           false, "vector<real>", "Layer inner scales [m]");
   config.add("atm.layer_z"      ,"", "atm.layer_z"      , argType::Required, "atm", "layer_z",       false, "vector<real>", "layer heights [m]");
   config.add("atm.h_obs"        ,"", "atm.h_obs"        , argType::Required, "atm", "h_obs",         false, "real"        , "height of observatory [m]");
   config.add("atm.H"            ,"", "atm.H"            , argType::Required, "atm", "H",             false, "real"        , "atmospheric scale heights [m]");
   config.add("atm.layer_Cn2"    ,"", "atm.layer_Cn2"    , argType::Required, "atm", "layer_Cn2",     false, "vector<real>", "Layer Cn^2.  Note that this does not set r_0.");  
   config.add("atm.layer_v_wind" ,"", "atm.layer_v_wind" , argType::Required, "atm", "layer_v_wind",  false, "vector<real>", "Layer wind speeds [m/s]");
   config.add("atm.layer_dir"    ,"", "atm.layer_dir"    , argType::Required, "atm", "layer_dir",     false, "vector<real>", "Layer wind directions [rad]");
   config.add("atm.v_wind"       ,"", "atm.v_wind"       , argType::Required, "atm", "v_wind",        false, "real"        , "Mean windspeed (5/3 momement), rescales layers [m/s]");
   config.add("atm.z_mean"       ,"", "atm.z_mean"       , argType::Required, "atm", "z_mean",        false, "real"        , "Mean layer height (5/3 momemnt), rescales layers [m/s]");   
   config.add("atm.nonKolmogorov","", "atm.nonKolmogorov", argType::Required, "atm", "nonKolmogorov", false, "bool"        , "Set to use a non-Kolmogorov PSD. See alpha and beta.");   
   config.add("atm.alpha"        ,"", "atm.alpha"        , argType::Required, "atm", "alpha"        , false, "real"        , "Non-kolmogorov PSD exponent.");   
   config.add("atm.beta"         ,"", "atm.beta"         , argType::Required, "atm", "beta"         , false, "real"        , "Non-kolmogorov PSD normalization constant.");   
}

template<typename realT>
void aoAtmosphere<realT>::loadConfig( app::appConfigurator & config )
{
   //Here "has side effecs" means that the set function does more than simply copy the value.
   
   //The order of lam_0, Cn2, and r_0 is so that r_0 overrides the value set with Cn2 if lam_0 != 0.
   //lam_0 comes first because it calibrates r0 and Cn2

   //lam_0
   config(m_lam_0, "atm.lam_0");

   //layer_Cn2
   std::vector<realT> lcn2 = m_layer_Cn2;
   config(lcn2, "atm.layer_Cn2"); 
   if(config.isSet("atm.layer_Cn2")) layer_Cn2(lcn2); 

   realT r0 = r_0();
   config(r0, "atm.r_0");
   if(config.isSet("atm.r_0")) r_0(r0, m_lam_0);

   config(m_L_0, "atm.L_0");
   
   config(m_l_0, "atm.l_0");
   
   //Has side effects:
   std::vector<realT> layz = m_layer_z;
   config(layz, "atm.layer_z"); //Do this no matter what to record source
   if(config.isSet("atm.layer_z")) layer_z(layz); //but only call this if changed

   config(m_h_obs, "atm.h_obs");
   config(m_H, "atm.H");

   //Has side effects:
   std::vector<realT> lvw = m_layer_v_wind;
   config(lvw, "atm.layer_v_wind"); //Do this no matter what to record source
   if(config.isSet("atm.layer_v_wind")) layer_v_wind(lvw); //but only call this if changed

   //Has side effects:
   std::vector<realT> ld = m_layer_dir;
   config(ld, "atm.layer_dir"); //Do this no matter what to record source
   if(config.isSet("atm.layer_dir")) layer_dir(ld); //but only call this if changed

   realT vw = m_v_wind;
   config(vw,"atm.v_wind"); //Do this no matter what to record source
   if(config.isSet("atm.v_wind")) v_wind(vw); //but only call this if changed

   realT zm = m_z_mean;
   config(zm,"atm.z_mean"); //Do this no matter what to record source
   if(config.isSet("atm.z_mean")) z_mean(zm); //but only call this if changed
   
   config(m_nonKolmogorov, "atm.nonKolmogorov");

   realT a = m_alpha;
   config(a, "atm.alpha");
   if(config.isSet("atm.alpha")) alpha(a); //this sets m_nonKolmogorov

   realT b = m_beta;
   config(b, "atm.beta");
   if(config.isSet("atm.beta")) alpha(b); //this sets m_nonKolmogorov

}


extern template
class aoAtmosphere<float>;

extern template
class aoAtmosphere<double>;

extern template
class aoAtmosphere<long double>;

#ifdef HASQUAD
extern template
class aoAtmosphere<__float128>;
#endif

}//namespace analysis 
}//namespace AO
}//namespace mx

#endif //aoAtmosphere_hpp
