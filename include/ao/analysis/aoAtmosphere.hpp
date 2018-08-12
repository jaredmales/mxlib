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

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <mx/mxlib.h>

#include "aoConstants.hpp"
using namespace mx::AO::constants;

namespace mx
{
namespace AO
{
namespace analysis 
{
   
///A class to specify atmosphere parameters and perform related calculations.
/** 
  * \todo Add layer outer scales.
  * \todo Add calculation of f_G and tau_0.
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
   realT _r_0; ///< Fried's parameter, in m
   
   realT _lam_0 ;///<Wavelength of Fried's paraemter, in m

   realT _L_0; ///<The outer scale, in m

   std::vector<realT> _layer_z; ///< Vector of layer heights, in m, above the observatory.
   
   realT _h_obs; ///< Height of the observatory above sea level, in m.  
   
   realT _H; ///< The atmospheric scale height, in m.
   
   std::vector<realT> _layer_Cn2; ///< Vector of layer strengths. 
   
   std::vector<realT> _layer_v_wind; ///< Vector of layer wind speeds, in m/s.
   
   std::vector<realT> _layer_dir; ///<Vector of layer wind directions, in radians.

   bool _v_wind_updated; ///< whether or not _v_wind has been updated after changes
   
   realT _v_wind; ///< \f$ C_n^2 \f$ averaged windspeed
   
   realT _dir_wind; ///< \f$ C_n^2 \f$ averaged direction
   
   bool _z_mean_updated; ///< whether or not _z_mean has been update after changes
   
   realT _z_mean; ///< \f$ C_n^2 \f$ averaged layer height
   
    
   bool m_nonKolmogorov {false};
   
   realT m_alpha{2};
   
   realT m_beta{1};
   
public:
   
   ///Get the value of Fried's parameter r_0 at the reference wavelength lam_0.
   /**
     * \returns the curret value of _r_0.
     */ 
   realT r_0();
   
   ///Get the value of Fried's parameter r_0 at the specified wavelength.
   /**
     * \note This does not set the value of r_0.
     * 
     * \returns the current value of _r_0 * pow(lam, 6/5).
     */
   realT r_0( const realT & lam /**< [in] the wavelength at which to calculate r_0*/ );
   
   ///Set the value of Fried's parameter and the reference wavelength.
   /** If the provided reference wavelength is \<=0, then 0.5 microns is used.
     * 
     */ 
   void r_0( const realT & r0, ///< [in] is the new value of _r_0  
             const realT & l0 ///< [in] is the new value of _lam_0, if 0 then 0.5 microns is the default.
           ); 
   
   ///Get the current value of the reference wavelength.
   /** This is the wavelength at which r_0 is specified.
     *
     * \returns the current value of _lam_0. 
     */
   realT lam_0();
   
   ///Get the value of the outer scale.
   /**
     * \returns the current value of _L_0.
     */ 
   realT L_0();
   
   ///Set the value of the outer scale.
   /**
     */ 
   void L_0( const realT & L0 /**< [in] is the new value of _L_0 */);

   ///Get the height of a single layer.
   /**
     * \returns the value of _layer_z[n].
     */
   realT layer_z( const int n /**< [in] specifies the layer. */);
   
   ///Get the vector layer heights.
   /** 
     * \returns a copy of the vector of layer heights: _layer_z.
     */ 
   std::vector<realT> layer_z();
   
   ///Set the vector of layer heights.
   /**
     */
   void layer_z(const std::vector<realT> & layz /**< [in] new vector of layer heights, which is copied to _layer_z. */ );
   
   /// Get the height of the observatory.
   /** 
     * \return the current value of _h_obs.
     */ 
   realT h_obs();
   
   /// Set the height of the observatory.
   /**
     */
   void h_obs( realT nh /**< [in] the new height of the observatory [m] */);
   
   /// Get the atmospheric scale height.
   /** 
     * \returns the current value of _H.
     */ 
   realT H();
   
   /// Set the atmospheric scale height.
   /**
     */
   void H( realT nH /**< [in] the new value of _H [m] */);
   
   ///Get the strength of a single layer.
   /**
     * \returns the value of _layer_Cn2[n].
     */
   realT layer_Cn2(const int n /**< [in] specifies the layer. */);
   
   ///Get the vector of layer strengths.
   /** 
     * \returns a copy of the vector of layer strengths: _layer_Cn2.
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
     * \returns the value of _layer_v_wind[n].
     */
   realT layer_v_wind(const int n /**< [in] specifies the layer. */);
   
   ///Get the vector of layer windspeeds
   /** 
     * \returns a copy of the vector of layer windspeeds: _layer_v_wind.
     */ 
   std::vector<realT>  layer_v_wind();
   
   ///Set the vector of layer windspeeds.
   /**
     * \param 
     */
   void layer_v_wind(const std::vector<realT> & spd /**< [in] the new vector, which is copied to _layer_v_wind. */);
   
   ///Get the wind direction of a single layer.
   /**
     * \returns the value of _layer_dir[n].
     */
   realT layer_dir(const int n /**< [in] specifies the layer. */);
   
   ///Get the vector of layer wind directions
   /** 
     * \returns a copy of the vector of layer wind directions: _layer_dir.
     */
   std::vector<realT> layer_dir();
   
   ///Set the vector of layer wind directions.
   /**
     * \param 
     */
   void layer_dir(const std::vector<realT> & d /**< [in] the new vector, which is copied to _layer_dir. */);
   
   ///Get the number of layers
   /**
     * \returns the size of the _layer_Cn2 vector.
     */ 
   size_t n_layers();
   
   ///Get the weighted mean wind speed
   /** Returns the weighted mean wind speed according to the 5/3's turbulence moment.  This is defined as
     * 
     \f[
      \bar{v} = \left[\sum_i C_N^2(z_i) v_i^{5/3} \right]^{3/5}
     \f]
     * See Hardy (1998) Section 3.3.6. \cite hardy_1998.
     * 
     * This is only re-calculated if either _layer_Cn2 or _layer_v_wind is changed, otherwise this just
     * returns the value of _v_wind.
     * 
     * \returns the current value of _v_wind.
     */ 
   realT v_wind();
   
   ///Get the weighted mean wind direction
   /** Returns the weighted mean wind direction according to the 5/3's turbulence moment.  This is defined as
     * 
     \f[
      \bar{\theta} = \left[\sum_i C_N^2(z_i) \theta_i^{5/3} \right]^{3/5}
     \f]
     * See Hardy (1998) Section 3.3.6. \cite hardy_1998.
     * 
     * This is only re-calculated if either _layer_Cn2 or _layer_v_wind is changed, otherwise this just
     * returns the value of _dir_wind.
     * 
     * \returns the current value of _dir_wind.
     */ 
   realT dir_wind();
   
protected:
   ///Recalculate _v_wind
   /** Called by v_wind() whenever _v_wind_updated is true.
     *
     * \todo handle dir_wind averaging across 0/2pi  
     */
   void update_v_wind();
  
public:   
   ///Set the weighted mean _v_wind and renormalize the layer wind speeds
   /** Calling this function changes the values of _layer_v_wind so that the 
     * Layer averaged 5/3 \f$C_n^2\f$ moment of wind speed is the new value specified by vw.
     * 
     */
   void v_wind(const realT & vw /**< [in] the new value of _v_wind. */);
   
   ///Get the weighted mean layer height
   /** Returns the weighted layer height according to the 5/3's turbulence moment.  This is defined as
     \f[
      \bar{z} = \left[\sum_i C_N^2(z_i) z_i^{5/3} \right]^{3/5}
     \f]
     * See Hardy (1998) Section 3.3.6 and 3.7.2. \cite hardy_1998.
     * 
     * This is only re-calculated if either _layer_Cn2 or _layer_z is changed, otherwise this just
     * returns the value of _z_mean.
     * 
     * \returns the current value of _z_mean.
     */ 
   realT z_mean();
   
protected:
   ///Recalculate _z_mean
   /** Called by z_mean() whenever _z_mean_updated is true.
     */
   void update_z_mean();
  
public:   
   ///Set the weighted mean _z_mean and renormalize the layer heights
   /** Calling this function changes the values of _layer_z so that the 
     * Layer averaged 5/3 \f$C_n^2\f$ moment of the height is the new value specified.
     * 
     */
   void z_mean(const realT & zm /**< [in] is the new value of _v_wind. */);
   
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
     * I m_nonKolmogorov is false, this returns
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
   
   realT X_Z( realT k, ///< [in] the spatial frequency, in inverse metes 
              realT lambda_sci, ///< [in] is the science observation wavelength.
              realT lambda_wfs, ///< [in] is the wavefront sensor wavelength.
              realT secZ ///< [in] is the secant of the zenith distance.
            );
   
   ///Calculate the full-width at half-maximum of a seeing limited image for this atmosphere.
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
   
   ///Load the default atmosphere model from Guyon (2005).
   /** Sets the parameters from Table 4 of Guyon (2005) \cite guyon_2005.
     */
   void loadGuyon2005();
   
   ///Load parameters corresponding to the median atmosphere of the GMT site survey at LCO.
   /**
    */
   void loadLCO();
   

   ///Set a single layer model.
   /** 
     * 
     */
   void setSingleLayer( realT lz,   ///< [in] the layer height
                        realT vw,   ///< [in] the layer wind-seed
                        realT dir   ///< [in] the layer wind direction.
                      )
   {      
      layer_Cn2(std::vector<realT>({1}));
      layer_z(std::vector<realT>({lz}));
      layer_v_wind(std::vector<realT>({vw}));
      layer_dir(std::vector<realT>({dir}));
   }
      
   ///Output current parameters to a stream
   /** Prints a formatted list of all current parameters.
     *
     * \tparam iosT is a std::ostream-like type.
     */ 
   template<typename iosT>
   iosT & dumpAtmosphere( iosT & ios /**< [in] a std::ostream-like stream. */);
};

template<typename realT>
aoAtmosphere<realT>::aoAtmosphere()
{
   _lam_0 = 0.5e-6;
   
   _L_0 = 0;
   
   _h_obs = 0;
   _H = 8000;
   
   _v_wind_updated = false;
   _z_mean_updated = false;
}

template<typename realT>
realT aoAtmosphere<realT>::r_0()
{
   return _r_0;
}

template<typename realT>
realT aoAtmosphere<realT>::r_0(const realT & lam)
{
   return _r_0*pow(lam/_lam_0, six_fifths<realT>());
}

template<typename realT>
void aoAtmosphere<realT>::r_0(const realT & r0, const realT  & l0)
{
   _r_0 = r0;
   
   if( l0 > 0)
   {
      _lam_0 = l0;
   }
   else
   {
      _lam_0 = 0.5e-6;
   }
}

template<typename realT>
realT aoAtmosphere<realT>::lam_0()
{
   return _lam_0;
}

template<typename realT>
realT aoAtmosphere<realT>::L_0()
{
   return _L_0;
}

template<typename realT>
void aoAtmosphere<realT>::L_0(const realT & l0)
{
   _L_0 = l0;
}

template<typename realT>
realT aoAtmosphere<realT>::layer_z(const int n)
{
   return _layer_z[n];
}

template<typename realT>
std::vector<realT> aoAtmosphere<realT>::layer_z()
{
   return _layer_z;
}

template<typename realT>
realT aoAtmosphere<realT>::h_obs()
{
   return _h_obs;
}
   
template<typename realT>
void aoAtmosphere<realT>::h_obs( realT nh )
{
   _h_obs = nh;
}

template<typename realT>
realT aoAtmosphere<realT>::H()
{
   return _H;
}
   
template<typename realT>
void aoAtmosphere<realT>::H( realT nH )
{
   _H = nH;
}
   
template<typename realT>
void aoAtmosphere<realT>::layer_z(const std::vector<realT> & z)
{
   _layer_z = z;
   _z_mean_updated = false;
}

template<typename realT>
realT aoAtmosphere<realT>::layer_Cn2(const int n)
{
   return _layer_Cn2[n];
}

template<typename realT>
std::vector<realT>  aoAtmosphere<realT>::layer_Cn2()
{
   return _layer_Cn2;
}
   
template<typename realT>
void aoAtmosphere<realT>::layer_Cn2(const std::vector<realT> & cn2, const realT l0)
{
   _layer_Cn2 = cn2; 
   
   
   realT layer_norm = 0;
   
   for(size_t i=0;i < _layer_Cn2.size();++i) 
   {
      layer_norm += cn2[i];
   }
      
   for(size_t i=0;i< _layer_Cn2.size();++i) _layer_Cn2[i] = _layer_Cn2[i]/layer_norm;
   
   if(l0 > 0)
   {   
      _r_0 = 1.0 / pow(layer_norm * 5.520e13, three_fifths<realT>() );  
      _lam_0 = l0;
   }

   _v_wind_updated = false;
   _z_mean_updated = false;
}

template<typename realT>
realT aoAtmosphere<realT>::layer_v_wind(const int n)
{
   return _layer_v_wind[n];
}


template<typename realT>
std::vector<realT> aoAtmosphere<realT>::layer_v_wind()
{
   return _layer_v_wind;
}

template<typename realT>
void aoAtmosphere<realT>::layer_v_wind(const std::vector<realT> & v)
{
   _layer_v_wind = v;
   _v_wind_updated = false;
}

template<typename realT>
realT  aoAtmosphere<realT>::layer_dir(const int n)
{
   return _layer_dir[n];
}

template<typename realT>
std::vector<realT>  aoAtmosphere<realT>::layer_dir()
{
   return _layer_dir;
}

template<typename realT>
void aoAtmosphere<realT>::layer_dir(const std::vector<realT> & dir)
{
   _layer_dir = dir;
   _v_wind_updated = false;
}

template<typename realT>
size_t aoAtmosphere<realT>::n_layers()
{
   return _layer_Cn2.size();
}

template<typename realT>
realT aoAtmosphere<realT>::v_wind()
{
   if(_v_wind_updated == false) update_v_wind();
   return _v_wind;
}

template<typename realT>
realT aoAtmosphere<realT>::dir_wind()
{
   if(_v_wind_updated == false) update_v_wind();
   return _dir_wind;
}

template<typename realT>
void aoAtmosphere<realT>::update_v_wind()
{

   if( _layer_v_wind.size() == 0 || _layer_Cn2.size() == 0)
   {
      _v_wind = 0;
      _dir_wind = 0;
      _v_wind_updated = true;
      return;
   }
   
   _v_wind = 0;
   
   realT s = 0;
   realT c = 0;
   
   for(size_t i=0;i<_layer_Cn2.size(); ++i)
   {
      _v_wind += _layer_Cn2[i] * pow(_layer_v_wind[i], five_thirds<realT>() );
      s += pow(_layer_v_wind[i], five_thirds<realT>())*sin(_layer_dir[i]) ;//pow( sin(_layer_dir[i]), five_thirds<realT>() );
      c += pow(_layer_v_wind[i], five_thirds<realT>())*cos(_layer_dir[i]) ;//pow( cos(_layer_dir[i]), five_thirds<realT>() );
   }
      
   _v_wind = pow(_v_wind, three_fifths<realT>());
   
   //_dir_wind = atan(pow(s, three_fifths<realT>()) / pow(c, three_fifths<realT>()));
   _dir_wind = atan( s / c);
   if(_dir_wind < 0) _dir_wind += pi<realT>();
   
   _v_wind_updated = true;
   
}

template<typename realT>
void aoAtmosphere<realT>::v_wind(const realT & vw)
{
   if(_v_wind_updated == false) update_v_wind();
   
   realT vw_old = _v_wind;
   
   _v_wind = vw;
   
   //Now update the layers if needed
   if( _layer_v_wind.size() > 0)
   {
      for(size_t i=0; i< _layer_v_wind.size(); ++i)
      {
         _layer_v_wind[i] = _layer_v_wind[i]*(vw/vw_old);
      }
   }
}

template<typename realT>
realT aoAtmosphere<realT>::z_mean()
{
   if(_z_mean_updated == false) update_z_mean();
   return _z_mean;
}

template<typename realT>
void aoAtmosphere<realT>::update_z_mean()
{
   if( _layer_z.size() == 0 || _layer_Cn2.size() == 0)
   {
      _z_mean = 0;
      _z_mean_updated = true;
      return;
   }
   
   _z_mean = 0;
   
   for(size_t i=0;i<_layer_Cn2.size(); ++i)
   {
      _z_mean += _layer_Cn2[i] * pow(_layer_z[i], five_thirds<realT>() );
   }
      
   _z_mean = pow(_z_mean, three_fifths<realT>());
   
   _z_mean_updated = true;
   
}

template<typename realT>
void aoAtmosphere<realT>::z_mean(const realT & zm)
{
   if(_z_mean_updated == false) update_z_mean();
   
   realT zh_old = _z_mean;
   
   _z_mean = zm;
   
   //Now update the layers if needed
   if( _layer_z.size() > 0)
   {
      for(size_t i=0; i< _layer_z.size(); ++i)
      {
         _layer_z[i] = _layer_z[i]*(zm/zh_old);
      }
   }
}

template<typename realT>
realT aoAtmosphere<realT>::alpha()
{
   if(!m_nonKolmogorov)
   {
      return eleven_thirds<realT>();
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
      return constants::a_PSD<realT>()*pow(_r_0, -five_thirds<realT>());
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
    
   for(size_t i=0;i<_layer_Cn2.size(); ++i)
   {
      c += _layer_Cn2[i] * pow( cos(pi<realT>()*k*k*lam_sci * _layer_z[i] * secZ), 2);
   }
   
   return c;
}
 
template<typename realT>
realT aoAtmosphere<realT>::dX(realT f, realT lam_sci, realT lam_wfs)
{
   realT c = 0;
 
   for(size_t i=0;i<_layer_Cn2.size(); ++i)
   {   
      c += _layer_Cn2[i]* pow((cos( pi<realT>()*f*f*lam_sci * _layer_z[i]) - cos( pi<realT>()*f*f*lam_wfs * _layer_z[i])), 2);
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
 
   for(size_t i=0;i<_layer_Cn2.size(); ++i)
   {   
      c += _layer_Cn2[i]*pow(sin( pi<realT>()*k*k*lam_sci * _layer_z[i] * secZ), 2);
   }
   return c;
}

template<typename realT>
realT aoAtmosphere<realT>::dY(realT f, realT lam_sci, realT lam_wfs)
{
   realT c = 0;

   for(size_t i=0;i<_layer_Cn2.size(); ++i)
   {   
      c += _layer_Cn2[i]*pow( (sin(pi<realT>()*f*f*lam_sci * layer_z[i]) - sin( pi<realT>()*f*f*lam_wfs * layer_z[i])), 2);
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
   realT x0 = (n_air(lambda_wfs) - n_air(lambda_i)) * _H*tanZ*secZ; 
   realT x;
   
   for(size_t i = 0; i < _layer_Cn2.size(); ++i)
   {
      x = x0*(1-exp((_layer_z[i]+_h_obs)/_H));
      c += _layer_Cn2[i] * pow( cos(pi<realT>()*k*k*lambda_i * _layer_z[i]*secZ), 2) * pow( sin(pi<realT>()*x*k*cos(0.*3.14/180.)), 2);
   }
   
   return 4*c;
}

template<typename realT>
realT aoAtmosphere<realT>::fwhm(realT lam_sci)
{
   realT r0lam = r_0(lam_sci);
   
   //std::cerr << r0lam << "\n";
   realT fwhm = 0.98*(lam_sci/r0lam);
   //std::cerr << fwhm << "\n";
   
   if( L_0() > 0) fwhm *= sqrt( 1 - 2.183*pow(r0lam/L_0(), 0.356));
   
   return fwhm;
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
   
   L_0(25.0);
   
   h_obs(2400);
}


template<typename realT>
template<typename iosT>
iosT & aoAtmosphere<realT>::dumpAtmosphere( iosT & ios)
{
   ios << "# Atmosphere Parameters:\n";
   ios << "#    r_0 = " << r_0() << '\n';
   ios << "#    lam_0 = " << lam_0() << '\n';
   ios << "#    L_0 = " << L_0() << '\n';
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
 
}//namespace analysis 
}//namespace AO
}//namespace mx

#endif //aoAtmosphere_hpp
