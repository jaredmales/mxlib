/** \file aoSystem.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines an analytical AO system
  * \ingroup mxAOAnalytic_files
  * 
  */

#ifndef __aoSystem_hpp__
#define __aoSystem_hpp__


#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

// #include <mx/eigenImage.hpp>
// #include <mx/airy.hpp>
#include <mx/roots.hpp>
#include <mx/mxError.hpp>

#include "aoConstants.hpp"
using namespace mx::AO::constants;

#include "aoAtmosphere.hpp"
#include "aoWFS.hpp"


namespace mx
{
namespace AO
{

/// Describes an analytic adaptive optics (AO) system.
/** 
  * Templatized by the turbulence power spectral density (PSD) and the wavefront sensor (WFS)
  * sensitivity function. 
  *
  * \tparam floatT the floating point type used for all calculations
  * \tparam inputSpecT specifies the turbulence spatial PSD type
  * \tparam wfsBetaT specifies the WFS sensitivity type.  
  */
template<typename floatT, class inputSpectT, class wfsBetaT>
class aoSystem
{

public:
      
   aoAtmosphere<floatT> atm;
   inputSpectT psd;
   wfsBetaT wfsBeta;

   
protected:
   
   floatT _F0; ///< 0 mag flux from star at WFS [photons/sec]
   floatT _D; ///< Telescope diameter [m]

   floatT _d_min; ///< Minimum AO system actuator pitch [m]
   floatT _d_opt; ///< Current optimum AO system actuator pitch [m]
   bool _optd; ///< Flag controlling whether actuator pitch is optimized (true) or just uses _d_min (false).  Default: true.
   
   floatT _lam_wfs; ///< WFS wavelength [m]
   floatT _npix_wfs; ///< Number of WFS pixels
   floatT _ron_wfs; ///< WFS readout noise [electrons/pix]
   floatT _Fbg; ///< Background flux, [photons/sec/pixel]
   
   
   floatT _minTauWFS; ///< Minimum WFS exposure time [sec]
   floatT _deltaTau; ///< Loop latency [sec]
   bool _optTau; ///< Flag controlling whether optimum integration time is calculated (true) or if _minTauWFS is used (false). Default: true.

   floatT _lam_sci; ///< Science wavelength [m]

   floatT _fit_mn_max; ///< Maximum spatial frequency index to use for fitting error calculation.
   
   floatT _ncp_wfe; ///<Static WFE [m rms]
   floatT _ncp_alpha; ///< Power-law exponent for the NCP aberations.  Default is 2.
   
   floatT _starMag; ///< The magnitude of the star.
   
   bool _specsChanged;///< Flag to indicate that a specification has changed.
   bool _dminChanged;///< Flag to indicate that d_min has changed.
   
   floatT _wfeMeasurement; ///< Total WFE due to measurement a error [rad^2 at _lam_sci]
   floatT _wfeTimeDelay; ///< Total WFE due to time delay [rad^2 at _lam_sci]
   floatT _wfeFitting; ///< Total WFE due to fitting error [rad^2 at _lam_sci]
   floatT _wfeNCP; ///< Total WFE due to NCP errors [rad^2 at _lam_sci]
   
   floatT _wfeVar; ///< The WFE variance, in meters^2.  Never use this directly, instead use wfeVar().
   floatT _strehl; ///<Strehl ratio, a calculated quantity.  Never use this directdy, instead use strehl().
   
   
public:
   
   ///Default c'tor
   aoSystem();

   ///Initialize all members.
   void initialize();

   ///Load the default parameters from Guyon, 2005.
   /**  
     *
     */ 
   void load_guyon_2005(); 
   
   void load_GMT_model();
   void load_MagAO_model();

   /// Set the value of the 0 magnitude photon rate
   /** This is the photon rate at the WFS, photons/sec.
     * 
     */
   void F0(floatT nF0 /**< [in] is the new value of _F0.*/);
   
   /// Get the value of the 0 magnitude photon rate
   /**
     * \returns the current value of _F0.
     */ 
   floatT F0();
   
   /// Set the value of the Star's magnitude
   /**
     */
   void starMag(floatT nmag /**< [in] is the new value of _starMag.*/);
   
   /// Get the value of the Star's magnitude
   /**
     * \returns the current value of _starMag
     */ 
   floatT starMag();
   
   ///The photon flux at a given star magnitude.
   /**
     * \returns the photon flux of the star in the bandpass assumed by F0.
     */ 
   floatT Fg( floatT mag /**< [in] is the magnitude of the star. */);
   
   /// Get the photon rate at the current Star magnitude.
   /** Calculates \f$ F_\gamma = F_0 10^{-0.4 m} \f$ where \f$ F_0 \f$ is _F0 and \f$ m \f$ is _starMag.
     * 
     * \returns the current value of the current photon rate.
     */ 
   floatT Fg();
   
   /// Set the value of the primary mirror diameter.
   /** 
     */
   void D( floatT nD /**< [in] is the new value of _D. */);
   
   /// Get the value of the primary mirror diamter
   /**
     * \returns the current value of _D.
     */ 
   floatT D();
   
   /// Set the value of the minimum subaperture sampling.
   /**
     */
   void d_min( floatT nd /**< [in] is the new value of _d_min */);
   
   /// Get the value of the minimum subaperture sampling.
   /**
     * \returns the new value of _d_min.
     */ 
   floatT d_min();
   
   /// Set whether or not the value of d is optimized or just set to _d_min.
   /**
     */
   void optd( bool od /**< [in] is the new value of _optd */);
   
   /// Get the value of _optd.
   /**
     * \returns the new value of _optd.
     */ 
   bool optd();
   
   /// Set the value of the WFS wavelength.
   /**
     */
   void lam_wfs( floatT nlam /**< [in] is the new value of _lam_wfs */);
   
   /// Get the value of the WFS wavelength.
   /**
     * \returns the current value of _lam_wfs.
     */ 
   floatT lam_wfs();
   
   /// Set the number of pixels in the WFS
   /**
     */
   void npix_wfs( floatT npix /**< [in] is the new value of _npix_wfs */);
   
   /// Get the number of pixels in the WFS
   /**
     * \returns the current value of _npix_wfs
     */ 
   floatT npix_wfs();
   
   /// Set the value of the WFS readout noise
   /**
     */
   void ron_wfs( floatT nron /**< [in] is the new value of _ron_wfs */);
   
   /// Get the value of the WFS readout noise
   /**
     * \returns the current value of _ron_wfs
     */ 
   floatT ron_wfs();
   
   /// Set the value of the background blux.
   /**
     */
   void Fbg(floatT fbg /**< [in] is the new value of _Fbg */);
   
   /// Get the value of the background flux.
   /**
     * \returns
     */ 
   floatT Fbg();
   
   /// Set the value of the minimum WFS exposure time.
   /**
     */
   void minTauWFS(floatT ntau /**< [in] is the new value of _minTauWFS */);
   
   /// Get the value of the minimum WFS exposure time.
   /**
     * \returns the current value of _minTauWFS.
     */ 
   floatT minTauWFS();
   
   /// Set the value of _deltaTau.
   /**
     */
   void deltaTau(floatT ndel /**< [in] is the new value of _deltaTau*/);
   
   /// Get the value of _deltaTau.
   /**
     * \returns the current value of _deltaTau.
     */ 
   floatT deltaTau();

   /// Set the value of _optTau.
   /**
     */   
   void optTau( bool ot /**< [in] is the new value of _optTau */);
   
   /// Get the value of _optTau.
   /**
     * \returns the current value of _optTau.
     */ 
   bool optTau();
   
   /// Set the science wavelength.
   /**
     */
   void lam_sci(floatT nlam /**< [in] is the new value of _lam_sci */);
   
   /// Get the science wavelength.
   /**
     * \returns the current value of _lam_sci.
     */ 
   floatT lam_sci();
   
   /// Set the value of the non-common path WFE.
   /**
     */
   void ncp_wfe(floatT nwfe /**< [in] is the new value of _ncp_wfe*/);
   
   /// Get the value of the non-common path WFE.
   /**
     * \returns the current value of _ncp_wfe.
     */ 
   floatT ncp_wfe();
   

   
   /** \name Measurement Error
     * Calculating the WFE due to WFS measurement noise.
     * @{
     */
      
   ///Calculate the signal to noise ratio squared (S/N)^2 for the WFS measurement
   /** The S/N squared is 
      \f[
      (S/N)^2 = \frac{ F_\gamma^2 \tau_{wfs}^2 }{ F_\gamma \tau_{wfs} + n_{pix} F_{bg} \tau_{wfs} + n_{pix} \sigma_{ron}^2 }
      \f]
    
     * 
     * \returns the S/N squared
     */
   floatT signal2Noise2( floatT & tau_wfs /**< [in/out] specifies the WFS exposure time.  If 0, then optimumTauWFS is used*/);
   
   ///Calculate the measurement noise at a spatial frequency
   /** Calculates the wavefront phase variance due measurement noise at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     *
     * \returns the measurement noise in rad^2 rms at the science wavelength
     */ 
   floatT measurementError( floatT m, ///< [in] specifies the u component of the spatial frequency  
                            floatT n  ///< [in] specifies the v component of the spatial frequency
                          );

   ///Calculate the total measurement error over all corrected spatial frequencies
   /**
     * \returns the total WFE due to measurement error.
     */ 
   floatT measurementError();
   
   ///@}
   
   /** \name Time Delay Error
     * Calculating the WFE due to time delay.
     * @{
     */
   
   ///Calculate the time delay at a spatial frequency at the optimum exposure time.
   /** Calculates the wavefront phase variance due to time delay at \f$ f = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \returns the measurement noise in rad^2 rms at the science wavelength
     */ 
   floatT timeDelayError( floatT m, ///< [in] specifies the u component of the spatial frequency
                          floatT n  ///< [in] specifies the v component of the spatial frequency
                        );

   /// Calculate the time delay error over all corrected spatial frequencies
   /**
     * 
     * \returns the total WFE due to time delay.
     */ 
   floatT timeDelayError();

   ///@}

   /** \name Fitting Error
     * Calculating the WFE due to uncorrected spatial frequencies.
     * @{
     */
   
   /// Calculate the fitting error at a specific spatial frequency.   
   /**
     * \returns the fitting error in rad^2 rms at the science wavelength at (m,n).
     */ 

   floatT fittingError( floatT m,  ///< [in] specifies the u component of the spatial frequency
                        floatT n   ///< [in] specifies the v component of the spatial frequency
                      );

   /// Calculate the total fitting error over all uncorrected spatial frequencies.
   /**
     *
     * \returns the total fitting error.
     */ 
   floatT fittingError();

   ///@}
   
   /** \name Optimum Parameters 
     * Functions to calculate the optimum integration time and actuator spacing.
     * 
     * @{
     */
   
   ///Calculate the optimum exposure time for a given spatial frequency.
   /** Finds the optimum exposure time at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \returns the optimum expsoure time.
     */
   floatT optimumTauWFS( floatT m, ///< [in] is the spatial frequency index in u
                         floatT n  ///< [in] is the spatial frequency index in v
                       );

   ///Calculate the optimum actuator spacing.
   /** Finds the value of _d_opt where the fitting error is less than than the combined time delay and measurement error.
     *
     * \returns the current value of _d_opt. 
     */
   floatT d_opt();
   
   /// @}


   

   /// Calculate the NCP variance at a spatial frequency.
   /** Finds the NCP variance at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \returns the NCP variance at k.
     */
   floatT ncpError( int m, ///< [in] is the spatial frequency index in u
                    int n  ///< [in] is the spatial frequency index in v
                  );
   
   /// Calculate the total NCP variance in rad^2.
   /**
     * \returns the total NCP variance.
     */ 
   floatT ncpError();
   
   /** \name Overall WFE and Strehl Ratio 
     * Functions to calculate the total WFE and Strehl ratio.
     * 
     * @{
     */
   
protected:
   ///Calculate the component WFE, total WFE, and Strehl ratio.
   /** This should only be called when something changes.
     */
   void calcStrehl();
   
public:
   ///Get the current value of the total WFE variance.
   /** If no changes, merely returns _wfeVar.  Calls calcStrehl if there are changes.
     *
     * \returns the current value of _wfeVar;
     */
   floatT wfeVar();
   
   ///Get the current Strehl ratio.
   /** If no changes, merely returns _strehl.  Calls calcStrehl if there are changes.
     * Strehl is calculated using the extended Marechal approximation:
     * 
     \f[
      S = e^{-\sigma_{wfe}^2}
      \f]
     * where \f$ \sigma_{wfe}^2 \f$ is the current value of _wfeVar.
     *
     * \returns the current value of _strehl. 
     */
   floatT strehl();

   /// @}
   
   
   /** \name Contrast
     * Calculating contrast
     * 
     * @{
     */
   
   floatT C0( floatT m,  ///< [in] is the spatial frequency index in u
              floatT n   ///< [in] is the spatial frequency index in v
            );
   
   template<typename imageT>
   void C0Map( imageT & im /**< */);
   
   floatT C1( floatT m,  ///< [in] is the spatial frequency index in u
              floatT n   ///< [in] is the spatial frequency index in v
            );
   
   template<typename imageT>
   void C1Map( imageT & im /**< */);
   
   ///Calculate the contrast due to measurement and time delay errors at a spatial frequency.
   /** Contrast C2 is just the total variance due to time delay and measurement errors, 
     * divided by the Strehl ratio.
     * 
     * \param m is the spatial frequency index in u.
     * \param n is the spatial frequency index in v.
     * 
     * \returns C2.
     */
   floatT C2( floatT m,  ///< [in] is the spatial frequency index in u
              floatT n   ///< [in]  is the spatial frequency index in v
            );
   
   template<typename imageT>
   void C2Map( imageT & im /**< */ );
      
   
   
   ///@}
   
   
   template<typename iosT>
   iosT & dumpAOSystem( iosT & ios /**< */);
   
};



template<typename floatT, class inputSpectT, class wfsBetaT>
aoSystem<floatT, inputSpectT, wfsBetaT>::aoSystem()
{
   initialize();
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::initialize()
{
   _F0 = 0;
   _D = 0;

   _d_min = 0;
   _d_opt = 0;
   _optd = true;
   
   _lam_wfs = 0;
   _npix_wfs = 0;
   _ron_wfs = 0;
   _Fbg = 0;
   
   _minTauWFS = 0;
   _deltaTau = 0;
   _optTau = true;
   
   _lam_sci = 0;

   _fit_mn_max = 100;
   
   _ncp_wfe = 0;
   _ncp_alpha = 2;
   
   _starMag = 0;
   
   _specsChanged = true;
   _dminChanged = true;
   
   _wfeMeasurement = 0;
   _wfeTimeDelay = 0;
   _wfeFitting = 0;
   _wfeNCP = 0;
   
   _wfeVar = 0;
   _strehl = 0;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::load_guyon_2005()
{
   atm.load_Guyon_2005();
   
   _F0 = 1.75e9;
   _lam_wfs = 0.55e-6;
   _lam_sci = 1.6e-6;
   D(8.);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::load_GMT_model()
{
   atm.load_GMT_model();
   
   
   _F0 = _F0 = (6.7e10)*0.25* 0.75 * 0.8* (0.25*pi<floatT>()*24.5*24.5*(1.-0.29*0.29));
   _lam_wfs = 0.675e-6;
   _lam_sci = 0.9e-6;
   D(25.4);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::load_MagAO_model()
{   
   atm.load_MagAO_model();   
   
   D(6.5);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::F0(floatT nF0)
{
   _F0 = nF0;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::F0()
{
   return _F0;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::starMag(floatT nmag)
{
   _starMag = nmag;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::starMag()
{
   return _starMag;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::Fg(floatT mag)
{
   return _F0*pow(10.0, -0.4*mag);
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::Fg()
{
   return Fg(_starMag);
}


template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::D(floatT nD)
{
   _D = nD;
   psd.D(_D);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::D()
{
   return _D;
}


template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::d_min(floatT nd)
{
   _d_min = nd;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::d_min()
{
   return _d_min;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::optd(bool od)
{
   _optd = od;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
bool aoSystem<floatT, inputSpectT, wfsBetaT>::optd()
{
   return _optd;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::lam_wfs(floatT nlam)
{
   _lam_wfs = nlam;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::lam_wfs()
{
   return _lam_wfs;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::npix_wfs(floatT npix)
{
   _npix_wfs = npix;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::npix_wfs()
{
   return _npix_wfs;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::ron_wfs(floatT nron)
{
   _ron_wfs = nron;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::ron_wfs()
{
   return _ron_wfs;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::Fbg(floatT fbg)
{
   _Fbg = fbg;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::Fbg()
{
   return _Fbg;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::minTauWFS(floatT ntau)
{
   _minTauWFS = ntau;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::minTauWFS()
{
   return _minTauWFS;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::deltaTau(floatT ndel)
{
   _deltaTau = ndel;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::deltaTau()
{
   return _deltaTau;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::optTau(bool ot)
{
   _optTau = ot;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
bool aoSystem<floatT, inputSpectT, wfsBetaT>::optTau()
{
   return _optTau;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::lam_sci(floatT nlam)
{
   _lam_sci = nlam;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::lam_sci()
{
   return _lam_sci;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::ncp_wfe(floatT nwfe)
{
   _ncp_wfe = nwfe;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::ncp_wfe()
{
   return _ncp_wfe;
}


template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::signal2Noise2( floatT & tau_wfs )
{      
   floatT F = Fg();
               
   return pow(F*tau_wfs,2)/((F+_npix_wfs*_Fbg)*tau_wfs + _npix_wfs*_ron_wfs*_ron_wfs);
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::measurementError( floatT m, 
                                                                  floatT n )
{
   if(m ==0 and n == 0) return 0;
   
   floatT tau_wfs;
   
   if(_optTau) tau_wfs = optimumTauWFS(m, n);
   else tau_wfs = _minTauWFS;
   
   floatT beta_p = wfsBeta.beta_p(m,n,_D);
            
   floatT snr2 = signal2Noise2( tau_wfs );
         
   return pow(beta_p,2)/snr2*pow(_lam_wfs/_lam_sci, 2);
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::measurementError()
{   
   int mn_max = floor(0.5*_D/d_opt());
   floatT sum = 0;

   for(int m=-mn_max; m <= mn_max; ++m)
   {
      for(int n=-mn_max; n<mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += measurementError(m,n);
      }
   }

   return sum;
}
         
template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::timeDelayError( floatT m, 
                                                                floatT n )
{
   if(m ==0 and n == 0) return 0;
   
   floatT k = sqrt(m*m + n*n)/_D;
   
   floatT tau_wfs;
   
   if(_optTau) tau_wfs = optimumTauWFS(m, n);
   else tau_wfs = _minTauWFS;
   
   floatT tau = tau_wfs + _deltaTau;
   
   //std::cout << m << " " << n << " " << tau << "\n";
         
   return psd(atm, k, _lam_sci)/pow(_D,2) * pow(two_pi<floatT>()*atm.v_wind()*k,2) * pow(tau,2);      
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::timeDelayError()
{   
   int mn_max = floor(0.5*_D/d_opt());
   
   floatT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += timeDelayError(m,n);
      }
   }

   return sum;
}



template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::fittingError( floatT m, 
                                                              floatT n )
{
   floatT k = sqrt(m*m+n*n)/_D;
      
   return psd(atm,k, _lam_sci)/(_D*_D);
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::fittingError()
{   
   int mn_max = _D/(2.0*d_opt());

   floatT sum = 0;
   
   for(int m = -_fit_mn_max; m <= _fit_mn_max; ++m)
   {
      for(int n = -_fit_mn_max; n <= _fit_mn_max; ++n)
      {
         if( abs(m) <= mn_max && abs(n) <= mn_max) continue;
         
         sum += fittingError(m,n);
      }
   }
            
   return sum;
}


template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::optimumTauWFS( floatT m, 
                                                               floatT n )
{
   if(_D == 0)
   {
      mxError("aoSystem::optimumTauWFS", MXE_PARAMNOTSET, "Diameter (D) not set.");
      return -1;
   }
   
   if(_F0 == 0)
   {
      mxError("aoSystem::optimumTauWFS", MXE_PARAMNOTSET, "0-mag photon flux (F0) not set.");
      return -1;
   }
   
   
   //if(m == 0) return 0; //Just a big number
   
   floatT k = sqrt(m*m + n*n)/_D;
   
   floatT F = Fg();

   floatT beta_p = wfsBeta.beta_p(m,n,_D);

   //Set up for root finding:
   floatT a, b, c, d, e;
   
   floatT Atmp = 2*pow(atm.lam_0(),2)*psd(atm,k)/pow(_D,2)*pow(2*pi<floatT>()*atm.v_wind()*k,2);
   floatT Dtmp = pow(_lam_wfs*beta_p/F,2);
   
   a = Atmp;
   b = Atmp *_deltaTau;
   c = 0;
   d = -Dtmp * (F+_npix_wfs*_Fbg);
   e = -Dtmp * 2*(_npix_wfs*pow(_ron_wfs,2));
   
   std::vector<std::complex<floatT> > x;
   
   //Get the roots
   mx::quarticRoots(x, a, b , c, d, e);
   
   //Now pick the largest positive real root
   floatT tauopt = 0.0;
   
   for(int i=0; i<4; i++)
   {
      if( real(x[i]) > 0 && imag(x[i]) == 0 && real(x[i]) > tauopt) tauopt = real(x[i]);
   }
   
   if(tauopt < _minTauWFS) tauopt = _minTauWFS;
   
   return tauopt;
   
   
}


template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::d_opt()
{
   if(!_dminChanged) return _d_opt;
   
   if(!_optd) 
   {
      _d_opt = _d_min;
      _dminChanged = false;
      
      return _d_opt;
   }
   
   floatT d = _d_min;
      
   int m = _D/(2*d);
   int n = 0;
   
   while( measurementError(m,n) + timeDelayError(m,n) > fittingError(m, n) && d < _D/2 ) 
   {
      d += _d_min/1000.0;
      m = _D/(2*d);
      n = 0;
   }
   
   _d_opt = d;
   _dminChanged = false;
   
   return _d_opt;
}


template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::ncpError( int m, 
                                                          int n )
{
   if(m ==0 and n == 0) return 0;
   
   floatT k = sqrt(m*m + n*n)/_D;
   
   return (_ncp_alpha - 2)/(two_pi<floatT>()) * pow(_D, -_ncp_alpha) * ncpError() * pow(k, -_ncp_alpha);
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::ncpError()
{
   return pow(_ncp_wfe,2)*pow(2.0*pi<floatT>()/_lam_sci,2);
}

template<typename floatT, class inputSpectT, class wfsBetaT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::calcStrehl()
{  
   _wfeMeasurement = measurementError();
   _wfeTimeDelay = timeDelayError();
   _wfeFitting = fittingError();
   _wfeNCP = ncpError();
   
   //std::cerr << _wfeMeasurement << " " << _wfeTimeDelay << " " << _wfeFitting << " " << _wfeNCP << "\n";
   
   _wfeVar = _wfeMeasurement + _wfeTimeDelay  + _wfeFitting  + _wfeNCP;
   
   _strehl = exp(-1 * _wfeVar);
   
   _specsChanged = false;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::wfeVar()
{
   if(_specsChanged || _dminChanged ) calcStrehl();
   
   return _wfeVar;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::strehl()
{
   if( _specsChanged || _dminChanged ) calcStrehl();
   
   return _strehl;
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::C0(  floatT m, 
                                                     floatT n )
{
   if(m ==0 && n == 0) return 0;

   floatT k = sqrt(m*m + n*n)/_D;
            
   floatT sig = psd(atm,k)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.X(k, _lam_sci);
   
   return sig/strehl();
}

template<typename floatT, class inputSpectT, class wfsBetaT>
template<typename imageT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::C0Map( imageT & im )
{
   int dim1=im.rows();
   int dim2=im.cols();

   int mc = 0.5*(dim1-1);
   int nc = 0.5*(dim2-1);
   
   int m, n;
   
   int mmax = _D/(2.*d_opt());
   int nmax = mmax;
   
   for(int i=0; i< dim1; ++i)
   {
      m = i - mc;
      
      for(int j=0; j< dim2; ++j)
      {
         n = j - nc;
                  
         im(i,j) = C0(m, n );
      }
   }
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::C1(  floatT m, 
                                                     floatT n )
{
   if(m ==0 && n == 0) return 0;

   floatT k = sqrt(m*m + n*n)/_D;
            
   floatT sig = psd(atm,k)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.Y(k, _lam_sci);
   
   return sig/strehl();
}

template<typename floatT, class inputSpectT, class wfsBetaT>
template<typename imageT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::C1Map( imageT & im )
{
   int dim1=im.rows();
   int dim2=im.cols();

   int mc = 0.5*(dim1-1);
   int nc = 0.5*(dim2-1);
   
   int m, n;
   
   int mmax = _D/(2.*d_opt());
   int nmax = mmax;
   
   for(int i=0; i< dim1; ++i)
   {
      m = i - mc;
      
      for(int j=0; j< dim2; ++j)
      {
         n = j - nc;
                  
         im(i,j) = C1(m, n );
      }
   }
}

template<typename floatT, class inputSpectT, class wfsBetaT>
floatT aoSystem<floatT, inputSpectT, wfsBetaT>::C2(  floatT m, 
                                                     floatT n )
{
   if(m ==0 && n == 0) return 0;

   int mn_max = _D/(2.*d_opt());
   
   if(abs(m) > mn_max || abs(n) > mn_max)
   {
      return fittingError( m, n )/strehl();
   }
         
   floatT sig = measurementError(m, n) + timeDelayError(m,n);
   
   return sig/strehl();
}

template<typename floatT, class inputSpectT, class wfsBetaT>
template<typename imageT>
void aoSystem<floatT, inputSpectT, wfsBetaT>::C2Map( imageT & im )
{
   int dim1=im.rows();
   int dim2=im.cols();

   int mc = 0.5*(dim1-1);
   int nc = 0.5*(dim2-1);
   
   int m, n;
   
   int mmax = _D/(2.*d_opt());
   int nmax = mmax;
   
   for(int i=0; i< dim1; ++i)
   {
      m = i - mc;
      
      for(int j=0; j< dim2; ++j)
      {
         
         n = j - nc;

         if( m==0 && n == 0)
         {
            im(i,j) = 0;
            continue;
         }
         if(m == 0) m = 1;
                  
         im(i,j) = C2(m, n );
      }
   }
}







   
   
   
   
template<typename floatT, class inputSpectT, class wfsBetaT>
template<typename iosT>
iosT & aoSystem<floatT, inputSpectT, wfsBetaT>::dumpAOSystem( iosT & ios)
{
   ios << "# AO Params:\n";
   ios << "#    D = " << D() << '\n';
   ios << "#    d_min = " << d_min() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    F0 = " << F0() << '\n';
   ios << "#    starMag = " << starMag() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   
   ios << "#    lam_wfs = " << lam_wfs() << '\n';
   ios << "#    npix_wfs = " << npix_wfs() << '\n';
   ios << "#    ron_wfs = " << ron_wfs() << '\n';
   ios << "#    Fbg = " << Fbg() << '\n';
   ios << "#    minTauWFS = " << minTauWFS() << '\n';
   ios << "#    deltaTau = " << deltaTau() << '\n';
   
   wfsBeta.dumpWFS(ios);
   psd.dumpPSD(ios);
   atm.dumpAtmosphere(ios);
   
   ios << "# Software versions: " << '\n';
   ios << "#    mxlib_comp sha1 = " << mxlib_compiled_git_sha1() << '\n';
   ios << "#    mxlib_comp modified = " << mxlib_compiled_git_repo_modified() << '\n';
   ios << "#    mxlib_uncomp sha1 = " << MXLIB_UNCOMP_CURRENT_SHA1 << '\n';
   ios << "#    mxlib_uncomp modified = " << MXLIB_UNCOMP_REPO_MODIFIED << '\n';
//   ios << "#    mx::AOAnalytic sha1 = " << MXAOANALYTIC_CURRENT_SHA1 << '\n';
//   ios << "#    mx::AOanalytic modified = " << MXAOANALYTIC_REPO_MODIFIED << '\n';
      
   return ios;
}



} //namespace AO
} //namespace mx

#endif //__aoSystem_hpp__
