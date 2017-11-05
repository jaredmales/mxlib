/** \file aoSystem.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines an analytical AO system
  * \ingroup mxAO_files
  * 
  */

#ifndef aoSystem_hpp
#define aoSystem_hpp


#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <mx/math/roots.hpp>
#include <mx/mxError.hpp>

#include "aoConstants.hpp"
using namespace mx::AO::constants;

#include "aoAtmosphere.hpp"
#include "aoWFS.hpp"

#define FITTING_ERROR_NO 0
#define FITTING_ERROR_ZERO 1
#define FITTING_ERROR_X 2
#define FITTING_ERROR_Y 3

namespace mx
{
namespace AO
{
namespace analysis 
{
   
/// Describes an analytic adaptive optics (AO) system.
/** 
  * Templatized by the turbulence power spectral density (PSD).
  *
  * \tparam realT the floating point type used for all calculations
  * \tparam inputSpecT specifies the turbulence spatial PSD type
  * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
  * 
  * \ingroup mxAOAnalytic
  */
template<typename realT, class inputSpectT, typename iosT = std::ostream>
class aoSystem
{

public:
      
   aoAtmosphere<realT> atm;
   inputSpectT psd;

protected:
   
   realT _F0; ///< 0 mag flux from star at WFS [photons/sec]
   realT _D; ///< Telescope diameter [m]

   realT _d_min; ///< Minimum AO system actuator pitch [m]
   realT _d_opt; ///< Current optimum AO system actuator pitch [m]
   bool _optd; ///< Flag controlling whether actuator pitch is optimized (true) or just uses _d_min (false).  Default: true.
   
   wfs<realT, iosT> * _wfsBeta; ///< The WFS beta_p class.
   
   realT _lam_wfs; ///< WFS wavelength [m]
   realT _npix_wfs; ///< Number of WFS pixels
   realT _ron_wfs; ///< WFS readout noise [electrons/pix]
   realT _Fbg; ///< Background flux, [photons/sec/pixel]
   
   
   realT _minTauWFS; ///< Minimum WFS exposure time [sec]
   realT _deltaTau; ///< Loop latency [sec]
   bool _optTau; ///< Flag controlling whether optimum integration time is calculated (true) or if _minTauWFS is used (false). Default: true.

   realT _lam_sci; ///< Science wavelength [m]

   realT _zeta; ///<  Zenith angle [radians]
   realT _secZeta; ///< Secant of the Zenith angle (calculated)
   
   int _fit_mn_max; ///< Maximum spatial frequency index to use for fitting error calculation.
   
   realT _ncp_wfe; ///<Static WFE [m rms]
   realT _ncp_alpha; ///< Power-law exponent for the NCP aberations.  Default is 2.
   
   realT _starMag; ///< The magnitude of the star.
   
   bool _specsChanged;///< Flag to indicate that a specification has changed.
   bool _dminChanged;///< Flag to indicate that d_min has changed.
   
   realT _wfeMeasurement; ///< Total WFE due to measurement a error [rad^2 at _lam_sci]
   realT _wfeTimeDelay; ///< Total WFE due to time delay [rad^2 at _lam_sci]
   realT _wfeFitting; ///< Total WFE due to fitting error [rad^2 at _lam_sci]
   realT _wfeChromScintOPD; ///< Total WFE due to the chromaticity of scintillation OPD [rad^2 at lam_sc]
   realT _wfeChromIndex; ///< Total WFE due to the chromaticity of the index of refraction [rad^2 at lam_Sci]
   realT _wfeAnisoOPD; ///< Total WFE due to dispersive anisoplanatism OPD.
   
   realT _wfeNCP; ///< Total WFE due to NCP errors [rad^2 at _lam_sci]
   
   realT _wfeVar; ///< The WFE variance, in meters^2.  Never use this directly, instead use wfeVar().
   realT _strehl; ///<Strehl ratio, a calculated quantity.  Never use this directdy, instead use strehl().
   
   
public:
   
   ///Default c'tor
   /** Calls initialize().
     */
   aoSystem();

   ///Destructor
   ~aoSystem();
   
protected:
   ///Initialize all members.
   void initialize();

public:
      
   ///Load the default parameters from Guyon, 2005 \cite guyon_2005.
   /**  
     *
     */ 
   void loadGuyon2005(); 
   
   ///Load parameters corresponding to the MagAO-X system.
   void loadMagAOX();
   
   ///Load parameters corresponding to the G-MagAO-X system.
   void loadGMagAOX();
   

   /// Set the value of the 0 magnitude photon rate
   /** This is the photon rate at the WFS, photons/sec.
     * 
     */
   void F0(realT nF0 /**< [in] is the new value of _F0.*/);
   
   /// Get the value of the 0 magnitude photon rate
   /**
     * \returns the current value of _F0.
     */ 
   realT F0();
   
   /// Set the value of the Star's magnitude
   /**
     */
   void starMag(realT nmag /**< [in] is the new value of _starMag.*/);
   
   /// Get the value of the Star's magnitude
   /**
     * \returns the current value of _starMag
     */ 
   realT starMag();
   
   ///The photon flux at a given star magnitude.
   /**
     * \returns the photon flux of the star in the bandpass assumed by F0.
     */ 
   realT Fg( realT mag /**< [in] is the magnitude of the star. */);
   
   /// Get the photon rate at the current Star magnitude.
   /** Calculates \f$ F_\gamma = F_0 10^{-0.4 m} \f$ where \f$ F_0 \f$ is _F0 and \f$ m \f$ is _starMag.
     * 
     * \returns the current value of the current photon rate.
     */ 
   realT Fg();
   
   /// Set the value of the primary mirror diameter.
   /** 
     */
   void D( realT nD /**< [in] is the new value of _D. */);
   
   /// Get the value of the primary mirror diamter
   /**
     * \returns the current value of _D.
     */ 
   realT D();
   
   /// Set the value of the minimum subaperture sampling.
   /**
     */
   void d_min( realT nd /**< [in] is the new value of _d_min */);
   
   /// Get the value of the minimum subaperture sampling.
   /**
     * \returns the new value of _d_min.
     */ 
   realT d_min();
   
   /// Set whether or not the value of d is optimized or just set to _d_min.
   /**
     */
   void optd( bool od /**< [in] is the new value of _optd */);
   
   /// Get the value of _optd.
   /**
     * \returns the new value of _optd.
     */ 
   bool optd();
   
   template<typename wfsT>
   void wfsBeta( const wfsT & w)
   {
      _wfsBeta = (wfs<realT,iosT> *) &w;
   }
   
   template<typename wfsT>
   void wfsBeta( const wfsT * w)
   {
      _wfsBeta = (wfs<realT,iosT> *) w;
   }
   
   realT beta_p( realT m, realT n)
   {
      if( _wfsBeta == 0) wfsBetaUnalloc();
      
      return _wfsBeta->beta_p(m, n, _D, d_opt(), atm.r_0(_lam_sci) );
   }
   
   ///Check for unassigned wfs pointer
   /** Prints error and exit()s.
     */
   int wfsBetaUnalloc()
   {
      mxError("aoSystem", MXE_PARAMNOTSET, "The WFS is not assigned.");
      exit(-1);
   }

   /// Set the value of the WFS wavelength.
   /**
     */
   void lam_wfs( realT nlam /**< [in] is the new value of _lam_wfs */);
   
   /// Get the value of the WFS wavelength.
   /**
     * \returns the current value of _lam_wfs.
     */ 
   realT lam_wfs();
   
   /// Set the number of pixels in the WFS
   /**
     */
   void npix_wfs( realT npix /**< [in] is the new value of _npix_wfs */);
   
   /// Get the number of pixels in the WFS
   /**
     * \returns the current value of _npix_wfs
     */ 
   realT npix_wfs();
   
   /// Set the value of the WFS readout noise
   /**
     */
   void ron_wfs( realT nron /**< [in] is the new value of _ron_wfs */);
   
   /// Get the value of the WFS readout noise
   /**
     * \returns the current value of _ron_wfs
     */ 
   realT ron_wfs();
   
   /// Set the value of the background blux.
   /**
     */
   void Fbg(realT fbg /**< [in] is the new value of _Fbg */);
   
   /// Get the value of the background flux.
   /**
     * \returns
     */ 
   realT Fbg();
   
   /// Set the value of the minimum WFS exposure time.
   /**
     */
   void minTauWFS(realT ntau /**< [in] is the new value of _minTauWFS */);
   
   /// Get the value of the minimum WFS exposure time.
   /**
     * \returns the current value of _minTauWFS.
     */ 
   realT minTauWFS();
   
   /// Set the value of _deltaTau.
   /**
     */
   void deltaTau(realT ndel /**< [in] is the new value of _deltaTau*/);
   
   /// Get the value of _deltaTau.
   /**
     * \returns the current value of _deltaTau.
     */ 
   realT deltaTau();

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
   void lam_sci(realT nlam /**< [in] is the new value of _lam_sci */);
   
   /// Get the science wavelength.
   /**
     * \returns the current value of _lam_sci.
     */ 
   realT lam_sci();
   
   
   /// Set the zenith angle, and its secant.
   /** 
     */ 
   void zeta( realT nz /**< [in] The new value of _zeta */ );
   
   /// Get the zenith angle
   /**
     * \return the current value of _zeta
     */ 
   realT zeta();
   
   /// Get the zecant of the zenith angle
   /**
     * \return the current value of _secZeta
     */ 
   realT secZeta();
   
   /// Set the value of _fit_mn_max
   /**
     */
   void fit_mn_max( int mnm /**< [in] is the new value of _fit_mn_max */ );
   
   /// Get the value of _fit_mn_max
   /**
     */
   int fit_mn_max();
   
   
   /// Set the value of the non-common path WFE.
   /**
     */
   void ncp_wfe(realT nwfe /**< [in] is the new value of _ncp_wfe*/);
   
   /// Get the value of the non-common path WFE.
   /**
     * \returns the current value of _ncp_wfe.
     */ 
   realT ncp_wfe();
   

   
   /// Set the value of the non-common path WFE PSD index.
   /**
     */
   void ncp_alpha(realT alpha /**< [in] is the new value of _ncp_alpha*/);
   
   /// Get the value of the non-common path WFE PSD index.
   /**
     * \returns the current value of _ncp_alpha.
     */ 
   realT ncp_alpha();
   
   
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
   realT signal2Noise2( realT & tau_wfs /**< [in/out] specifies the WFS exposure time.  If 0, then optimumTauWFS is used*/);
   
   ///Calculate the measurement noise at a spatial frequency
   /** Calculates the wavefront phase variance due measurement noise at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     *
     * \returns the measurement noise in rad^2 rms at the science wavelength
     */ 
   realT measurementError( realT m, ///< [in] specifies the u component of the spatial frequency  
                            realT n  ///< [in] specifies the v component of the spatial frequency
                          );

   ///Calculate the total measurement error over all corrected spatial frequencies
   /**
     * \returns the total WFE due to measurement error.
     */ 
   realT measurementError();
   
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
   realT timeDelayError( realT m, ///< [in] specifies the u component of the spatial frequency
                          realT n  ///< [in] specifies the v component of the spatial frequency
                        );

   /// Calculate the time delay error over all corrected spatial frequencies
   /**
     * 
     * \returns the total WFE due to time delay.
     */ 
   realT timeDelayError();

   ///@}

   /** \name Fitting Error
     * Calculating the WFE due to uncorrected spatial frequencies.
     * @{
     */
   
   /// Calculate the fitting error at a specific spatial frequency.   
   /**
     * \returns the fitting error in rad^2 rms at the science wavelength at (m,n).
     */ 
   realT fittingError( realT m,  ///< [in] specifies the u component of the spatial frequency
                        realT n   ///< [in] specifies the v component of the spatial frequency
                      );

   /// Calculate the total fitting error over all uncorrected spatial frequencies.
   /**
     *
     * \returns the total fitting error.
     */ 
   realT fittingError();

   ///@}
   
   
   /** \name Chromatic Errors
     * Calculating the WFE due to chromaticity in scintillaion, index of refraction, and dispersion.
     * @{
     */
   
   /// Calculate the wavefront error due to scintillation chromaticity in the OPD at a specific spatial frequency.   
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */    
   realT chromScintOPDError();

   /// Calculate the wavefront error due to scintillation chromaticity in amplitude at a specific spatial frequency.   
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */       
   realT chromScintAmpError();
   
   /// Calculate the wavefront error due to chromaticity in the index of refraction at a specific spatial frequency.   
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */       
   realT chromIndexError();

   /// Calculate the wavefront error due to dispersive anisoplanatism in the OPD at a specific spatial frequency.   
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */       
   realT dispAnisoOPDError();

   /// Calculate the wavefront error due to dispersive anisoplanatism in the amplitude at a specific spatial frequency.   
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */         
   realT dispAnisoAmpError();
   
   
   ///@}
   
   /** \name Optimum Parameters 
     * Functions to calculate the optimum integration time and actuator spacing.
     * 
     * @{
     */
   
   ///Calculate the optimum exposure time for a given spatial frequency.
   /** Finds the optimum exposure time at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \todo Check inclusion of X in parameters
     * 
     * \returns the optimum expsoure time.
     */
   realT optimumTauWFS( realT m, ///< [in] is the spatial frequency index in u
                         realT n  ///< [in] is the spatial frequency index in v
                       );

   ///Calculate the optimum actuator spacing.
   /** Finds the value of _d_opt where the fitting error is less than than the combined time delay and measurement error.
     *
     * \returns the current value of _d_opt. 
     */
   realT d_opt();
   
   /// @}


   

   /// Calculate the NCP variance at a spatial frequency.
   /** Finds the NCP variance at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \returns the NCP variance at k.
     */
   realT ncpError( int m, ///< [in] is the spatial frequency index in u
                   int n  ///< [in] is the spatial frequency index in v
                 );
   
   /// Calculate the total NCP variance in rad^2.
   /**
     * \returns the total NCP variance.
     */ 
   realT ncpError();
   
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
   realT wfeVar();
   
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
   realT strehl();

   /// @}
   
   
   /** \name Contrast
     * Calculating contrast
     * 
     * @{
     */
   
   ///Worker function for raw contrast fuctions
   /** Most of the logic in this calculation is the same, except for the calculation of variance.
     * The varFunc function pointer is used to calculate the variance.  This handles the other details such
     * as Strehl normalization, fitting error (if appropriate), and bounds checking.
     *
     * \tparam varFuncT  is a function-pointer-to-member (of this class) with signature realT(*f)(realT, realT)
     */
   template<typename varFuncT>
   realT C_( realT m, ///< [in] is the spatial frequency index in u/
             realT n, ///< [in] is the spatial frequency index in v.
             bool normStrehl, ///< [in] flag controls whether the contrast is normalized by Strehl ratio/
             varFuncT varFunc, ///< [in] the variance function to use.
             int doFittingError ///< [in] flag to describe how fitting error should be considered for this term: FITTING_ERROR_NO, FITTING_ERROR_ZERO, FITTING_ERROR_X, or FITTING_ERROR_Y.
           );

   ///Worker function for the contrast-map functions.
   /** The map calculation is the same for all terms, except for the calculation of contrast at each pixel.
     * The Cfunc function pointer is used to actually get the contrast.
     * 
     * \tparam imageT is an Eigen-like image
     * \tparam CfuncT is a function-pointer-to-member (of this class) with signature realT (*f)(realT, realT, bool)
     */ 
   template<typename imageT, typename CfuncT>
   void C_Map( imageT & im, ///< [out] the map image to be filled in.
               CfuncT Cfunc ///< [in] the raw contrast function to use for filling in the map. 
             );
   
   ///Calculate the residual variance due to uncorrected phase at a spatial frequency.
   /** Used to calculate contrast \ref C0().
     * 
     * \returns variance at (m,n).
     */
   realT C0var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to uncorrected phase, C0.
   /** Contrast C0 is the uncorrected phase, with the effects of scintillation included.  See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C0.
     */
   realT C0( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in] is the spatial frequency index in v
             bool normStrehl = true ///< [in] flag controls whether the contrast is normalized by Strehl ratio
           );
   
   ///Calculate a 2D map of contrast C0
   /**
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C0Map( imageT & map /**< [in] the map image to be filled in with contrast */ );
   
   ///Calculate the residual variance due to uncorrected amplitude at a spatial frequency.
   /** Used to calculate contrast \ref C1().
     * 
     * \returns variance at (m,n).
     */
   realT C1var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to uncorrected amplitude, C1.
   /** Contrast C1 is the uncorrected amplitude due to scintillation.  See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C0.
     */
   realT C1( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in] is the spatial frequency index in v
             bool normStrehl = true ///< [in] flag controls whether the contrast is normalized by Strehl ratio
           );
   
   ///Calculate a 2D map of contrast C1.
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C1Map( imageT & map /**< [in] the map image to be filled in with contrast */ );
   
   
   ///Calculate the residual variance due to measurement and time delay errors in phase/OPD at a spatial frequency.
   /** Used to calculate contrast \ref C2().
     * 
     * \returns variance at (m,n).
     */
   realT C2var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to measurement and time delay errors in phase/OPD at a spatial frequency.
   /** Contrast C2 is just the total variance due to time delay and measurement errors, 
     * divided by the Strehl ratio.   See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C2.
     */
   realT C2( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in]  is the spatial frequency index in v
             bool normStrehl = true ///< [in] flag controls whether the contrast is normalized by Strehl ratio.
           );

   ///Calculate a 2D map of contrast \ref C2().
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C2Map( imageT & map /**< [in] the map image to be filled in with contrast */ );
   
   ///Calculate the residual variance due to measurement and time delay errors in amplitude at a spatial frequency.
   /** Used to calculate contrast \ref C3().
     * 
     * \returns variance at (m,n).
     */
   realT C3var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to measurement and time delay errors in amplitude at a spatial frequency.
   /** Contrast C3 is just the total variance due to time delay and measurement errors, 
     * divided by the Strehl ratio.    See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C3.
     */
   realT C3( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in]  is the spatial frequency index in v
             bool normStrehl = true ///< [in] flag controls whether the contrast is normalized by Strehl ratio.
           );
   
   ///Calculate a 2D map of contrast C3.
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C3Map( imageT & map /**< [in] the map image to be filled in with contrast */ );

   
   ///Calculate the residual variance due to scintilation-OPD chromaticity.
   /** Used to calculate contrast \ref C4().
     * 
     * \returns variance at (m,n).
     */
   realT C4var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to scintilation-OPD chromaticity.
   /** Contrast C4 is due to the chromaticity of scintillation, causing the ODP measurement to be slightly incorrect at the science wavelength.
     * See Guyon (2005) \cite guyon_2005, and the updated derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C4.
     */
   realT C4( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in]  is the spatial frequency index in v
             bool normStrehl = true ///< [in] flag controls whether the contrast is normalized by Strehl ratio.
           );
   
   ///Calculate a 2D map of contrast C4
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C4Map( imageT & map /**< [in] the map image to be filled in with contrast */ );


   ///Calculate the residual variance due to to scintilation-amplitude chromaticity.
   /** Used to calculate contrast \ref C5().
     * 
     * \returns variance at (m,n).
     */
   realT C5var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to scintilation-amplitude chromaticity.
   /** Contrast C5 is due to the chromaticity of scintillation, causing the amplitude measurement to be slightly incorrect at
     * the science wavelength. See Guyon (2005) \cite guyon_2005, and the updated derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C4.
     */
   realT C5( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in]  is the spatial frequency index in v
             bool normStrehl = true ///< [in] flag controls whether the contrast is normalized by Strehl ratio.
           );
   
   ///Calculate a 2D map of contrast C5
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C5Map( imageT & map /**< [in] the map image to be filled in with contrast */ );

   
   ///Calculate the residual variance due to chromaticity of the index of refraction of air.
   /** Used to calculate contrast \ref C6().
     * 
     * \returns variance at (m,n).
     */
   realT C6var( realT m,  ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to chromaticity of the index of refraction of air.
   /** Contrast C6 is due to the index of refraction of air being wavelength dependent,  causing the ODP measurement to be slightly incorrect
     * at the science wavelength.   See Guyon (2005) \cite guyon_2005, and the updated derivation in Males \& Guyon (2017) \cite males_guyon_2017. 
     * 
     * \returns C6.
     */
   realT C6( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in] is the spatial frequency index in v
             bool normStrehl = true ///< flag controls whether the contrast is normalized by Strehl ratio.
            );
   
   ///Calculate a 2D map of contrast C6
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C6Map( imageT & map /**< [in] the map image to be filled in with contrast */ );
   
   
   ///Calculate the residual variance due to dispersive anisoplanatism.
   /** Used to calculate contrast \ref C7().
     * 
     * \returns variance at (m,n).
     */
   realT C7var( realT m, ///< [in] is the spatial frequency index in u
                realT n ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to dispersive anisoplanatism.
   /** Contrast C7 is due to atmospheric dispersion causing light at different wavelengths to take different paths through atmospheric turbulence.
     * See Fitzgerald (2017, in prep) \cite fitzgerald_2017
     * 
     * \returns C7.
     */
   realT C7( realT m,  ///< [in] is the spatial frequency index in u
             realT n,   ///< [in] is the spatial frequency index in v
             bool normStrehl = true ///< flag controls whether the contrast is normalized by Strehl ratio.
            );
   
   ///Calculate a 2D map of contrast C7
   /** The contrast is Strehl-normalized here.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C7Map( imageT & map /**< [in] the map image to be filled in with contrast */ );

   
   ///@}
   
   ///Output current parameters to a stream
   /** Outputs a formatted list of all current parameters.
     *
     */ 
   iosT & dumpAOSystem( iosT & ios /**< [in] a std::ostream-like stream. */);
   
};



template<typename realT, class inputSpectT, typename iosT>
aoSystem<realT, inputSpectT, iosT>::aoSystem()
{
   initialize();
}

template<typename realT, class inputSpectT, typename iosT>
aoSystem<realT, inputSpectT, iosT>::~aoSystem()
{
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::initialize()
{
   _wfsBeta = 0;
   
   _F0 = 0;
   _D = 0;

   _d_min = 0;
   _d_opt = 0;
   _optd = false;//true;
   
   _lam_wfs = 0;
   _npix_wfs = 0;
   _ron_wfs = 0;
   _Fbg = 0;
   
   _minTauWFS = 0;
   _deltaTau = 0;
   _optTau = true;
   
   _lam_sci = 0;

   _zeta = 0;
   _secZeta = 1;
   
   _fit_mn_max = 100;
   
   _ncp_wfe = 0;
   _ncp_alpha = 2;
   
   _starMag = 0;
   
   _specsChanged = true;
   _dminChanged = true;
   
   _wfeMeasurement = 0;
   _wfeTimeDelay = 0;
   _wfeFitting = 0;
   _wfeChromScintOPD = 0;
   _wfeChromIndex = 0;
   _wfeAnisoOPD = 0;
   
   _wfeNCP = 0;
   
   _wfeVar = 0;
   _strehl = 0;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadGuyon2005()
{
   atm.loadGuyon2005();
   
   F0(1.75e9*0.25*pi<realT>()*64.); //Converting to photons/sec
   lam_wfs(0.55e-6);
   lam_sci(1.6e-6);
   D(8.);
   starMag(5);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadMagAOX()
{   
   atm.loadLCO();   
   
   F0(7.6e10);
   lam_wfs(0.851e-6);
   lam_sci(0.656e-6);
   
   d_min( 6.5/48.0 );
   minTauWFS( 1./3630. );
   
   D(6.5);
   
   _specsChanged = true;
   _dminChanged = true;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadGMagAOX()
{

   loadMagAOX();
   
   
   F0(7.6e10 * 368.0 / (0.25*pi<realT>()*6.5*6.5*(1.-0.29*0.29))); //Scaled up.
   
   D(25.4);
   
   _specsChanged = true;
   _dminChanged = true;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::F0(realT nF0)
{
   _F0 = nF0;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::F0()
{
   return _F0;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::starMag(realT nmag)
{
   _starMag = nmag;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::starMag()
{
   return _starMag;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::Fg(realT mag)
{
   return _F0*pow(10.0, -0.4*mag);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::Fg()
{
   return Fg(_starMag);
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::D(realT nD)
{
   _D = nD;
   psd.D(_D);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::D()
{
   return _D;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::d_min(realT nd)
{
   _d_min = nd;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::d_min()
{
   return _d_min;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::optd(bool od)
{
   _optd = od;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::optd()
{
   return _optd;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::lam_wfs(realT nlam)
{
   _lam_wfs = nlam;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::lam_wfs()
{
   return _lam_wfs;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::npix_wfs(realT npix)
{
   _npix_wfs = npix;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::npix_wfs()
{
   return _npix_wfs;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ron_wfs(realT nron)
{
   _ron_wfs = nron;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ron_wfs()
{
   return _ron_wfs;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::Fbg(realT fbg)
{
   _Fbg = fbg;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::Fbg()
{
   return _Fbg;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::minTauWFS(realT ntau)
{
   _minTauWFS = ntau;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::minTauWFS()
{
   return _minTauWFS;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::deltaTau(realT ndel)
{
   _deltaTau = ndel;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::deltaTau()
{
   return _deltaTau;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::optTau(bool ot)
{
   _optTau = ot;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::optTau()
{
   return _optTau;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::lam_sci(realT nlam)
{
   _lam_sci = nlam;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::lam_sci()
{
   return _lam_sci;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::zeta(realT nz)
{
   _zeta = nz;
   _secZeta = 1/cos(_zeta);
   
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::zeta()
{
   return _zeta;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::secZeta()
{
   return _secZeta;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::fit_mn_max( int mnm )
{
   if(mnm < 0) mnm = 0;
   _fit_mn_max = mnm;
}

template<typename realT, class inputSpectT, typename iosT>
int aoSystem<realT, inputSpectT, iosT>::fit_mn_max()
{
   return _fit_mn_max;
}
   
template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ncp_wfe(realT nwfe)
{
   _ncp_wfe = nwfe;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncp_wfe()
{
   return _ncp_wfe;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ncp_alpha(realT alpha)
{
   _ncp_alpha = alpha;
   _specsChanged = true;
   _dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncp_alpha()
{
   return _ncp_alpha;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::signal2Noise2( realT & tau_wfs )
{      
   realT F = Fg();
               
   return pow(F*tau_wfs,2)/((F+_npix_wfs*_Fbg)*tau_wfs + _npix_wfs*_ron_wfs*_ron_wfs);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementError( realT m, 
                                                                  realT n )
{
   if(m ==0 and n == 0) return 0;
   
   realT tau_wfs;
   
   if(_optTau) tau_wfs = optimumTauWFS(m, n);
   else tau_wfs = _minTauWFS;
   
   if (_wfsBeta == 0) wfsBetaUnalloc();
   
   realT beta_p = _wfsBeta->beta_p(m,n,_D, d_opt(), atm.r_0(_lam_sci));
            
   realT snr2 = signal2Noise2( tau_wfs );
         
  
   return pow(beta_p,2)/snr2*pow(_lam_wfs/_lam_sci, 2);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementError()
{   
   int mn_max = floor(0.5*_D/d_opt());
   realT sum = 0;

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
         
template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayError( realT m, 
                                                    realT n )
{
   if(m ==0 and n == 0) return 0;
   
   realT k = sqrt(m*m + n*n)/_D;
   
   realT tau_wfs;
   
   if(_optTau) tau_wfs = optimumTauWFS(m, n);
   else tau_wfs = _minTauWFS;
   
   realT tau = tau_wfs + _deltaTau;
   
   //std::cout << m << " " << n << " " << tau << "\n";
         
   return psd(atm, k, 0, _secZeta)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.X(k, _lam_sci, _secZeta) * pow(two_pi<realT>()*atm.v_wind()*k,2) * pow(tau,2);      
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayError()
{   
   int mn_max = floor(0.5*_D/d_opt());
   
   realT sum = 0;

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



template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::fittingError( realT m, 
                                                  realT n )
{
   realT k = sqrt(m*m+n*n)/_D;
      
   return psd(atm, k, 0, _secZeta)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::fittingError()
{   
   int mn_max = _D/(2.0*d_opt());

   realT sum = 0;
   
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

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::chromScintOPDError()
{
   int mn_max = floor(0.5*_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += C4var(m,n);
      }
   }

   return sum;
   
}

template<typename realT, class inputSpectT, typename iosT>   
realT aoSystem<realT, inputSpectT, iosT>::chromScintAmpError()
{
   return 0;
#if 0
   int mn_max = floor(0.5*_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += C5var(m,n);
      }
   }

   return sum;
#endif   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::chromIndexError()
{
   int mn_max = floor(0.5*_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += C6var(m,n);
      }
   }

   return sum;
   
}  

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::dispAnisoOPDError()
{
   int mn_max = floor(0.5*_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += C7var(m,n);
      }
   }

   return sum;
   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::dispAnisoAmpError()
{
   return 0;
   
#if 0
   int mn_max = floor(0.5*_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         sum += C8var(m,n);
      }
   }

   return sum;
#endif   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::optimumTauWFS( realT m, 
                                                               realT n )
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
   
   realT k = sqrt(m*m + n*n)/_D;
   
   realT F = Fg();

   if (_wfsBeta == 0) wfsBetaUnalloc();
      
   realT beta_p = _wfsBeta->beta_p(m,n,_D, d_opt(), atm.r_0(_lam_sci));

   //Set up for root finding:
   realT a, b, c, d, e;
   
   realT Atmp = 2*pow(atm.lam_0(),2)*psd(atm, k, 0, _secZeta)/pow(_D,2)*atm.X(k, _lam_sci, _secZeta)*pow(2*pi<realT>()*atm.v_wind()*k,2);
   realT Dtmp = pow(_lam_wfs*beta_p/F,2);
   
   a = Atmp;
   b = Atmp *_deltaTau;
   c = 0;
   d = -Dtmp * (F+_npix_wfs*_Fbg);
   e = -Dtmp * 2*(_npix_wfs*pow(_ron_wfs,2));
   
   std::vector<std::complex<realT> > x;
   
   //Get the roots
   mx::quarticRoots(x, a, b , c, d, e);
   
   //Now pick the largest positive real root
   realT tauopt = 0.0;
   
   for(int i=0; i<4; i++)
   {
      if( real(x[i]) > 0 && imag(x[i]) == 0 && real(x[i]) > tauopt) tauopt = real(x[i]);
   }
   
   if(tauopt < _minTauWFS) tauopt = _minTauWFS;
   
   return tauopt;
   
   
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::d_opt()
{
   if(!_dminChanged) return _d_opt;
   
   if( _d_min == 0 )
   {
      _d_opt = 1e-50;
      _dminChanged = false;
      return _d_opt;
   }
   
   if(!_optd) 
   {
      _d_opt = _d_min;
      _dminChanged = false;
      
      return _d_opt;
   }
   
   realT d = _d_min;
      
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


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncpError( int m, 
                                                          int n )
{
   if(m ==0 and n == 0) return 0;
   
   realT k = sqrt(m*m + n*n)/_D;
   
   return (_ncp_alpha - 2)/(two_pi<realT>()) * pow(_D, -_ncp_alpha) * ncpError() * pow(k, -_ncp_alpha);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncpError()
{
   return pow(_ncp_wfe,2)*pow(2.0*pi<realT>()/_lam_sci,2);
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::calcStrehl()
{  
   _wfeMeasurement = measurementError();
   _wfeTimeDelay = timeDelayError();
   _wfeFitting = fittingError();
   
   _wfeChromScintOPD = chromScintOPDError();
   _wfeChromIndex = chromIndexError();
   _wfeAnisoOPD = dispAnisoOPDError();
   
   _wfeNCP = ncpError();
   
   _wfeVar = _wfeMeasurement + _wfeTimeDelay  + _wfeFitting  + _wfeChromScintOPD +_wfeChromIndex + _wfeAnisoOPD + _wfeNCP;
   
   _strehl = exp(-1 * _wfeVar);
   
   _specsChanged = false;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::wfeVar()
{
   if(_specsChanged || _dminChanged ) calcStrehl();
   
   return _wfeVar;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::strehl()
{
   if( _specsChanged || _dminChanged ) calcStrehl();
   
   return _strehl;
}

template<typename realT, class inputSpectT, typename iosT>
template<typename varFuncT>
realT aoSystem<realT, inputSpectT, iosT>::C_(  realT m, 
                                         realT n,
                                         bool normStrehl,
                                         varFuncT varFunc,
                                         int doFittingError
                                      )
{
   if(m ==0 && n == 0) return 0;
   
   realT S = 1;
   
   if(normStrehl) S = strehl();


   if( doFittingError != FITTING_ERROR_NO)
   {
      int mn_max = _D/(2.*d_opt());
   
      if(mn_max > 0 && (abs(m) > mn_max || abs(n) > mn_max))
      {
         if(doFittingError == FITTING_ERROR_ZERO) return 0;
         
         realT fe = fittingError(m,n);
         
         realT k = sqrt(m*m + n*n)/_D;
         
         if(doFittingError == FITTING_ERROR_X) fe *= atm.X(k, _lam_sci, _secZeta);
         else if(doFittingError == FITTING_ERROR_Y) fe *= atm.Y(k, _lam_sci, _secZeta);
         else
         {
            std::cerr << "Unknown doFittingError\n";
            exit(-1);
         }
         
         return fe / S;
      }
   }
   
   realT var = (this->*varFunc)(m, n);
   
   return var/S;
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT, typename CfuncT>
void aoSystem<realT, inputSpectT, iosT>::C_Map( imageT & im,
                                          CfuncT Cfunc
                                        )
{
   int dim1=im.rows();
   int dim2=im.cols();

   int mc = 0.5*(dim1-1);
   int nc = 0.5*(dim2-1);
   
   int m, n;
   
   int mmax = _D/(2.*d_opt());
   int nmax = mmax;
   
   //std::cerr << dim1 << " " << dim2 << "\n";
   for(int i=0; i< dim1; ++i)
   {
      m = i - mc;
      
      for(int j=0; j< dim2; ++j)
      {
         n = j - nc;
                  
         im(i,j) = (this->*Cfunc)(m, n, true);
      }
   }
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C0var( realT m, 
                                           realT n
                                         )
{
   realT k = sqrt(m*m + n*n)/_D;
   return psd(atm, k, 0, _secZeta)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.X(k, _lam_sci, _secZeta);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C0( realT m, 
                                        realT n,
                                        bool normStrehl
                                      )
{
   
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C0var, FITTING_ERROR_NO);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C0Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C0);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C1var( realT m, 
                                                     realT n
                                                   )
{
   realT k = sqrt(m*m + n*n)/_D;
            
   return psd(atm, k, 0, _secZeta)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.Y(k, _lam_sci, _secZeta);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C1( realT m, 
                                                  realT n,
                                                  bool normStrehl
                                                )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C1var, FITTING_ERROR_NO);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C1Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C1);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C2var(  realT m, 
                                                      realT n
                                                   )
{
   return measurementError(m, n) + timeDelayError(m,n);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C2(  realT m, 
                                                   realT n,
                                                   bool normStrehl
                                                )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C2var, FITTING_ERROR_X);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C2Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C2);   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C4var( realT m, 
                                           realT n
                                         )
{
   realT k = sqrt(m*m + n*n)/_D;
   
   return psd(atm, k, 0, _secZeta)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.dX(k, _lam_sci, _lam_wfs);
}

                                                  
template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C4( realT m, 
                                                  realT n,
                                                  bool normStrehl
                                                )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C4var, FITTING_ERROR_ZERO);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C4Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C4);   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C6var( realT m, 
                                           realT n
                                         )
{
   realT ni = atm.n_air(_lam_sci);
   realT nw = atm.n_air(_lam_wfs);
   
   return C0var(m, n) * pow( (ni-nw)/ni, 2);   
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C6( realT m, 
                                        realT n,
                                        bool normStrehl
                                      )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C6var, FITTING_ERROR_ZERO);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C6Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C6);   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C7var( realT m, 
                                                     realT n
                                                   )
{
   realT k = sqrt(m*m + n*n)/_D;
     
   return psd(atm, k, 0, _secZeta)/pow(_D,2) * pow(atm.lam_0()/_lam_sci, 2) * atm.X_Z(k, _lam_wfs, _lam_sci, _secZeta);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C7( realT m, 
                                                  realT n,
                                                  bool normStrehl
                                                )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C7var, FITTING_ERROR_ZERO);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C7Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C7);   
}

template<typename realT, class inputSpectT, typename iosT>
//template<typename iosT>
iosT & aoSystem<realT, inputSpectT, iosT>::dumpAOSystem( iosT & ios)
{
   ios << "# AO Params:\n";
   ios << "#    D = " << D() << '\n';
   ios << "#    d_min = " << d_min() << '\n';
   ios << "#    optd = " << std::boolalpha << _optd << '\n';
   ios << "#    d_opt = " << d_opt() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    F0 = " << F0() << '\n';
   ios << "#    starMag = " << starMag() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    zeta    = " << zeta() << '\n';
   
   ios << "#    lam_wfs = " << lam_wfs() << '\n';
   ios << "#    npix_wfs = " << npix_wfs() << '\n';
   ios << "#    ron_wfs = " << ron_wfs() << '\n';
   ios << "#    Fbg = " << Fbg() << '\n';
   ios << "#    minTauWFS = " << minTauWFS() << '\n';
   ios << "#    deltaTau = " << deltaTau() << '\n';
   ios << "#    optTau = " << std::boolalpha << _optTau << '\n';
   ios << "#    fit_mn_max = " << _fit_mn_max << '\n';
   ios << "#    ncp_wfe = " << _ncp_wfe << '\n';
   ios << "#    ncp_alpha = " << _ncp_alpha << '\n';
   
   
   if (_wfsBeta == 0) wfsBetaUnalloc();
   _wfsBeta->dumpWFS(ios);
   psd.dumpPSD(ios);
   atm.dumpAtmosphere(ios);
   
   ios << "# Software versions: " << '\n';
   ios << "#    mxlib_comp sha1 = " << mxlib_compiled_git_sha1() << '\n';
   ios << "#    mxlib_comp modified = " << mxlib_compiled_git_repo_modified() << '\n';
   ios << "#    mxlib_uncomp sha1 = " << MXLIB_UNCOMP_CURRENT_SHA1 << '\n';
   ios << "#    mxlib_uncomp modified = " << MXLIB_UNCOMP_REPO_MODIFIED << '\n';
      
   return ios;
}


} //namespace analysis
} //namespace AO
} //namespace mx

#endif //aoSystem_hpp
