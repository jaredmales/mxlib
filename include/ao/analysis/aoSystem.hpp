/** \file aoSystem.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines an analytical AO system
  * \ingroup mxAO_files
  * 
  */

//***********************************************************************//
// Copyright 2015-2018 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef aoSystem_hpp
#define aoSystem_hpp


#include "../../math/constants.hpp"
#include "../../math/roots.hpp"
#include "../../mxError.hpp"
#include "../../mxException.hpp"

#include "aoConstants.hpp"


#include "aoAtmosphere.hpp"
#include "aoPSDs.hpp"
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
   
   realT m_F0 {0}; ///< 0 mag flux from star at WFS [photons/sec]
   realT m_D {0}; ///< Telescope diameter [m]

   realT m_d_min {0}; ///< Minimum AO system actuator pitch [m]
   realT m_d_opt {0}; ///< Current optimum AO system actuator pitch [m]
   bool m_optd {false}; ///< Flag controlling whether actuator pitch is optimized (true) or just uses m_d_min (false).  Default: true.
   realT m_optd_delta {1.0}; ///< The fractional change from d_min used in optimization.  Set to 1 for integer binnings, > 1 for finer sampling.
   
   wfs<realT, iosT> * m_wfsBeta {nullptr}; ///< The WFS beta_p class.
   
   realT m_lam_wfs {0}; ///< WFS wavelength [m]

   //WFS detector parameters, which may be a function of binning
   std::vector<realT> m_npix_wfs;  ///< Number of WFS pixels
   std::vector<realT> m_ron_wfs;   ///< WFS readout noise [electrons/pix]
   std::vector<realT> m_Fbg;       ///< Background flux, [photons/sec/pixel]
   std::vector<realT> m_minTauWFS; ///< Minimum WFS exposure time [sec]
   
   bool m_bin_npix {true}; ///< Flag controlling whether or not to bin WFS pixels according to the actuator spacing.
      
   realT m_tauWFS {0}; ///< Actual WFS exposure time [sec]
   
   realT m_deltaTau {0}; ///< Loop latency [sec]
   
   bool m_optTau {true}; ///< Flag controlling whether optimum integration time is calculated (true) enforcing m_minTauWFS, or if m_tauWFS is used (false). Default: true.

   realT m_lam_sci {0}; ///< Science wavelength [m]

   realT m_zeta {0}; ///<  Zenith angle [radians]
   realT m_secZeta {1}; ///< Secant of the Zenith angle (calculated)
   
   int m_fit_mn_max {100}; ///< Maximum spatial frequency index to use for fitting error calculation.
   
   realT m_spatialFilter_ku {std::numeric_limits<realT>::max()}; ///< The spatial filter cutoff in u, [m^-1]
   realT m_spatialFilter_kv {std::numeric_limits<realT>::max()}; ///< The spatial filter cutoff in v, [m^-1]
   
   realT m_ncp_wfe {0}; ///<Static WFE [m rms]
   realT m_ncp_alpha {2}; ///< Power-law exponent for the NCP aberations.  Default is 2.
   
   realT m_starMag {0}; ///< The magnitude of the star.
   
   bool m_specsChanged {true};///< Flag to indicate that a specification has changed.
   bool m_dminChanged {true};///< Flag to indicate that d_min has changed.
   
   bool m_circularLimit {false}; ///< Flag to indicate that the spatial frequency limit is circular, not square.
   
   realT m_wfeMeasurement {0}; ///< Total WFE due to measurement a error [rad^2 at m_lam_sci]
   realT m_wfeTimeDelay {0}; ///< Total WFE due to time delay [rad^2 at m_lam_sci]
   realT m_wfeFitting {0}; ///< Total WFE due to fitting error [rad^2 at m_lam_sci]
   realT m_wfeChromScintOPD {0}; ///< Total WFE due to the chromaticity of scintillation OPD [rad^2 at lam_sc]
   realT m_wfeChromIndex {0}; ///< Total WFE due to the chromaticity of the index of refraction [rad^2 at lam_Sci]
   realT m_wfeAnisoOPD {0}; ///< Total WFE due to dispersive anisoplanatism OPD.
   
   realT m_wfeNCP  {0}; ///< Total WFE due to NCP errors [rad^2 at m_lam_sci]
   
   realT m_wfeVar {0}; ///< The WFE variance, in meters^2.  Never use this directly, instead use wfeVar().
   realT m_strehl {0}; ///<Strehl ratio, a calculated quantity.  Never use this directdy, instead use strehl().
   
   
public:
   
   ///Default c'tor
   /** Calls initialize().
     */
   aoSystem();

   ///Destructor
   ~aoSystem();
   

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
   void F0(realT nF0 /**< [in] is the new value of m_F0.*/);
   
   /// Get the value of the 0 magnitude photon rate
   /**
     * \returns the current value of m_F0.
     */ 
   realT F0();
   
   /// Set the value of the Star's magnitude
   /**
     */
   void starMag(realT nmag /**< [in] is the new value of m_starMag.*/);
   
   /// Get the value of the Star's magnitude
   /**
     * \returns the current value of m_starMag
     */ 
   realT starMag();
   
   int circularLimit( bool cl );
   
   bool circularLimit();
   
   ///The photon flux at a given star magnitude.
   /**
     * \returns the photon flux of the star in the bandpass assumed by F0.
     */ 
   realT Fg( realT mag /**< [in] is the magnitude of the star. */);
   
   /// Get the photon rate at the current Star magnitude.
   /** Calculates \f$ F_\gamma = F_0 10^{-0.4 m} \f$ where \f$ F_0 \f$ is m_F0 and \f$ m \f$ is m_starMag.
     * 
     * \returns the current value of the current photon rate.
     */ 
   realT Fg();
   
   /// Set the value of the primary mirror diameter.
   /** 
     */
   void D( realT nD /**< [in] is the new value of m_D. */);
   
   /// Get the value of the primary mirror diamter
   /**
     * \returns the current value of m_D.
     */ 
   realT D();
   
   /// Set the value of the minimum subaperture sampling.
   /**
     */
   void d_min( realT nd /**< [in] is the new value of m_d_min */);
   
   /// Get the value of the minimum subaperture sampling.
   /**
     * \returns the new value of m_d_min.
     */ 
   realT d_min();
   
   /// Set whether or not the value of d is optimized or just set to m_d_min.
   /**
     */
   void optd( bool od /**< [in] is the new value of m_optd */);
   
   /// Get the value of m_optd.
   /**
     * \returns the new value of m_optd_delta.
     */ 
   bool optd();
   
   /// Set the fractional change in actuator spacing for optimization.
   /** Sets the fraction of m_d_min by which the optimizer changes actautor spacing.
     */
   void optd_delta( realT odd /**< [in] is the new value of m_optd_delta */);
   
   /// Get the value of the fractional change in actuator spacing for optimization..
   /**
     * \returns the value of m_optd_delta.
     */ 
   realT optd_delta();
   
   template<typename wfsT>
   void wfsBeta( const wfsT & w)
   {
      m_wfsBeta = (wfs<realT,iosT> *) &w;
   }
   
   template<typename wfsT>
   void wfsBeta( const wfsT * w)
   {
      m_wfsBeta = (wfs<realT,iosT> *) w;
   }
   
   realT beta_p( realT m, realT n)
   {
      if( m_wfsBeta == 0) wfsBetaUnalloc();
      
      return m_wfsBeta->beta_p(m, n, m_D, d_opt(), atm.r_0(m_lam_sci) );
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
   void lam_wfs( realT nlam /**< [in] is the new value of m_lam_wfs */);
   
   /// Get the value of the WFS wavelength.
   /**
     * \returns the current value of m_lam_wfs.
     */ 
   realT lam_wfs();
   
   /// Set the number of pixels in the WFS
   /** This sets the vector to length 1 with a single value
     */
   void npix_wfs( realT npix /**< [in] is the new single value of m_npix_wfs */);
   
   /// Set the number of pixels in the WFS for each binning mode
   /** This sets the vector.
     */
   void npix_wfs( const std::vector<realT> & npix /**< [in] is the new values of m_npix_wfs */);

   /// Set the number of pixels in the WFS for a given binning mode
   /** This sets the entry in the vector, if it has the proper length.
     */
   void npix_wfs( int idx,   ///< [in] the index of the binning mode
                  realT npix ///< [in] is the new value of m_npix_wfs[idx] 
                );

   /// Get the number of pixels in the WFS for a given binning mode
   /**
     * \returns the current value of m_npix_wfs[idx]
     */ 
   realT npix_wfs(size_t idx /**< [in] the index of the binning mode.*/);
   
   /// Get the number of pixels in the WFS for all binning modes
   /**
     * \returns the current value of m_npix_wfs
     */ 
   std::vector<realT> npix_wfs();

   /// Set the value of the WFS readout noise
   /** This sets the vector to length 1 with a single value
     */
   void ron_wfs( realT nron /**< [in] is the new single value of m_ron_wfs */);
   
   /// Set the value of the WFS readout noise for each binning mode
   /** This sets the vector.
     */
   void ron_wfs( const std::vector<realT> & nron /**< [in] is the new value of m_ron_wfs */);

   /// Set the value of the WFS readout noise for a given bining mode
   /** This sets the entry in the vector if the vector has the proper length
     */
   void ron_wfs( int idx,   ///< [in] the index of the binning mode
                 realT nron ///< [in] is the new value of m_ron_wfs [idx]
               );

   /// Get the value of the WFS readout noise for a given binning mode
   /**
     * \returns the current value of m_ron_wfs[idx]
     */ 
   realT ron_wfs( size_t idx /**< [in] the index of the binning mode.*/);
   
   /// Get the value of the WFS readout noise for all binning modes
   /**
     * \returns the current value of m_ron_wfs
     */ 
   std::vector<realT> ron_wfs();

   /// Set a single value of the background flux
   /** This sets the vector to length 1 with a single value
     */
   void Fbg(realT fbg /**< [in] is the new value of m_Fbg */);
   
   /// Set the value of the background fluxes.
   /** This sets the value of the background flux for each binning mode
     */
   void Fbg(const std::vector<realT> & fbg /**< [in] is the new value of m_Fbg */);

   /// Set a single value of the background flux for a given binning mode
   /** This sets the entry in the vector if the vector has the proper length
     */
   void Fbg( int idx,  ///< [in] the index of the binning mode
             realT fbg ///< [in] is the new value of m_Fbg[idx]
            );

   /// Get the value of the background flux for a given binning mode
   /**
     * \returns m_Fbg[idx]
     */ 
   realT Fbg( size_t idx /**< [in] the index of the binning mode.*/);
   
   /// Get the value of the background flux for all binning modes
   /**
     * \returns the current value of m_Fbg
     */ 
   std::vector<realT> Fbg();
   
   /// Set a single value of the minimum WFS exposure time.
   /** This sets the vector to length 1 with a single value
     */
   void minTauWFS(realT ntau /**< [in] is the new value of m_minTauWFS */);
   
   /// Set the value of the minimum WFS exposure times.
   /** This sets the vector of the minimum WFS exposure times
     */
   void minTauWFS(const std::vector<realT> &  ntau /**< [in] is the new value of m_minTauWFS */);

   /// Set a single value of the minimum WFS exposure time for a given binning mode.
   /** This sets the entry in the vector if the vector has sufficient length.
     */
   void minTauWFS(size_t idx, ///< [in] the index of the binning mode
                  realT ntau  ///< [in] is the new value of m_minTauWFS
                 );

   /// Get the value of the minimum WFS exposure time for a given binning.
   /**
     * \returns the current value of m_minTauWFS[idx].
     */ 
   realT minTauWFS(size_t idx /**< [in] the index of the binning mode.*/);
   
   /// Get the values of the minimum WFS exposure time.
   /**
     * \returns the current value of m_minTauWFS.
     */ 
   std::vector<realT> minTauWFS();

   /// Set the value of the pixel binning flag
   /**
     */
   void bin_npix(bool bnp /**< [in] is the new value of m_bin_npix */);
   
   /// Get the value of the pixel binngin flag
   /**
     * \returns
     */ 
   bool bin_npix();
   
   
   
   /// Set the value of the WFS exposure time.
   /**
     */
   void tauWFS(realT ntau /**< [in] is the new value of m_tauWFS */);
   
   /// Get the value of the minimum WFS exposure time.
   /**
     * \returns the current value of m_tauWFS.
     */ 
   realT tauWFS();
   
   /// Set the value of m_deltaTau.
   /**
     */
   void deltaTau(realT ndel /**< [in] is the new value of m_deltaTau*/);
   
   /// Get the value of m_deltaTau.
   /**
     * \returns the current value of m_deltaTau.
     */ 
   realT deltaTau();

   /// Set the value of m_optTau.
   /**
     */   
   void optTau( bool ot /**< [in] is the new value of m_optTau */);
   
   /// Get the value of m_optTau.
   /**
     * \returns the current value of m_optTau.
     */ 
   bool optTau();
   
   /// Set the science wavelength.
   /**
     */
   void lam_sci(realT nlam /**< [in] is the new value of m_lam_sci */);
   
   /// Get the science wavelength.
   /**
     * \returns the current value of m_lam_sci.
     */ 
   realT lam_sci();
   
   
   /// Set the zenith angle, and its secant.
   /** 
     */ 
   void zeta( realT nz /**< [in] The new value of m_zeta */ );
   
   /// Get the zenith angle
   /**
     * \return the current value of m_zeta
     */ 
   realT zeta();
   
   /// Get the zecant of the zenith angle
   /**
     * \return the current value of m_secZeta
     */ 
   realT secZeta();
   
   /// Set the value of m_fit_mn_max
   /**
     */
   void fit_mn_max( int mnm /**< [in] is the new value of m_fit_mn_max */ );
   
   /// Get the value of m_fit_mn_max
   /**
     */
   int fit_mn_max();
   
   /// Set the value of spatialFilter_ku
   /**
     */
   void spatialFilter_ku( realT ku /**< [in] is the new value of spatialFilter_ku*/ );
   
   /// Get the value of spatialFilter_ku
   /** \returns the current value of m_spatialFilter_ku
     */
   realT spatialFilter_ku();
   
   /// Set the value of spatialFilter_kv
   /**
     */
   void spatialFilter_kv( realT kv /**< [in] is the new value of spatialFilter_kv*/ );
   
   /// Get the value of spatialFilter_kv
   /** \returns the current value of m_spatialFilter_kv
     */
   realT spatialFilter_kv();
   
   /// Set the value of the non-common path WFE.
   /**
     */
   void ncp_wfe(realT nwfe /**< [in] is the new value of m_ncp_wfe*/);
   
   /// Get the value of the non-common path WFE.
   /**
     * \returns the current value of m_ncp_wfe.
     */ 
   realT ncp_wfe();
   

   
   /// Set the value of the non-common path WFE PSD index.
   /**
     */
   void ncp_alpha(realT alpha /**< [in] is the new value of m_ncp_alpha*/);
   
   /// Get the value of the non-common path WFE PSD index.
   /**
     * \returns the current value of m_ncp_alpha.
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
   realT signal2Noise2( realT & tau_wfs, ///< [in/out] specifies the WFS exposure time.  If 0, then optimumTauWFS is used
                        realT d          ///< [in] the actuator spacing in meters, used if binning WFS pixels
                      );
   
   ///Calculate the measurement noise at a spatial frequency and specified actuator spacing
   /** Calculates the wavefront phase variance due measurement noise at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     *
     * \returns the measurement noise in rad^2 rms at the science wavelength
     */ 
   realT measurementError( realT m, ///< [in] specifies the u component of the spatial frequency  
                           realT n, ///< [in] specifies the v component of the spatial frequency
                           realT d  ///< [in] the actuator spacing in meters
                         );
   
   ///Calculate the measurement noise at a spatial frequency using the optimum actuator spacing
   /** 
     * Calculates the wavefront phase variance due measurement noise at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * If optd == false this uses the minimum actuator spacing.
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
   
   ///Calculate the time delay at a spatial frequency at the optimum exposure time and the specified actuator spacing
   /** Calculates the wavefront phase variance due to time delay at \f$ f = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \returns the measurement noise in rad^2 rms at the science wavelength
     */ 
   realT timeDelayError( realT m, ///< [in] specifies the u component of the spatial frequency
                         realT n, ///< [in] specifies the v component of the spatial frequency
                         realT d  ///< [in] the actuator spacing, in meters
                       );
   
   ///Calculate the time delay at a spatial frequency at the optimum exposure time at the optimum actuator spacing.
   /** Calculates the wavefront phase variance due to time delay at \f$ f = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * If optd == false this uses the minimum actuator spacing.
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
   
   ///Calculate the optimum exposure time for a given spatial frequency at a specified actuator spacing
   /** Finds the optimum exposure time at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * 
     * \todo Check inclusion of X in parameters
     * 
     * \returns the optimum expsoure time.
     */
   realT optimumTauWFS( realT m, ///< [in] is the spatial frequency index in u
                        realT n,  ///< [in] is the spatial frequency index in v
                        realT d
                      );
   
   ///Calculate the optimum exposure time for a given spatial frequency at the optimum actuator spacing.
   /** Finds the optimum exposure time at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     * If optd == false this uses the minimum actuator spacing.
     * 
     * \todo Check inclusion of X in parameters
     * 
     * \returns the optimum expsoure time.
     */
   realT optimumTauWFS( realT m, ///< [in] is the spatial frequency index in u
                        realT n  ///< [in] is the spatial frequency index in v
                      );

   ///Calculate the optimum actuator spacing.
   /** Finds the value of m_d_opt where the fitting error is less than than the combined time delay and measurement error.
     *
     * \returns the current value of m_d_opt. 
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
   /** If no changes, merely returns m_wfeVar.  Calls calcStrehl if there are changes.
     *
     * \returns the current value of m_wfeVar;
     */
   realT wfeVar();
   
   ///Get the current Strehl ratio.
   /** If no changes, merely returns m_strehl.  Calls calcStrehl if there are changes.
     * Strehl is calculated using the extended Marechal approximation:
     * 
     \f[
      S = e^{-\sigma_{wfe}^2}
      \f]
     * where \f$ \sigma_{wfe}^2 \f$ is the current value of m_wfeVar.
     *
     * \returns the current value of m_strehl. 
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
     * \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     * \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     * \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     * \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     *  \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     *  \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     *  \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
     * 
     * \tparam imageT is an Eigen-like image type.
     */ 
   template<typename imageT>
   void C4Map( imageT & map /**< [in] the map image to be filled in with contrast */ );


   ///Calculate the residual variance due to to scintilation-amplitude chromaticity.
   /** Used to calculate contrast \ref C5().
     * 
     *  \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     *  \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     * \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
     *  \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
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
}

template<typename realT, class inputSpectT, typename iosT>
aoSystem<realT, inputSpectT, iosT>::~aoSystem()
{
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadGuyon2005()
{
   atm.loadGuyon2005();
   
   F0(1.75e9*0.25*math::pi<realT>()*64.); //Converting to photons/sec
   lam_wfs(0.55e-6);
   lam_sci(1.6e-6);
   D(8.);
   starMag(5);
   
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadMagAOX()
{   
   atm.loadLCO();   
   F0(7.6e10);
   lam_wfs(0.851e-6);
   lam_sci(0.656e-6);
   
   d_min( 6.5/48.0 );
   minTauWFS( (realT) (1./3622.) );
   tauWFS(1./3622.);
   
   D(6.5);
   
   m_specsChanged = true;
   m_dminChanged = true;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadGMagAOX()
{

   loadMagAOX();
   
   
   F0(7.6e10 * 368.0 / (0.25*math::pi<realT>()*6.5*6.5*(1.-0.29*0.29))); //Scaled up.
   
   D(25.4);
   
   m_specsChanged = true;
   m_dminChanged = true;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::F0(realT nF0)
{
   m_F0 = nF0;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::F0()
{
   return m_F0;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::starMag(realT nmag)
{
   m_starMag = nmag;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::starMag()
{
   return m_starMag;
}

template<typename realT, class inputSpectT, typename iosT>
int aoSystem<realT, inputSpectT, iosT>::circularLimit(bool cl)
{
   m_circularLimit = cl;
   m_specsChanged = true;
   m_dminChanged = true;
   
   return 0;
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::circularLimit()
{
   return m_circularLimit;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::Fg(realT mag)
{
   return m_F0*pow(10.0, -0.4*mag);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::Fg()
{
   return Fg(m_starMag);
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::D(realT nD)
{
   m_D = nD;
   psd.D(m_D);
   
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::D()
{
   return m_D;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::d_min(realT nd)
{
   m_d_min = nd;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::d_min()
{
   return m_d_min;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::optd(bool od)
{
   m_optd = od;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::optd()
{
   return m_optd;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::optd_delta(realT odd)
{
   m_optd_delta = odd;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::optd_delta()
{
   return m_optd_delta;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::lam_wfs(realT nlam)
{
   m_lam_wfs = nlam;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::lam_wfs()
{
   return m_lam_wfs;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::npix_wfs(realT npix)
{
   m_npix_wfs.resize(1);
   m_npix_wfs[0] = npix;

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::npix_wfs(const std::vector<realT> & npix)
{
   m_npix_wfs.resize(npix.size());
   for(size_t n=0; n < npix.size(); ++n)
   {
      m_npix_wfs[n] = npix[n];
   }

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::npix_wfs( int idx,
                                                   realT npix
                                                 )
{
   if(m_npix_wfs.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_npix_wfs"));
   }

   m_npix_wfs[idx] = npix;

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::npix_wfs(size_t idx)
{
   if(m_npix_wfs.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_npix_wfs"));
   }

   return m_npix_wfs[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::npix_wfs()
{
   return m_npix_wfs;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ron_wfs(realT nron)
{
   m_ron_wfs.resize(1);
   m_ron_wfs[0] = nron;

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ron_wfs(const std::vector<realT> & nron)
{
   m_ron_wfs.resize(nron.size());
   for(size_t n=0; n < nron.size(); ++n)
   {
      m_ron_wfs[n] = nron[n];
   }

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ron_wfs( int idx,
                                                  realT nron
                                                )
{
   if(m_ron_wfs.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_ron_wfs"));
   }

   m_ron_wfs[idx] = nron;

   m_specsChanged = true;
   m_dminChanged = true;
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ron_wfs(size_t idx)
{
   if(m_ron_wfs.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_ron_wfs"));
   }

   return m_ron_wfs[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::ron_wfs()
{
   return m_ron_wfs;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::Fbg(realT fbg)
{
   m_Fbg.resize(1);
   m_Fbg[0] = fbg;

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::Fbg(const std::vector<realT> & fbg)
{
   m_Fbg.resize(fbg.size());
   for(size_t n=0; n < fbg.size(); ++n)
   {
      m_Fbg[n] = fbg[n];
   }
   
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::Fbg( int idx,
                                              realT fbg
                                             )
{
   if(m_Fbg.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_Fbg"));
   }

   m_Fbg[idx] = fbg;

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::Fbg(size_t idx)
{
   if(m_Fbg.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_Fbg"));
   }

   return m_Fbg[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::Fbg()
{
   return m_Fbg;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::minTauWFS(realT ntau)
{
   m_minTauWFS.resize(1);
   m_minTauWFS[0] = ntau;
 
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::minTauWFS(const std::vector<realT> & ntau)
{
   m_minTauWFS.resize(ntau.size());
   for(size_t n=0; n < ntau.size(); ++n)
   {
      m_minTauWFS[n] = ntau[n];
   }
   
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::minTauWFS( size_t idx,
                                                    realT ntau
                                                  )
{
   if(m_minTauWFS.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_ntau_wfs"));
   }

   m_minTauWFS[idx] = ntau;
 
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::minTauWFS(size_t idx)
{
   if(m_minTauWFS.size() < idx+1)
   {
      throw(mxException("mxlib", MXE_SIZEERR, "A size was calculated incorrectly.", __FILE__, __LINE__, "idx larger than m_ntau_wfs"));
   }

   return m_minTauWFS[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::minTauWFS()
{
   return m_minTauWFS;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::bin_npix(bool bnp)
{
   if(bnp != m_bin_npix)
   {
      m_bin_npix = bnp;
      m_specsChanged = true;
      m_dminChanged = true;
   }
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::bin_npix()
{
   return m_bin_npix;
}


template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::tauWFS(realT ntau)
{
   m_tauWFS = ntau;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::tauWFS()
{
   return m_tauWFS;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::deltaTau(realT ndel)
{
   m_deltaTau = ndel;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::deltaTau()
{
   return m_deltaTau;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::optTau(bool ot)
{
   m_optTau = ot;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::optTau()
{
   return m_optTau;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::lam_sci(realT nlam)
{
   m_lam_sci = nlam;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::lam_sci()
{
   return m_lam_sci;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::zeta(realT nz)
{
   m_zeta = nz;
   m_secZeta = 1/cos(m_zeta);
   
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::zeta()
{
   return m_zeta;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::secZeta()
{
   return m_secZeta;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::fit_mn_max( int mnm )
{
   if(mnm < 0) mnm = 0;
   m_fit_mn_max = mnm;
}

template<typename realT, class inputSpectT, typename iosT>
int aoSystem<realT, inputSpectT, iosT>::fit_mn_max()
{
   return m_fit_mn_max;
}
 
template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::spatialFilter_ku( realT kx )
{
   m_spatialFilter_ku = fabs(kx);
   m_specsChanged = true; //not sure if needed
   m_dminChanged = true; //not sure if needed
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::spatialFilter_ku()
{
   return m_spatialFilter_ku;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::spatialFilter_kv( realT ky )
{
   m_spatialFilter_kv = fabs(ky);
   m_specsChanged = true; //not sure if needed
   m_dminChanged = true; //not sure if needed
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::spatialFilter_kv()
{
   return m_spatialFilter_kv;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ncp_wfe(realT nwfe)
{
   m_ncp_wfe = nwfe;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncp_wfe()
{
   return m_ncp_wfe;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::ncp_alpha(realT alpha)
{
   m_ncp_alpha = alpha;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncp_alpha()
{
   return m_ncp_alpha;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::signal2Noise2( realT & tau_wfs,
                                                         realT d
                                                       )
{      
   realT F = Fg();
               
   double binfact = 1.0;
   int binidx = 0; //index into WFS configurations
   if( m_bin_npix )
   {
      if(m_npix_wfs.size() > 1)
      {
         binidx = d/m_d_min-1; //d_opt should be >= d_min.  This will switch WFS binnings when possible due to number of DOF.
         
         if(binidx < 0) binidx = 0;
         if(binidx >= m_npix_wfs.size()) binidx = m_npix_wfs.size()-1;
         
      }
      else
      {
         int intbin = d/m_d_min;

         binfact = 1./pow((realT) intbin,2);
      }
   }

   return pow(F*tau_wfs,2)/((F+m_npix_wfs[binidx]*binfact*m_Fbg[binidx])*tau_wfs + m_npix_wfs[binidx]*binfact*pow(m_ron_wfs[binidx],2));
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementError( realT m, 
                                                            realT n )
{
   return measurementError(m, n, d_opt());
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementError( realT m, 
                                                            realT n,
                                                            realT d
                                                          )
{
   if(m ==0 and n == 0) return 0;
   
   realT tau_wfs;
   
   if(m_optTau) tau_wfs = optimumTauWFS(m, n, d);
   else tau_wfs = m_tauWFS;
   
   if (m_wfsBeta == 0) wfsBetaUnalloc();
   
   realT beta_p = m_wfsBeta->beta_p(m,n,m_D, d, atm.r_0(m_lam_wfs));
            
   realT snr2 = signal2Noise2( tau_wfs, d );
   
   return pow(beta_p,2)/snr2*pow(m_lam_wfs/m_lam_sci, 2);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementError()
{   
   int mn_max = floor(0.5*m_D/d_opt());
   realT sum = 0;

   for(int m=-mn_max; m <= mn_max; ++m)
   {
      for(int n=-mn_max; n<mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         if( m_circularLimit )
         {
            if( m*m + n*n > mn_max*mn_max ) continue;
         }
         
         sum += measurementError(m,n);
      }
   }

   return sum;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayError( realT m, 
                                                          realT n,
                                                          realT d
                                                        )
{
   if(m ==0 and n == 0) return 0;
   
   realT k = sqrt(m*m + n*n)/m_D;
   
   realT tau_wfs;
   
   if(m_optTau) tau_wfs = optimumTauWFS(m, n, d);
   else tau_wfs = m_tauWFS;
   
   realT tau = tau_wfs + m_deltaTau;
         
   return psd(atm, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * sqrt(atm.X(k, m_lam_sci, m_secZeta)) * pow(math::two_pi<realT>()*atm.v_wind()*k,2) * pow(tau,2);      
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayError( realT m, 
                                                          realT n )
{
   return timeDelayError(m,n,d_opt());
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayError()
{   
   int mn_max = floor(0.5*m_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         if( m_circularLimit )
         {
            if( m*m + n*n > mn_max*mn_max ) continue;
         }
         
         sum += timeDelayError(m,n);
      }
   }

   return sum;
}



template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::fittingError( realT m, 
                                                        realT n )
{
   realT k = sqrt(m*m+n*n)/m_D;
      
   return psd(atm, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::fittingError()
{   
   int mn_max = m_D/(2.0*d_opt());

   realT sum = 0;
   
   for(int m = -m_fit_mn_max; m <= m_fit_mn_max; ++m)
   {
      for(int n = -m_fit_mn_max; n <= m_fit_mn_max; ++n)
      {
         if(m_circularLimit)
         {
            if( m*m + n*n <= mn_max*mn_max) continue;
         }
         else 
         {
            if( abs(m) <= mn_max && abs(n) <= mn_max) continue;
         }
         
         sum += fittingError(m,n);
      }
   }
            
   return sum;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::chromScintOPDError()
{
   int mn_max = floor(0.5*m_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         if( m_circularLimit )
         {
            if( m*m + n*n > mn_max*mn_max ) continue;
         }
         
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
   int mn_max = floor(0.5*m_D/d_opt());
   
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
   int mn_max = floor(0.5*m_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         if( m_circularLimit )
         {
            if( m*m + n*n > mn_max*mn_max ) continue;
         }
         
         sum += C6var(m,n);
      }
   }

   return sum;
   
}  

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::dispAnisoOPDError()
{
   int mn_max = floor(0.5*m_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         if( m_circularLimit )
         {
            if( m*m + n*n > mn_max*mn_max ) continue;
         }
         
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
   int mn_max = floor(0.5*m_D/d_opt());
   
   realT sum = 0;

   for(int m = -mn_max; m <= mn_max; ++m)
   {
      for(int n = -mn_max; n <= mn_max; ++n)
      {
         if(n == 0 && m == 0) continue;
         
         if( m_circularLimit )
         {
            if( m*m + n*n > mn_max*mn_max ) continue;
         }
         
         sum += C8var(m,n);
      }
   }

   return sum;
#endif   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::optimumTauWFS( realT m, 
                                                         realT n,
                                                         realT dact //here called dact due to parameter collision with root-finding
                                                       )
{
   if(m_D == 0)
   {
      mxError("aoSystem::optimumTauWFS", MXE_PARAMNOTSET, "Diameter (D) not set.");
      return -1;
   }
   
   if(m_F0 == 0)
   {
      mxError("aoSystem::optimumTauWFS", MXE_PARAMNOTSET, "0-mag photon flux (F0) not set.");
      return -1;
   }

   double binfact = 1.0;
   int binidx = 0; //index into WFS configurations
   if( m_bin_npix )
   {
      if(m_npix_wfs.size() > 1)
      {
         binidx = m_d_opt/m_d_min - 1; //d_opt should be >= d_min.  This will switch WFS binnings when possible due to number of DOF.
         if(binidx < 0) binidx = 0;
      }
      else
      {
         int intbin = m_d_opt/m_d_min;
         binfact = 1./pow((realT) intbin,2);
      }
   }
   
   realT k = sqrt(m*m + n*n)/m_D;
   
   realT F = Fg();

   if (m_wfsBeta == 0) wfsBetaUnalloc();
      
   realT beta_p = m_wfsBeta->beta_p(m,n,m_D, dact, atm.r_0(m_lam_wfs));

   //Set up for root finding:
   realT a, b, c, d, e;
   
   realT Atmp = 2*pow(atm.lam_0(),2)*psd(atm, k,  m_secZeta)/pow(m_D,2)*(atm.X(k, m_lam_wfs, m_secZeta))*pow(math::two_pi<realT>()*atm.v_wind()*k,2);
   realT Dtmp = pow(m_lam_wfs*beta_p/F,2);
   
   a = Atmp;
   b = Atmp *m_deltaTau;
   c = 0;
   d = -Dtmp * (F+m_npix_wfs[binidx]*binfact*m_Fbg[binidx]);
   e = -Dtmp * 2*(m_npix_wfs[binidx]*binfact*pow(m_ron_wfs[binidx],2));
   
   std::vector<std::complex<realT> > x;
   
   //Get the roots
   math::quarticRoots(x, a, b , c, d, e);
   
   //Now pick the largest positive real root
   realT tauopt = 0.0;
   
   for(int i=0; i<4; i++)
   {
      if( real(x[i]) > 0 && imag(x[i]) == 0 && real(x[i]) > tauopt) tauopt = real(x[i]);
   }
   
   if(tauopt < m_minTauWFS[binidx]) tauopt = m_minTauWFS[binidx];
   
   return tauopt;
   
   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::optimumTauWFS( realT m, 
                                                         realT n )
{
   return optimumTauWFS(m, n, d_opt());
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::d_opt()
{
   if(!m_dminChanged) return m_d_opt;
   
   if( m_d_min == 0 )
   {
      m_d_opt = 1e-50;
      m_dminChanged = false;
      return m_d_opt;
   }
   
   if(!m_optd) 
   {
      m_d_opt = m_d_min;
      m_dminChanged = false;
      
      return m_d_opt;
   }
   
   realT d = m_d_min;
      
   int m = m_D/(2*d);
   int n = 0;

   while( measurementError(m,n, d) + timeDelayError(m,n,d) > fittingError(m, n) && d < m_D/2 ) 
   {
      d += m_d_min/m_optd_delta;
      m = m_D/(2*d);
      n = 0;
   }
   
   m_d_opt = d;
   m_dminChanged = false;
   
   return m_d_opt;
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncpError( int m, 
                                                          int n )
{
   if(m ==0 and n == 0) return 0;
   
   realT k = sqrt(m*m + n*n)/m_D;
   
   return (m_ncp_alpha - 2)/(math::two_pi<realT>()) * pow(m_D, -m_ncp_alpha) * ncpError() * pow(k, -m_ncp_alpha);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::ncpError()
{
   return pow(m_ncp_wfe,2)*pow(math::two_pi<realT>()/m_lam_sci,2);
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::calcStrehl()
{  
   m_wfeMeasurement = measurementError();
   m_wfeTimeDelay = timeDelayError();
   m_wfeFitting = fittingError();
   
   m_wfeChromScintOPD = chromScintOPDError();
   m_wfeChromIndex = chromIndexError();
   m_wfeAnisoOPD = dispAnisoOPDError();
   
   m_wfeNCP = ncpError();
   
   m_wfeVar = m_wfeMeasurement + m_wfeTimeDelay  + m_wfeFitting  + m_wfeChromScintOPD +m_wfeChromIndex + m_wfeAnisoOPD + m_wfeNCP;
   
   m_strehl = exp(-1 * m_wfeVar);
   
   m_specsChanged = false;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::wfeVar()
{
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   
   return m_wfeVar;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::strehl()
{
   if( m_specsChanged || m_dminChanged ) calcStrehl();
   
   return m_strehl;
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
      int mn_max = m_D/(2.*d_opt());
   
      bool outside = false;
      if(m_circularLimit)
      {
         if(m*m + n*n > mn_max*mn_max) outside = true;
      }
      else
      {
         if(abs(m) > mn_max || abs(n) > mn_max) outside = true;
      }
      
      if(mn_max > 0 && outside )
      {
         if(doFittingError == FITTING_ERROR_ZERO) return 0;
         
         realT fe = fittingError(m,n);
         
         realT k = sqrt(m*m + n*n)/m_D;
         
         if(doFittingError == FITTING_ERROR_X) fe *= (atm.X(k, m_lam_sci, m_secZeta));
         else if(doFittingError == FITTING_ERROR_Y) fe *= (atm.Y(k, m_lam_sci, m_secZeta));
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
      
   for(int i=0; i< dim1; ++i)
   {
      int m = i - mc;
      
      for(int j=0; j< dim2; ++j)
      {
         int n = j - nc;
                  
         im(i,j) = (this->*Cfunc)(m, n, true);
      }
   }
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C0var( realT m, 
                                           realT n
                                         )
{
   realT k = sqrt(m*m + n*n)/m_D;
   return psd(atm, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * (atm.X(k, m_lam_sci, m_secZeta));
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
   realT k = sqrt(m*m + n*n)/m_D;
            
   return psd(atm, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * (atm.Y(k, m_lam_sci, m_secZeta));
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
   realT k = sqrt(m*m + n*n)/m_D;
   
   return psd(atm, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * atm.dX(k, m_lam_sci, m_lam_wfs);
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
   realT ni = atm.n_air(m_lam_sci);
   realT nw = atm.n_air(m_lam_wfs);
   
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
   realT k = sqrt(m*m + n*n)/m_D;
     
   return psd(atm, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * atm.X_Z(k, m_lam_wfs, m_lam_sci, m_secZeta);
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
iosT & aoSystem<realT, inputSpectT, iosT>::dumpAOSystem( iosT & ios)
{
   ios << "# AO Params:\n";
   ios << "#    D = " << D() << '\n';
   ios << "#    d_min = " << d_min() << '\n';
   ios << "#    optd = " << std::boolalpha << m_optd << '\n';
   ios << "#    d_opt_delta = " << optd_delta() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    F0 = " << F0() << '\n';
   ios << "#    starMag = " << starMag() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    zeta    = " << zeta() << '\n';   
   ios << "#    lam_wfs = " << lam_wfs() << '\n';
/*   ios << "#    npix_wfs = " << npix_wfs() << '\n';
   ios << "#    ron_wfs = " << ron_wfs() << '\n';
   ios << "#    Fbg = " << Fbg() << '\n';
   ios << "#    minTauWFS = " << minTauWFS() << '\n';
   */ios << "#    bin_npix = " << std::boolalpha << m_bin_npix << '\n';
   ios << "#    tauWFS = " << tauWFS() << '\n';
   ios << "#    optTau = " << std::boolalpha << m_optTau << '\n';
   ios << "#    deltaTau = " << deltaTau() << '\n';
   ios << "#    fit_mn_max = " << m_fit_mn_max << '\n';
   ios << "#    spatialFilter_ku = " << m_spatialFilter_ku << '\n';
   ios << "#    spatialFilter_kv = " << m_spatialFilter_kv << '\n';
   ios << "#    ncp_wfe = " << m_ncp_wfe << '\n';
   ios << "#    ncp_alpha = " << m_ncp_alpha << '\n';
   
   
   if (m_wfsBeta == 0) wfsBetaUnalloc();
   m_wfsBeta->dumpWFS(ios);
   psd.dumpPSD(ios);
   atm.dumpAtmosphere(ios);
   
   ios << "#    Software version: " << '\n';
   ios << "#       mxlib sha1 = " << MXLIB_UNCOMP_CURRENT_SHA1 << '\n';
   ios << "#       mxlib modified = " << MXLIB_UNCOMP_REPO_MODIFIED << '\n';
      
   return ios;
}

extern template
class aoSystem<float, vonKarmanSpectrum<float>, std::ostream>; 

extern template
class aoSystem<double, vonKarmanSpectrum<double>, std::ostream>; 

extern template
class aoSystem<float, vonKarmanSpectrum<float>, std::ostream>; 

extern template
class aoSystem<long double, vonKarmanSpectrum<long double>, std::ostream>; 

#ifdef HASQUAD
extern template
class aoSystem<__float128, vonKarmanSpectrum<__float128>, std::ostream>; 
#endif

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //aoSystem_hpp
