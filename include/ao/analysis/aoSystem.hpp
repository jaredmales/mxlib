/** \file aoSystem.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines an analytical AO system
  * \ingroup mxAO_files
  * 
  */

//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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
   
   realT m_D {0}; ///< Telescope diameter [m]

   std::vector<realT> m_d_min {0}; ///< Minimum AO system actuator pitch [m].  One per WFS mode.

   realT m_d_opt {0}; ///< Current optimum AO system actuator pitch [m]
   bool m_optd {false}; ///< Flag controlling whether actuator pitch is optimized (true) or just uses m_d_min (false).  Default: true.
   realT m_optd_delta {1.0}; ///< The fractional change from d_min used in optimization.  Set to 1 for integer binnings, > 1 for finer sampling.
   
   wfs<realT, iosT> * m_wfsBeta {nullptr}; ///< The WFS beta_p class.
   bool m_ownWfsBeta {false}; ///< Flag indicating if the WFS beta_p pointer is owned by this instance.
   
   realT m_lam_wfs {0}; ///< WFS wavelength [m]

   //WFS detector parameters, which may be a function of binning
   std::vector<realT> m_npix_wfs {0};  ///< Number of WFS pixels.  One per WFS mode.
   std::vector<realT> m_ron_wfs {0};   ///< WFS readout noise [electrons/pix].  One per WFS mode.
   std::vector<realT> m_Fbg {0};       ///< Background flux, [photons/sec/pixel].  One per WFS mode.
   std::vector<realT> m_minTauWFS {0}; ///< Minimum WFS exposure time [sec].  One per WFS mode.
   
   bool m_bin_npix {true}; ///< Flag controlling whether or not to bin WFS pixels according to the actuator spacing.
      
   int m_bin_opt {0}; ///< The optimum binning factor.  If WFS modes are used, this is the mode index (0 to N-1).  If not, it is 1 minus the pixel binning factor. It is always 1 minus the actuator binning factor.
   
   realT m_tauWFS {0}; ///< Actual WFS exposure time [sec]
   
   realT m_deltaTau {0}; ///< Loop latency [sec]
   
   bool m_optTau {true}; ///< Flag controlling whether optimum integration time is calculated (true) enforcing m_minTauWFS, or if m_tauWFS is used (false). Default: true.

   realT m_lam_sci {0}; ///< Science wavelength [m]

   realT m_zeta {0}; ///<  Zenith angle [radians]
   realT m_secZeta {1}; ///< Secant of the Zenith angle (calculated)
   
   int m_fit_mn_max {100}; ///< Maximum spatial frequency index to use for fitting error calculation.

   bool m_circularLimit {false}; ///< Flag to indicate that the spatial frequency limit is circular, not square.

   realT m_spatialFilter_ku {std::numeric_limits<realT>::max()}; ///< The spatial filter cutoff in u, [m^-1]
   realT m_spatialFilter_kv {std::numeric_limits<realT>::max()}; ///< The spatial filter cutoff in v, [m^-1]
   
   realT m_ncp_wfe {0}; ///<Static WFE [m rms]
   realT m_ncp_alpha {2}; ///< Power-law exponent for the NCP aberations.  Default is 2.

   realT m_F0 {0}; ///< 0 mag flux from star at WFS [photons/sec]
      
   realT m_starMag {0}; ///< The magnitude of the star.
   
   bool m_specsChanged {true};///< Flag to indicate that a specification has changed.
   bool m_dminChanged {true};///< Flag to indicate that d_min has changed.

   Eigen::Array<bool,-1,-1> m_controlledModes; ///< Map of which modes are under control.  Set by calcStrehl.

   realT m_wfeMeasurement {0}; ///< Total WFE due to measurement a error [rad^2 at m_lam_sci]
   realT m_wfeTimeDelay {0}; ///< Total WFE due to time delay [rad^2 at m_lam_sci]
   realT m_wfeFitting {0}; ///< Total WFE due to fitting error [rad^2 at m_lam_sci]
   realT m_wfeChromScintOPD {0}; ///< Total WFE due to the chromaticity of scintillation OPD [rad^2 at lam_sc]
   realT m_wfeChromIndex {0}; ///< Total WFE due to the chromaticity of the index of refraction [rad^2 at lam_Sci]
   realT m_wfeAnisoOPD {0}; ///< Total WFE due to dispersive anisoplanatism OPD.
   
   realT m_wfeNCP  {0}; ///< Total WFE due to NCP errors [rad^2 at m_lam_sci]
   
   realT m_wfeVar {0}; ///< The WFE variance, in meters^2.  Never use this directly, instead use wfeVar().
   realT m_strehl {0}; ///< Strehl ratio, a calculated quantity.  Never use this directdy, instead use strehl().
   
   
public:
   
   ///Default c'tor
   /** Calls initialize().
     */
   aoSystem();

   ///Destructor
   ~aoSystem();

   ///Load the default parameters from Guyon, 2005 \cite guyon_2005.
   /**  
     *
     */ 
   void loadGuyon2005(); 
   
   ///Load parameters corresponding to the MagAO-X system.
   void loadMagAOX();
   
   ///Load parameters corresponding to the G-MagAO-X system.
   void loadGMagAOX();
   
   /// Set the value of the primary mirror diameter.
   /** 
     */
   void D( realT nD /**< [in] is the new value of m_D. */);
   
   /// Get the value of the primary mirror diamter
   /**
     * \returns the current value of m_D.
     */ 
   realT D();
   
   /// Set the minimum subaperture sampling for each WFS mode.
   /** This sets the vector
     */
   void d_min( const std::vector<realT> & nd /**< [in] is the new values of m_d_min */);
   
   /// Set the minimum subaperture sampling for a given WFS mode
   /** This sets the entry in the vector, if it has the proper length.
     */
   void d_min( int idx, ///< [in] the index of the WFS mode
               realT nd ///< [in] is the new value of m_d_min[idx] 
             );

   /// Get the minimum subaperture sampling for a given WFS mode
   /**
     * \returns the current value of m_d_min[idx]
     */ 
   realT d_min( size_t idx /**< [in] the index of the WFS mode.*/);

   /// Get the minimum subaperture sampling for all WFS modes.
   /**
     * \returns the current value of m_d_min.
     */ 
   std::vector<realT> d_min();
   
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
      
   /// Set the WFS Beta pointer
   /** The WFS beta parameter determines the photon noise sensitivity.
     *
     * \tparam wfsT is a type derived from ao::analysis::wfs
     */
   template<typename wfsT>
   void wfsBeta( const wfsT & w /**< [in] an object derived from ao::analysis::wfs base class*/);
   
   /// Set the WFS Beta pointer
   /** The WFS beta parameter determines the photon noise sensitivity.
     *
     * \tparam wfsT is a type derived from ao::analysis::wfs
     */
   template<typename wfsT>
   void wfsBeta( const wfsT * w /**< [in] pointer to an object derived from ao::analysis::wfs base class. If nullptr then a new object is allocated and managed.*/);
   
   /// Get the WFS Beta pointer
   /**
     * \returns a pointer to the current m_wfsBeta 
     */ 
   wfs<realT, iosT> * wfsBeta();

   /// Get the value of beta_p for a spatial frequency
   /** beta_p is the photon noise sensitivity of the WFS.
     *
     * \returns beta_p as calculated by the WFS.
     */
   realT beta_p( realT m, ///< [in] the spatial frequency index 
                 realT n  ///< [in] the spatial frequency index
               );
   
   /// Get the value of beta_r for a spatial frequency
   /** beta_r is the read noise sensitivity of the WFS.
     *
     * \returns beta_r as calculated by the WFS.
     */
   realT beta_r( realT m, ///< [in] the spatial frequency index 
                 realT n  ///< [in] the spatial frequency index
               );

   /// Set the value of the WFS wavelength.
   /**
     */
   void lam_wfs( realT nlam /**< [in] is the new value of m_lam_wfs */);
   
   /// Get the value of the WFS wavelength.
   /**
     * \returns the current value of m_lam_wfs.
     */ 
   realT lam_wfs();
   
   /// Set the number of pixels in the WFS for each WFS mode
   /** This sets the vector.
     */
   void npix_wfs( const std::vector<realT> & npix /**< [in] is the new values of m_npix_wfs */);

   /// Set the number of pixels in the WFS for a given WFS mode
   /** This sets the entry in the vector, if it has the proper length.
     */
   void npix_wfs( int idx,   ///< [in] the index of the WFS mode
                  realT npix ///< [in] is the new value of m_npix_wfs[idx] 
                );

   /// Get the number of pixels in the WFS for a given WFS mode
   /**
     * \returns the current value of m_npix_wfs[idx]
     */ 
   realT npix_wfs(size_t idx /**< [in] the index of the WFS mode.*/);
   
   /// Get the number of pixels in the WFS for all WFS modes
   /**
     * \returns the current value of m_npix_wfs
     */ 
   std::vector<realT> npix_wfs();

   /// Set the value of the WFS readout noise for each WFS mode
   /** This sets the vector.
     */
   void ron_wfs( const std::vector<realT> & nron /**< [in] is the new value of m_ron_wfs */);

   /// Set the value of the WFS readout noise for a given bining mode
   /** This sets the entry in the vector if the vector has the proper length
     */
   void ron_wfs( int idx,   ///< [in] the index of the WFS mode
                 realT nron ///< [in] is the new value of m_ron_wfs [idx]
               );

   /// Get the value of the WFS readout noise for a given WFS mode
   /**
     * \returns the current value of m_ron_wfs[idx]
     */ 
   realT ron_wfs( size_t idx /**< [in] the index of the WFS mode.*/);
   
   /// Get the value of the WFS readout noise for all WFS modes
   /**
     * \returns the current value of m_ron_wfs
     */ 
   std::vector<realT> ron_wfs();
   
   /// Set the value of the background fluxes.
   /** This sets the value of the background flux for each WFS mode
     */
   void Fbg(const std::vector<realT> & fbg /**< [in] is the new value of m_Fbg */);

   /// Set a single value of the background flux for a given WFS mode
   /** This sets the entry in the vector if the vector has the proper length
     */
   void Fbg( int idx,  ///< [in] the index of the WFS mode
             realT fbg ///< [in] is the new value of m_Fbg[idx]
            );

   /// Get the value of the background flux for a given WFS mode
   /**
     * \returns m_Fbg[idx]
     */ 
   realT Fbg( size_t idx /**< [in] the index of the WFS mode.*/);
   
   /// Get the value of the background flux for all WFS modes
   /**
     * \returns the current value of m_Fbg
     */ 
   std::vector<realT> Fbg();
   
   /// Set the value of the minimum WFS exposure times.
   /** This sets the vector of the minimum WFS exposure times
     */
   void minTauWFS(const std::vector<realT> &  ntau /**< [in] is the new value of m_minTauWFS */);

   /// Set a single value of the minimum WFS exposure time for a given WFS mode.
   /** This sets the entry in the vector if the vector has sufficient length.
     */
   void minTauWFS(size_t idx, ///< [in] the index of the WFS mode
                  realT ntau  ///< [in] is the new value of m_minTauWFS
                 );

   /// Get the value of the minimum WFS exposure time for a given binning.
   /**
     * \returns the current value of m_minTauWFS[idx].
     */ 
   realT minTauWFS(size_t idx /**< [in] the index of the WFS mode.*/);
   
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
   
   /// Get the value of the optimum binning factor/index.
   /** If WFS modes are used, this is the mode index (0 to N-1).  If not, it is 1 minus the pixel binning factor.
     * It is always 1 minus the actuator binning factor.
     * This calls opt_d() to perform the optimization if needed.
     * 
     * \returns the current value of m_bin_opt.
     */ 
   int bin_opt();

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
   
   /// Set the value of the circularLimit flag
   /** If this is true, then the spatial frequency limits are treated as a circle
     * rather than a square.  This will create a circular dark hole.
     */
   void circularLimit( bool cl /**< [in] the new value of the circularLimit flag*/);
   
   /// Get the value of the circularLimit flag
   /** 
     * \returns the current value of m_circularLimit
     */
   bool circularLimit();

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

   /** \name Measurement Error
     * Calculating the WFE due to WFS measurement noise.
     * @{
     */
      
   /// Calculate the terms of the signal to noise ratio squared (S/N)^2 for the WFS measurement
   /** The S/N squared is 
      \f[
      (S/N)^2 = \frac{ F_\gamma^2 \tau_{wfs}^2 }{ F_\gamma \tau_{wfs} + n_{pix} F_{bg} \tau_{wfs} + n_{pix} \sigma_{ron}^2 }
      \f]
    
     * 
     * \returns the S/N squared
     */
   realT signal2Noise2( realT & Nph,   ///< [out] the number of photons   
                        realT tau_wfs, ///< [in] specifies the WFS exposure time.  If 0, then optimumTauWFS is used
                        realT d,       ///< [in] the actuator spacing in meters, used if binning WFS pixels
                        int b          ///< [in] the binning parameter.  Either the WFS mode index, or the binning factor minus 1.
                      );
   
   ///Calculate the measurement noise at a spatial frequency and specified actuator spacing
   /** Calculates the wavefront phase variance due measurement noise at \f$ k = (m/D)\hat{u} + (n/D)\hat{v} \f$.
     *
     * \returns the measurement noise in rad^2 rms at the science wavelength
     */ 
   realT measurementError( realT m, ///< [in] specifies the u component of the spatial frequency  
                           realT n, ///< [in] specifies the v component of the spatial frequency
                           realT d, ///< [in] the actuator spacing in meters
                           int b    ///< [in] the binning parameter.  Either the WFS mode index, or the binning factor minus 1.
                         );
   
   ///Calculate the total measurement error over all corrected spatial frequencies
   /** This totals the measurement WFE at the optimum actuator spacing and binning.
     *
     * \overload
     *  
     * \returns the total WFE due to measurement error.
     */ 
   realT measurementErrorTotal();
   
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
                         realT d, ///< [in] the actuator spacing, in meters
                         int b    ///< [in] the binning parameter.  Either the WFS mode index, or the binning factor minus 1.
                       );

   /// Calculate the time delay error over all corrected spatial frequencies
   /** This totals the time delay WFE at the optimum actuator spacing and binning.
     * 
     * \overload
     * 
     * \returns the total WFE due to time delay.
     */ 
   realT timeDelayErrorTotal();

   ///@}

   /** \name Fitting Error
     * Calculating the WFE due to uncorrected spatial frequencies.
     * @{
     */
   
   /// Calculate the fitting error at a specific spatial frequency.   
   /**
     * \returns the fitting error in rad^2 rms at the science wavelength at (m,n).
     */ 
   realT fittingError( realT m, ///< [in] specifies the u component of the spatial frequency
                       realT n  ///< [in] specifies the v component of the spatial frequency
                     );

   /// Calculate the total fitting error over all uncorrected spatial frequencies.
   /** This totals the fitting WFE for the optimum actuator spacing and binning.
     *
     * \overload
     * 
     * \returns the total fitting error.
     */ 
   realT fittingErrorTotal();

   ///@}
   
   
   /** \name Chromatic Errors
     * Calculating the WFE due to chromaticity in scintillaion, index of refraction, and dispersion.
     * @{
     */
   
   /// Calculate the wavefront error due to scintillation chromaticity in the OPD at a spatial frequency.
   /** This totals the WFE from scintillation chromaticity in the OPD.
     *
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */
   realT chromScintOPDError( int m, ///< [in] specifies the u component of the spatial frequency
                             int n  ///< [in] specifies the v component of the spatial frequency
                           );

   /// Calculate the wavefront error due to scintillation chromaticity in the OPD over all spatial frequencies.   
   /** This totals the WFE from scintillation chromaticity in the OPD for the optimum actuator spacing.
     *
     * \overload
     *  
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */    
   realT chromScintOPDErrorTotal();

   /// Calculate the wavefront error due to scintillation chromaticity in amplitude for a given spatial frequency.
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */
   realT chromScintAmpError( int m, ///< [in] specifies the u component of the spatial frequency
                             int n  ///< [in] specifies the v component of the spatial frequency
                           );

   /// Calculate the wavefront error due to scintillation chromaticity in amplitude over all spatial frequencies.   
   /**
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */       
   realT chromScintAmpError();

   /// Calculate the wavefront error due to chromaticity in the index of refraction at a given spatial frequency.
   /** This totals the WFE from chromaticity in the index of refraction.
     *
     * \overload
     *
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */
   realT chromIndexError( int m, ///< [in] specifies the u component of the spatial frequency
                          int n  ///< [in] specifies the v component of the spatial frequency
                        );

   /// Calculate the wavefront error due to chromaticity in the index of refraction at a specific spatial frequency.
   /** This totals the WFE from chromaticity in the index of refraction for the optimum actuator spacing.
     * 
     * \overload
     * 
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */       
   realT chromIndexErrorTotal();

   /// Calculate the wavefront error due to dispersive anisoplanatism in the OPD at a given spatial frequency
   /** This calculates the WFE from dispersive anisoplanatism in the OPD.
     *
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */
   realT dispAnisoOPDError( int m, ///< [in] specifies the u component of the spatial frequency
                            int n  ///< [in] specifies the v component of the spatial frequency
                          );

   /// Calculate the wavefront error due to dispersive anisoplanatism in the OPD over all specific spatial frequencies.   
   /** This totals the WFE from dispersive anisoplanatism in the OPD for the optimum actuator spacing.
     *
     * \overload
     *  
     * \returns the WFE in rad^2 rms at the science wavelength at (m,n).
     */       
   realT dispAnisoOPDErrorTotal();

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
                        realT n, ///< [in] is the spatial frequency index in v
                        realT d, ///< [in] the actuator spacing, in meters
                        int b    ///< [in] the binning parameter.  Either the WFS mode index, or the binning factor minus 1.
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
   
   ///Get the Strehl ratio for a given actuator pitch and WFS WFS mode.
   /** 
     * Strehl is calculated using the extended Marechal approximation:
     * 
     \f[
      S = e^{-\sigma_{wfe}^2}
      \f]
     * where \f$ \sigma_{wfe}^2 \f$ is the WFE for the spacing and binning
     *
     * \returns the strehl for these values
     */
   realT strehl( realT d, ///< [in] the actuator spacing, in meters
                 int b    ///< [in] the binning parameter.  Either the WFS mode index, or the binning factor minus 1.
               );

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
   
   /// Worker function for raw contrast fuctions
   /** The logic in this calculation is the same regardless of component, except for the calculation of variance.
     * The varFunc function pointer is used to calculate the variance, otherwise  This handles the other details such
     * as Strehl normalization, fitting error (if appropriate), and bounds checking.
     *
     * \note this gives the raw PSD, it must be convolved with the PSF for the true contrast.
     * 
     * \tparam varFuncT  is a function-pointer-to-member (of this class) with signature realT(*f)(realT, realT)
     */
   template<typename varFuncT>
   realT C_( int m,             ///< [in] is the spatial frequency index in u/
             int n,             ///< [in] is the spatial frequency index in v.
             bool normStrehl,   ///< [in] flag controls whether the contrast is normalized by Strehl ratio/
             varFuncT varFunc,  ///< [in] the variance function to use.
             int doFittingError ///< [in] flag to describe how fitting error should be considered for this term: FITTING_ERROR_NO, FITTING_ERROR_ZERO, FITTING_ERROR_X, or FITTING_ERROR_Y.
           );

   ///Worker function for the contrast-map functions.
   /** The map calculation is the same for all terms, except for the calculation of variance at each pixel.
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to uncorrected phase, C0.
   /** Contrast C0 is the uncorrected phase, with the effects of scintillation included.  See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \note this is the raw PSD, it must be convolved with the PSF for the true contrast.
     * 
     * \returns C0.
     */
   realT C0( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in] is the spatial frequency index in v
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to uncorrected amplitude, C1.
   /** Contrast C1 is the uncorrected amplitude due to scintillation.  See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \returns C0.
     */
   realT C1( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in] is the spatial frequency index in v
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to measurement and time delay errors in phase/OPD at a spatial frequency.
   /** Contrast C2 is just the total variance due to time delay and measurement errors, 
     * divided by the Strehl ratio.   See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \returns C2.
     */
   realT C2( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in]  is the spatial frequency index in v
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to measurement and time delay errors in amplitude at a spatial frequency.
   /** Contrast C3 is just the total variance due to time delay and measurement errors, 
     * divided by the Strehl ratio.    See Guyon (2005) \cite guyon_2005, and the updated
     * derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \returns C3.
     */
   realT C3( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in]  is the spatial frequency index in v
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to scintilation-OPD chromaticity.
   /** Contrast C4 is due to the chromaticity of scintillation, causing the ODP measurement to be slightly incorrect at the science wavelength.
     * See Guyon (2005) \cite guyon_2005, and the updated derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \returns C4.
     */
   realT C4( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in]  is the spatial frequency index in v
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to scintilation-amplitude chromaticity.
   /** Contrast C5 is due to the chromaticity of scintillation, causing the amplitude measurement to be slightly incorrect at
     * the science wavelength. See Guyon (2005) \cite guyon_2005, and the updated derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \returns C4.
     */
   realT C5( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in]  is the spatial frequency index in v
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
                realT n   ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to chromaticity of the index of refraction of air.
   /** Contrast C6 is due to the index of refraction of air being wavelength dependent,  causing the ODP measurement to be slightly incorrect
     * at the science wavelength.   See Guyon (2005) \cite guyon_2005, and the updated derivation in Males \& Guyon (2018) \cite males_guyon_2018. 
     * 
     * \returns C6.
     */
   realT C6( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in] is the spatial frequency index in v
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
                realT n  ///< [in] is the spatial frequency index in v
              );
   
   ///Calculate the contrast due to dispersive anisoplanatism.
   /** Contrast C7 is due to atmospheric dispersion causing light at different wavelengths to take different paths through atmospheric turbulence.
     * See Fitzgerald (2017, in prep) \cite fitzgerald_2017
     * 
     * \returns C7.
     */
   realT C7( realT m,               ///< [in] is the spatial frequency index in u
             realT n,               ///< [in] is the spatial frequency index in v
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
   

    /// Setup the configurator to configure this class
   /**
     * todo: "\test Loading aoAtmosphere config settings \ref tests_ao_analysis_aoAtmosphere_config "[test doc]" 
     */
   void setupConfig( app::appConfigurator & config /**< [in] the app::configurator object*/);

   /// Load the configuration of this class from a configurator
   /**
     * \todo: "\test Loading aoAtmosphere config settings \ref tests_ao_analysis_aoAtmosphere_config "[test doc]""
     */
   void loadConfig( app::appConfigurator & config /**< [in] the app::configurator object*/);
};



template<typename realT, class inputSpectT, typename iosT>
aoSystem<realT, inputSpectT, iosT>::aoSystem()
{
   wfsBeta<wfs<realT,iosT>>(nullptr); //allocate a default ideal WFS.
}

template<typename realT, class inputSpectT, typename iosT>
aoSystem<realT, inputSpectT, iosT>::~aoSystem()
{
   if(m_wfsBeta && m_ownWfsBeta)
   {
      delete m_wfsBeta;
   }
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadGuyon2005()
{
   atm.loadGuyon2005();
   
   F0(1.75e9*0.25*math::pi<realT>()*64.*0.18); //Converting to photons/sec
   lam_wfs(0.55e-6);
   lam_sci(1.6e-6);
   D(8.);
   starMag(5);
   
   //The rest of Guyon 05 is very ideal
   npix_wfs(std::vector<realT>({12868}));
   ron_wfs(std::vector<realT>({0.0}));
   Fbg(std::vector<realT>({0.0}));
   
   d_min( 8.0/1e3); //Allow super fine sampling
   minTauWFS( (realT) (1./1e9) ); // Just be super fast.
   tauWFS(1./1e9); //Just be super fast
   
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadMagAOX()
{   
   atm.loadLCO();   
   F0(4.2e10);
   lam_wfs(0.791e-6);
   lam_sci(0.656e-6);
   
   npix_wfs(std::vector<realT>({9024}));
   ron_wfs(std::vector<realT>({0.57}));
   Fbg(std::vector<realT>({0.22}));
   
   d_min(std::vector<realT>({6.5/48.0}));
   minTauWFS(std::vector<realT>({1./3622.}));
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
void aoSystem<realT, inputSpectT, iosT>::d_min(const std::vector<realT> & nd)
{
   m_d_min.assign(nd.begin(), nd.end());
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::d_min( int idx,
                                                realT nd
                                              )
{
   if(m_d_min.size() < idx+1)
   {
      mxThrowException(err::sizeerr, "aoSystem::d_min", "idx larger than m_d_min");
   }

   m_d_min[idx] = nd;

   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::d_min(size_t idx)
{
   if(m_d_min.size() < idx+1)
   {
      mxThrowException(err::sizeerr, "aoSystem::d_min", "idx larger than m_d_min");
   }

   return m_d_min[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::d_min()
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
template<typename wfsT>
void aoSystem<realT, inputSpectT, iosT>::wfsBeta( const wfsT & w )
{
   if(m_wfsBeta && m_ownWfsBeta)
   {
      delete m_wfsBeta;
      m_wfsBeta = nullptr;
   }

   m_wfsBeta = (wfs<realT,iosT> *) &w;
   m_ownWfsBeta = false;
}

template<typename realT, class inputSpectT, typename iosT>
template<typename wfsT>
void aoSystem<realT, inputSpectT, iosT>::wfsBeta( const wfsT * w )
{
   if(m_wfsBeta && m_ownWfsBeta)
   {
      delete m_wfsBeta;
      m_wfsBeta = nullptr;
   }

   if(w)
   {
      m_wfsBeta = (wfs<realT,iosT> *) w;
      m_ownWfsBeta = false;
   }
   else
   {
      m_wfsBeta = new wfsT;
      m_ownWfsBeta = true;
   }
}

template<typename realT, class inputSpectT, typename iosT>
wfs<realT, iosT> * aoSystem<realT, inputSpectT, iosT>::wfsBeta()
{
   return m_wfsBeta;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::beta_p( realT m, realT n)
{
   if( m_wfsBeta == 0) mxThrowException(err::paramnotset, "aoSystem::beta_p", "The WFS is not assigned."); 
   
   return m_wfsBeta->beta_p(m, n, m_D, d_opt(), atm.r_0(m_lam_sci) );
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::beta_r( realT m, realT n)
{
   if( m_wfsBeta == 0) mxThrowException(err::paramnotset, "aoSystem::beta_r", "The WFS is not assigned."); 
   
   return m_wfsBeta->beta_r(m, n, m_D, d_opt(), atm.r_0(m_lam_sci) );
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
      mxThrowException(err::sizeerr, "aoSystem::npix_wfs", "idx larger than m_npix_wfs");
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
      mxThrowException(err::sizeerr, "aoSystem::npix_wfs", "idx larger than m_npix_wfs");
   }

   return m_npix_wfs[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::npix_wfs()
{
   return m_npix_wfs;
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
      mxThrowException(err::sizeerr, "aoSystem::ron_wfs", "idx larger than m_ron_wfs");
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
      mxThrowException(err::sizeerr, "aoSystem::ron_wfs", "idx larger than m_ron_wfs");
   }

   return m_ron_wfs[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::ron_wfs()
{
   return m_ron_wfs;
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
      mxThrowException(err::sizeerr, "aoSyste,::Fbg", "idx larger than m_Fbg");
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
      mxThrowException(err::sizeerr, "aoSyste,::Fbg", "idx larger than m_Fbg");
   }

   return m_Fbg[idx];
}

template<typename realT, class inputSpectT, typename iosT>
std::vector<realT> aoSystem<realT, inputSpectT, iosT>::Fbg()
{
   return m_Fbg;
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
      mxThrowException(err::sizeerr, "aoSystem::minTauWFS", "idx larger than m_ntau_wfs");
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
      mxThrowException(err::sizeerr, "aoSystem::minTauWFS", "idx larger than m_ntau_wfs");
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
int aoSystem<realT, inputSpectT, iosT>::bin_opt()
{
   d_opt();
   return m_bin_opt;
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
void aoSystem<realT, inputSpectT, iosT>::circularLimit(bool cl)
{
   m_circularLimit = cl;
   m_specsChanged = true;
   m_dminChanged = true;
}

template<typename realT, class inputSpectT, typename iosT>
bool aoSystem<realT, inputSpectT, iosT>::circularLimit()
{
   return m_circularLimit;
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
realT aoSystem<realT, inputSpectT, iosT>::signal2Noise2( realT & Nph,
                                                         realT tau_wfs,
                                                         realT d,
                                                         int b
                                                       )
{      
   Nph = Fg()*tau_wfs;
   
   double binfact = 1.0;
   int binidx = 0;
   if( m_bin_npix )
   {
      if(m_npix_wfs.size() == 1) //Only if "true binning", otherwise the WFS mode vectors handle it.
      {
         binfact = 1./pow((realT) b+1,2);
      }
      else binidx = b;
   }

   return pow(Nph,2)/(  m_npix_wfs[binidx]*binfact * ( m_Fbg[binidx]*tau_wfs + pow(m_ron_wfs[binidx],2)));
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementError( realT m, 
                                                            realT n,
                                                            realT d,
                                                            int b
                                                          )
{
   if(m ==0 and n == 0) return 0;
   
   realT tau_wfs;
   
   if(m_optTau) tau_wfs = optimumTauWFS(m, n, d, b);
   else tau_wfs = m_tauWFS;

   if (m_wfsBeta == 0) mxThrowException(err::paramnotset, "aoSystem::measurementError", "The WFS is not assigned."); 
   
   realT beta_p = m_wfsBeta->beta_p(m,n,m_D, d, atm.r_0(m_lam_wfs));
   realT beta_r = m_wfsBeta->beta_r(m,n,m_D, d, atm.r_0(m_lam_wfs));
   realT Nph = 0;
   realT snr2 = signal2Noise2( Nph, tau_wfs, d, b );
   
   return (pow(beta_r,2)/snr2  + pow(beta_p,2)/Nph)  *pow(m_lam_wfs/m_lam_sci, 2);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::measurementErrorTotal()
{   
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   return m_wfeMeasurement;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayError( realT m, 
                                                          realT n,
                                                          realT d,
                                                          int b
                                                        )
{
   if(m ==0 and n == 0) return 0;
   
   realT k = sqrt(m*m + n*n)/m_D;
   
   realT tau_wfs;

   if(m_optTau) tau_wfs = optimumTauWFS(m, n, d, b);
   else tau_wfs = m_tauWFS;
   
   realT tau = tau_wfs + m_deltaTau;
         
   ///\todo handle multiple layers
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * sqrt(atm.X(k, m_lam_sci, m_secZeta)) * pow(math::two_pi<realT>()*atm.v_wind()*k,2) * pow(tau,2);      
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::timeDelayErrorTotal()
{
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   return m_wfeTimeDelay;
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::fittingError( realT m, 
                                                        realT n )
{
   realT k = sqrt(m*m+n*n)/m_D;
      
   ///\todo handle multiple layers
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::fittingErrorTotal()
{   
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   return m_wfeFitting;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::chromScintOPDError( int m,
                                                              int n
                                                            )
{
   return C4var(m,n);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::chromScintOPDErrorTotal()
{
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   return m_wfeChromScintOPD;
}

template<typename realT, class inputSpectT, typename iosT>   
realT aoSystem<realT, inputSpectT, iosT>::chromScintAmpError( int m,
                                                              int n
                                                            )
{
   return C5var(m,n);

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
realT aoSystem<realT, inputSpectT, iosT>::chromIndexError( int m,
                                                           int n
                                                         )
{
   return C6var(m,n);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::chromIndexErrorTotal()
{
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   return m_wfeChromIndex;
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::dispAnisoOPDError( int m,
                                                             int n
                                                           )
{
   return C7var(m,n);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::dispAnisoOPDErrorTotal()
{
   if(m_specsChanged || m_dminChanged ) calcStrehl();
   return m_wfeAnisoOPD;
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
                                                         realT dact, //here called dact due to parameter collision with root-finding
                                                         int bbin
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
         binidx = bbin;
      }
      else
      {
         binfact = 1.0/pow((realT) (bbin+1),2);
      }
   }
   
   realT k = sqrt(m*m + n*n)/m_D;
   
   realT F = Fg();

   if (m_wfsBeta == 0) mxThrowException(err::paramnotset, "aoSystem::beta_p", "The WFS is not assigned."); 
      
   realT beta_p = m_wfsBeta->beta_p(m,n,m_D, dact, atm.r_0(m_lam_wfs));
   realT beta_r = m_wfsBeta->beta_r(m,n,m_D, dact, atm.r_0(m_lam_wfs));

   //Set up for root finding:
   realT a, b, c, d, e;
   
   ///\todo handle multiple layers
   realT Atmp = 2*pow(atm.lam_0(),2)*psd(atm, 0, k,  m_secZeta)/pow(m_D,2)*(atm.X(k, m_lam_wfs, m_secZeta))*pow(math::two_pi<realT>()*atm.v_wind()*k,2);
   
   a = Atmp;
   b = Atmp *m_deltaTau;
   c = 0;
   d = -1*(pow(m_lam_wfs,2) / F) * ( (m_npix_wfs[binidx]*binfact*m_Fbg[binidx] / F)*pow(beta_r,2) + pow(beta_p,2)); 
   e = -2*pow(m_lam_wfs,2) * (m_npix_wfs[binidx]*binfact)*pow(m_ron_wfs[binidx],2) * pow(beta_r,2) / pow(F,2);
   
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
   d_opt(); //also gets m_bin_opt;
   return optimumTauWFS(m, n, m_d_opt, m_bin_opt);
}


template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::d_opt()
{
   if(!m_dminChanged) return m_d_opt;
   
   //Validate the modes
   ///\todo this should be in a separate valideModes function
   if(m_npix_wfs.size() < 1)
   {
      mxThrowException(err::sizeerr, "aoSystem::d_opt", "npix_wfs must have at least one entry");
   }

   if(m_d_min.size() != m_npix_wfs.size())
   {
      mxThrowException(err::sizeerr, "aoSystem::d_opt", "d_min must be the same size as npix_wfs");
   }

   if(m_ron_wfs.size() != m_npix_wfs.size())
   {
      mxThrowException(err::sizeerr, "aoSystem::d_opt", "ron_wfs must be the same size as npix_wfs");
   }

   if(m_Fbg.size() != m_npix_wfs.size())
   {
      mxThrowException(err::sizeerr, "aoSystem::d_opt", "F_bg must be the same size as npix_wfs");
   }

   if(m_minTauWFS.size() != m_npix_wfs.size())
   {
      mxThrowException(err::sizeerr, "aoSystem::d_opt", "minTauWFS must be the same size as npix_wfs");
   }

   ///\todo investigate why we tolerate m_d_min = 0 here.  Is this just a config error?
   if( m_d_min.size() == 1 && m_d_min[0] == 0 )
   {
      m_d_opt = 1e-50;
      m_bin_opt = 0;
      m_dminChanged = false;
      return m_d_opt;
   }
   
   if(!m_optd) 
   {
      m_d_opt = m_d_min[0];
      m_bin_opt = 0;
      m_dminChanged = false;
      
      return m_d_opt;
   }
   
   
   realT best_d = 0;

   if(m_bin_npix)
   {
      if(m_npix_wfs.size() > 1) //Optimize over the WFS modes
      {
         int best_idx = 0;
         best_d = m_d_min[0];

         realT bestStrehlOverall = 0;
         
         for(int b=0; b < m_npix_wfs.size(); ++b)
         {
            realT d = m_d_min[b];
            realT s = strehl(d,b); 
            realT bestStrehl = s;

            //Find the actuator pitch at which AO error is less than fitting error at Nyquist
            while( d < m_D/2 ) 
            {
               d += m_d_min[b]/m_optd_delta;
               s = strehl(d,b);
             
               if(s < bestStrehl) //Strehl got worse.  Can't be <= otherwise small m_optd_deltas will immediately out
               {
                  d -= m_d_min[b]/m_optd_delta; //go back to previous d
                  break;
               }
               else bestStrehl = s;
            }

            //Check if this is optimum so far
            if(bestStrehl >= bestStrehlOverall) // >= to favor fewer actuators
            {
               bestStrehlOverall = bestStrehl;
               best_idx = b;
               best_d = d;

            }
         }

         m_bin_opt = best_idx;
      }
      else
      {
         realT d = m_d_min[0];
         m_bin_opt = 0; 

         realT s = strehl(d,m_bin_opt); 
         realT bestStrehl = s;

         while( d < m_D/2 ) 
         {
            d += m_d_min[0]/m_optd_delta;
            m_bin_opt = d/m_d_min[0] - 1;
            if(m_bin_opt < 0) m_bin_opt = 0;

            s = strehl(d,m_bin_opt);

            if(s < bestStrehl) //Strehl got worse.  Can't be <= otherwise small m_optd_deltas will immediately out
            {
               d -= m_d_min[0]/m_optd_delta; //go back to previous d
               m_bin_opt = d/m_d_min[0] - 1;
               if(m_bin_opt < 0) m_bin_opt = 0;
               break;
            }
            else bestStrehl = s;
         }
         best_d = d;
      }
   }
   else
   {
      realT d = m_d_min[0];
      m_bin_opt = 0;

      realT s = strehl(d,m_bin_opt); 
      realT bestStrehl = s;

      //Find the actuator pitch which minimizes total error (AO error is less than fitting error at Nyquist)
      while( d < m_D/2 ) 
      {
         d += m_d_min[0]/m_optd_delta;
         
         s = strehl(d,m_bin_opt);
         if(s < bestStrehl) //Strehl got worse.  Can't be <= otherwise small m_optd_deltas will immediately out
         {
            d -= m_d_min[0]/m_optd_delta; //go back to previous d
            break;
         }
         else bestStrehl = s;
      }
      best_d = d;
   }
   
   m_d_opt = best_d;
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
realT aoSystem<realT, inputSpectT, iosT>::strehl( realT d,
                                                  int b
                                                )
{  

   int mn_max = floor(0.5*m_D/d);

   realT wfeVar = 0;

   bool controlled;
   for(int m = -m_fit_mn_max; m <= m_fit_mn_max; ++m)
   {
      for(int n = -m_fit_mn_max; n <= m_fit_mn_max; ++n)
      {
         controlled = false;
         if(m_circularLimit)
         {
            if( m*m + n*n <= mn_max*mn_max) controlled = true;
         }
         else
         {
            if( abs(m) <= mn_max && abs(n) <= mn_max) controlled = true;
         }

         if(controlled)
         {
            realT wfeMeasurement = measurementError(m,n, d, b);
            realT wfeTimeDelay = timeDelayError(m,n,d, b);
            realT wfeChromScintOPD = chromScintOPDError(m,n);
            realT wfeChromIndex = chromIndexError(m,n);
            realT wfeAnisoOPD = dispAnisoOPDError(m,n);

            realT wfeFitting = fittingError(m,n);

            if(wfeFitting < wfeMeasurement + wfeTimeDelay + wfeChromScintOPD + wfeChromIndex + wfeAnisoOPD)
            {
               wfeVar += wfeFitting;
            }
            else
            {
               wfeVar +=  wfeMeasurement + wfeTimeDelay + wfeChromScintOPD + wfeChromIndex + wfeAnisoOPD;
            }
         }
         else
         {
            wfeVar += fittingError(m,n);
         }

         realT wfeNCP = ncpError(m, n);
         wfeVar += wfeNCP;
      }
   }

   return exp(-1 * wfeVar);
   
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::calcStrehl()
{  
   realT d = d_opt(); //have to run this to get m_bin_opt set.
   int b = m_bin_opt;

   int mn_max = floor(0.5*m_D/d);

   m_wfeMeasurement = 0;
   m_wfeTimeDelay = 0;
   m_wfeChromScintOPD = 0;
   m_wfeChromIndex = 0;
   m_wfeAnisoOPD = 0;

   m_wfeFitting = 0;

   m_wfeNCP = 0;


   m_controlledModes.resize(2*m_fit_mn_max+1, 2*m_fit_mn_max+1);
   m_controlledModes.setZero();

   bool controlled;
   for(int m = -m_fit_mn_max; m <= m_fit_mn_max; ++m)
   {
      for(int n = -m_fit_mn_max; n <= m_fit_mn_max; ++n)
      {
         controlled = false;
         if(m_circularLimit)
         {
            if( m*m + n*n <= mn_max*mn_max) controlled = true;
         }
         else
         {
            if( abs(m) <= mn_max && abs(n) <= mn_max) controlled = true;
         }

         if(controlled)
         {
            realT wfeMeasurement = measurementError(m,n, d, b);
            realT wfeTimeDelay = timeDelayError(m,n,d, b);
            realT wfeChromScintOPD = chromScintOPDError(m,n);
            realT wfeChromIndex = chromIndexError(m,n);
            realT wfeAnisoOPD = dispAnisoOPDError(m,n);

            realT wfeFitting = fittingError(m,n);

            if(wfeFitting < wfeMeasurement + wfeTimeDelay + wfeChromScintOPD + wfeChromIndex + wfeAnisoOPD)
            {
               m_wfeFitting += wfeFitting;
            }
            else //it's worth controlling this mode.
            {
               m_wfeMeasurement += wfeMeasurement;
               m_wfeTimeDelay += wfeTimeDelay;
               m_wfeChromScintOPD += wfeChromScintOPD;
               m_wfeChromIndex += wfeChromIndex;
               m_wfeAnisoOPD += wfeAnisoOPD;
               m_controlledModes(m_fit_mn_max+m, m_fit_mn_max+n) = 1;
            }
         }
         else
         {
            m_wfeFitting += fittingError(m,n);
         }

         m_wfeNCP += ncpError(m, n);
      }
   }
   
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
realT aoSystem<realT, inputSpectT, iosT>::C_( int m,
                                              int n,
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
   
      /*bool outside = false;
      if(m_circularLimit)
      {
         if(m*m + n*n > mn_max*mn_max) outside = true;
      }
      else
      {
         if(abs(m) > mn_max || abs(n) > mn_max) outside = true;
      }*/
//      bool outside =;
      if(  m_controlledModes(m_fit_mn_max+m, m_fit_mn_max+n)  == false) //  mn_max > 0 && outside )
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
   //get if not doing fitting error or if inside control region:

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
   ///\todo handle multiple layers
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * (atm.X(k, m_lam_sci, m_secZeta));
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
            
   ///\todo handle multiple layers
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * (atm.Y(k, m_lam_sci, m_secZeta));
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
   d_opt();
   return measurementError(m, n, d_opt(), m_bin_opt) + timeDelayError(m,n, d_opt(), m_bin_opt);
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
realT aoSystem<realT, inputSpectT, iosT>::C3var( realT m, 
                                                 realT n
                                               )
{
   return 0;//measurementError(m, n) + timeDelayError(m,n);
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C3( realT m, 
                                              realT n,
                                              bool normStrehl
                                            )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C3var, FITTING_ERROR_ZERO);//FITTING_ERROR_Y);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C3Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C3);   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C4var( realT m, 
                                           realT n
                                         )
{
   realT k = sqrt(m*m + n*n)/m_D;
   
   ///\todo handle multiple layers
   //This does not need to be divided by X b/c we haven't multiplied by it, this is is C0/X.
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * atm.dX(k, m_lam_sci, m_lam_wfs);
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
realT aoSystem<realT, inputSpectT, iosT>::C5var( realT m, 
                                           realT n
                                         )
{
   realT k = sqrt(m*m + n*n)/m_D;

   ///\todo handle multiple layers
   //This does not need to be divided by Y b/c we haven't multiplied by it, this is is C1/Y.
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * atm.dY(k, m_lam_sci, m_lam_wfs);
}

                                                  
template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C5( realT m, 
                                              realT n,
                                              bool normStrehl
                                            )
{
   return C_(m,n,normStrehl,&aoSystem<realT, inputSpectT, iosT>::C5var, FITTING_ERROR_ZERO);
}

template<typename realT, class inputSpectT, typename iosT>
template<typename imageT>
void aoSystem<realT, inputSpectT, iosT>::C5Map( imageT & im )
{
   C_Map(im, &aoSystem<realT, inputSpectT, iosT>::C5);   
}

template<typename realT, class inputSpectT, typename iosT>
realT aoSystem<realT, inputSpectT, iosT>::C6var( realT m, 
                                           realT n
                                         )
{
   realT ni = atm.n_air(m_lam_sci);
   realT nw = atm.n_air(m_lam_wfs);
   
   return C0var(m, n) * pow( (ni-nw)/(ni-1), 2);  //Fixes error in Eqn 29 
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
     
   ///\todo this needs to handle multiple layers -- do we violate an assumption of L_0 isn't constant?
   return psd(atm, 0, k, m_secZeta)/pow(m_D,2) * pow(atm.lam_0()/m_lam_sci, 2) * atm.X_Z(k, m_lam_wfs, m_lam_sci, m_secZeta);
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
   
   if(m_d_min.size() > 0)
   {
      ios << "#    d_min = " << d_min((size_t) 0);
      for(size_t n=1; n < m_d_min.size(); ++n) ios << ',' << d_min(n);
      ios << '\n';
   }
   else ios << "#    d_min = null\n";

   ios << "#    optd = " << std::boolalpha << m_optd << '\n';
   ios << "#    d_opt_delta = " << optd_delta() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    F0 = " << F0() << '\n';
   ios << "#    starMag = " << starMag() << '\n';
   ios << "#    lam_sci = " << lam_sci() << '\n';
   ios << "#    zeta    = " << zeta() << '\n';   
   ios << "#    lam_wfs = " << lam_wfs() << '\n';
   
   if(npix_wfs().size() > 0)
   {
      ios << "#    npix_wfs = " << npix_wfs((size_t) 0);
      for(size_t n=1; n < npix_wfs().size(); ++n) ios << ',' << npix_wfs(n);
      ios << '\n';
   }
   else ios << "#    npix_wfs = null\n";

   if(ron_wfs().size() > 0)
   {
      ios << "#    ron_wfs = " << ron_wfs((size_t) 0);
      for(size_t n=1; n < ron_wfs().size(); ++n) ios << ',' << ron_wfs(n);
      ios << '\n';
   }
   else ios << "#    ron_wfs = null\n";

   if(Fbg().size() > 0)
   {
      ios << "#    Fbg = " << Fbg((size_t) 0);
      for(size_t n=1; n < Fbg().size(); ++n) ios << ',' << Fbg(n);
      ios << '\n';
   }
   else ios << "#    Fbg = null\n";

   if(minTauWFS().size() > 0)
   {
      ios << "#    minTauWFS = " << minTauWFS((size_t) 0);
      for(size_t n=1; n < minTauWFS().size(); ++n) ios << ',' << minTauWFS(n);
      ios << '\n';
   }
   else ios << "#    minTauWFS = null\n";

   ios << "#    bin_npix = " << std::boolalpha << m_bin_npix << '\n';
   ios << "#    tauWFS = " << tauWFS() << '\n';
   ios << "#    optTau = " << std::boolalpha << m_optTau << '\n';
   ios << "#    deltaTau = " << deltaTau() << '\n';
   ios << "#    fit_mn_max = " << m_fit_mn_max << '\n';
   ios << "#    spatialFilter_ku = " << m_spatialFilter_ku << '\n';
   ios << "#    spatialFilter_kv = " << m_spatialFilter_kv << '\n';
   ios << "#    ncp_wfe = " << m_ncp_wfe << '\n';
   ios << "#    ncp_alpha = " << m_ncp_alpha << '\n';
   
   
   if (m_wfsBeta == 0) mxThrowException(err::paramnotset, "aoSystem::beta_p", "The WFS is not assigned."); 
   m_wfsBeta->dumpWFS(ios);
   psd.dumpPSD(ios);
   atm.dumpAtmosphere(ios);
   
   dumpGitStatus(ios);

   return ios;
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::setupConfig( app::appConfigurator & config )
{
   using namespace mx::app;

   //AO System configuration
   config.add("aosys.wfs"               ,"", "aosys.wfs"             , argType::Required, "aosys", "wfs",              false, "string", "The WFS type: idealWFS, unmodPyWFS, asympModPyWFS, shwfs, calculatedWFS");
   config.add("aosys.wfs_beta_p"        ,"", "aosys.wfs_beta_p"      , argType::Required, "aosys", "wfs_beta_p",       false, "string", "The beta_p file path for calcualtedWFS");
   config.add("aosys.wfs_beta_r"        ,"", "aosys.wfs_beta_r"      , argType::Required, "aosys", "wfs_beta_r",       false, "string", "The beta_r file path for calcualtedWFS");
   config.add("aosys.wfs_sensitivity"   ,"", "aosys.wfs_sensitivity" , argType::Required,     "aosys", "wfs_sensitivity",  false, "bool", "Flag indicating that beta_p/beta_r are sensitivities (inverse) [default false]");
   config.add("aosys.D"                 ,"", "aosys.D"               , argType::Required, "aosys", "D",                false, "real", "The telescope diameter [m]");
   config.add("aosys.d_min"             ,"", "aosys.d_min"           , argType::Required, "aosys", "d_min",            false, "real", "The minimum actuator spacing [m]");
   config.add("aosys.optd"              ,"", "aosys.optd"            , argType::Optional, "aosys", "optd",             false, "bool", "Whether or not the actuator spacing is optimized");
   config.add("aosys.optd_delta"        ,"", "aosys.optd_delta"      , argType::Required, "aosys", "optd_delta",       false, "bool", "The fractional change from d_min used in optimization.  Set to 1 (default) for integer binnings, > 1 for finer sampling.");
   config.add("aosys.F0"                ,"", "aosys.F0"              , argType::Required, "aosys", "F0",               false, "real", "Zero-mag photon flux, [photons/sec]");     
   config.add("aosys.lam_wfs"           ,"", "aosys.lam_wfs"         , argType::Required, "aosys", "lam_wfs",          false, "real", "WFS wavelength [m]" );
   config.add("aosys.npix_wfs"          ,"", "aosys.npix_wfs"        , argType::Required, "aosys", "npix_wfs",         false, "vector<real>", "The number of pixels in the WFS");
   config.add("aosys.ron_wfs"           ,"", "aosys.ron_wfs"         , argType::Required, "aosys", "ron_wfs",          false, "vector<real>", "WFS readout noise [photons/read]");
   config.add("aosys.bin_npix"          ,"", "aosys.bin_npix"        , argType::Required, "aosys", "bin_npix",         false, "bool", "Whether or not WFS pixels are re-binned along with actuator spacing optimization");
   config.add("aosys.Fbg"               ,"", "aosys.Fbg"             , argType::Required, "aosys", "Fbg",              false, "vector<real>", "Background counts, [counts/pix/sec]");\
   config.add("aosys.tauWFS"            ,"", "aosys.tauWFS"          , argType::Required, "aosys", "tauWFS",           false, "real", "WFS integration time [s]");
   config.add("aosys.minTauWFS"         ,"", "aosys.minTauWFS"       , argType::Required, "aosys", "minTauWFS",        false, "vector<real>", "Minimum WFS integration time [s]");
   config.add("aosys.deltaTau"          ,"", "aosys.deltaTau"        , argType::Required, "aosys", "deltaTau",         false, "real", "Loop delay [s]");
   config.add("aosys.optTau"            ,"", "aosys.optTau"          , argType::Optional, "aosys", "optTau",           false, "bool", "Whether or not the integration time is optimized");
   config.add("aosys.lam_sci"           ,"", "aosys.lam_sci"         , argType::Required, "aosys", "lam_sci",          false, "real", "Science wavelength [m]");
   config.add("aosys.zeta"              ,"", "aosys.zeta"            , argType::Required, "aosys", "zeta",             false, "real", "Zenith distance [rad]");
   config.add("aosys.fit_mn_max"        ,"", "aosys.fit_mn_max"      , argType::Required, "aosys", "fit_mn_max",       false, "real", "Maximum spatial frequency index to use for analysis");
   config.add("aosys.circularLimit"     ,"", "aosys.circularLimit"   , argType::Optional, "aosys", "circularLimit",    false, "bool", " Flag to indicate that the spatial frequency limit is circular, not square.");
   config.add("aosys.spatialFilter_ku", "", "aosys.spatialFilter_ku", argType::Required, "aosys", "spatialFilter_ku", false, "real", "Spatial filter cutoff frequency in u [m^-1]");
   config.add("aosys.spatialFilter_kv", "", "aosys.spatialFilter_kv", argType::Required, "aosys", "spatialFilter_kv", false, "real", "Spatial filter cutoff frequency in v [m^-1]");
   config.add("aosys.ncp_wfe"           ,"", "aosys.ncp_wfe"         , argType::Required, "aosys", "ncp_wfe",          false, "real", "NCP WFE between 1 lambda/D and fit_mn_max [rad^2]");
   config.add("aosys.ncp_alpha"         ,"", "aosys.ncp_alpha"       , argType::Required, "aosys", "ncp_alpha",        false, "real", "PSD index for NCP WFE");
   config.add("aosys.starMag"           ,"", "aosys.starMag"         , argType::Required, "aosys", "starMag",          false, "real", "Star magnitude");
   config.add("aosys.starMags"          ,"", "aosys.starMags"        , argType::Required, "aosys", "starMags",         false, "real vector", "A vector of star magnitudes");
   
   atm.setupConfig(config);
   psd.setupConfig(config);
   
}

template<typename realT, class inputSpectT, typename iosT>
void aoSystem<realT, inputSpectT, iosT>::loadConfig( app::appConfigurator & config )
{
   //WFS 
   if(config.isSet("aosys.wfs"))
   {
      std::string wfsStr;
      config(wfsStr, "aosys.wfs");
      
      if(wfsStr == "ideal")
      {
         wfsBeta<wfs<realT>>(nullptr);
      }
      else if(wfsStr == "unmodPyWFS")
      {
         wfsBeta<pywfsUnmod<realT>>(nullptr);
      }
      else if(wfsStr == "asympModPyWFS")
      {
         wfsBeta<pywfsModAsymptotic<realT>>(nullptr);
      }
      else if(wfsStr == "SHWFS")
      {
         wfsBeta<shwfs<realT>>(nullptr);
      }
      else if(wfsStr == "calculatedWFS")
      {
         wfsBeta<calculatedWFS<realT>>(nullptr);

         calculatedWFS<realT> * cwfs = static_cast<calculatedWFS<realT> *>(m_wfsBeta);
         config(cwfs->m_beta_p_file, "aosys.wfs_beta_p");
         config(cwfs->m_beta_r_file, "aosys.wfs_beta_r");
         bool sens = cwfs->m_sensitivity;
         config(sens, "aosys.wfs_sensitivity");
         if(config.isSet("aosys.wfs_sensitivity")) 
         {
            cwfs->m_sensitivity = sens;
         }
      }
      else
      {
         mxThrowException(err::invalidconfig, "aoSystem::loadConfig", "unknown WFS " + wfsStr + " specified");  
      }
   }

   //We load the default value, get the config value (which will be the default if not set), and
   //then if it was set, call the setter function so any side effects are captured.

   //diameter
   realT nD = D();
   config(nD, "aosys.D");
   if(config.isSet("aosys.D") ) D(nD);
   
   //d_min
   std::vector<realT> nd_min = d_min();
   config(nd_min, "aosys.d_min");
   if(config.isSet("aosys.d_min")) d_min(nd_min);

   bool noptd = optd();
   config(noptd, "aosys.optd");
   if(config.isSet("aosys.optd")) optd(noptd);
   
   realT noptd_delta =optd_delta();
   config(noptd_delta, "aosys.optd_delta");
   if(config.isSet("aosys.optd_delta")) optd_delta(noptd_delta);

   realT nlam_wfs = lam_wfs();
   config(nlam_wfs, "aosys.lam_wfs");
   if(config.isSet("aosys.lam_wfs") ) lam_wfs(nlam_wfs);

   //npix_wfs   
   std::vector<realT> nnpix_wfs = npix_wfs();
   config(nnpix_wfs, "aosys.npix_wfs");   
   if(config.isSet("aosys.npix_wfs")) npix_wfs(nnpix_wfs);

   //ron_wfs   
   std::vector<realT> nron_wfs = ron_wfs();
   config(nron_wfs, "aosys.ron_wfs");   
   if(config.isSet("aosys.ron_wfs")) ron_wfs(nron_wfs);

   //Fbg   
   std::vector<realT> nFbg = Fbg();
   config(nFbg, "aosys.Fbg");   
   if(config.isSet("aosys.Fbg")) Fbg(nFbg);

   //minTauWFS   
   std::vector<realT> nminTauWFS = minTauWFS();
   config(nminTauWFS, "aosys.minTauWFS");   
   if(config.isSet("aosys.minTauWFS")) minTauWFS(nminTauWFS);

   //bin_npix
   bool nbin_npix = bin_npix();
   config(nbin_npix, "aosys.bin_npix");
   if(config.isSet("aosys.bin_npix")) bin_npix(nbin_npix);

   //tauWFS
   realT ntauWFS = tauWFS();
   config( ntauWFS, "aosys.tauWFS");
   if(config.isSet("aosys.tauWFS")) tauWFS(ntauWFS);
   
   //deltaTau
   realT ndeltaTau = deltaTau();
   config( ndeltaTau, "aosys.deltaTau");
   if(config.isSet("aosys.deltaTau")) deltaTau(ndeltaTau);

   //optTau   
   bool noptTau = optTau();
   config(noptTau, "aosys.optTau");
   if(config.isSet("aosys.optTau")) optTau(noptTau);

   //lam_sci
   realT nlam_sci = lam_sci();
   config( nlam_sci, "aosys.lam_sci");
   if(config.isSet("aosys.lam_sci") ) lam_sci( nlam_sci );
 
   //zeta
   realT nzeta = zeta();
   config(nzeta , "aosys.zeta");
   if(config.isSet("aosys.zeta") ) zeta(nzeta);

   //fit_mn_max
   realT fmnm = fit_mn_max();
   config( fmnm, "aosys.fit_mn_max");
   if(config.isSet("aosys.fit_mn_max") ) fit_mn_max(fmnm);
   
   //circularLimit
   bool cl = circularLimit();
   config(cl, "aosys.circularLimit");
   if(config.isSet("aosys.circularLimit") ) circularLimit(cl);

   //spatialFilter_ku
   realT ku = spatialFilter_ku();
   config( ku, "aosys.spatialFilter_ku");
   if(config.isSet("aosys.spatialFilter_ku") ) spatialFilter_ku(ku);
     
   //spatialFilter_kv
   realT kv = spatialFilter_kv();
   config( kv, "aosys.spatialFilter_kv");
   if(config.isSet("aosys.spatialFilter_kv") ) spatialFilter_kv(kv);

   //ncp_wfe
   realT nwfe = ncp_wfe();
   config( nwfe, "aosys.ncp_wfe");
   if(config.isSet("aosys.ncp_wfe") ) ncp_wfe(nwfe);
    
   //ncp_alpha
   realT na = ncp_alpha();
   config( na, "aosys.ncp_alpha");
   if(config.isSet("aosys.ncp_alpha") ) ncp_alpha( na );

   //F0
   realT nF0 = F0();
   config(nF0, "aosys.F0");
   if(config.isSet("aosys.F0") ) F0(nF0);

   //star_mag
   realT smag = starMag();
   config( smag, "aosys.starMag");
   if(config.isSet("aosys.starMag") ) starMag( smag );

   atm.loadConfig(config);
   psd.loadConfig(config);
}

extern template
class aoSystem<float, vonKarmanSpectrum<float>, std::ostream>; 

extern template
class aoSystem<double, vonKarmanSpectrum<double>, std::ostream>; 

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
