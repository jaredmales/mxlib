/** \file fourierTemporalPSD.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Calculation of the temporal PSD of Fourier modes.
 * \ingroup mxAOm_files
 *
 */

//***********************************************************************//
// Copyright 2016-2022 Jared R. Males (jaredmales@gmail.com)
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

#ifndef fourierTemporalPSD_hpp
#define fourierTemporalPSD_hpp

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>

#include <sys/stat.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Eigen/Dense>

#include "../../mxlib.hpp"
#include "../../math/constants.hpp"
#include "../../math/floatUtils.hpp"
#include "../../math/func/jinc.hpp"
#include "../../math/func/airyPattern.hpp"
#include "../../math/vectorUtils.hpp"
#include "../../ioutils/fits/fitsFile.hpp"
#include "../../sigproc/fourierModes.hpp"
#include "../../sigproc/psdVarMean.hpp"
#include "../../ioutils/stringUtils.hpp"
#include "../../ioutils/readColumns.hpp"
#include "../../ioutils/binVector.hpp"
#include "../../ioutils/fileUtils.hpp"

#include "../../ipc/ompLoopWatcher.hpp"

#include "aoSystem.hpp"
#include "aoPSDs.hpp"
#include "wfsNoisePSD.hpp"
#include "clAOLinearPredictor.hpp"
#include "clGainOpt.hpp"
#include "varmapToImage.hpp"
#include "speckleAmpPSD.hpp"

#include "aoConstants.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

#ifndef WSZ

/** \def WFZ
 * \brief Size of the GSL integration workspace
 */
#define WSZ 100000

#endif

enum basis : unsigned int
{
    basic,   ///< The basic sine and cosine Fourier modes
    modified ///< The modified Fourier basis from \cite males_guyon_2017
};

/// Policy for handling GSL quadrature non-convergence statuses.
enum class fourierTemporalPSDPolicy
{
    permissive, ///< Retain the best finite approximation and record the status.
    strict      ///< Record every status and return an error if any integration does not converge.
};

/// Aggregated GSL quadrature diagnostics for a Fourier temporal PSD calculation.
template <typename realT>
struct fourierTemporalPSDReport
{
    /// Summary of one GSL status code.
    struct statusSummary
    {
        size_t count{ 0 };                     ///< Number of occurrences of this status.
        std::map<size_t, size_t> countByLayer; ///< Number of occurrences in each atmospheric layer.
        realT maximumAbsoluteError{ 0 };       ///< Largest GSL absolute-error estimate.
        realT maximumToleranceRatio{ 0 };      ///< Largest error estimate relative to the requested tolerance.
        size_t worstLayer{ 0 };                ///< Layer containing the largest tolerance ratio.
        realT worstFrequency{ 0 };             ///< Frequency containing the largest tolerance ratio.
    };

    size_t integrationsAttempted{ 0 };      ///< Total number of quadrature calls.
    size_t integrationsConverged{ 0 };      ///< Number of quadrature calls returning `GSL_SUCCESS`.
    std::map<int, statusSummary> gslStatus; ///< Summaries keyed by the raw GSL status code.

    /// Reset all accumulated diagnostics.
    void clear();

    /// Record one quadrature result.
    void record( int status,              /**< [in] raw GSL status code */
                 size_t layer,            /**< [in] atmospheric layer index */
                 realT frequency,         /**< [in] temporal frequency */
                 realT result,            /**< [in] quadrature result */
                 realT absoluteError,     /**< [in] GSL absolute-error estimate */
                 realT absoluteTolerance, /**< [in] requested absolute tolerance */
                 realT relativeTolerance /**< [in] requested relative tolerance */ );

    /// Merge another report into this report.
    void merge( const fourierTemporalPSDReport &other /**< [in] report to merge */ );

    /// Return the total number of non-successful integrations.
    [[nodiscard]] size_t failureCount() const;

    /// Write a human-readable summary of the accumulated quadrature diagnostics.
    void write( std::ostream &output /**< [out] stream receiving the summary */ ) const;
};

/// \cond fourierTemporalPSD_detail
namespace fourierTemporalPSD_detail
{

/// Function type used to allocate a GSL integration workspace.
using gslWorkspaceAllocator = gsl_integration_workspace *(*)( size_t );

/// Deleter providing RAII ownership for a GSL integration workspace.
struct gslWorkspaceDeleter
{
    /// Free an allocated GSL integration workspace.
    void operator()( gsl_integration_workspace *workspace /**< [in] workspace to free */ ) const noexcept
    {
        if( workspace != nullptr )
        {
            gsl_integration_workspace_free( workspace );
        }
    }
};

/// Unique ownership handle for a GSL integration workspace.
using gslWorkspacePtr = std::unique_ptr<gsl_integration_workspace, gslWorkspaceDeleter>;

/// Return whether a GSL status represents a potentially usable non-converged approximation.
inline bool isConvergenceStatus( int status )
{
    return status == GSL_EMAXITER || status == GSL_EROUND || status == GSL_ESING || status == GSL_EDIVERGE;
}

/// Convert a fatal GSL status to an mxlib status.
inline error_t gslStatusToError( int status )
{
    if( status == GSL_ENOMEM )
    {
        return error_t::allocerr;
    }

    if( status == GSL_EDOM || status == GSL_EINVAL )
    {
        return error_t::invalidconfig;
    }

    return error_t::liberr;
}

/// Apply the configured non-convergence policy to one GSL status.
inline error_t applyPolicy( int status, fourierTemporalPSDPolicy policy )
{
    if( status == GSL_SUCCESS )
    {
        return error_t::noerror;
    }

    if( isConvergenceStatus( status ) )
    {
        return policy == fourierTemporalPSDPolicy::permissive ? error_t::noerror : error_t::liberr;
    }

    return gslStatusToError( status );
}

/// Mutex serializing scoped changes to GSL's process-global error handler.
inline std::mutex &gslErrorHandlerMutex()
{
    static std::mutex mutex;
    return mutex;
}

/// Disable the GSL error handler for a complete top-level PSD calculation and restore it on exit.
/** The mutex serializes handler changes made by this implementation. Unrelated code cannot be protected from GSL's
 * process-global handler state unless it coordinates with the same mutex.
 */
class scopedGslErrorHandlerOff
{
  public:
    /// Lock handler management and retain the previously installed handler.
    scopedGslErrorHandlerOff() : m_lock( gslErrorHandlerMutex() ), m_previous( gsl_set_error_handler_off() )
    {
    }

    /// Disallow copying ownership of the saved handler.
    scopedGslErrorHandlerOff( const scopedGslErrorHandlerOff & ) = delete;

    /// Disallow copy assignment of the handler guard.
    scopedGslErrorHandlerOff &operator=( const scopedGslErrorHandlerOff & ) = delete;

    /// Restore the previously installed handler before releasing the lock.
    ~scopedGslErrorHandlerOff()
    {
        static_cast<void>( gsl_set_error_handler( m_previous ) );
    }

  private:
    std::unique_lock<std::mutex> m_lock;        ///< Lock held while the handler is disabled.
    gsl_error_handler_t *m_previous{ nullptr }; ///< Handler restored on destruction.
};

} // namespace fourierTemporalPSD_detail
/// \endcond

template <typename realT>
void fourierTemporalPSDReport<realT>::clear()
{
    integrationsAttempted = 0;
    integrationsConverged = 0;
    gslStatus.clear();
}

template <typename realT>
void fourierTemporalPSDReport<realT>::record( int status,
                                              size_t layer,
                                              realT frequency,
                                              realT result,
                                              realT absoluteError,
                                              realT absoluteTolerance,
                                              realT relativeTolerance )
{
    ++integrationsAttempted;
    if( status == GSL_SUCCESS )
    {
        ++integrationsConverged;
        return;
    }

    statusSummary &summary = gslStatus[status];
    ++summary.count;
    ++summary.countByLayer[layer];
    summary.maximumAbsoluteError = std::max( summary.maximumAbsoluteError, std::abs( absoluteError ) );

    const realT requestedTolerance = std::max( std::abs( absoluteTolerance ), std::abs( relativeTolerance * result ) );
    const realT toleranceRatio = requestedTolerance > 0 ? std::abs( absoluteError ) / requestedTolerance
                                                        : std::numeric_limits<realT>::infinity();
    if( toleranceRatio >= summary.maximumToleranceRatio )
    {
        summary.maximumToleranceRatio = toleranceRatio;
        summary.worstLayer = layer;
        summary.worstFrequency = frequency;
    }
}

template <typename realT>
void fourierTemporalPSDReport<realT>::merge( const fourierTemporalPSDReport &other )
{
    integrationsAttempted += other.integrationsAttempted;
    integrationsConverged += other.integrationsConverged;

    for( const auto &[status, otherSummary] : other.gslStatus )
    {
        statusSummary &summary = gslStatus[status];
        summary.count += otherSummary.count;
        for( const auto &[layer, count] : otherSummary.countByLayer )
        {
            summary.countByLayer[layer] += count;
        }
        summary.maximumAbsoluteError = std::max( summary.maximumAbsoluteError, otherSummary.maximumAbsoluteError );
        if( otherSummary.maximumToleranceRatio >= summary.maximumToleranceRatio )
        {
            summary.maximumToleranceRatio = otherSummary.maximumToleranceRatio;
            summary.worstLayer = otherSummary.worstLayer;
            summary.worstFrequency = otherSummary.worstFrequency;
        }
    }
}

template <typename realT>
size_t fourierTemporalPSDReport<realT>::failureCount() const
{
    return integrationsAttempted - integrationsConverged;
}

template <typename realT>
void fourierTemporalPSDReport<realT>::write( std::ostream &output ) const
{
    output << "GSL quadrature: " << integrationsConverged << '/' << integrationsAttempted << " converged\n";
    for( const auto &[status, summary] : gslStatus )
    {
        output << "  " << gsl_strerror( status ) << " (" << status << "): " << summary.count << ", max absolute error "
               << summary.maximumAbsoluteError << ", max tolerance ratio " << summary.maximumToleranceRatio
               << " at layer " << summary.worstLayer << ", frequency " << summary.worstFrequency << ", layers {";
        bool firstLayer = true;
        for( const auto &[layer, count] : summary.countByLayer )
        {
            if( !firstLayer )
            {
                output << ", ";
            }
            output << layer << ": " << count;
            firstLayer = false;
        }
        output << "}\n";
    }
}

// Forward declaration
template <typename realT, typename aosysT>
realT F_basic( realT kv, void *params );

// Forward declaration
template <typename realT, typename aosysT>
realT F_mod( realT kv, void *params );

/// Class to manage the calculation of temporal PSDs of the Fourier modes in atmospheric turbulence.
/** Works with both basic (sines/cosines) and modified Fourier modes.
 *
 * \tparam realT is a real floating point type for calculations.  Currently must be double due to gsl_integration.
 * \tparam aosysT is an AO system type, usually of type ao_system.
 *
 * \todo Split off the integration parameters in a separate structure.
 * \todo once integration parameters are in a separate structure, make this a class with protected members.
 * \ingroup mxAOAnalytic
 */
template <typename _realT, typename aosysT>
struct fourierTemporalPSD
{
    /// The type for arithmetic
    typedef _realT realT;

    /// The complex type for arithmetic
    typedef std::complex<realT> complexT;

    /// Quadrature report type used by this specialization.
    typedef fourierTemporalPSDReport<realT> reportT;

    /// Pointer to an AO system structure.
    aosysT *m_aosys{ nullptr };

    realT m_f{ 0 };                 ///< the current temporal frequency
    realT m_m{ 0 };                 ///< the spatial frequency m index
    realT m_n{ 0 };                 ///< the spatial frequency n index
    realT m_cq{ 0 };                ///< The cosine of the wind direction
    realT m_sq{ 0 };                ///< The sine of the wind direction
    realT m_spatialFilter{ false }; ///< Flag indicating if a spatial filter is applied

    bool m_strehlOG{ false };
    bool m_uncorrectedOG{ false };

    realT m_f0{ 0 }; ///< the Berdja boiling parameter

    int m_p{ 1 }; ///< The parity of the mode, +/- 1.  If _useBasis==MXAO_FTPSD_BASIS_BASIC then +1 indicates cosine, -1
                  ///< indicates sine.
    int _layer_i; ///< The index of the current layer.

    int _useBasis; ///< Set to  MXAO_FTPSD_BASIS_BASIC/MODIFIED/PROJECTED_* to use the basic sin/cos modes, the modified
                   ///< Fourier modes, or a projection of them.

  protected:
    /// Unique ownership of the GSL integration workspace used by worker instances.
    fourierTemporalPSD_detail::gslWorkspacePtr m_workspace;

    /// Allocation function used when a worker lazily creates its GSL workspace.
    fourierTemporalPSD_detail::gslWorkspaceAllocator m_workspaceAllocator{ gsl_integration_workspace_alloc };

  public:
    realT _absTol;                            ///< The absolute tolerance to use in the GSL integrator
    realT _relTol;                            ///< The relative tolerance to use in the GSL integrator

    int m_mode_i;                             ///< Projected basis mode index

    Eigen::Array<realT, -1, -1> m_modeCoeffs; ///< Coeeficients of the projection onto the Fourier modes
    realT m_minCoeffVal;

    std::vector<realT> Jps;
    std::vector<realT> Jms;
    std::vector<int> ps;
    std::vector<realT> ms;
    std::vector<realT> ns;

    void initProjection()
    {
        Jps.resize( m_modeCoeffs.cols() );
        Jms.resize( m_modeCoeffs.cols() );
        ps.resize( m_modeCoeffs.cols() );
        ms.resize( m_modeCoeffs.cols() );
        ns.resize( m_modeCoeffs.cols() );

        for( int i = 0; i < m_modeCoeffs.cols(); ++i )
        {
            int m, n, p;
            sigproc::fourierModeCoordinates( m, n, p, i );
            ps[i] = p;
            ms[i] = m;
            ns[i] = n;
        }
    }

  public:
    /// Default c'tor
    fourierTemporalPSD();

    /// Disallow copying unique workspace ownership.
    fourierTemporalPSD( const fourierTemporalPSD & ) = delete;

    /// Disallow copy assignment of unique workspace ownership.
    fourierTemporalPSD &operator=( const fourierTemporalPSD & ) = delete;

    /// Move workspace ownership and evaluator state.
    fourierTemporalPSD( fourierTemporalPSD && ) noexcept = default;

    /// Move-assign workspace ownership and evaluator state.
    fourierTemporalPSD &operator=( fourierTemporalPSD && ) noexcept = default;

    /// Release owned resources.
    ~fourierTemporalPSD() = default;

  protected:
    /// Construct with a custom workspace allocator.
    explicit fourierTemporalPSD(
        fourierTemporalPSD_detail::gslWorkspaceAllocator allocator /**< [in] workspace allocation function */ );

    /// Initialize parameters to default values.
    void initialize();

    /// Allocate the worker workspace if it is not already available.
    error_t allocateWorkspace();

    /// Validate state and arguments shared by single- and multilayer calculations.
    error_t validatePsdInputs( const std::vector<realT> &PSD,  /**< [in] output storage to validate */
                               const std::vector<realT> &freq, /**< [in] temporal-frequency grid */
                               realT m,                        /**< [in] first spatial-frequency index */
                               realT n,                        /**< [in] second spatial-frequency index */
                               int p,                          /**< [in] Fourier-mode parity */
                               realT fmax,                     /**< [in] maximum exactly integrated frequency */
                               int layer_i,                    /**< [in] layer index, or -1 to validate all layers */
                               fourierTemporalPSDPolicy policy /**< [in] non-convergence policy */ );

    /// Validate the configured atmosphere and optionally a requested layer.
    error_t validateAtmosphere( int layer_i /**< [in] layer index, or -1 to validate all layers */ );

  public:
    /** \name GSL Integration Tolerances
     * For good results it seems that absolute tolerance (absTol) needs to be 1e-10.  Lower tolerances cause some
     * frequencies to drop out, etc. Relative tolerance (relTol) seems to be less sensitive, and 1e-4 works on cases
     * tested as of 1 Jan, 2017.
     *
     * See the documentation for the GSL Library integrators at
     * (https://www.gnu.org/software/gsl/manual/htmlm_node/QAGI-adaptive-integration-on-infinite-intervals.html)
     * @{
     */

    /// Set absolute tolerance
    /**
     * \param at is the new absolute tolerance.
     */
    void absTol( realT at );

    /// Get the current absolute tolerance
    /**
     * \returns _absTol
     */
    realT absTol();

    /// Set relative tolerance
    /**
     * \param rt is the new relative tolerance.
     */
    void relTol( realT rt );

    /// Get the current relative tolerance
    /**
     * \returns _relTol
     */
    realT relTol();

    ///@}

    /// Determine the frequency of the highest V-dot-k peak
    /**
     * \param m the spatial frequency u index
     * \param n the spatial frequency v index
     *
     * \return the frequency of the fastest peak
     */
    realT fastestPeak( int m, int n );

  protected:
    /// Calculate a single-layer temporal PSD while the caller manages the GSL error handler.
    error_t singleLayerPSDImpl( std::vector<realT> &PSD,  /**< [out] calculated PSD */
                                std::vector<realT> &freq, /**< [in] temporal-frequency grid */
                                realT m,                  /**< [in] first spatial-frequency index */
                                realT n,                  /**< [in] second spatial-frequency index */
                                int layer_i,              /**< [in] atmospheric-layer index */
                                int p,                    /**< [in] Fourier-mode parity */
                                realT fmax,               /**< [in] maximum exactly integrated frequency */
                                reportT &report,          /**< [out] accumulated quadrature report */
                                fourierTemporalPSDPolicy policy /**< [in] non-convergence policy */ );

  public:
    /// Calculate the temporal PSD for a Fourier mode for a single layer.
    /** `PSD` and `freq` must have the same nonzero size. Frequencies must be finite, nonnegative, and strictly
     * increasing. The AO system, integration controls, requested layer, and atmosphere are validated before
     * calculation. A precondition or allocation failure leaves `PSD` unchanged and clears `report` when supplied.
     *
     * When extending beyond `fmax`, up to the last 50 exactly integrated bins are averaged after projection to the
     * first tail frequency. If fewer than 50 exact bins are available, all available exact bins are used. At least one
     * exact bin is required to initialize the tail.
     *
     * In permissive mode, finite best approximations returned with `GSL_EMAXITER`, `GSL_EROUND`, `GSL_ESING`, or
     * `GSL_EDIVERGE` are retained and summarized in `report`. In strict mode the calculation continues to characterize
     * all such failures but returns `error_t::liberr` and the output PSD must be discarded.
     *
     * \returns `error_t::noerror` on success, an argument/configuration error for invalid inputs, or
     * `error_t::liberr` when strict quadrature handling detects non-convergence.
     */
    error_t singleLayerPSD(
        std::vector<realT> &PSD,   /**< [out] calculated PSD */
        std::vector<realT> &freq,  /**< [in] temporal-frequency grid */
        realT m,                   /**< [in] first spatial-frequency index */
        realT n,                   /**< [in] second spatial-frequency index */
        int layer_i,               /**< [in] atmospheric-layer index */
        int p,                     /**< [in] Fourier-mode parity */
        realT fmax = 0,            /**< [in] maximum exactly integrated frequency, or 0 for the grid maximum */
        reportT *report = nullptr, /**< [out] optional quadrature report */
        fourierTemporalPSDPolicy policy = fourierTemporalPSDPolicy::permissive /**< [in] non-convergence policy */ );

    ///\cond multilayerm_parallel
    // Conditional to exclude from Doxygen.

  protected:
    // Type to allow overloading of the multiLayerPSD workers based on whether they are parallelized or not.
    template <bool m_parallel>
    struct isParallel
    {
    };

    // Parallelized version of multiLayerPSD, with OMP directives.
    error_t multiLayerPSDImpl( std::vector<realT> &PSD,         /**< [out] calculated PSD */
                               std::vector<realT> &freq,        /**< [in] temporal-frequency grid */
                               realT m,                         /**< [in] first spatial-frequency index */
                               realT n,                         /**< [in] second spatial-frequency index */
                               int p,                           /**< [in] Fourier-mode parity */
                               realT fmax,                      /**< [in] maximum exactly integrated frequency */
                               reportT &report,                 /**< [out] accumulated quadrature report */
                               fourierTemporalPSDPolicy policy, /**< [in] non-convergence policy */
                               isParallel<true> parallel /**< [in] parallel dispatch tag */ );

    // Non-Parallelized version of multiLayerPSD, without OMP directives.
    error_t multiLayerPSDImpl( std::vector<realT> &PSD,         /**< [out] calculated PSD */
                               std::vector<realT> &freq,        /**< [in] temporal-frequency grid */
                               realT m,                         /**< [in] first spatial-frequency index */
                               realT n,                         /**< [in] second spatial-frequency index */
                               int p,                           /**< [in] Fourier-mode parity */
                               realT fmax,                      /**< [in] maximum exactly integrated frequency */
                               reportT &report,                 /**< [out] accumulated quadrature report */
                               fourierTemporalPSDPolicy policy, /**< [in] non-convergence policy */
                               isParallel<false> parallel /**< [in] sequential dispatch tag */ );

    ///\endcond

  public:
    /// Calculate the temporal PSD for a Fourier mode in a multi-layer model.
    /** `PSD` and `freq` must have the same nonzero size. Frequencies must be finite, nonnegative, and strictly
     * increasing. The AO system, integration controls, and complete atmosphere are validated before calculation. A
     * precondition failure leaves `PSD` unchanged and clears `report` when supplied.
     *
     * \tparam parallel controls whether layers are calculated in parallel.  Default is true.  Set to false if this is
     * called inside a parallelized loop, as in \ref makePSDGrid.
     *
     * In permissive mode, recognized convergence failures are retained and summarized in `report`. Strict mode returns
     * `error_t::liberr` if any layer has such a failure; the output PSD is incomplete and must be discarded whenever
     * this function returns an error.
     *
     * \returns `error_t::noerror` on success, or the first layer error in atmospheric-layer order.
     */
    template <bool parallel = true>
    error_t multiLayerPSD(
        std::vector<realT> &PSD,   /**< [out] calculated PSD */
        std::vector<realT> &freq,  /**< [in] temporal-frequency grid */
        realT m,                   /**< [in] first spatial-frequency index */
        realT n,                   /**< [in] second spatial-frequency index */
        int p,                     /**< [in] Fourier-mode parity */
        realT fmax = 0,            /**< [in] maximum exactly integrated frequency, or 0 for the default cutoff */
        reportT *report = nullptr, /**< [out] optional quadrature report */
        fourierTemporalPSDPolicy policy = fourierTemporalPSDPolicy::permissive /**< [in] non-convergence policy */ );

    /// Calculate PSDs over a grid of spatial frequencies.
    /** The grid of spatial frequencies is square, set by the maximum value of m and n.
     *
     * The PSDs are written as mx::binVector binary files to a directory.  We do not use FITS since
     * this adds overhead and cfitisio handles parallelization poorly due to the limitation on number of created file
     * pointers.
     *
     */
    void makePSDGrid( const std::string &dir, ///< [in] the directory for output of the PSDs
                      int mnMax,              ///< [in] the maximum value of m and n in the grid.
                      realT dFreq,            ///< [in] the temporal frequency spacing.
                      realT maxFreq,          ///< [in] the maximum temporal frequency to calculate
                      realT fmax = 0 ///< [in] [optional] set the maximum temporal frequency for the calculation. The
                                     ///< PSD is filled in with a -17/3 power law past
                                     /// this frequency.  If 0, then it is taken to be 150 Hz + 2*fastestPeak(m,n).
    );

    /// Analyze a PSD grid under closed-loop control.
    /** This always analyzes the simple integrator, and can also analyze the linear predictor controller. Outputs maps
     * of optimum gains, predictor coefficients, variances, and contrasts for a range of guide star magnitudes.
     * Optionally calculates speckle lifetimes.  Optionally writes the closed-loop PSDs and transfer functions.
     */
    int analyzePSDGrid( const std::string &subDir, /**< [out] the sub-directory of psdDir where to write the
                                                        results. Is created. */
                        const std::string &psdDir, ///< [in]  the directory containing the grid of PSDs.
                        int mnMax,                 ///< [in]  the maximum value of m and n in the grid.
                        int mnCon,                 ///< [in]  the maximum value of m and n which can be controlled.
                        realT gfixed,              ///< [in]  if \> 0 then this fixed gain is used in the SI.
                        int lpNc,                  /**< [in]  the number of linear predictor coefficients to analyze.
                                                            If 0 then LP is not analyzed.*/
                        realT lpRegPrecision,      /**< [in]  the initial precision for the LP regularization
                                                              algorithm. Normal value is 2. Higher is faster. Decrease
                                                              if getting stuck in local minima.*/
                        std::vector<realT> &mags,  ///< [in]  the guide star magnitudes to analyze for.
                        int lifetimeTrials = 0,    /**< [in] [optional] number of trials used for calculating speckle
                                                    lifetimes.  If 0,lifetimes are not calculated. */
                        bool ucLifeTs = false,     /**< [in] [optional] flag controlling whether lifetimes are
                                                                              calculated for uncontrolled modes.*/
                        bool writePSDs = false,    ///< [in] [optional] flag controlling if resultant PSDs are saved
                        bool writeXfer = false     /**< [in] [optional] flag controlling if resultant
                                                                        transfer functions are saved*/
    );

    int intensityPSD( const std::string &subDir,  // sub-directory of psdDir which contains the controlled system
                                                  // results, and where the lifetimes will be written.
                      const std::string &psdDir,  // directory containing the grid of PSDs
                      const std::string &CvdPath, // path to the covariance decomposition
                      int mnMax,                  ///< [in]  the maximum value of m and n in the grid.
                      int mnCon,                  ///< [in]  the maximum value of m and n which can be controlled.
                      std::vector<realT> &mags,   ///< [in]  the guide star magnitudes
                      int lifetimeTrials,         /**< [in]  [optional] number of trials used for calculating
                                                                        speckle lifetimes.  If 0, lifetimes are not
                                                                        calculated.*/
                      bool writePSDs              /**< [in]  [optional] flag controlling if resultant
                                                                        PSDs are saved*/ );

    /** \name Disk Storage
     * These methods handle writing to and reading from disk.  The calculated PSDs are store in the mx::BinVector binary
     * format.
     *
     * A grid of PSDs is specified by its directory name.  The directory contains one frequency file (freq.binv), and a
     * set of PSD files, named according to psd_\<m\>_\<n\>_.binv.
     *
     *
     * @{
     */
    /// Get the frequency scale for a PSD grid.
    /**
     */
    int getGridFreq( std::vector<realT> &freq, ///< [out] the vector to populate with the frequency scale.
                     const std::string &dir    ///< [in] specifies the directory containing the grid.
    );

    /// Get a single PSD from a PSD grid.
    /**
     */
    int getGridPSD( std::vector<realT> &psd, ///< [out] the vector to populate with the PSD.
                    const std::string &dir,  ///< [in] specifies the directory containing the grid.
                    int m,                   ///< [in] specifies the u component of spatial frequency.
                    int n                    ///< [in] specifies the v component of spatial frequency.
    );

    /// Get both the frequency scale and a single PSD from a PSD grid.
    /**
     */
    int getGridPSD( std::vector<realT> &freq, ///< [out] the vector to populate with the frequency scale.
                    std::vector<realT> &psd,  ///< [out] the vector to populate with the PSD.
                    const std::string &dir,   ///< [in] specifies the directory containing the grid.
                    int m,                    ///< [in] specifies the u component of spatial frequency.
                    int n                     ///< [in] specifies the v component of spatial frequency.
    );

    ///@}
};

template <typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::fourierTemporalPSD()
{
    m_aosys = nullptr;
    initialize();
}

template <typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::fourierTemporalPSD( fourierTemporalPSD_detail::gslWorkspaceAllocator allocator )
    : m_workspaceAllocator( allocator )
{
    m_aosys = nullptr;
    initialize();
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::allocateWorkspace()
{
    if( m_workspace != nullptr )
    {
        return error_t::noerror;
    }

    if( m_workspaceAllocator == nullptr )
    {
        return internal::mxlib_error_report( error_t::invalidconfig, "GSL workspace allocator is null" );
    }

    m_workspace.reset( m_workspaceAllocator( WSZ ) );
    if( m_workspace == nullptr )
    {
        return internal::mxlib_error_report( error_t::allocerr, "could not allocate GSL integration workspace" );
    }

    return error_t::noerror;
}

template <typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::initialize()
{
    _useBasis = basis::modified;

    _absTol = 1e-10;
    _relTol = 1e-4;
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::validateAtmosphere( int layer_i )
{
    if( m_aosys == nullptr )
    {
        return internal::mxlib_error_report( error_t::invalidconfig, "AO system pointer is null" );
    }

    if( !math::isFinite( m_aosys->D() ) || m_aosys->D() <= 0 )
    {
        return internal::mxlib_error_report( error_t::invalidconfig,
                                             "AO system aperture diameter must be finite and positive" );
    }

    auto &atmosphere = m_aosys->atm;
    const size_t layerCount = atmosphere.n_layers();
    if( layerCount == 0 )
    {
        return internal::mxlib_error_report( error_t::sizeerr, "atmosphere must contain at least one layer" );
    }

    const std::vector<realT> outerScale = atmosphere.L_0();
    const std::vector<realT> innerScale = atmosphere.l_0();
    const std::vector<realT> altitude = atmosphere.layer_z();
    const std::vector<realT> strength = atmosphere.layer_Cn2();
    const std::vector<realT> windSpeed = atmosphere.layer_v_wind();
    const std::vector<realT> windDirection = atmosphere.layer_dir();
    if( outerScale.size() != layerCount || innerScale.size() != layerCount || altitude.size() != layerCount ||
        strength.size() != layerCount || windSpeed.size() != layerCount || windDirection.size() != layerCount )
    {
        return internal::mxlib_error_report( error_t::sizeerr, "atmosphere layer-vector sizes do not match" );
    }

    if( atmosphere.nonKolmogorov() )
    {
        const std::vector<realT> alpha = atmosphere.alpha();
        const std::vector<realT> beta = atmosphere.beta();
        const std::vector<realT> beta0 = atmosphere.beta_0();
        if( alpha.size() != layerCount || beta.size() != layerCount || beta0.size() != layerCount )
        {
            return internal::mxlib_error_report( error_t::sizeerr,
                                                 "non-Kolmogorov atmosphere vector sizes do not match" );
        }

        for( size_t index = 0; index < layerCount; ++index )
        {
            if( !math::isFinite( alpha[index] ) || !math::isFinite( beta[index] ) || beta[index] <= 0 ||
                !math::isFinite( beta0[index] ) || beta0[index] < 0 )
            {
                return internal::mxlib_error_report( error_t::invalidconfig,
                                                     "non-Kolmogorov atmosphere parameters are invalid" );
            }
        }
    }
    else if( !math::isFinite( atmosphere.r_0() ) || atmosphere.r_0() <= 0 || !math::isFinite( atmosphere.lam_0() ) ||
             atmosphere.lam_0() <= 0 )
    {
        return internal::mxlib_error_report( error_t::invalidconfig,
                                             "atmosphere Fried parameter and wavelength must be finite and positive" );
    }

    bool hasPositiveStrength = false;
    realT totalStrength = 0;
    for( size_t index = 0; index < layerCount; ++index )
    {
        if( !math::isFinite( outerScale[index] ) || outerScale[index] <= 0 || !math::isFinite( innerScale[index] ) ||
            innerScale[index] < 0 || !math::isFinite( altitude[index] ) || altitude[index] < 0 ||
            !math::isFinite( strength[index] ) || strength[index] < 0 || !math::isFinite( windSpeed[index] ) ||
            windSpeed[index] <= 0 || !math::isFinite( windDirection[index] ) )
        {
            return internal::mxlib_error_report( error_t::invalidconfig, "atmosphere layer parameters are invalid" );
        }
        hasPositiveStrength = hasPositiveStrength || strength[index] > 0;
        totalStrength += strength[index];
    }

    if( !hasPositiveStrength || !math::isFinite( totalStrength ) || totalStrength <= 0 )
    {
        return internal::mxlib_error_report( error_t::invalidconfig,
                                             "atmosphere must contain a positive layer strength" );
    }

    if( layer_i < -1 || ( layer_i >= 0 && static_cast<size_t>( layer_i ) >= layerCount ) )
    {
        return internal::mxlib_error_report( error_t::invalidarg, "atmosphere layer index is out of range" );
    }

    return error_t::noerror;
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::validatePsdInputs( const std::vector<realT> &PSD,
                                                              const std::vector<realT> &freq,
                                                              realT m,
                                                              realT n,
                                                              int p,
                                                              realT fmax,
                                                              int layer_i,
                                                              fourierTemporalPSDPolicy policy )
{
    if( freq.empty() || PSD.size() != freq.size() )
    {
        return internal::mxlib_error_report( error_t::sizeerr,
                                             "PSD and frequency vectors must have the same nonzero size" );
    }

    for( size_t index = 0; index < freq.size(); ++index )
    {
        if( !math::isFinite( freq[index] ) || freq[index] < 0 || ( index > 0 && freq[index] <= freq[index - 1] ) )
        {
            return internal::mxlib_error_report(
                error_t::invalidarg,
                "frequency grid must be finite, nonnegative, and strictly increasing" );
        }
    }

    if( !math::isFinite( m ) || !math::isFinite( n ) || !math::isFinite( fmax ) || fmax < 0 )
    {
        return internal::mxlib_error_report(
            error_t::invalidarg,
            "mode coordinates must be finite and frequency cutoff must be finite and nonnegative" );
    }

    if( p != -1 && p != 1 )
    {
        return internal::mxlib_error_report( error_t::invalidarg, "Fourier-mode parity must be -1 or +1" );
    }

    if( _useBasis != basis::basic && _useBasis != basis::modified )
    {
        return internal::mxlib_error_report( error_t::invalidarg, "value of _useBasis is not valid" );
    }

    if( policy != fourierTemporalPSDPolicy::permissive && policy != fourierTemporalPSDPolicy::strict )
    {
        return internal::mxlib_error_report( error_t::invalidarg, "quadrature policy is not valid" );
    }

    if( !math::isFinite( _absTol ) || _absTol <= 0 || !math::isFinite( _relTol ) || _relTol <= 0 || _relTol >= 1 )
    {
        return internal::mxlib_error_report(
            error_t::invalidconfig,
            "GSL absolute tolerance must be positive and relative tolerance must be between zero and one" );
    }

    if( !math::isFinite( m_f0 ) || m_f0 < 0 )
    {
        return internal::mxlib_error_report( error_t::invalidconfig,
                                             "turbulence boiling parameter must be finite and nonnegative" );
    }

    const error_t atmosphereStatus = validateAtmosphere( layer_i );
    if( atmosphereStatus != error_t::noerror )
    {
        return atmosphereStatus;
    }

    if( _useBasis == basis::modified )
    {
        if( !math::isFinite( m_aosys->lam_sci() ) || m_aosys->lam_sci() <= 0 || !math::isFinite( m_aosys->lam_wfs() ) ||
            m_aosys->lam_wfs() <= 0 || !math::isFinite( m_aosys->zeta() ) ||
            std::abs( m_aosys->zeta() ) >= math::half_pi<realT>() )
        {
            return internal::mxlib_error_report(
                error_t::invalidconfig,
                "AO wavelengths must be positive and zenith angle must lie strictly between -pi/2 and pi/2" );
        }
    }

    if( !math::isFinite( m_aosys->spatialFilter_ku() ) || m_aosys->spatialFilter_ku() <= 0 ||
        !math::isFinite( m_aosys->spatialFilter_kv() ) || m_aosys->spatialFilter_kv() <= 0 )
    {
        return internal::mxlib_error_report( error_t::invalidconfig,
                                             "AO spatial-filter limits must be finite and positive" );
    }

    return error_t::noerror;
}

template <typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::absTol( realT at )
{
    _absTol = at;
}

template <typename realT, typename aosysT>
realT fourierTemporalPSD<realT, aosysT>::absTol()
{
    return _absTol;
}

template <typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::relTol( realT rt )
{
    _relTol = rt;
}

template <typename realT, typename aosysT>
realT fourierTemporalPSD<realT, aosysT>::relTol()
{
    return _relTol;
}

template <typename realT, typename aosysT>
realT fourierTemporalPSD<realT, aosysT>::fastestPeak( int m, int n )
{
    realT ku, kv, vu, vv;

    ku = ( (realT)m / m_aosys->D() );
    kv = ( (realT)n / m_aosys->D() );

    realT f, fmax = 0;

    for( size_t i = 0; i < m_aosys->atm.n_layers(); ++i )
    {
        vu = m_aosys->atm.layer_v_wind( i ) * cos( m_aosys->atm.layer_dir( i ) );
        vv = m_aosys->atm.layer_v_wind( i ) * sin( m_aosys->atm.layer_dir( i ) );

        f = fabs( ku * vu + kv * vv );
        if( f > fmax )
            fmax = f;
    }

    return fmax;
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::singleLayerPSD( std::vector<realT> &PSD,
                                                           std::vector<realT> &freq,
                                                           realT m,
                                                           realT n,
                                                           int layer_i,
                                                           int p,
                                                           realT fmax,
                                                           reportT *report,
                                                           fourierTemporalPSDPolicy policy )
{
    reportT localReport;
    reportT &activeReport = report == nullptr ? localReport : *report;
    activeReport.clear();

    const error_t status = validatePsdInputs( PSD, freq, m, n, p, fmax, layer_i, policy );
    if( status != error_t::noerror )
    {
        return status;
    }

    fourierTemporalPSD_detail::scopedGslErrorHandlerOff handlerGuard;
    return singleLayerPSDImpl( PSD, freq, m, n, layer_i, p, fmax, activeReport, policy );
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::singleLayerPSDImpl( std::vector<realT> &PSD,
                                                               std::vector<realT> &freq,
                                                               realT m,
                                                               realT n,
                                                               int layer_i,
                                                               int p,
                                                               realT fmax,
                                                               reportT &report,
                                                               fourierTemporalPSDPolicy policy )
{
    if( fmax == 0 )
        fmax = freq[freq.size() - 1];

    if( freq[0] > fmax )
    {
        return internal::mxlib_error_report(
            error_t::invalidarg,
            "at least one exact frequency bin is required to initialize the PSD tail" );
    }

    realT v_wind = m_aosys->atm.layer_v_wind( layer_i );
    realT q_wind = m_aosys->atm.layer_dir( layer_i );

    // Rotate the basis
    realT cq = cos( q_wind );
    realT sq = sin( q_wind );

    realT scale = 2 * ( 1 / v_wind ); // Factor of 2 for negative frequencies.

    // Create a local instance so that we're reentrant
    fourierTemporalPSD<realT, aosysT> params( m_workspaceAllocator );
    const error_t workspaceStatus = params.allocateWorkspace();
    if( workspaceStatus != error_t::noerror )
    {
        return workspaceStatus;
    }

    params.m_aosys = m_aosys;
    params._layer_i = layer_i;
    params.m_m = m * cq + n * sq;
    params.m_n = -m * sq + n * cq;
    params.m_cq = cq; // for de-rotating ku and kv for spatial filtering
    params.m_sq = sq; // for de-rotation ku and kv for spatial filtering
    if( m_aosys->spatialFilter_ku() < std::numeric_limits<realT>::max() ||
        m_aosys->spatialFilter_kv() < std::numeric_limits<realT>::max() )
        params.m_spatialFilter = true;

    params.m_p = p;
    params.m_f0 = m_f0;

    params.m_mode_i = m_mode_i;
    params.m_modeCoeffs = m_modeCoeffs;
    params.m_minCoeffVal = m_minCoeffVal;

    realT result{ 0 };
    realT error{ 0 };
    error_t returnStatus = error_t::noerror;

    // Setup the GSL calculation
    gsl_function func;
    switch( _useBasis )
    {
        case basis::basic: // MXAO_FTPSD_BASIS_BASIC:
            func.function = &F_basic<realT, aosysT>;
            break;
        case basis::modified: // MXAO_FTPSD_BASIS_MODIFIED:
            func.function = &F_mod<realT, aosysT>;
            break;
        default:
            return internal::mxlib_error_report( error_t::invalidarg, "value of _useBasis is not valid." );
    }

    func.params = &params;

    // Here we only calculate up to fmax.
    size_t i = 0;
    while( freq[i] <= fmax )
    {
        params.m_f = freq[i];

        const int ec = gsl_integration_qagi( &func, _absTol, _relTol, WSZ, params.m_workspace.get(), &result, &error );
        report.record( ec, static_cast<size_t>( layer_i ), freq[i], result, error, _absTol, _relTol );

        const error_t integrationStatus = fourierTemporalPSD_detail::applyPolicy( ec, policy );
        if( integrationStatus != error_t::noerror )
        {
            if( !fourierTemporalPSD_detail::isConvergenceStatus( ec ) )
            {
                return internal::mxlib_error_report( integrationStatus,
                                                     std::string( "gsl_integration_qagi failed: " ) +
                                                         gsl_strerror( ec ) );
            }

            returnStatus = integrationStatus;
        }

        if( !math::isFinite( result ) || !math::isFinite( error ) )
        {
            return internal::mxlib_error_report( error_t::liberr,
                                                 "gsl_integration_qagi returned a nonfinite result or error estimate" );
        }

        PSD[i] = scale * result;

        ++i;
        if( i >= freq.size() )
            break;
    }

    // Now fill in from fmax to the actual max frequency with a -(alpha+2) power law.
    size_t j = i;

    if( j == freq.size() )
        return returnStatus;

    // First average up to the last 50 exactly integrated bins after projecting them to the first tail frequency.
    constexpr size_t maximumTailAverageCount = 50;
    const size_t tailAverageCount = std::min( i, maximumTailAverageCount );
    PSD[j] = 0;
    for( size_t k = tailAverageCount; k > 0; --k )
    {
        PSD[j] +=
            PSD[i - k] * pow( freq[i - k] / freq[j], m_aosys->atm.alpha( layer_i ) + 2 ); // seventeen_thirds<realT>());
    }
    PSD[j] /= static_cast<realT>( tailAverageCount );
    ++j;
    ++i;
    if( j == freq.size() )
        return returnStatus;
    while( j < freq.size() )
    {
        PSD[j] =
            PSD[i - 1] * pow( freq[i - 1] / freq[j], m_aosys->atm.alpha( layer_i ) + 2 ); // seventeen_thirds<realT>());
        ++j;
    }

    return returnStatus;
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::multiLayerPSDImpl( std::vector<realT> &PSD,
                                                              std::vector<realT> &freq,
                                                              realT m,
                                                              realT n,
                                                              int p,
                                                              realT fmax,
                                                              reportT &report,
                                                              fourierTemporalPSDPolicy policy,
                                                              isParallel<true> parallel )
{
    static_cast<void>( parallel );

    const size_t layerCount = m_aosys->atm.n_layers();
    std::vector<error_t> layerStatus( layerCount, error_t::noerror );
    std::vector<reportT> layerReport( layerCount );

#pragma omp parallel
    {
        // Records each layer PSD
        std::vector<realT> single_PSD( freq.size() );

#pragma omp for
        for( size_t i = 0; i < m_aosys->atm.n_layers(); ++i )
        {
            std::fill( single_PSD.begin(), single_PSD.end(), 0 );
            layerStatus[i] =
                singleLayerPSDImpl( single_PSD, freq, m, n, static_cast<int>( i ), p, fmax, layerReport[i], policy );

// Now add the single layer PSD to the overall PSD, weighted by Cn2
#pragma omp critical
            if( layerStatus[i] == error_t::noerror )
            {
                for( size_t j = 0; j < freq.size(); ++j )
                {
                    PSD[j] += m_aosys->atm.layer_Cn2( i ) * single_PSD[j];
                }
            }
        }
    }

    error_t returnStatus = error_t::noerror;
    for( size_t i = 0; i < layerCount; ++i )
    {
        report.merge( layerReport[i] );
        if( returnStatus == error_t::noerror && layerStatus[i] != error_t::noerror )
        {
            returnStatus = layerStatus[i];
        }
    }

    return returnStatus;
}

template <typename realT, typename aosysT>
error_t fourierTemporalPSD<realT, aosysT>::multiLayerPSDImpl( std::vector<realT> &PSD,
                                                              std::vector<realT> &freq,
                                                              realT m,
                                                              realT n,
                                                              int p,
                                                              realT fmax,
                                                              reportT &report,
                                                              fourierTemporalPSDPolicy policy,
                                                              isParallel<false> parallel )
{
    static_cast<void>( parallel );

    // Records each layer PSD
    std::vector<realT> single_PSD( freq.size() );
    error_t returnStatus = error_t::noerror;

    for( size_t i = 0; i < m_aosys->atm.n_layers(); ++i )
    {
        std::fill( single_PSD.begin(), single_PSD.end(), 0 );
        reportT layerReport;
        const error_t layerStatus =
            singleLayerPSDImpl( single_PSD, freq, m, n, static_cast<int>( i ), p, fmax, layerReport, policy );
        report.merge( layerReport );
        if( returnStatus == error_t::noerror && layerStatus != error_t::noerror )
        {
            returnStatus = layerStatus;
        }

        // Now add the single layer PSD to the overall PSD, weighted by Cn2
        if( layerStatus == error_t::noerror )
        {
            for( size_t j = 0; j < freq.size(); ++j )
            {
                PSD[j] += m_aosys->atm.layer_Cn2( i ) * single_PSD[j];
            }
        }
    }

    return returnStatus;
}

template <typename realT, typename aosysT>
template <bool parallel>
error_t fourierTemporalPSD<realT, aosysT>::multiLayerPSD( std::vector<realT> &PSD,
                                                          std::vector<realT> &freq,
                                                          realT m,
                                                          realT n,
                                                          int p,
                                                          realT fmax,
                                                          reportT *report,
                                                          fourierTemporalPSDPolicy policy )
{
    reportT localReport;
    reportT &activeReport = report == nullptr ? localReport : *report;
    activeReport.clear();

    const error_t validationStatus = validatePsdInputs( PSD, freq, m, n, p, fmax, -1, policy );
    if( validationStatus != error_t::noerror )
    {
        return validationStatus;
    }

    // PSD is zeroed every time to make sure we don't accumulate on repeated calls
    for( size_t j = 0; j < PSD.size(); ++j )
        PSD[j] = 0;

    if( fmax == 0 )
    {
        fmax = 150 + 2 * fastestPeak( m, n );
    }

    fourierTemporalPSD_detail::scopedGslErrorHandlerOff handlerGuard;
    return multiLayerPSDImpl( PSD, freq, m, n, p, fmax, activeReport, policy, isParallel<parallel>() );
}

template <typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::makePSDGrid(
    const std::string &dir, int mnMax, realT dFreq, realT maxFreq, realT fmax )
{
    std::vector<realT> freq;

    std::vector<sigproc::fourierModeDef> spf;

    std::string fn;

    sigproc::makeFourierModeFreqs_Rect( spf, 2 * mnMax );

    // Calculate number of samples, and make sure we get to at least maxFreq
    int N = (int)( maxFreq / dFreq );
    if( N * dFreq < maxFreq )
        N += 1;

    /*** Dump Params to file ***/
    mkdir( dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

    std::ofstream fout;
    fn = dir + '/' + "params.txt";
    fout.open( fn );

    fout << "#---------------------------\n";
    m_aosys->dumpAOSystem( fout );
    fout << "#---------------------------\n";
    fout << "# PSD Grid Parameters\n";
    fout << "#    absTol " << _absTol << '\n';
    fout << "#    relTol " << _relTol << '\n';
    fout << "#    useBasis " << _useBasis << '\n';
    fout << "#    makePSDGrid call:\n";
    fout << "#       mnMax = " << mnMax << '\n';
    fout << "#       dFreq = " << dFreq << '\n';
    fout << "#       maxFreq = " << maxFreq << '\n';
    fout << "#       fmax = " << fmax << '\n';
    fout << "#---------------------------\n";

    fout.close();

    // Make directory
    std::string psddir = dir + "/psds";
    mkdir( psddir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

    // Create frequency scale.
    math::vectorScale( freq, N, dFreq, 0 ); // dFreq); //offset from 0 by dFreq, so f=0 not included

    fn = psddir + '/' + "freq.binv";

    ioutils::writeBinVector( fn, freq );

    size_t nLoops = 0.5 * spf.size();

    ipc::ompLoopWatcher<> watcher( nLoops, std::cout );

#pragma omp parallel
    {
        std::vector<realT> PSD;
        PSD.resize( N );
        std::string fname;

        int m, n;

#pragma omp for
        for( size_t i = 0; i < nLoops; ++i )
        {
            m = spf[i * 2].m;
            n = spf[i * 2].n;

            if( fabs( (realT)m / m_aosys->D() ) >= m_aosys->spatialFilter_ku() ||
                fabs( (realT)n / m_aosys->D() ) >= m_aosys->spatialFilter_kv() )
            {
                watcher.incrementAndOutputStatus();
                continue;
            }

            multiLayerPSD<false>( PSD, freq, m, n, 1, fmax );

            fname = std::format( "{}/psd_{}_{}.binv", psddir, m, n );
            //    psddir + '/' + "psd_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) +
            //    ".binv";

            ioutils::writeBinVector( fname, PSD );

            watcher.incrementAndOutputStatus();
        }
    }
}

template <typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::analyzePSDGrid( const std::string &subDir,
                                                       const std::string &psdDir,
                                                       int mnMax,
                                                       int mnCon,
                                                       realT gfixed,
                                                       int lpNc,
                                                       realT lpRegPrecision,
                                                       std::vector<realT> &mags,
                                                       int lifetimeTrials,
                                                       bool uncontrolledLifetimes,
                                                       bool writePSDs,
                                                       bool writeXfer )
{

    std::string dir = psdDir + "/" + subDir;

    /*** Dump Params to file ***/
    mkdir( dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

    std::ofstream fout;
    std::string fn = dir + '/' + "params.txt";
    fout.open( fn );

    fout << "#---------------------------\n";
    m_aosys->dumpAOSystem( fout );
    fout << "#---------------------------\n";
    fout << "# Analysis Parameters\n";
    fout << "#    mnMax    = " << mnMax << '\n';
    fout << "#    mnCon    = " << mnCon << '\n';
    fout << "#    lpNc     = " << lpNc << '\n';
    fout << "#    mags     = ";
    for( size_t i = 0; i < mags.size() - 1; ++i )
        fout << mags[i] << ",";
    fout << mags[mags.size() - 1] << '\n';
    fout << "#    lifetimeTrials = " << lifetimeTrials << '\n';
    fout << "#    uncontrolledLifetimes = " << uncontrolledLifetimes << '\n';
    fout << "#    writePSDs = " << std::boolalpha << writePSDs << '\n';
    fout << "#    writeXfer = " << std::boolalpha << writeXfer << '\n';

    fout.close();

    //**** Calculating A Variance Map ****//

    realT fs = 1.0 / m_aosys->tauWFS();
    realT tauWFS = m_aosys->tauWFS();
    realT deltaTau = m_aosys->deltaTau();

    std::vector<sigproc::fourierModeDef> fms;

    sigproc::makeFourierModeFreqs_Rect( fms, 2 * mnMax );
    size_t nModes = 0.5 * fms.size();

    Eigen::Array<realT, -1, -1> gains, vars, speckleLifetimes, gains_lp, vars_lp, speckleLifetimes_lp;

    gains.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
    vars.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
    speckleLifetimes.resize( 2 * mnMax + 1, 2 * mnMax + 1 );

    gains( mnMax, mnMax ) = 0;
    vars( mnMax, mnMax ) = 0;
    speckleLifetimes( mnMax, mnMax ) = 0;

    gains_lp.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
    vars_lp.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
    speckleLifetimes_lp.resize( 2 * mnMax + 1, 2 * mnMax + 1 );

    gains_lp( mnMax, mnMax ) = 0;
    vars_lp( mnMax, mnMax ) = 0;
    speckleLifetimes_lp( mnMax, mnMax ) = 0;

    bool doLP = false;
    if( lpNc > 1 )
        doLP = true;
    Eigen::Array<realT, -1, -1> lpC;

    if( doLP )
    {
        lpC.resize( nModes, lpNc );
        lpC.setZero();
    }

    std::vector<realT> S_si, S_lp;

    if( writePSDs )
    {
        for( size_t s = 0; s < mags.size(); ++s )
        {
            std::string psdOutDir = std::format( "{}/outputPSDS_{}_si", dir, mags[s] );
            // dir + "/" + "outputPSDs_" + ioutils::convert ToString( mags[s] ) + "_si";
            mkdir( psdOutDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

            if( doLP )
            {
                std::string psdOutDir = std::format( "{}/outputPSDS_{}_lp", dir, mags[s] );
                // dir + "/" + "outputPSDs_" + ioutils::convert ToString( mags[s] ) + "_lp";
                mkdir( psdOutDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
            }
        }
    }

    m_aosys->beta_p( 1, 1 );

    ipc::ompLoopWatcher<> watcher( nModes * mags.size(), std::cout );
    std::atomic<int> analysisStatus{ static_cast<int>( error_t::noerror ) };

    for( size_t s = 0; s < mags.size(); ++s )
    {
        m_aosys->starMag( mags[s] );

        // In non-parallel space calculate OG=Strehl

        realT opticalGain{ 1.0 };

        if( m_strehlOG )
        {
            // Iterative optical gain
            /// \todo need upstream NCP and CP NCP and NCP NCP
            realT ncp = m_aosys->ncp_wfe(); // save ncp and then set it to zero for this part.
            m_aosys->ncp_wfe( 0 );

            realT lam_sci = m_aosys->lam_sci();
            m_aosys->lam_sci( m_aosys->lam_wfs() );

            realT S = m_aosys->strehl();
            std::cerr << S << "\n";

            for( int s = 0; s < 4; ++s )
            {
                m_aosys->opticalGain( S );
                m_aosys->optd( m_aosys->optd() ); // just trigger a re-calc
                S = m_aosys->strehl();
                std::cerr << S << "\n";
            }

            opticalGain = S;

            m_aosys->lam_sci( lam_sci );
            m_aosys->ncp_wfe( ncp );
        }

        m_aosys->optd( m_aosys->optd() ); // just trigger a re-calc
        realT strehl = m_aosys->strehl();

#pragma omp parallel
        {
            realT localMag = mags[s];

            realT var0;

            realT gopt, var;

            realT gopt_lp, var_lp;

            std::vector<realT> tfreq; // The frequency scale of the PSDs
            std::vector<realT> tPSDp; // The open-loop turbulence PSD for a Fourier mode
            std::vector<realT>
                tPSDpPOL; // The pseudo-open-loop turbulence PSD for a Fourier mode, with optical gain effects included

            std::vector<realT> tfreqHF; // The above-Nyquist frequencies, saved if outputing the PSDS.
            std::vector<realT> tPSDpHF; // The above-Nyquist component of the open-loop PSD, saved if outputing the
                                        // PSDs.

            //**< Get the frequency grid, and nyquist limit it to f_s/2
            getGridPSD( tfreq, tPSDp, psdDir, 0, 1 ); // To get the freq grid

            size_t imax = 0;
            while( tfreq[imax] <= 0.5 * fs )
            {
                ++imax;
                if( imax > tfreq.size() - 1 )
                    break;
            }

            if( imax < tfreq.size() - 1 && tfreq[imax] <= 0.5 * fs * ( 1.0 + 1e-7 ) )
            {
                ++imax;
            }

            if( writePSDs )
            {
                tfreqHF.assign( tfreq.begin(), tfreq.end() );
            }

            tfreq.erase( tfreq.begin() + imax, tfreq.end() );

            tPSDpPOL.resize( tfreq.size() ); // pre=allocate
            //**>

            std::vector<realT> tPSDn; // The open-loop WFS noise PSD
            tPSDn.resize( tfreq.size() );

            //**< Setup the controllers
            mx::AO::analysis::clAOLinearPredictor<realT> tflp;
            tflp.m_precision0 = lpRegPrecision;
            /*tflp.m_min_sc0 = 0;
            tflp.m_max_sc0 = 1000;*/

            mx::AO::analysis::clGainOpt<realT> go_si( tauWFS, deltaTau );
            mx::AO::analysis::clGainOpt<realT> go_lp( tauWFS, deltaTau );

            go_si.f( tfreq );
            go_lp.f( tfreq );

            realT gmax = 0;
            realT gmax_lp = 0;
            //**>

            int m, n;

            //**< For use in lifetime calculations
            sigproc::psdVarMean<sigproc::psdVarMeanParams<realT>> pvm;
            std::vector<std::complex<realT>> ETFxn;
            std::vector<std::complex<realT>> NTFxn;

            if( lifetimeTrials > 0 )
            {
                ETFxn.resize( tfreq.size() );
                NTFxn.resize( tfreq.size() );

                if( writeXfer )
                {
                    std::string tfOutFile = std::format( "{}/outputTF_{}_si", dir, mags[s] );
                    // dir + "/" + "outputTF_" + ioutils::convert ToString( mags[s] ) + "_si/";
                    ioutils::createDirectories( tfOutFile );
                }

                if( doLP )
                {
                    if( writeXfer )
                    {
                        std::string tfOutFile = std::format( "{}/outputTF_{}_lp", dir, mags[s] );
                        // dir + "/" + "outputTF_" + ioutils::convert ToString( mags[s] ) + "_lp/";
                        ioutils::createDirectories( tfOutFile );
                    }
                }
            }

//**>

/*#pragma omp critical
std::cerr << __FILE__ << " " << __LINE__ << "\n";
*/
// want to schedule dynamic with small chunks so maximal processor usage,
// otherwise we can end up with a small number of cores being used at the end
#pragma omp for schedule( dynamic, 5 )
            for( size_t i = 0; i < nModes; ++i )
            {
                // Determine the spatial frequency at this step
                m = fms[2 * i].m;
                n = fms[2 * i].n;

                if( fabs( (realT)m / m_aosys->D() ) >= m_aosys->spatialFilter_ku() ||
                    fabs( (realT)n / m_aosys->D() ) >= m_aosys->spatialFilter_kv() )
                {
                    gains( mnMax + m, mnMax + n ) = 0;
                    gains( mnMax - m, mnMax - n ) = 0;

                    gains_lp( mnMax + m, mnMax + n ) = 0;
                    gains_lp( mnMax - m, mnMax - n ) = 0;

                    vars( mnMax + m, mnMax + n ) = 0;
                    vars( mnMax - m, mnMax - n ) = 0;

                    vars_lp( mnMax + m, mnMax + n ) = 0;
                    vars_lp( mnMax - m, mnMax - n ) = 0;
                    speckleLifetimes( mnMax + m, mnMax + n ) = 0;
                    speckleLifetimes( mnMax - m, mnMax - n ) = 0;
                    speckleLifetimes_lp( mnMax + m, mnMax + n ) = 0;
                    speckleLifetimes_lp( mnMax - m, mnMax - n ) = 0;
                }
                else
                {

                    realT k = sqrt( m * m + n * n ) / m_aosys->D();

                    //**< Get the open-loop turb. PSD
                    getGridPSD( tPSDp, psdDir, m, n );

                    // Get integral of entire open-loop PSD
                    var0 = sigproc::psdVar( tfreq, tPSDp );

                    if( writePSDs )
                    {
                        tPSDpHF.assign( tPSDp.begin() + imax, tPSDp.end() );
                    }

                    // erase points above Nyquist limit
                    tPSDp.erase( tPSDp.begin() + imax, tPSDp.end() );

                    // And now determine the variance which has been erased.
                    // limVar is the out-of-band variance, which we add back in for completeness
                    realT limVar = 0; // var0 - sigproc::psdVar( tfreq, tPSDp);

                    // And construct the POL psd
                    if( m_uncorrectedOG )
                    {
                        for( size_t n = 0; n < tPSDp.size(); ++n )
                        {
                            tPSDpPOL[n] = tPSDp[n] * pow( opticalGain, 2 );
                        }
                    }
                    else
                    {
                        for( size_t n = 0; n < tPSDp.size(); ++n )
                        {
                            tPSDpPOL[n] = tPSDp[n];
                        }
                    }
                    //**>

                    //**< Determine if we're inside the hardwarecontrol limit
                    bool inside = false;

                    if( m_aosys->circularLimit() )
                    {
                        if( m * m + n * n <= mnCon * mnCon )
                            inside = true;
                    }
                    else
                    {
                        if( fabs( m ) <= mnCon && fabs( n ) <= mnCon )
                            inside = true;
                    }
                    //**>

                    /* This is to select specific points for troubleshooting*/
                    // if( !( ( m ==-2 && n == 84 ) || (m==-2 && n == 2800)) ) inside = false;

                    if( inside )
                    {
                        // Get the WFS noise PSD (which is already resized to match tfreq)
                        wfsNoisePSD<realT>( tPSDn,
                                            m_aosys->beta_p( m, n ) / sqrt( opticalGain ),
                                            m_aosys->Fg( localMag ),
                                            tauWFS,
                                            m_aosys->npix_wfs( (size_t)0 ),
                                            m_aosys->Fbg( (size_t)0 ),
                                            m_aosys->ron_wfs( (size_t)0 ) );

                        gmax = 0;
                        if( gfixed > 0 )
                        {
                            gopt = gfixed;
                            var = go_si.clVariance( tPSDp, tPSDn, gopt );
                        }
                        else
                        {
                            // Calculate gain using the POL PSD
                            error_t gainStatus = go_si.optGainOpenLoop( gopt, var, tPSDpPOL, tPSDn, true );
                            if( gainStatus != error_t::noerror )
                            {
                                int expected = static_cast<int>( error_t::noerror );
                                analysisStatus.compare_exchange_strong( expected, static_cast<int>( gainStatus ) );
                                gopt = 0;
                                var = go_si.clVariance( tPSDp, tPSDn, gopt );
                            }

                            if( m_uncorrectedOG )
                            {
                                gopt *= opticalGain;
                            }

                            // But use the variance from the actual POL
                            var = go_si.clVariance( tPSDp, tPSDn, gopt );

                            // Output gain curve for this mode (for troubleshooting)
                            /*#pragma omp critical
                            {
                                std::string foutn = "gcurve_";
                                foutn += std::to_string(m) + "_" + std::to_string(n) + ".dat";

                                std::ofstream foutf(foutn);

                                for(size_t n = 0; n < 10000; ++n)
                                {
                                    realT gg = (1.0*n)/10000.;
                                    foutf << gg << " " << go_si.clVariance(tPSDp, tPSDn, gg) << "\n";
                                }

                                foutf.close();

                                std::cerr << "\n" << gmax << " " << gopt << " " << var << " " << go_si.clVariance(tPSDp,
                            tPSDn, 0.64) << "\n";
                            }*/
                        }

                        var += limVar;

                        if( doLP )
                        {
                            realT min_sc;
                            error_t rv = tflp.regularizeCoefficients( gmax_lp,
                                                                      gopt_lp,
                                                                      var_lp,
                                                                      min_sc,
                                                                      go_lp,
                                                                      tPSDpPOL,
                                                                      tPSDn,
                                                                      lpNc );

                            if( rv != error_t::noerror )
                            {
                                int expected = static_cast<int>( error_t::noerror );
                                analysisStatus.compare_exchange_strong( expected, static_cast<int>( rv ) );
                            }

                            for( int n = 0; n < lpNc; ++n )
                            {
                                lpC( i, n ) = go_lp.a( n );
                            }

                            if( m_uncorrectedOG )
                            {
                                go_lp.aScale( opticalGain );
                                go_lp.bScale( opticalGain );
                                gopt_lp *= opticalGain;
                            }

                            var_lp = go_lp.clVariance( tPSDp, tPSDn, gopt_lp );
                            var_lp += limVar;
                        }
                        else
                        {
                            gopt_lp = 0;
                        }
                    }
                    else
                    {
                        // Zero the noise PSD
                        tPSDn.assign( tPSDn.size(), 0.0 );

                        gopt = 0;
                        var = var0;
                        var_lp = var0;
                        gopt_lp = 0;
                        go_lp.a( std::vector<realT>( { 1 } ) );
                        go_lp.b( std::vector<realT>( { 1 } ) );
                    }

                    //**< Determine if closed-loop is making a difference:

                    if( gopt > 0 && var > var0 )
                    {
                        gopt = 0;
                        var = var0;
                    }

                    if( gopt_lp > gopt && var_lp > var )
                    {
                        // Set LP to SI (or off if SI is off)
                        gopt_lp = gopt;
                        var_lp = var;
                        go_lp.a( std::vector<realT>( { 1 } ) );
                        go_lp.b( std::vector<realT>( { 1 } ) );
                    }
                    //**>

                    //**< Fill in the gain and variance maps
                    gains( mnMax + m, mnMax + n ) = gopt;
                    gains( mnMax - m, mnMax - n ) = gopt;

                    gains_lp( mnMax + m, mnMax + n ) = gopt_lp;
                    gains_lp( mnMax - m, mnMax - n ) = gopt_lp;

                    vars( mnMax + m, mnMax + n ) = var;
                    vars( mnMax - m, mnMax - n ) = var;

                    vars_lp( mnMax + m, mnMax + n ) = var_lp;
                    vars_lp( mnMax - m, mnMax - n ) = var_lp;
                    //**>

                    //**< Calculate Speckle Lifetimes
                    if( ( lifetimeTrials > 0 || writeXfer ) && ( uncontrolledLifetimes || inside ) )
                    {
                        std::vector<realT> spfreq, sppsd;

                        if( gopt > 0 )
                        {
                            for( size_t i = 0; i < tfreq.size(); ++i )
                            {
                                ETFxn[i] = go_si.clETF( i, gopt );
                                NTFxn[i] = go_si.clNTF( i, gopt );
                            }
                        }
                        else
                        {
                            for( size_t i = 0; i < tfreq.size(); ++i )
                            {
                                ETFxn[i] = 1;
                                NTFxn[i] = 0;
                            }
                        }

                        if( writeXfer )
                        {
                            std::string tfOutFile = std::format( "{}/outputTF_{}_si", dir, mags[s] );
                            //    dir + "/" + "outputTF_" + ioutils::convert ToString( mags[s] ) + "_si/";

                            std::string etfOutFile = std::format( "{}/etf_{}_{}.binv", tfOutFile, m, n );
                            // tfOutFile + "etf_" + ioutils::convert ToString( m ) + '_' +
                            //                          ioutils::convert ToString( n ) + ".binv";
                            ioutils::writeBinVector( etfOutFile, ETFxn );

                            std::string ntfOutFile = std::format( "{}/ntf_{}_{}.binv", tfOutFile, m, n );
                            // tfOutFile + "ntf_" + ioutils::convert ToString( m ) + '_' +
                            //                          ioutils::convert ToString( n ) + ".binv";
                            ioutils::writeBinVector( ntfOutFile, NTFxn );

                            if( i == 0 ) // Write freq on the first one
                            {
                                std::string fOutFile = tfOutFile + "freq.binv";
                                ioutils::writeBinVector( fOutFile, tfreq );
                            }
                        }

                        if( lifetimeTrials > 0 )
                        {
                            speckleAmpPSD( spfreq, sppsd, tfreq, tPSDp, ETFxn, tPSDn, NTFxn, lifetimeTrials );
                            realT spvar = sigproc::psdVar( spfreq, sppsd );

                            realT splifeT = 100.0;
                            realT error;

                            realT tau = pvm( error, spfreq, sppsd, splifeT ) * ( splifeT ) / spvar;

                            speckleLifetimes( mnMax + m, mnMax + n ) = tau;
                            speckleLifetimes( mnMax - m, mnMax - n ) = tau;
                        }

                        if( doLP )
                        {
                            if( gopt_lp > 0 )
                            {
                                for( size_t i = 0; i < tfreq.size(); ++i )
                                {
                                    ETFxn[i] = go_lp.clETF( i, gopt_lp );
                                    NTFxn[i] = go_lp.clNTF( i, gopt_lp );
                                }
                            }
                            else
                            {
                                for( size_t i = 0; i < tfreq.size(); ++i )
                                {
                                    ETFxn[i] = 1;
                                    NTFxn[i] = 0;
                                }
                            }

                            if( writeXfer )
                            {
                                std::string tfOutFile = std::format( "{}/outputTF_{}_lp", dir, mags[s] );
                                // dir + "/" + "outputTF_" + ioutils::convert ToString( mags[s] ) + "_lp/";

                                std::string etfOutFile = std::format( "{}/etf_{}_{}.binv", tfOutFile, m, n );
                                // tfOutFile + "etf_" + ioutils::convert ToString( m ) + '_' +
                                //         ioutils::convert ToString( n ) + ".binv";
                                ioutils::writeBinVector( etfOutFile, ETFxn );

                                std::string ntfOutFile = std::format( "{}/ntf_{}_{}.binv", tfOutFile, m, n );
                                // tfOutFile + "ntf_" + ioutils::convert ToString( m ) + '_' +
                                //                        ioutils::convert ToString( n ) + ".binv";
                                ioutils::writeBinVector( ntfOutFile, NTFxn );

                                if( i == 0 ) // Write freq on the first one
                                {
                                    std::string fOutFile = tfOutFile + "freq.binv";
                                    ioutils::writeBinVector( fOutFile, tfreq );
                                }
                            }

                            if( lifetimeTrials > 0 )
                            {
                                speckleAmpPSD( spfreq, sppsd, tfreq, tPSDp, ETFxn, tPSDn, NTFxn, lifetimeTrials );
                                realT spvar = sigproc::psdVar( spfreq, sppsd );

                                realT splifeT = 100.0;
                                realT error;

                                realT tau = pvm( error, spfreq, sppsd, splifeT ) * ( splifeT ) / spvar;

                                speckleLifetimes_lp( mnMax + m, mnMax + n ) = tau;
                                speckleLifetimes_lp( mnMax - m, mnMax - n ) = tau;
                            }
                        } // if(doLP)
                    } // if( (lifetimeTrials > 0 || writeXfer) && ( uncontrolledLifetimes || inside ))
                    //**>

                    // Calculate the controlled PSDs and output
                    if( writePSDs )
                    {
                        std::string psdOutFile =
                            std::format( "{}/outputPSDs_{}_si/psd_{}_{}.binv", dir, mags[s], m, n );

                        //   dir + "/" + "outputPSDs_" + ioutils::convert ToString( mags[s] ) + "_si/";
                        // psdOutFile +=
                        //    "psd_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) + ".binv";

                        std::vector<realT> psdOut( tPSDp.size() + tPSDpHF.size() );

                        // Calculate the output PSD if gains are applied
                        if( gopt > 0 )
                        {
                            realT ETF, NTF;

                            for( size_t i = 0; i < tfreq.size(); ++i )
                            {
                                go_si.clTF2( ETF, NTF, i, gopt );
                                psdOut[i] = tPSDp[i] * ETF + tPSDn[i] * NTF;
                            }

                            for( size_t i = 0; i < tPSDpHF.size(); ++i )
                            {
                                psdOut[tfreq.size() + i] = tPSDpHF[i];
                            }
                        }
                        else // otherwise just copy
                        {
                            for( size_t i = 0; i < tfreq.size(); ++i )
                            {
                                psdOut[i] = tPSDp[i];
                            }

                            for( size_t i = 0; i < tPSDpHF.size(); ++i )
                            {
                                psdOut[tfreq.size() + i] = tPSDpHF[i];
                            }
                        }

                        ioutils::writeBinVector( psdOutFile, psdOut );

                        if( i == 0 ) // Write freq on the first one
                        {
                            psdOutFile = std::format( "{}/outputPSDs_{}_si/freq.binv", dir, mags[s] );
                            //    dir + "/" + "outputPSDs_" + ioutils::convert ToString( mags[s] ) + "_si/freq.binv";
                            // ioutils::writeBinVector( psdOutFile, tfreqHF );
                        }

                        if( doLP )
                        {
                            std::string psdOutFile =
                                std::format( "{}/outputPSDs_{}_lp/psd_{}_{}.binv", dir, mags[s], m, n );

                            // psdOutFile = dir + "/" + "outputPSDs_" + ioutils::convert ToString( mags[s] ) + "_lp/";
                            // psdOutFile +=
                            //     "psd_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) +
                            //     ".binv";

                            // Calculate the output PSD if gains are applied
                            if( gopt_lp > 0 )
                            {
                                realT ETF, NTF;

                                for( size_t i = 0; i < tfreq.size(); ++i )
                                {
                                    go_lp.clTF2( ETF, NTF, i, gopt_lp );
                                    psdOut[i] = tPSDp[i] * ETF + tPSDn[i] * NTF;
                                }
                                for( size_t i = 0; i < tPSDpHF.size(); ++i )
                                {
                                    psdOut[tfreq.size() + i] = tPSDpHF[i];
                                }
                            }
                            else // otherwise just copy
                            {
                                for( size_t i = 0; i < tfreq.size(); ++i )
                                {
                                    psdOut[i] = tPSDp[i];
                                }

                                for( size_t i = 0; i < tPSDpHF.size(); ++i )
                                {
                                    psdOut[tfreq.size() + i] = tPSDpHF[i];
                                }
                            }

                            ioutils::writeBinVector( psdOutFile, psdOut );

                            if( i == 0 )
                            {
                                psdOutFile = std::format( "{}/outputPSDs_{}_lp/freq.binv", dir, mags[s] );
                                // psdOutFile =
                                //     dir + "/" + "outputPSDs_" + ioutils::convert ToString( mags[s] ) +
                                //     "_lp/freq.binv";
                                // ioutils::writeBinVector( psdOutFile, tfreq );
                            }
                        }
                    }
                }
                watcher.incrementAndOutputStatus();

            } // omp for i..nModes
        } // omp Parallel

        if( analysisStatus.load() != static_cast<int>( error_t::noerror ) )
        {
            return analysisStatus.load();
        }

        Eigen::Array<realT, -1, -1> cim;

        fits::fitsFile<realT> ff;
        std::string fn = std::format( "{}/gainmap_{}_si.fits", dir, mags[s] );
        // dir + "/gainmap_" + ioutils::convert ToString( mags[s] ) + "_si.fits";
        ff.write( fn, gains );

        fn = std::format( "{}/varmap_{}_si.fits", dir, mags[s] );
        // dir + "/varmap_" + ioutils::convert ToString( mags[s] ) + "_si.fits";
        ff.write( fn, vars );

        cim = vars;

        realT Ssi = exp( -1 * cim.sum() );
        S_si.push_back( strehl );
        cim /= strehl;

        fn = std::format( "{}/contrast_{}_si.fits", dir, mags[s] );
        // dir + "/contrast_" + ioutils::convert ToString( mags[s] ) + "_si.fits";
        ff.write( fn, cim );

        if( lifetimeTrials > 0 )
        {
            fn = std::format( "{}/speckleLifetimes_{}_si.fits", dir, mags[s] );
            // dir + "/speckleLifetimes_" + ioutils::convert ToString( mags[s] ) + "_si.fits";
            ff.write( fn, speckleLifetimes );
        }

        if( doLP )
        {
            fn = std::format( "{}/gainmap_{}_lp.fits", dir, mags[s] );
            // dir + "/gainmap_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
            ff.write( fn, gains_lp );

            fn = std::format( "{}/lpcmap_{}_lp.fits", dir, mags[s] );
            // dir + "/lpcmap_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
            ff.write( fn, lpC );

            fn = std::format( "{}/varmap_{}_lp.fits", dir, mags[s] );
            // dir + "/varmap_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
            ff.write( fn, vars_lp );

            cim = vars_lp;

            // Scale Strehl by the ratio of total variance
            realT Slp = strehl * exp( -1 * cim.sum() ) /
                        Ssi; // This is a hack until we do a real fitting error calculation or something
            S_lp.push_back( Slp );
            cim /= Slp;

            fn = std::format( "{}/contrast_{}_lp.fits", dir, mags[s] );
            // dir + "/contrast_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
            ff.write( fn, cim );

            if( lifetimeTrials > 0 )
            {
                fn = std::format( "{}/speckleLifetimes_{}_lp.fits", dir, mags[s] );
                // dir + "/speckleLifetimes_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
                ff.write( fn, speckleLifetimes_lp );
            }
        }

    } // s (mag)

    fn = dir + "/strehl_si.txt";
    fout.open( fn );
    for( size_t i = 0; i < mags.size(); ++i )
    {
        fout << mags[i] << " " << S_si[i] << "\n";
    }

    fout.close();

    if( doLP )
    {
        fn = dir + "/strehl_lp.txt";
        fout.open( fn );
        for( size_t i = 0; i < mags.size(); ++i )
        {
            fout << mags[i] << " " << S_lp[i] << "\n";
        }

        fout.close();
    }

    return 0;
}

template <typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::intensityPSD(
    const std::string &subDir,  // sub-directory of psdDir which contains the controlled system results,
                                // and where the lifetimes will be written.
    const std::string &psdDir,  // directory containing the PSDS
    const std::string &CvdPath, // path to the covariance decomposition
    int mnMax,
    int mnCon,
    std::vector<realT> &mags,
    int lifetimeTrials,
    bool writePSDs )
{

    std::string dir = psdDir + "/" + subDir;

    /*** Dump Params to file ***/
    mkdir( dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

    std::ofstream fout;
    std::string fn = dir + '/' + "splife_params.txt";
    fout.open( fn );

    fout << "#---------------------------\n";
    m_aosys->dumpAOSystem( fout );
    fout << "#---------------------------\n";
    fout << "# Analysis Parameters\n";
    fout << "#    mnMax    = " << mnMax << '\n';
    fout << "#    mnCon    = " << mnCon << '\n';
    fout << "#    mags     = ";
    for( size_t i = 0; i < mags.size() - 1; ++i )
        fout << mags[i] << ",";
    fout << mags[mags.size() - 1] << '\n';
    fout << "#    lifetimeTrials = " << lifetimeTrials << '\n';
    // fout << "#    uncontrolledLifetimes = " << uncontrolledLifetimes << '\n';
    fout << "#    writePSDs = " << std::boolalpha << writePSDs << '\n';

    fout.close();

    realT fs = 1.0 / m_aosys->tauWFS();
    realT tauWFS = m_aosys->tauWFS();
    realT deltaTau = m_aosys->deltaTau();

    //** Get the Fourier Grid
    std::vector<sigproc::fourierModeDef> fms;

    sigproc::makeFourierModeFreqs_Rect( fms, 2 * mnMax );
    size_t nModes = 0.5 * fms.size();
    std::cerr << "nModes: " << nModes << " (" << fms.size() << ")\n";

    Eigen::Array<realT, -1, -1> speckleLifetimes;
    Eigen::Array<realT, -1, -1> speckleLifetimes_lp;

    speckleLifetimes.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
    speckleLifetimes( mnMax, mnMax ) = 0;

    speckleLifetimes_lp.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
    speckleLifetimes_lp( mnMax, mnMax ) = 0;

    /*********************************************************************/
    // 0) Get the frequency grid, and nyquist limit it to f_s/2
    /*********************************************************************/

    std::vector<realT> tfreq;
    std::vector<realT> tPSDp; // The open-loop OPD PSD
    std::vector<realT> tPSDn; // The open-loop WFS noise PSD
    std::vector<complexT> tETF;
    std::vector<complexT> tNTF;

    if( getGridFreq( tfreq, psdDir ) < 0 )
        return -1;

    size_t imax = 0;
    while( tfreq[imax] <= 0.5 * fs )
    {
        ++imax;
        if( imax > tfreq.size() - 1 )
            break;
    }

    if( imax < tfreq.size() - 1 && tfreq[imax] <= 0.5 * fs * ( 1.0 + 1e-7 ) )
        ++imax;

    tfreq.erase( tfreq.begin() + imax, tfreq.end() );

    // Now allocate memory
    tPSDn.resize( tfreq.size() );
    std::vector<std::vector<realT>> sqrtOPDPSD;
    sqrtOPDPSD.resize( nModes );

    std::vector<std::vector<realT>> opdPSD;
    opdPSD.resize( nModes );

    std::vector<realT> psd2sided;

    // Store the mode variance for later normalization
    std::vector<realT> modeVar;
    modeVar.resize( nModes );

    /*********************************************************************/
    // 1)  Read in each PSD, and load it into the array in FFT order
    /*********************************************************************/

    for( size_t i = 0; i < nModes; ++i )
    {
        // Determine the spatial frequency at this step
        int m = fms[2 * i].m;
        int n = fms[2 * i].n;

        //**< Get the open-loop turb. PSD
        if( getGridPSD( tPSDp, psdDir, m, n ) < 0 )
            return -1;
        tPSDp.erase( tPSDp.begin() + imax, tPSDp.end() ); // Nyquist limit
        modeVar[i] = sigproc::psdVar( tfreq, tPSDp );

        // And now normalize
        sigproc::normPSD( tPSDp, tfreq, 1.0 );                             // Normalize
        sigproc::augment1SidedPSD( psd2sided, tPSDp, !( tfreq[0] == 0 ) ); // Convert to FFT storage order

        opdPSD[i].resize( psd2sided.size() );
        sqrtOPDPSD[i].resize( psd2sided.size() );

        for( size_t j = 0; j < psd2sided.size(); ++j )
        {
            opdPSD[i][j] = psd2sided[j] * modeVar[i];
            sqrtOPDPSD[i][j] = sqrt( psd2sided[j] ); // Store the square-root for later
        }
    }

    size_t sz2Sided = psd2sided.size();

    std::vector<realT> freq2sided;
    freq2sided.resize( sz2Sided );
    sigproc::augment1SidedPSDFreq( freq2sided, tfreq );

    tPSDp.resize( tfreq.size() );
    tETF.resize( tfreq.size() );
    tNTF.resize( tfreq.size() );

    std::vector<std::vector<realT>> sqrtNPSD;
    sqrtNPSD.resize( nModes );

    std::vector<realT> noiseVar;
    noiseVar.resize( nModes );

    std::vector<std::vector<complexT>> ETFsi;
    std::vector<std::vector<complexT>> ETFlp;
    ETFsi.resize( nModes );
    ETFlp.resize( nModes );

    std::vector<std::vector<complexT>> NTFsi;
    std::vector<std::vector<complexT>> NTFlp;
    NTFsi.resize( nModes );
    NTFlp.resize( nModes );

    std::string tfInFile;
    std::string etfInFile;
    std::string ntfInFile;

    improc::eigenImage<realT> Cvd;
    fits::fitsFile<realT> ff;
    ff.read( Cvd, CvdPath );

    std::vector<std::complex<realT>> tPSDc, psd2sidedc;

    /*********************************************************************/
    // 2)  Analyze each star magnitude
    /*********************************************************************/
    ipc::ompLoopWatcher<> watcher( lifetimeTrials * mags.size(), std::cout );
    for( size_t s = 0; s < mags.size(); ++s )
    {
        /*********************************************************************/
        // 2.0)  Read in the transfer functions for each mode
        /*********************************************************************/
        for( size_t i = 0; i < nModes; ++i )
        {
            // Determine the spatial frequency at this step
            int m = fms[2 * i].m;
            int n = fms[2 * i].n;

            //**< Determine if we're inside the hardwarecontrol limit
            bool inside = false;

            if( m_aosys->circularLimit() )
            {
                if( m * m + n * n <= mnCon * mnCon )
                    inside = true;
            }
            else
            {
                if( fabs( m ) <= mnCon && fabs( n ) <= mnCon )
                    inside = true;
            }

            // Get the WFS noise PSD (which is already resized to match tfreq)
            wfsNoisePSD<realT>( tPSDn,
                                m_aosys->beta_p( m, n ),
                                m_aosys->Fg( mags[s] ),
                                tauWFS,
                                m_aosys->npix_wfs( (size_t)0 ),
                                m_aosys->Fbg( (size_t)0 ),
                                m_aosys->ron_wfs( (size_t)0 ) );
            sigproc::augment1SidedPSD( psd2sided, tPSDn, !( tfreq[0] == 0 ) ); // Convert to FFT storage order

            // Pre-calculate the variance of the noise for later use
            noiseVar[i] = sigproc::psdVar( tfreq, tPSDn );

            sqrtNPSD[i].resize( psd2sided.size() );
            for( size_t j = 0; j < psd2sided.size(); ++j )
                sqrtNPSD[i][j] = sqrt( psd2sided[j] );

            ETFsi[i].resize( sz2Sided );
            ETFlp[i].resize( sz2Sided );
            NTFsi[i].resize( sz2Sided );
            NTFlp[i].resize( sz2Sided );

            if( inside )
            {
                tfInFile = std::format( "{}/outputTF_{}_si", dir, mags[s] );
                // dir + "/" + "outputTF_" + ioutils::convert ToString( mags[s] ) + "_si/";

                etfInFile = std::format( "{}/etf_{}_{}.binv", tfInFile, m, n );
                // tfInFile + "etf_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) + ".binv";
                ioutils::readBinVector( tPSDc, etfInFile );
                sigproc::augment1SidedPSD( psd2sidedc, tPSDc, !( tfreq[0] == 0 ), 1 ); // Convert to FFT storage order
                for( size_t j = 0; j < psd2sidedc.size(); ++j )
                    ETFsi[i][j] = psd2sidedc[j];

                ntfInFile = std::format( "{}/ntf_{}_{}.binv", tfInFile, m, n );
                // tfInFile + "ntf_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) + ".binv";
                ioutils::readBinVector( tPSDc, ntfInFile );
                sigproc::augment1SidedPSD( psd2sidedc, tPSDc, !( tfreq[0] == 0 ), 1 ); // Convert to FFT storage order
                for( size_t j = 0; j < psd2sidedc.size(); ++j )
                    NTFsi[i][j] = psd2sidedc[j];

                tfInFile = std::format( "{}/outputTF_{}_lp", dir, mags[s] );
                // dir + "/" + "outputTF_" + ioutils::convert ToString( mags[s] ) + "_lp/";

                etfInFile = std::format( "{}/etf_{}_{}.binv", tfInFile, m, n );
                // tfInFile + "etf_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) + ".binv";
                ioutils::readBinVector( tPSDc, etfInFile );
                sigproc::augment1SidedPSD( psd2sidedc, tPSDc, !( tfreq[0] == 0 ), 1 ); // Convert to FFT storage order
                for( size_t j = 0; j < psd2sidedc.size(); ++j )
                    ETFlp[i][j] = psd2sidedc[j];

                ntfInFile = std::format( "{}/ntf_{}_{}.binv", tfInFile, m, n );
                // tfInFile + "ntf_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) + ".binv";
                ioutils::readBinVector( tPSDc, ntfInFile );
                sigproc::augment1SidedPSD( psd2sidedc, tPSDc, !( tfreq[0] == 0 ), 1 ); // Convert to FFT storage order
                for( size_t j = 0; j < psd2sidedc.size(); ++j )
                    NTFlp[i][j] = psd2sidedc[j];
            }
            else
            {
                for( int q = 0; q < ETFsi.size(); ++q )
                {
                    ETFsi[i][q] = 1;
                    NTFsi[i][q] = 0;
                    ETFlp[i][q] = 1;
                    NTFlp[i][q] = 0;
                }
            }
        }

        sigproc::averagePeriodogram<realT> tavgPgram(
            sz2Sided / 1.,
            1 / fs ); // this is just to get the size right, per-thread instances below
        std::vector<std::vector<realT>> spPSDs;
        spPSDs.resize( nModes );
        for( size_t pp = 0; pp < spPSDs.size(); ++pp )
        {
            spPSDs[pp].resize( tavgPgram.size() );
            for( size_t nn = 0; nn < spPSDs[pp].size(); ++nn )
                spPSDs[pp][nn] = 0;
        }

        std::vector<std::vector<realT>> spPSDslp;
        spPSDslp.resize( nModes );
        for( size_t pp = 0; pp < spPSDslp.size(); ++pp )
        {
            spPSDslp[pp].resize( tavgPgram.size() );
            for( size_t nn = 0; nn < spPSDslp[pp].size(); ++nn )
                spPSDslp[pp][nn] = 0;
        }

#pragma omp parallel
        {
            // Normally distributed random numbers
            math::normDistT<realT> normVar;
            normVar.seed();

            // FFTs for going to Fourier domain and back to time domain.
            math::ft::fftT<realT, std::complex<realT>, 1, 0> fftF( sqrtOPDPSD[0].size() );
            math::ft::fftT<std::complex<realT>, realT, 1, 0> fftB( sqrtOPDPSD[0].size(), math::ft::dir::backward );

            // Working memory
            std::vector<std::complex<realT>> tform1( sqrtOPDPSD[0].size() );
            std::vector<std::complex<realT>> tform2( sqrtOPDPSD[0].size() );
            std::vector<std::complex<realT>> Ntform1( sqrtOPDPSD[0].size() );
            std::vector<std::complex<realT>> Ntform2( sqrtOPDPSD[0].size() );

            std::vector<std::complex<realT>> tform1lp( sqrtOPDPSD[0].size() );
            std::vector<std::complex<realT>> tform2lp( sqrtOPDPSD[0].size() );
            std::vector<std::complex<realT>> Ntform1lp( sqrtOPDPSD[0].size() );
            std::vector<std::complex<realT>> Ntform2lp( sqrtOPDPSD[0].size() );

            // OPD-PSD filter
            sigproc::psdFilter<realT, 1> pfilt;
            pfilt.psdSqrt( &sqrtOPDPSD[0], tfreq[1] - tfreq[0] ); // Pre-configure

            // Noise-PSD filter
            sigproc::psdFilter<realT, 1> nfilt;
            nfilt.psdSqrt( &sqrtNPSD[0], tfreq[1] - tfreq[0] ); // Pre-configure

            // The h time-series
            std::vector<std::vector<realT>> hts;
            hts.resize( 2 * nModes );

            // The correlated h time-series
            std::vector<std::vector<realT>> htsCorr;
            htsCorr.resize( 2 * nModes );

            for( size_t pp = 0; pp < hts.size(); ++pp )
            {
                hts[pp].resize( sqrtOPDPSD[0].size() );
                htsCorr[pp].resize( sqrtOPDPSD[0].size() );
            }

            // The noise time-serieses
            std::vector<realT> N_n;
            N_n.resize( sz2Sided );

            std::vector<realT> N_nm;
            N_nm.resize( sz2Sided );

            // Periodogram averager
            sigproc::averagePeriodogram<realT> avgPgram( sz2Sided / 1., 1 / fs ); //, 0, 1);
            avgPgram.win( sigproc::window::hann );

            // The periodogram output
            std::vector<realT> tpgram( avgPgram.size() );

            // Holds the speckle time-series
            improc::eigenCube<realT> spTS;
            spTS.resize( 2 * mnMax + 1, 2 * mnMax + 1, tform1.size() );

            improc::eigenCube<realT> spTSlp;
            spTSlp.resize( 2 * mnMax + 1, 2 * mnMax + 1, tform1.size() );

// Here's where the big loop of n-trials should start
#pragma omp for
            for( int zz = 0; zz < lifetimeTrials; ++zz )
            {
                std::complex<realT> scale = exp( std::complex<realT>( 0, math::half_pi<realT>() ) ) /
                                            std::complex<realT>( ( tform1.size() ), 0 );

                /*********************************************************************/
                // 2.1)  Generate filtered noise for each mode, with temporal phase shifts at each spatial frequency
                /*********************************************************************/
                for( size_t pp = 0; pp < nModes; ++pp )
                {
                    // Fill in standard normal noise
                    for( size_t nn = 0; nn < hts[2 * pp].size(); ++nn )
                    {
                        hts[2 * pp][nn] = normVar;
                    }

                    // Set sqrt(PSD), just a pointer switch
                    pfilt.psdSqrt( &sqrtOPDPSD[pp], tfreq[1] - tfreq[0] );

                    // And now filter the noise to a time-series of h
                    pfilt( hts[2 * pp] );

                    /**/
                    // Then construct 2nd mode with temporal shift
                    fftF( tform1.data(), hts[2 * pp].data() );

                    // Apply the phase shift to form the 2nd time series
                    for( size_t nn = 0; nn < hts[2 * pp].size(); ++nn )
                        tform1[nn] = tform1[nn] * scale;

                    fftB( hts[2 * pp + 1].data(), tform1.data() );
                    /**/
                }
                //** At this point we have correlated time-series, with the correct temporal PSD, but not yet spatially
                // correlated

                /*********************************************************************/
                // 2.2)  Correlate the time-series for each mode
                /*********************************************************************/
                // #pragma omp parallel for
                for( size_t pp = 0; pp < hts.size(); ++pp )
                {
                    for( size_t nn = 0; nn < hts[0].size(); ++nn )
                    {
                        htsCorr[pp][nn] = 0;
                    }

                    for( size_t qq = 0; qq <= pp; ++qq )
                    {
                        realT cvd = Cvd( qq, pp );
                        realT *d1 = htsCorr[pp].data();
                        realT *d2 = hts[qq].data();
                        for( size_t nn = 0; nn < hts[0].size(); ++nn )
                        {
                            d1[nn] += d2[nn] * cvd;
                        }
                    }
                }

                /*
                for(size_t pp=0; pp < hts.size(); ++pp)
                {
                   for(size_t nn=0; nn< hts[0].size(); ++nn)
                   {
                      htsCorr[pp][nn] = hts[pp][nn];
                   }
                }*/

                /*********************************************************************/
                // 2.2.a)  Re-normalize b/c the correlation step above does not result in correct variances
                ///\todo should be able to scale the covar by r0, and possibly D
                /*********************************************************************/
                for( size_t pp = 0; pp < nModes; ++pp )
                {
                    math::vectorMeanSub( htsCorr[2 * pp] );
                    math::vectorMeanSub( htsCorr[2 * pp + 1] );

                    realT var = math::vectorVariance( htsCorr[2 * pp] );
                    realT norm = sqrt( modeVar[pp] / var );
                    for( size_t nn = 0; nn < htsCorr[2 * pp].size(); ++nn )
                        htsCorr[2 * pp][nn] *= norm;

                    var = math::vectorVariance( htsCorr[2 * pp + 1] );
                    norm = sqrt( modeVar[pp] / var );
                    for( size_t nn = 0; nn < htsCorr[2 * pp + 1].size(); ++nn )
                        htsCorr[2 * pp + 1][nn] *= norm;
                }

                scale = std::complex<realT>( tform1.size(), 0 );

                /*********************************************************************/
                // 2.3)  Generate speckle intensity time-series
                /*********************************************************************/
                for( size_t pp = 0; pp < nModes; ++pp )
                {
                    // Now we take them back to the FD and apply the xfer
                    // and add the noise

                    fftF( tform1.data(), htsCorr[2 * pp].data() );
                    fftF( tform2.data(), htsCorr[2 * pp + 1].data() );

                    // Make some noise
                    for( int nn = 0; nn < sz2Sided; ++nn )
                    {
                        N_n[nn] = normVar;
                        N_nm[nn] = normVar;
                    }

                    // Filter it
                    // Set sqrt(PSD), just a pointer switch
                    pfilt.psdSqrt( &sqrtNPSD[pp], tfreq[1] - tfreq[0] );
                    nfilt.filter( N_n );
                    nfilt.filter( N_nm );

                    // Normalize it
                    realT Nactvar = 0.5 * ( math::vectorVariance( N_n ) + math::vectorVariance( N_nm ) );
                    realT norm = sqrt( noiseVar[pp] / Nactvar );
                    for( size_t q = 0; q < N_n.size(); ++q )
                        N_n[q] *= norm;
                    for( size_t q = 0; q < N_nm.size(); ++q )
                        N_nm[q] *= norm;

                    // And move them to the Fourier domain
                    fftF( Ntform1.data(), N_n.data() );
                    fftF( Ntform2.data(), N_nm.data() );

                    // Apply the closed loop transfers
                    for( size_t mm = 0; mm < tform1.size(); ++mm )
                    {
                        // Apply the augmented ETF to two time-series
                        tform1lp[mm] = tform1[mm] * ETFlp[pp][mm] / scale;
                        tform2lp[mm] = tform2[mm] * ETFlp[pp][mm] / scale;

                        Ntform1lp[mm] = Ntform1[mm] * NTFlp[pp][mm] / scale;
                        Ntform2lp[mm] = Ntform2[mm] * NTFlp[pp][mm] / scale;

                        tform1[mm] *= ETFsi[pp][mm] / scale;
                        tform2[mm] *= ETFsi[pp][mm] / scale;

                        Ntform1[mm] *= NTFsi[pp][mm] / scale;
                        Ntform2[mm] *= NTFsi[pp][mm] / scale;
                    }

                    // And make the speckle TS
                    int m = fms[2 * pp].m;
                    int n = fms[2 * pp].n;

                    //<<<<<<<<****** Transform back to the time domain.
                    fftB( htsCorr[2 * pp].data(), tform1.data() );
                    fftB( htsCorr[2 * pp + 1].data(), tform2.data() );
                    fftB( N_n.data(), Ntform1.data() );
                    fftB( N_nm.data(), Ntform2.data() );

                    for( int i = 0; i < hts[2 * pp].size(); ++i )
                    {
                        realT h1 = htsCorr[2 * pp][i] + N_n[i];
                        realT h2 = htsCorr[2 * pp + 1][i] + N_nm[i];

                        spTS.image( i )( mnMax + m, mnMax + n ) = ( pow( h1, 2 ) + pow( h2, 2 ) );
                        spTS.image( i )( mnMax - m, mnMax - n ) = spTS.image( i )( mnMax + m, mnMax + n );
                    }

                    fftB( htsCorr[2 * pp].data(), tform1lp.data() );
                    fftB( htsCorr[2 * pp + 1].data(), tform2lp.data() );
                    fftB( N_n.data(), Ntform1lp.data() );
                    fftB( N_nm.data(), Ntform2lp.data() );

                    for( int i = 0; i < hts[2 * pp].size(); ++i )
                    {
                        realT h1 = htsCorr[2 * pp][i] + N_n[i];
                        realT h2 = htsCorr[2 * pp + 1][i] + N_nm[i];

                        spTSlp.image( i )( mnMax + m, mnMax + n ) = ( pow( h1, 2 ) + pow( h2, 2 ) );
                        spTSlp.image( i )( mnMax - m, mnMax - n ) = spTSlp.image( i )( mnMax + m, mnMax + n );
                    }
                }

                /*********************************************************************/
                // 2.5)  Calculate speckle PSD for each mode
                /*********************************************************************/
                std::vector<realT> speckAmp( spTS.planes() );
                std::vector<realT> speckAmplp( spTS.planes() );

                for( size_t pp = 0; pp < nModes; ++pp )
                {
                    int m = fms[2 * pp].m;
                    int n = fms[2 * pp].n;

                    realT mn = 0;
                    realT mnlp = 0;
                    for( int i = 0; i < spTS.planes(); ++i )
                    {
                        speckAmp[i] = spTS.image( i )( mnMax + m, mnMax + n );
                        speckAmplp[i] = spTSlp.image( i )( mnMax + m, mnMax + n );

                        mn += speckAmp[i];
                        mnlp += speckAmplp[i];
                    }
                    mn /= speckAmp.size();
                    mnlp /= speckAmplp.size();

                    // mean subtract
                    for( int i = 0; i < speckAmp.size(); ++i )
                        speckAmp[i] -= mn;
                    for( int i = 0; i < speckAmplp.size(); ++i )
                        speckAmplp[i] -= mnlp;

                    // Calculate PSD of the speckle amplitude
                    avgPgram( tpgram, speckAmp );
                    for( size_t nn = 0; nn < spPSDs[pp].size(); ++nn )
                        spPSDs[pp][nn] += tpgram[nn];

                    avgPgram( tpgram, speckAmplp );
                    for( size_t nn = 0; nn < spPSDslp[pp].size(); ++nn )
                        spPSDslp[pp][nn] += tpgram[nn];
                }

                watcher.incrementAndOutputStatus();
            } // for(int zz=0; zz<lifetimeTrials; ++zz)
        } // #pragma omp parallel

        std::vector<realT> spFreq( spPSDs[0].size() );
        for( size_t nn = 0; nn < spFreq.size(); ++nn )
            spFreq[nn] = tavgPgram[nn];

        improc::eigenImage<realT> taus, tauslp;
        taus.resize( 2 * mnMax + 1, 2 * mnMax + 1 );
        tauslp.resize( 2 * mnMax + 1, 2 * mnMax + 1 );

        improc::eigenCube<realT> imc, imclp;
        std::vector<realT> clPSD;

        if( writePSDs )
        {
            imc.resize( 2 * mnMax + 1, 2 * mnMax + 1, spPSDs[0].size() );
            imclp.resize( 2 * mnMax + 1, 2 * mnMax + 1, spPSDs[0].size() );
            clPSD.resize( sz2Sided );
        }

        sigproc::psdVarMean<sigproc::psdVarMeanParams<realT>> pvm;
        /*********************************************************************/
        // 3.0)  Calculate lifetimes from the PSDs
        /*********************************************************************/
        for( size_t pp = 0; pp < nModes; ++pp )
        {
            spPSDs[pp][0] = spPSDs[pp][1];     // deal with under-estimated mean.
            spPSDslp[pp][0] = spPSDslp[pp][1]; // deal with under-estimated mean.

            int m = fms[2 * pp].m;
            int n = fms[2 * pp].n;

            realT var;
            if( writePSDs ) // Have to normalize the intensity for some reason if we want to use the PSDs
            {
                for( size_t nn = 0; nn < spPSDs[pp].size(); ++nn )
                {
                    spPSDs[pp][nn] /= lifetimeTrials;
                }

                for( size_t nn = 0; nn < sz2Sided; ++nn )
                {
                    clPSD[nn] =
                        opdPSD[pp][nn] * norm( ETFsi[pp][nn] ) + pow( sqrtNPSD[pp][nn], 2 ) * norm( NTFsi[pp][nn] );
                }

                var = sigproc::psdVar( freq2sided, clPSD );

                realT pvar = sigproc::psdVar( spFreq, spPSDs[pp] );

                for( size_t nn = 0; nn < spPSDs[pp].size(); ++nn )
                {
                    spPSDs[pp][nn] *= var / pvar;
                    imc.image( nn )( mnMax + m, mnMax + n ) = spPSDs[pp][nn];
                    imc.image( nn )( mnMax - m, mnMax - n ) = spPSDs[pp][nn];
                }

                // lp
                for( size_t nn = 0; nn < spPSDslp[pp].size(); ++nn )
                {
                    spPSDslp[pp][nn] /= lifetimeTrials;
                }

                for( size_t nn = 0; nn < sz2Sided; ++nn )
                {
                    clPSD[nn] =
                        opdPSD[pp][nn] * norm( ETFlp[pp][nn] ) + pow( sqrtNPSD[pp][nn], 2 ) * norm( NTFlp[pp][nn] );
                }

                var = sigproc::psdVar( freq2sided, clPSD );

                pvar = sigproc::psdVar( spFreq, spPSDslp[pp] );

                for( size_t nn = 0; nn < spPSDslp[pp].size(); ++nn )
                {
                    spPSDslp[pp][nn] *= var / pvar;
                    imclp.image( nn )( mnMax + m, mnMax + n ) = spPSDslp[pp][nn];
                    imclp.image( nn )( mnMax - m, mnMax - n ) = spPSDslp[pp][nn];
                }
            }

            var = sigproc::psdVar( spFreq, spPSDs[pp] );

            realT T = ( 1.0 / ( spFreq[1] - spFreq[0] ) ) * 10;
            realT error;
            realT tau = pvm( error, spFreq, spPSDs[pp], T ) * ( T ) / var;
            taus( mnMax + m, mnMax + n ) = tau;
            taus( mnMax - m, mnMax - n ) = tau;

            var = sigproc::psdVar( spFreq, spPSDslp[pp] );

            tau = pvm( error, spFreq, spPSDslp[pp], T ) * ( T ) / var;
            tauslp( mnMax + m, mnMax + n ) = tau;
            tauslp( mnMax - m, mnMax - n ) = tau;
        }

        /*********************************************************************/
        // 4.0)  Write the results to disk
        /*********************************************************************/
        fn = std::format( "{}/speckleLifetimes_{}_si.fits", dir, mags[s] );
        // dir + "/speckleLifetimes_" + ioutils::convert ToString( mags[s] ) + "_si.fits";
        ff.write( fn, taus );

        fn = std::format( "{}/speckleLifetimes_{}_lp.fits", dir, mags[s] );
        // dir + "/speckleLifetimes_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
        ff.write( fn, tauslp );

        if( writePSDs )
        {
            fn = std::format( "{}/specklePSDs_{}_si.fits", dir, mags[s] );
            // dir + "/specklePSDs_" + ioutils::convert ToString( mags[s] ) + "_si.fits";
            ff.write( fn, imc );

            fn = std::format( "{}/speckleLifetimes_{}_lp.fits", dir, mags[s] );
            // dir + "/specklePSDs_" + ioutils::convert ToString( mags[s] ) + "_lp.fits";
            ff.write( fn, imclp );
        }

    } // mags

    return 0;
}

template <typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::getGridFreq( std::vector<realT> &freq, const std::string &dir )
{
    std::string fn;
    fn = dir + "/psds/freq.binv";
    return ioutils::readBinVector( freq, fn );
}

template <typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::getGridPSD( std::vector<realT> &psd, const std::string &dir, int m, int n )
{
    std::string fn;
    fn = std::format( "{}/psds/psd_{}_{}.binv", dir, m, n );
    // dir + "/psds/psd_" + ioutils::convert ToString( m ) + '_' + ioutils::convert ToString( n ) + ".binv";
    return ioutils::readBinVector( psd, fn );
}

template <typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::getGridPSD(
    std::vector<realT> &freq, std::vector<realT> &psd, const std::string &dir, int m, int n )
{
    int rv = getGridFreq( freq, dir );
    if( rv < 0 )
        return rv;

    return getGridPSD( psd, dir, m, n );
}

/// Worker function for GSL Integration for the basic sin/cos Fourier modes.
/** \ingroup mxAOAnalytic
 */
template <typename realT, typename aosysT>
realT F_basic( realT kv, void *params )
{
    fourierTemporalPSD<realT, aosysT> *Fp = (fourierTemporalPSD<realT, aosysT> *)params;

    realT f = Fp->m_f;
    realT v_wind = Fp->m_aosys->atm.layer_v_wind( Fp->_layer_i );

    realT D = Fp->m_aosys->D();
    realT m = Fp->m_m;
    realT n = Fp->m_n;
    int p = Fp->m_p;

    realT ku = f / v_wind;

    realT kp = sqrt( pow( ku + m / D, 2 ) + pow( kv + n / D, 2 ) );
    realT kpp = sqrt( pow( ku - m / D, 2 ) + pow( kv - n / D, 2 ) );

    realT Q1 = math::func::jinc( math::pi<realT>() * D * kp );

    realT Q2 = math::func::jinc( math::pi<realT>() * D * kpp );

    realT Q = ( Q1 + p * Q2 );

    realT P =
        Fp->m_aosys->psd( Fp->m_aosys->atm, Fp->_layer_i, sqrt( pow( ku, 2 ) + pow( kv, 2 ) ), Fp->m_aosys->secZeta() );

    return P * Q * Q;
}

template <typename realT>
void turbBoilCubic(
    realT &a, realT &b, realT &c, realT &d, const realT &kv, const realT &f, const realT &Vu, const realT &f0, int pm )
{
    a = Vu * Vu * Vu;
    b = -( 3 * Vu * Vu * f + pm * f0 * f0 * f0 );
    c = 3 * f * f * Vu;
    d = -( f * f * f + pm * f0 * f0 * f0 * kv * kv );
}

/// Worker function for GSL Integration for the modified Fourier modes.
/** \ingroup mxAOAnalytic
 */
template <typename realT, typename aosysT>
realT F_mod( realT kv, void *params )
{
    fourierTemporalPSD<realT, aosysT> *Fp = (fourierTemporalPSD<realT, aosysT> *)params;

    realT f = Fp->m_f;
    realT v_wind = Fp->m_aosys->atm.layer_v_wind( Fp->_layer_i );

    realT D = Fp->m_aosys->D();
    realT m = Fp->m_m;
    realT n = Fp->m_n;

    realT f0 = Fp->m_f0;

    realT ku;
    if( f0 == 0 )
    {
        ku = f / v_wind;

        if( Fp->m_spatialFilter )
        {
            // de-rotate the spatial frequency vector back to pupil coordinates
            realT dku = ku * Fp->m_cq - kv * Fp->m_sq;
            realT dkv = ku * Fp->m_sq + kv * Fp->m_cq;
            // Return if spatially filtered
            if( fabs( dku ) >= Fp->m_aosys->spatialFilter_ku() )
                return 0;

            if( fabs( dkv ) >= Fp->m_aosys->spatialFilter_kv() )
                return 0;
        }

        realT kp = sqrt( pow( ku + m / D, 2 ) + pow( kv + n / D, 2 ) );
        realT kpp = sqrt( pow( ku - m / D, 2 ) + pow( kv - n / D, 2 ) );

        realT Jp = math::func::jinc( math::pi<realT>() * D * kp );

        realT Jm = math::func::jinc( math::pi<realT>() * D * kpp );

        realT QQ = 2 * ( Jp * Jp + Jm * Jm );

        realT P = Fp->m_aosys->psd( Fp->m_aosys->atm,
                                    Fp->_layer_i,
                                    sqrt( pow( ku, 2 ) + pow( kv, 2 ) ),
                                    Fp->m_aosys->lam_sci(),
                                    Fp->m_aosys->lam_wfs(),
                                    Fp->m_aosys->secZeta() );

        return P * QQ;
    }
    else
    {
        realT a, b, c, d, p, q;

        turbBoilCubic( a, b, c, d, kv, f, v_wind, f0, 1 );
        mx::math::cubicDepressed( p, q, a, b, c, d );
        realT t = mx::math::cubicRealRoot( p, q );

        ku = t - b / ( 3 * a );

        if( Fp->m_spatialFilter )
        {
            // de-rotate the spatial frequency vector back to pupil coordinates
            realT dku = ku * Fp->m_cq - kv * Fp->m_sq;
            realT dkv = ku * Fp->m_sq + kv * Fp->m_cq;
            // Return if spatially filtered
            if( fabs( dku ) >= Fp->m_aosys->spatialFilter_ku() )
                return 0;

            if( fabs( dkv ) >= Fp->m_aosys->spatialFilter_kv() )
                return 0;
        }

        realT kp = sqrt( pow( ku + m / D, 2 ) + pow( kv + n / D, 2 ) );
        realT kpp = sqrt( pow( ku - m / D, 2 ) + pow( kv - n / D, 2 ) );

        realT Jp = math::func::jinc( math::pi<realT>() * D * kp );

        realT Jm = math::func::jinc( math::pi<realT>() * D * kpp );

        realT QQ = 2 * ( Jp * Jp + Jm * Jm );

        realT P1 = Fp->m_aosys->psd( Fp->m_aosys->atm,
                                     Fp->_layer_i,
                                     sqrt( pow( ku, 2 ) + pow( kv, 2 ) ),
                                     Fp->m_aosys->lam_sci(),
                                     Fp->m_aosys->lam_wfs(),
                                     Fp->m_aosys->secZeta() );

        P1 *= QQ;

        turbBoilCubic( a, b, c, d, kv, f, v_wind, f0, -1 );
        mx::math::cubicDepressed( p, q, a, b, c, d );
        t = mx::math::cubicRealRoot( p, q );

        ku = t - b / ( 3 * a );

        if( Fp->m_spatialFilter )
        {
            // de-rotate the spatial frequency vector back to pupil coordinates
            realT dku = ku * Fp->m_cq - kv * Fp->m_sq;
            realT dkv = ku * Fp->m_sq + kv * Fp->m_cq;
            // Return if spatially filtered
            if( fabs( dku ) >= Fp->m_aosys->spatialFilter_ku() )
                return 0;

            if( fabs( dkv ) >= Fp->m_aosys->spatialFilter_kv() )
                return 0;
        }

        kp = sqrt( pow( ku + m / D, 2 ) + pow( kv + n / D, 2 ) );
        kpp = sqrt( pow( ku - m / D, 2 ) + pow( kv - n / D, 2 ) );

        Jp = math::func::jinc( math::pi<realT>() * D * kp );

        Jm = math::func::jinc( math::pi<realT>() * D * kpp );

        QQ = 2 * ( Jp * Jp + Jm * Jm );

        realT P2 = Fp->m_aosys->psd( Fp->m_aosys->atm,
                                     Fp->_layer_i,
                                     sqrt( pow( ku, 2 ) + pow( kv, 2 ) ),
                                     Fp->m_aosys->lam_sci(),
                                     Fp->m_aosys->lam_wfs(),
                                     Fp->m_aosys->secZeta() );

        P2 *= QQ;

        return 0.5 * ( P1 + P2 );
    }
}

/*extern template
struct fourierTemporalPSD<float, aoSystem<float, vonKarmanSpectrum<float>, std::ostream>>;*/

extern template struct fourierTemporalPSD<double, aoSystem<double, vonKarmanSpectrum<double>, std::ostream>>;

/*
extern template
struct fourierTemporalPSD<long double, aoSystem<long double, vonKarmanSpectrum<long double>, std::ostream>>;

#ifdef HASQUAD
extern template
struct fourierTemporalPSD<__float128, aoSystem<__float128, vonKarmanSpectrum<__float128>, std::ostream>>;
#endif
*/

} // namespace analysis
} // namespace AO
} // namespace mx

#endif // fourierTemporalPSD_hpp
