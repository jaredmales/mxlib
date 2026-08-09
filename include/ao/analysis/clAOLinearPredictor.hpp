/** \file clAOLinearPredictor.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Provides a class to manage closed loop gain linear predictor determination.
 * \ingroup mxAO_files
 *
 */

#ifndef clAOLinearPredictor_hpp
#define clAOLinearPredictor_hpp

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "../../mxlib.hpp"

#include "../../math/floatUtils.hpp"
#include "../../math/geo.hpp"

#include "../../sigproc/psdUtils.hpp"

#include "../../sigproc/autocorrelation.hpp"
#include "../../sigproc/linearPredictor.hpp"

#include "clGainOpt.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

#define CLAOLP_BREADCRUMB

// #define CLAOLP_BREADCRUMB std::cerr << __FILE__ << ' ' << __LINE__ << '\n';

/// Class to manage the calculation of linear predictor coefficients for a closed-loop AO system.
/**
 * \tparam _realT the real floating point type in which to do all arithmetic.
 *
 * \ingroup mxAOAnalytic
 */
template <typename _realT>
struct clAOLinearPredictor
{
    typedef _realT realT; ///< Floating-point type used for predictor calculations.

  public:
    /// Result from one evaluated regularization scale.
    struct regResult
    {
        realT sc;   ///< Regularization scale in dB.
        realT gopt; ///< Optimum gain at this scale.
        realT gmax; ///< Maximum stable gain at this scale.
        realT var;  ///< Closed-loop variance at the optimum gain.
    };

    /// Termination state of the most recent regularization search.
    enum class regularizationStatus
    {
        notRun,             ///< No regularization search has been attempted.
        converged,          ///< The requested precision was reached.
        boundaryLimited,    ///< The optimum remained on the expanded search boundary.
        invalidControls,    ///< The configured search controls were invalid.
        iterationLimit,     ///< The search exhausted its iteration limit.
        calculationFailure, ///< Coefficient or gain calculation failed.
    };

    /// Diagnostic summary of the most recent regularization search.
    struct regularizationReport
    {
        regularizationStatus status{ regularizationStatus::notRun }; ///< Search termination state.
        int iterations{ 0 };                                         ///< Refinement iterations attempted.
        std::size_t evaluations{ 0 };                                ///< Regularization scales evaluated.
    };

    std::vector<realT> m_PSDtn;                     ///< Working memory for the regularized PSD

    std::vector<realT> m_psd2s;                     ///< Working memory for the 2-sided regularized PSD

    std::vector<realT> m_ac;                        ///< Working memory to hold the autocorrelation.

    sigproc::autocorrelationFromPSD<realT> m_acpsd; ///< Converts the working PSD to an autocorrelation.

    sigproc::linearPredictor<realT> m_lp;           ///< Linear predictor used to calculate coefficients.

    realT m_min_var0{ 0 };                          ///< Initial minimum variance, with zero requesting initialization.
    realT m_min_sc0{ 10 };                          ///< Initial minimum regularization scale in dB.
    realT m_precision0{ 2 };                        ///< Initial regularization scale spacing in dB.
    realT m_max_sc0{ 100 };                         ///< Initial maximum regularization scale in dB.
    realT m_dPrecision{ 3 };                        ///< Divisor applied to the spacing during refinement.

    realT m_gmax_lp{ 5 };                           ///< The maximum allowable gain for LP.

    // Stopping conditions:
    realT m_minPrecision{ 0.001 };               ///< Minimum requested regularization spacing in dB.
    int m_maxIts{ 100 };                         ///< Maximum number of search refinement iterations.

    int m_extrap{ 1 };                           ///< The LP extrapolation length in loop steps. Normally it is 1 step.

    std::vector<regResult> m_regResults;         ///< Per-scale telemetry collected when requested.

    regularizationReport m_regularizationReport; ///< Diagnostic summary of the latest search.

  public:
    /// Construct a closed-loop linear-predictor calculator with default search controls.
    clAOLinearPredictor() = default;

    /// Calculate the LP coefficients for a turbulence PSD and a noise PSD.
    /** This combines the two PSDs, augments to two-sided, and calls the linearPredictor.calcCoefficients method.
     *
     * A regularization constant can be added to the PSD as well.
     *
     * \returns `error_t::noerror` on success, otherwise `error_t::liberr`.
     */
    mx::error_t calcCoefficients( std::vector<realT> &PSDt, /**< [in] the turbulence PSD */
                                  std::vector<realT> &PSDn, /**< [in] the WFS noise PSD */
                                  realT PSDreg,             /**< [in] the regularizing constant. Set to 0 to not use. */
                                  int Nc,                   /**< [in] the number of LP coefficients */
                                  realT condition = 0       /**< [in] the condition number for the SVD.  If 0 then
                                                                      levinson recursion is used. */
    )
    {
        CLAOLP_BREADCRUMB;
        m_PSDtn.resize( PSDt.size() );

        CLAOLP_BREADCRUMB;
        for( size_t i = 0; i < PSDt.size(); ++i )
        {
            m_PSDtn[i] = PSDt[i] + PSDn[i] + PSDreg;
        }

        CLAOLP_BREADCRUMB;
        sigproc::augment1SidedPSD( m_psd2s, m_PSDtn, 1 );

        CLAOLP_BREADCRUMB;
        m_ac.resize( m_psd2s.size() );

        CLAOLP_BREADCRUMB;
        m_acpsd( m_ac, m_psd2s );
        CLAOLP_BREADCRUMB;

        if( m_lp.calcCoefficients( m_ac, Nc, m_extrap, condition ) != 0 )
        {
            return internal::mxlib_error_report( error_t::liberr, "linearPredictor::calcCoefficients failed" );
        }

        return error_t::noerror;
    }

    /// Worker function for regularizing the PSD for coefficient calculation.
    /**
     * \tparam telem if true then the results are collected in m_regResults.
     *
     * On first call (min_var = 0):
     *     loop over scale factors from min_sc to max_sc (<=) in steps of precision.
     *
     * On subsequent calls, when min_var and min_sc are passed back in
     *     loop over scale factors from min_sc-precision to max_sc in steps of
     *
     * \returns `error_t::noerror` on success, otherwise the coefficient-calculation error.
     */
    template <bool telem>
    mx::error_t
    _regularizeCoefficients( realT &min_var,  /**< [in,out] the minimum variance found; set to 0 on initial call */
                             realT &min_sc,   /**< [in,out] the scale factor at the minimum variance */
                             realT precision, /**< [in] the step size for the scale factor */
                             realT max_sc,    /**< [in] the maximum scale factor to test */
                             clGainOpt<realT> &go_lp,  /**< [in] the gain optimization object */
                             std::vector<realT> &PSDt, /**< [in] the turbulence PSD */
                             std::vector<realT> &PSDn, /**< [in] the WFS noise PSD */
                             int Nc                    /**< [in] the number of coefficients */
    )
    {
        CLAOLP_BREADCRUMB;

        realT gmax_lp;
        realT gopt_lp;
        realT var_lp;

        realT sc0;

        if( min_var == 0 ) // first call
        {
            sc0 = min_sc;
            min_var = std::numeric_limits<realT>::max();
        }
        else
        {
            sc0 = min_sc - precision * m_dPrecision;
        }

        CLAOLP_BREADCRUMB;

        // auto it = std::max_element(std::begin(PSDt), std::end(PSDt));
        realT psdReg = PSDt[0]; //*it/10;

        CLAOLP_BREADCRUMB;
        // Test from sc0 to max_sc in steps of precision
        // for( realT sc = sc0; sc <= max_sc; sc += precision )
        for( realT sc = max_sc; sc >= sc0; sc -= precision )
        {
            CLAOLP_BREADCRUMB;
            ++m_regularizationReport.evaluations;
            error_t rv = calcCoefficients( PSDt, PSDn, psdReg * pow( 10, -sc / 10 ), Nc );
            if( rv != error_t::noerror )
            {
                m_regularizationReport.status = regularizationStatus::calculationFailure;
                return rv;
            }

            CLAOLP_BREADCRUMB;

            CLAOLP_BREADCRUMB;
            go_lp.a( m_lp.m_c );
            go_lp.b( m_lp.m_c );

            CLAOLP_BREADCRUMB;
            realT ll = 0, ul = 0;
            gmax_lp = go_lp.maxStableGain( ll, ul );
            if( gmax_lp > m_gmax_lp )
            {
                gmax_lp = m_gmax_lp;
            }

            CLAOLP_BREADCRUMB;
            gopt_lp = go_lp.optGainOpenLoop( var_lp, PSDt, PSDn, gmax_lp, false );

            if( telem )
            {
                m_regResults.push_back( { sc, gopt_lp, gmax_lp, var_lp } );
            }

            CLAOLP_BREADCRUMB;
            if( var_lp < min_var )
            {
                min_var = var_lp;
                min_sc = sc;
            }

            // A jump by a factor of 10 indicates the wall
            if( var_lp > 10 * min_var )
            {
                return error_t::noerror;
            }

            CLAOLP_BREADCRUMB;
        }

        CLAOLP_BREADCRUMB;
        return error_t::noerror;
    }

    /// Regularize the PSD and calculate the associated LP coefficients.
    /** The PSD is regularized by adding a constant to it.  This constant is found by minimizing the variance of the
     * residual PSD.
     *
     * \tparam telem if true then the results are collected in m_regResults
     *
     * \returns `error_t::noerror` for a converged or boundary-limited search, `error_t::invalidconfig` for invalid
     * controls, `error_t::timeout` on iteration exhaustion, or the coefficient-calculation error.
     */
    template <bool telem = false>
    mx::error_t
    regularizeCoefficients( realT &gmax_lp,           /**< [out] the maximum gain calculated for the regularized PSD */
                            realT &gopt_lp,           /**< [out] the optimum gain calculated for the regularized PSD */
                            realT &var_lp,            /**< [out] the variance at the optimum gain */
                            realT &min_sc,            /**< [out] the optimum regularization scale factor */
                            clGainOpt<realT> &go_lp,  /**< [in] the gain optimization object */
                            std::vector<realT> &PSDt, /**< [in] the turbulence PSD */
                            std::vector<realT> &PSDn, /**< [in] the WFS noise PSD */
                            int Nc                    /**< [in] the number of coefficients */
    )
    {
        CLAOLP_BREADCRUMB;

        m_regularizationReport = {};

        const realT intervalWidth = m_max_sc0 - m_min_sc0;
        if( !math::isFinite( m_min_sc0 ) || !math::isFinite( m_max_sc0 ) || !math::isFinite( m_precision0 ) ||
            !math::isFinite( m_minPrecision ) || !math::isFinite( m_dPrecision ) || m_minPrecision <= 0 ||
            m_precision0 <= m_minPrecision || intervalWidth <= 0 || m_precision0 > intervalWidth || m_dPrecision <= 1 ||
            m_maxIts <= 0 )
        {
            m_regularizationReport.status = regularizationStatus::invalidControls;
            return internal::mxlib_error_report( error_t::invalidconfig,
                                                 "invalid linear-predictor regularization search controls" );
        }

        realT min_var = m_min_var0;
        min_sc = m_min_sc0;
        realT precision = m_precision0;
        realT max_sc = m_max_sc0;

        if( telem )
        {
            m_regResults.reserve( m_maxIts * 50 );
        }

        CLAOLP_BREADCRUMB;
        int its = 0;
        while( precision > m_minPrecision && its < m_maxIts )
        {
            CLAOLP_BREADCRUMB;
            const bool firstIteration = its == 0;
            error_t rv = _regularizeCoefficients<telem>( min_var, min_sc, precision, max_sc, go_lp, PSDt, PSDn, Nc );
            ++its;
            m_regularizationReport.iterations = its;
            if( rv != error_t::noerror )
            {
                return rv;
            }

            CLAOLP_BREADCRUMB;
            if( min_sc == max_sc )
            {
                if( firstIteration )
                {
                    min_sc -= precision;
                    max_sc = 200;
                }
                else
                {
                    m_regularizationReport.status = regularizationStatus::boundaryLimited;
                    break;
                }
            }
            else
            {
                max_sc = min_sc + precision;
                precision /= m_dPrecision;
            }
        }

        if( precision > m_minPrecision && its >= m_maxIts &&
            m_regularizationReport.status != regularizationStatus::boundaryLimited )
        {
            m_regularizationReport.status = regularizationStatus::iterationLimit;
            return internal::mxlib_error_report( error_t::timeout,
                                                 "linear-predictor regularization reached its iteration limit" );
        }

        if( m_regularizationReport.status != regularizationStatus::boundaryLimited )
        {
            m_regularizationReport.status = regularizationStatus::converged;
        }

        CLAOLP_BREADCRUMB;
        // Now record final values
        error_t rv = calcCoefficients( PSDt, PSDn, PSDt[0] * pow( 10, -min_sc / 10 ), Nc );
        if( rv != error_t::noerror )
        {
            m_regularizationReport.status = regularizationStatus::calculationFailure;
            return rv;
        }

        CLAOLP_BREADCRUMB;
        go_lp.a( m_lp.m_c );
        go_lp.b( m_lp.m_c );

        CLAOLP_BREADCRUMB;
        realT ll = 0, ul = 0;
        gmax_lp = go_lp.maxStableGain( ll, ul );
        gopt_lp = go_lp.optGainOpenLoop( var_lp, PSDt, PSDn, gmax_lp, false );

        CLAOLP_BREADCRUMB;
        return error_t::noerror;
    }

    /// Regularize the PSD and calculate the associated LP coefficients.
    /** The PSD is regularized by adding a constant to it.  This constant is found by minimizing the variance of the
     * residual PSD.
     *
     * \tparam printout if true then per-scale results are collected in m_regResults.
     *
     * \returns `error_t::noerror` on success, otherwise the regularization error.
     */
    template <bool printout = false>
    mx::error_t optimizeNc( realT &gmax_lp,           /**< [out] maximum gain for the selected predictor */
                            realT &gopt_lp,           /**< [out] optimum gain for the selected predictor */
                            int &Nc,                  /**< [out] selected number of coefficients */
                            realT &var_lp,            /**< [out] variance at the optimum gain */
                            clGainOpt<realT> &go_lp,  /**< [in] the gain optimization object */
                            std::vector<realT> &PSDt, /**< [in] the turbulence PSD */
                            std::vector<realT> &PSDn, /**< [in] the WFS noise PSD */
                            int minNc,                /**< [in] minimum number of coefficients */
                            int maxNc /**< [in] maximum number of coefficients */ )
    {
        realT minVar = std::numeric_limits<realT>::max();

        for( int n = minNc; n <= maxNc; ++n )
        {
            realT _gmax_lp;
            realT _gopt_lp;
            realT _var_lp;
            realT min_sc;
            error_t rv = regularizeCoefficients<printout>( _gmax_lp, _gopt_lp, _var_lp, min_sc, go_lp, PSDt, PSDn, n );
            if( rv != error_t::noerror )
            {
                return rv;
            }

            if( _var_lp < minVar )
            {
                gmax_lp = _gmax_lp;
                gopt_lp = _gopt_lp;
                var_lp = _var_lp;
                Nc = n;

                minVar = var_lp;
            }
        }

        return error_t::noerror;
    }
};

} // namespace analysis
} // namespace AO
} // namespace mx

#endif // clAOLinearPredictor_hpp
