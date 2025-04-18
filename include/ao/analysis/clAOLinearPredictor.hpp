/** \file clAOLinearPredictor.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Provides a class to manage closed loop gain linear predictor determination.
 * \ingroup mxAO_files
 *
 */

#ifndef clAOLinearPredictor_hpp
#define clAOLinearPredictor_hpp

#include <vector>

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

//#define CLAOLP_BREADCRUMB std::cerr << __FILE__ << ' ' << __LINE__ << '\n';


/// Class to manage the calculation of linear predictor coefficients for a closed-loop AO system.
/**
 * \tparam _realT the real floating point type in which to do all arithmetic.
 *
 * \ingroup mxAOAnalytic
 */
template <typename _realT>
struct clAOLinearPredictor
{
    typedef _realT realT;

public:

    struct regResult
    {
        realT sc;
        realT gopt;
        realT gmax;
        realT var;
    };

    std::vector<realT> m_PSDtn; ///< Working memory for the regularized PSD

    std::vector<realT> m_psd2s; ///< Working memory for the 2-sided regularized PSD

    std::vector<realT> m_ac; ///< Working memory to hold the autocorrelation.

    sigproc::autocorrelationFromPSD<realT> m_acpsd;

    sigproc::linearPredictor<realT> m_lp;

    realT m_min_var0{ 0 };
    realT m_min_sc0{ 10 };
    realT m_precision0{ 2 };
    realT m_max_sc0{ 100 };
    realT m_dPrecision{ 3 };

    realT m_gmax_lp{ 5 }; ///< The maximum allowable gain for LP.

    // Stopping conditions:
    realT m_minPrecision{ 0.001 };
    int m_maxIts{ 100 };

    int m_extrap {1}; ///< The LP extrapolation length in loop steps.  Normally it is 1 step.

    std::vector<regResult> m_regResults;
public:

    clAOLinearPredictor()
    {}

    /// Calculate the LP coefficients for a turbulence PSD and a noise PSD.
    /** This combines the two PSDs, augments to two-sided, and calls the linearPredictor.calcCoefficients method.
     *
     * A regularization constant can be added to the PSD as well.
     *
     */
    int calcCoefficients( std::vector<realT> &PSDt, ///< [in] the turbulence PSD
                          std::vector<realT> &PSDn, ///< [in] the WFS noise PSD
                          realT PSDreg,             ///< [in] the regularizing constant.  Set to 0 to not use.
                          int Nc,                   ///< [in] the number of LP coefficients.
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

        return m_lp.calcCoefficients( m_ac, Nc, m_extrap , condition );

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
     */
    template <bool telem>
    int _regularizeCoefficients( realT &min_var,  ///< [in.out] the minimum variance found.  Set to 0 on initial call
                                 realT &min_sc,   ///< [in.out] the scale factor at the minimum variance.
                                 realT precision, ///< [in] the step-size for the scale factor
                                 realT max_sc,    ///< [in] the maximum scale factor to test
                                 clGainOpt<realT> &go_lp,  ///< [in] the gain optimization object
                                 std::vector<realT> &PSDt, ///< [in] the turbulence PSD
                                 std::vector<realT> &PSDn, ///< [in] the WFS noise PSD
                                 int Nc                    ///< [in] the number of coefficients
    )
    {
        CLAOLP_BREADCRUMB;

        realT gmax_lp;
        realT gopt_lp;
        realT var_lp;

        realT sc0;

        if( min_var == 0 ) //first call
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
        //for( realT sc = sc0; sc <= max_sc; sc += precision )
        for( realT sc = max_sc; sc >= sc0; sc -= precision )
        {
            CLAOLP_BREADCRUMB;
            if( calcCoefficients( PSDt, PSDn, psdReg * pow( 10, -sc / 10 ), Nc ) < 0 )
            {
                return -1;
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
                m_regResults.push_back({sc, gopt_lp, gmax_lp, var_lp});
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
                return 0;
            }

            CLAOLP_BREADCRUMB;
        }

        CLAOLP_BREADCRUMB;
        return -1;
    }

    /// Regularize the PSD and calculate the associated LP coefficients.
    /** The PSD is regularized by adding a constant to it.  This constant is found by minimizing the variance of the
     * residual PSD.
     *
     * \tparam telem if true then the results are collected in m_regResults
     */
    template <bool telem = false>
    int regularizeCoefficients( realT &gmax_lp,           ///< [out] the maximum gain calculated for the regularized PSD
                                realT &gopt_lp,           ///< [out] the optimum gain calculated for the regularized PSD
                                realT &var_lp,            ///< [out] the variance at the optimum gain.
                                realT &min_sc,            ///< [out] the optimum regularization scale factor
                                clGainOpt<realT> &go_lp,  ///< [in] the gain optimization object
                                std::vector<realT> &PSDt, ///< [in] the turbulence PSD
                                std::vector<realT> &PSDn, ///< [in] the WFS noise PSD
                                int Nc                    ///< [in] the number of coefficients
    )
    {

        CLAOLP_BREADCRUMB;

        realT min_var = m_min_var0;
        min_sc = m_min_sc0;
        realT precision = m_precision0;
        realT max_sc = m_max_sc0;

        if(telem)
        {
            m_regResults.reserve(m_maxIts * 50);
        }

        CLAOLP_BREADCRUMB;
        int its = 0;
        while( precision > m_minPrecision && its < m_maxIts )
        {
            CLAOLP_BREADCRUMB;
            if(_regularizeCoefficients<telem>( min_var, min_sc, precision, max_sc, go_lp, PSDt, PSDn, Nc ) < 0)
            {
                return -1;
            }

            CLAOLP_BREADCRUMB;
            if( min_sc == max_sc )
            {
                if( its == 0 )
                {
                    min_sc -= precision;
                    max_sc = 200;
                }
                else
                {
                    // std::cerr << "Error in regularizeCoefficients.\n";
                    // return -1;
                    break;
                }
            }
            else
            {
                max_sc = min_sc + precision;
                precision /= m_dPrecision;
            }

            ++its;
        }

        CLAOLP_BREADCRUMB;
        // Now record final values
        if( calcCoefficients( PSDt, PSDn, PSDt[0] * pow( 10, -min_sc / 10 ), Nc ) < 0 )
        {
            return -1;
        }

        CLAOLP_BREADCRUMB;
        go_lp.a( m_lp.m_c );
        go_lp.b( m_lp.m_c );

        CLAOLP_BREADCRUMB;
        realT ll = 0, ul = 0;
        gmax_lp = go_lp.maxStableGain( ll, ul );
        gopt_lp = go_lp.optGainOpenLoop( var_lp, PSDt, PSDn, gmax_lp, false );

        CLAOLP_BREADCRUMB;
        return 0;
    }

    /// Regularize the PSD and calculate the associated LP coefficients.
    /** The PSD is regularized by adding a constant to it.  This constant is found by minimizing the variance of the
     * residual PSD.
     *
     * \tparam printout if true then the results are printed to stdout as they are calculated.
     */
    template <bool printout = false>
    int optimizeNc( realT &gmax_lp, ///< [out] the maximum gain calculated for the regularized PSD
                    realT &gopt_lp, ///< [out] the optimum gain calculated for the regularized PSD
                    int &Nc,
                    realT &var_lp,            ///< [out] the variance at the optimum gain.
                    clGainOpt<realT> &go_lp,  ///< [in] the gain optimization object
                    std::vector<realT> &PSDt, ///< [in] the turbulence PSD
                    std::vector<realT> &PSDn, ///< [in] the WFS noise PSD
                    int minNc,                ///< [in] the number of coefficients
                    int maxNc )
    {
        realT minVar = std::numeric_limits<realT>::max();

        for( int n = minNc; n <= maxNc; ++n )
        {
            realT _gmax_lp;
            realT _gopt_lp;
            realT _var_lp;
            regularizeCoefficients<printout>( _gmax_lp, _gopt_lp, _var_lp, go_lp, PSDt, PSDn, n );

            if( _var_lp < minVar )
            {
                gmax_lp = _gmax_lp;
                gopt_lp = _gopt_lp;
                var_lp = _var_lp;
                Nc = n;

                minVar = var_lp;
            }
        }

        return 0;
    }
};

} // namespace analysis
} // namespace AO
} // namespace mx

#endif // clAOLinearPredictor_hpp
