/** \file clGainOpt.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Provides a class to manage closed loop gain optimization.
 * \ingroup mxAO_files
 *
 */

//***********************************************************************//
// Copyright 2016-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef clGainOpt_hpp
#define clGainOpt_hpp

#ifdef MX_INCLUDE_BOOST
    #include <boost/math/tools/minima.hpp>
#endif

#include <type_traits>

#include <Eigen/Dense>

#include "../../sys/timeUtils.hpp"

#include "../../math/constants.hpp"

// #define ALLOW_F_ZERO

namespace mx
{
namespace AO
{
namespace analysis
{

// forward declaration of worker functor
template <typename realT>
struct clGainOptOptGain_OL;

/// A class to manage optimizing closed-loop gains
/**
 * \tparam _realT the real floating point type in which to do all arithmetic.  Is used to define the complex type as
 * well.
 *
 * \ingroup mxAOAnalytic
 */
template <typename _realT>
struct clGainOpt
{
    typedef _realT realT;                  ///< The real data type
    typedef std::complex<_realT> complexT; ///< The complex data type

  protected:
    int m_N;                 ///< Number of integrations in the (optional) moving average.  Default is 1.
    realT m_Ti;              ///< The loop sampling interval
    realT m_tau;             ///< The loop delay

    std::vector<realT> m_b;  ///< Vector of FIR coefficients
    std::vector<realT> m_a;  ///< Vector of IIR coefficients

    std::vector<realT> m_f;  ///< Vector of frequencies

    bool m_fChanged{ true }; ///< True if frequency or max size of m_a and m_b changes

    bool m_changed{ true };  ///< True if any of the members which make up the basic transfer functions are changed

    Eigen::Array<realT, -1, -1> m_cs;
    Eigen::Array<realT, -1, -1> m_ss;

    std::vector<std::complex<realT>> m_H_dm;
    std::vector<std::complex<realT>> m_H_wfs;
    std::vector<std::complex<realT>> m_H_ma;
    std::vector<std::complex<realT>> m_H_del;
    std::vector<std::complex<realT>> m_H_con;

  public:
    /** Parameters for stability analysis
     * @{
     */

    realT m_maxFindMin; ///< The Minimum value for the maximum stable gain finding algorithm.

    ///@}

    /** Parameters for minimization finding
     * @{
     */

    realT m_minFindMin;         ///< The Minimum value for the minimum finding algorithm.
    realT m_minFindMaxFact;     ///< The maximum value, as a multiplicative factor of maximum gain
    int m_minFindBits;          ///< The bits of precision to use for minimum finding. Defaults to
                                ///< std::numeric_limits<realT>::digits.
    uintmax_t m_minFindMaxIter; ///< The maximum iterations allowed for minimization.

    ///@}

    /// Default c'tor.
    clGainOpt();

    /// C'tor setting the loop timings.
    /**
     */
    clGainOpt( realT Ti, ///< [in] the desired loop sampling interval.
               realT tau ///< [in] the desired loop delay.
    );

    /// Initialize this instance.
    void init();

    /// Get the number of integrations in the (optional) moving average
    /**
     * \returns the current value of m_N.
     */
    int N();

    /// Set the number of integrations in the moving average
    /**
     */
    void N( int newN /**< [in] the value of m_N. */ );

    /// Get the loop sampling interval
    /**
     * \returns the current value of m_Ti.
     */
    realT Ti();

    /// Set the loop sampling interval
    /**
     */
    void Ti( realT newTi /**< [in]  the new value of m_Ti. */ );

    /// Get the loop delay
    /**
     * \returns the current value of m_tau.
     */
    realT tau();

    /// Set the loop delay
    /**
     */
    void tau( realT newTau /**< [in] the new value of m_tau.*/ );

    /// Set the vector of FIR coefficients
    /**
     */
    void b( const std::vector<realT> &newB /**< [in] a vector of coefficients, which is copied to m_b.*/ );

    /// Set the vector of FIR coefficients
    /**
     */
    void b( const Eigen::Array<realT, -1, -1>
                &newB /**< [in] a column-vector Eigen::Array of coefficients,
                                which is copied to m_b.*/ );

    const std::vector<realT> &b()
    {
        return m_b;
    }

    void bScale( realT scale )
    {
        for( size_t n = 0; n < m_b.size(); ++n )
        {
            m_b[n] *= scale;
        }
    }

    /// Set the vector of IIR coefficients
    /**
     */
    void a( const std::vector<realT> &newA /**< [in] a vector of coefficients, which is copied to m_a. */ );

    /// Set the vector of IIR coefficients
    /**
     */
    void a( const Eigen::Array<realT, -1, -1> &newA /**< [in] a column-vector Eigen::Array of
                                                              coefficients, which is copied to m_a.*/ );

    realT a( size_t i )
    {
        return m_a[i];
    }

    void aScale( realT scale )
    {
        for( size_t n = 0; n < m_a.size(); ++n )
        {
            m_a[n] *= scale;
        }
    }

    /// Set the FIR and IIR coefficients so that the control law is a leaky integrator.
    /** Set remember to 1.0 for a pure integrator control law.
     *
     */
    void setLeakyIntegrator( realT remember  /**< [in] a number usually close to 1 setting the amount "remembered"
                                                       from previous iterations.*/);

    /// Set the vector of frequencies
    /**
     */
    void f( realT *newF, ///< [in] a pointer to an array containing the new frequencies
            size_t nF    ///< [in] the number of elements of size(realT) in newF.
    );

    /// Set the vector of frequencies
    /**
     */
    void f( std::vector<realT> &newF /**< [in] a vector containing the new frequencies */ );

    /// Get the size of the frequency vector
    /**
     * \returns m_f.size()
     */
    size_t f_size()
    {
        return m_f.size();
    }

    /// Get the i-th value of frequency.
    /** No range checks are conducted.
     *
     * \returns the value of m_f[i]
     *
     */
    realT f( size_t i /**< [in] the index of the frequency to return*/ );

    /// Calculate the open-loop transfer function
    /**
     * \return the complex value of the open-loop transfer function at f.
     */
    complexT olXfer( int fi,          ///< [in] the index of the frequency
                     complexT &H_dm,  ///< [out] the transfer function of the DM
                     complexT &H_del, ///< [out] the delay transfer function
                     complexT &H_con  ///< [out] the controller transfer function.
    );

    /// Calculate the open-loop transfer function
    /**
     * \overload
     *
     * \returns the complex value of the open-loop transfer function at f[fi].
     */
    complexT olXfer( int fi /**< [in] the index of the frequency */ );

    /// Return the closed loop error transfer function (ETF) at frequency f for gain g.
    /**
     * \returns the closed loop ETF at f and g.
     */
    complexT clETF( int fi, ///< [in] the index of the frequency at which to calculate the ETF
                    realT g ///< [in] the loop gain.
    );

    /// Return the closed loop error transfer function (ETF) phase at frequency f for gain g.
    /**
     * \returns the phase of the closed loop ETF at f and g.
     */
    realT clETFPhase( int fi, ///< [in] the index of the frequency at which to calculate the ETF
                      realT g ///< [in] the loop gain.
    );

    /// Return the norm of the closed loop error transfer function (ETF) at frequency f for gain g.
    /**
     * \returns the norm of the closed loop ETF at f and g.
     */
    realT clETF2( int fi, ///< [in] the index of the frequency at which to calculate the ETF.
                  realT g ///< [in] the loop gain.
    );

    /// Return the closed loop noise transfer function (NTF) at frequency f for gain g.
    /**
     * \returns the closed loop NTF at f and g.
     */
    complexT clNTF( int fi, ///< [in] the index of the frequency at which to calculate the NTF
                    realT g ///< [in] the loop gain.
    );

    /// Return the norm of the closed loop noise transfer function (NTF) at frequency f for gain g.
    /**
     * \returns the value of the closed loop NTF at f and g.
     */
    realT clNTF2( int fi, ///< [in] the index of the frequency at which to calculate the NTF
                  realT g ///< [in] the loop gain.
    );

    /// Return the norm of the closed loop transfer functions at frequency f for gain g.
    /** Calculates both the error transfer function (ETF) and the noise transfer function (NTF).
     * This minimizes the various complex number operations, compared to calling both clETF2 and clNTF2.
     *
     */
    void clTF2( realT &ETF, ///< [out] is set to the ETF at f and g
                realT &NTF, ///< [out] is set to the NTF at f and g
                int fi,     ///< [in] is the index of the frequency at which to calculate the ETF
                realT g     ///< [in] is the loop gain.
    );

    /// Calculate the closed loop variance given open-loop PSDs and gain
    /** Calculates the following quantities.
      \f[
      \sigma_{err}^2 = \sum_i \left| ETF(f_i) \right|^2 PSD_{err}(fi) \Delta f\\
      \sigma_{noise}^2 = \sum_i \left| NTF(f_i) \right|^2 PSD_{noise}(fi) \Delta f\\
      \sigma^2 = \sigma_{err}^2 + \sigma_{noise}^2
      \f]
      * \f$ \sigma^2 \f$ is returned, and \f$ \sigma_{err}^2 \f$ and \f$ \sigma_{noise}^2 \f$ are available as the
      optional
      * arguments varErr and varNoise.
      *
      * \returns the total variance (error + noise) in closed loop
      */
    realT clVariance( realT &varErr,                      ///< [out] the variance in the residual process error.
                      realT &varNoise,                    ///< [out] the variance in the residual measurement noise.
                      const std::vector<realT> &PSDerr,   ///< [in] the open-loop process error PSD.
                      const std::vector<realT> &PSDnoise, ///< [in] the open-loop measurement noise PSD.
                      realT g                             ///< [in] the gain.
    );

    /// Calculate the closed loop variance given open-loop PSDs and gain
    /** Overload of clVariance without the varErr and varNoise output parameters.
     *
     * \overload
     *
     * \returns the total variance (error + noise) in closed loop
     */
    realT clVariance( const std::vector<realT> &PSDerr,   ///< [in] the open-loop process error PSD.
                      const std::vector<realT> &PSDnoise, ///< [in] the open-loop measurement noise PSD.
                      realT g                             ///< [in] the gain.
    );

    /// Find the maximum stable gain for the loop parameters
    /**
     *
     * \returns the maximum stable gain for the loop parameters
     */

    /// Find the maximum stable gain for the loop parameters
    /** Conducts a search along the Nyquist contour of the open-loop transfer function to find
     * the most-negative crossing of the real axis.
     *
     * \returns the maximum stable gain for the loop parameters
     */
    realT maxStableGain( realT &ll, ///< [in.out] the lower limit used for the search
                         realT &ul  ///< [in.out] the upper limit used for hte search
    );

    /// Find the maximum stable gain for the loop parameters
    /** Conducts a search along the Nyquist contour of the open-loop transfer function to find
     * the most-negative crossing of the real axis.
     *
     * This version allows constant arguments.
     * \overload
     *
     * \returns the maximum stable gain for the loop parameters
     */
    realT maxStableGain( const realT &ll, ///< [in] the lower limit used for the search
                         const realT &ul  ///< [in] the upper limit used for hte search
    );

    /// Find the maximum stable gain for the loop parameters
    /** Conducts a search along the Nyquist contour of the open-loop transfer function to find
     * the most-negative crossing of the real axis.
     *
     * This version uses m_maxFindMin for the lower limit and no upper limit.
     *
     * \overload
     *
     * \returns the maximum stable gain for the loop parameters
     */
    realT maxStableGain();

    /// Return the optimum closed loop gain given an open loop PSD
    /** Uses _gmax for the upper limit.
     * \returns the optimum gain
     */
    realT optGainOpenLoop( realT &var,                         ///< [out] the variance at the optimum gain
                           const std::vector<realT> &PSDerr,   ///< [in] open loop error PSD
                           const std::vector<realT> &PSDnoise, ///< [in] open loop measurement noise PSD
                           bool gridSearch                     /**< [in] flag controlling whether an intial grid
                                                                         search is done to find the global minimum*/
    );

    /// Return the optimum closed loop gain given an open loop PSD
    /**
     * \returns the optimum gain.
     */
    realT optGainOpenLoop( realT &var,                         ///< [out] the variance at the optimum gain
                           const std::vector<realT> &PSDerr,   ///< [in] open loop error PSD
                           const std::vector<realT> &PSDnoise, ///< [in] open loop measurement noise PSD
                           realT &gmax,                        /**< [in] maximum gain to consider.
                                                                         If 0, then _gmax is used.*/
                           bool gridSearch                     /**< [in] flag controlling whether an intial grid
                                                                         search is done to find the global minimum*/
    );

    /// Calculate the pseudo open-loop PSD given a closed loop PSD
    /**
     * \returns 0 on success
     */
    int pseudoOpenLoop( std::vector<realT> &PSD, /**< [in.out] input closed loop PSD, on output contains the pseudo open
                                                               loop error PSD */
                        realT g                  ///< [in] the loop gain when PSD was measured.
    );

    int nyquist( std::vector<realT> &re, std::vector<realT> &im, realT g );
};

template <typename realT>
clGainOpt<realT>::clGainOpt()
{
    init();
}

template <typename realT>
clGainOpt<realT>::clGainOpt( realT Ti, realT tau )
{
    init();

    m_Ti = Ti;
    m_tau = tau;
}

template <typename realT>
void clGainOpt<realT>::init()
{
    m_N = 1;

    setLeakyIntegrator( 1.0 );

    m_Ti = 1. / 1000.;
    m_tau = 2.5 * m_Ti;

    m_maxFindMin = 0.0;

    m_minFindMin = 1e-9;
    m_minFindMaxFact = 0.999;
    m_minFindBits = std::numeric_limits<realT>::digits;
    m_minFindMaxIter = 10000;

    m_fChanged = true;
    m_changed = true;
}

template <typename realT>
int clGainOpt<realT>::N()
{
    return m_N;
}

template <typename realT>
void clGainOpt<realT>::N( int newN )
{
    if( m_N == newN )
    {
        return;
    }

    m_N = newN;
    m_changed = true;
}

template <typename realT>
realT clGainOpt<realT>::Ti()
{
    return m_Ti;
}

template <typename realT>
void clGainOpt<realT>::Ti( realT newTi )
{
    if( m_Ti == newTi )
    {
        return;
    }

    m_Ti = newTi;
    m_changed = true;
}

template <typename realT>
realT clGainOpt<realT>::tau()
{
    return m_tau;
}

template <typename realT>
void clGainOpt<realT>::tau( realT newTau )
{
    if( m_tau == newTau )
    {
        return;
    }

    m_tau = newTau;
    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::b( const std::vector<realT> &newB )
{
    if( newB.size() > (size_t)m_cs.cols() )
    {
        m_fChanged = true;
    }

    m_b = newB;
    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::b( const Eigen::Array<realT, -1, -1> &newB )
{
    if( newB.cols() > m_cs.cols() )
    {
        m_fChanged = true;
    }

    m_b.resize( newB.cols() );

    for( size_t i = 0; i < m_b.size(); ++i )
    {
        m_b[i] = newB( 0, i );
    }

    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::a( const std::vector<realT> &newA )
{
    if( newA.size() + 1 > (size_t)m_cs.cols() )
    {
        m_fChanged = true;
    }

    m_a = newA;
    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::a( const Eigen::Array<realT, -1, -1> &newA )
{
    if( newA.cols() + 1 > m_cs.cols() )
    {
        m_fChanged = true;
    }

    m_a.resize( newA.cols() );

    for( size_t i = 0; i < m_a.size(); ++i )
    {
        m_a[i] = newA( 0, i );
    }

    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::setLeakyIntegrator( realT remember )
{
    m_b.resize( 1 );
    m_b[0] = 1.0;

    m_a.resize( 1 );
    m_a[0] = remember;

    m_fChanged = true;
    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::f( realT *newF, size_t nF )
{
    m_f.resize( nF );
    for( int i = 0; i < nF; ++i )
        m_f[i] = newF[i];

    m_fChanged = true;
    m_changed = true;
}

template <typename realT>
void clGainOpt<realT>::f( std::vector<realT> &newF )
{
    m_f = newF;
    m_fChanged = true;

    m_changed = true;
}

template <typename realT>
realT clGainOpt<realT>::f( size_t i )
{

    return m_f[i];
}

template <typename realT>
std::complex<realT> clGainOpt<realT>::olXfer( int fi )
{
    complexT H_dm;
    complexT H_del;
    complexT H_con;

    return olXfer( fi, H_dm, H_del, H_con );
}

// If PRECALC_TRIG is defined, then the cosine and sine tables are pre-calculated and used instead of repeated exp(-i)
// calls. This is much much faster, though uses more memory.  In general, only undefine this for testing or debugging.
#define PRECALC_TRIG

template <typename realT>
std::complex<realT> clGainOpt<realT>::olXfer( int fi, complexT &H_dm, complexT &H_del, complexT &H_con )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
    {
#else
    if( m_f[fi] < 0 )
    {
#endif

        H_dm = 0;
        H_del = 0;
        H_con = 0;
        return 0;
    }

#ifdef PRECALC_TRIG
    if( m_fChanged )
    {
        size_t jmax = std::max( m_a.size() + 1, m_b.size() );

        m_cs.resize( m_f.size(), jmax );
        m_ss.resize( m_f.size(), jmax );

        for( size_t i = 0; i < m_f.size(); ++i )
        {
            m_cs( i, 0 ) = 1.0;
            m_ss( i, 0 ) = 0.0;

            for( size_t j = 1; j < jmax; ++j )
            {
                m_cs( i, j ) = cos( math::two_pi<realT>() * m_f[i] * m_Ti * realT( j ) );
                m_ss( i, j ) = sin( math::two_pi<realT>() * m_f[i] * m_Ti * realT( j ) );
            }
        }

        m_fChanged = false;
    }
#endif

    if( m_changed )
    {
        m_H_dm.resize( m_f.size(), 0 );
        m_H_wfs.resize( m_f.size(), 0 );
        m_H_ma.resize( m_f.size(), 0 );
        m_H_del.resize( m_f.size(), 0 );
        m_H_con.resize( m_f.size(), 0 );

        size_t jmax = std::min( m_a.size(), m_b.size() );

        // #pragma omp parallel for
        for( size_t i = 0; i < m_f.size(); ++i )
        {
#ifndef ALLOW_F_ZERO
            if( m_f[i] <= 0 )
                continue;
#else
            if( m_f[i] < 0 )
                continue;
#endif

            complexT s = complexT( 0.0, math::two_pi<realT>() * m_f[i] );

            complexT expsT = exp( -s * m_Ti );

            if( m_f[i] == 0 )
            {
                m_H_dm[i] = std::complex<realT>( 1, 0 );
            }
            else
            {
                m_H_dm[i] = ( realT( 1 ) - expsT ) / ( s * m_Ti );
            }

            m_H_wfs[i] = m_H_dm[i];

            m_H_ma[i] = 1; // realT(1./m_N)*(realT(1) - pow(expsT,m_N))/(realT(1) - expsT);

            m_H_del[i] = exp( -s * m_tau );

            complexT FIR = complexT( m_b[0], 0 );

            complexT IIR = complexT( 0.0, 0.0 );
            for( size_t j = 1; j < jmax; ++j )
            {
#ifdef PRECALC_TRIG
                realT cs = m_cs( i, j );
                realT ss = m_ss( i, j );
                FIR += m_b[j] * complexT( cs, -ss );
                IIR += m_a[j - 1] * complexT( cs, -ss );
#else
                complexT expZ = exp( -s * m_Ti * realT( j ) );
                FIR += m_b[j] * expZ;
                IIR += m_a[j - 1] * expZ;
#endif
            }

            for( size_t jj = jmax; jj < m_a.size() + 1; ++jj )
            {
#ifdef PRECALC_TRIG
                realT cs = m_cs( i, jj );
                realT ss = m_ss( i, jj );
                IIR += m_a[jj - 1] * complexT( cs, -ss );
#else
                complexT expZ = exp( -s * m_Ti * realT( jj ) );
                IIR += m_a[jj - 1] * expZ;
#endif
            }

            for( size_t jj = jmax; jj < m_b.size(); ++jj )
            {
#ifdef PRECALC_TRIG
                realT cs = m_cs( i, jj );
                realT ss = m_ss( i, jj );
                FIR += m_b[jj] * complexT( cs, -ss );
#else
                complexT expZ = exp( -s * m_Ti * realT( jj ) );
                FIR += m_b[jj] * expZ;
#endif
            }

            m_H_con[i] = FIR / ( realT( 1.0 ) - IIR );

            /*if( i == 0 || i == 1)
            {
               std::cerr << i << " " << m_f[fi] << " " << s << " " << expsT << " " << m_H_wfs[i] << " " << m_H_dm[i] <<
            " "
            << m_H_con[i] << " " << m_H_del[i] << "\n";
               //exit(0);
            }/**/
        }

        m_changed = false;
    }

    H_dm = m_H_dm[fi];
    H_del = m_H_del[fi]; //*m_H_ma[fi];
    H_con = m_H_con[fi];

    return ( m_H_dm[fi] * m_H_wfs[fi] * m_H_del[fi] * m_H_con[fi] );
}

template <typename realT>
std::complex<realT> clGainOpt<realT>::clETF( int fi, realT g )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
        return 0;
#else
    if( m_f[fi] < 0 )
        return 0;
#endif

    return ( realT( 1 ) / ( realT( 1 ) + g * olXfer( fi ) ) );
}

template <typename realT>
realT clGainOpt<realT>::clETFPhase( int fi, realT g )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
        return 0;
#else
    if( m_f[fi] < 0 )
        return 0;
#endif

    return std::arg( ( realT( 1 ) / ( realT( 1 ) + g * olXfer( fi ) ) ) );
}

template <typename realT>
realT clGainOpt<realT>::clETF2( int fi, realT g )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
        return 0;
#else
    if( m_f[fi] < 0 )
        return 0;
#endif

    return norm( realT( 1 ) / ( realT( 1 ) + g * olXfer( fi ) ) );
}

template <typename realT>
std::complex<realT> clGainOpt<realT>::clNTF( int fi, realT g )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
        return 0;
#else
    if( m_f[fi] < 0 )
        return 0;
#endif

    complexT H_dm, H_del, H_con;

    complexT olX = olXfer( fi, H_dm, H_del, H_con ); // H_dm*H_wfs*H_ma*H_del*H_con;

    return -( H_dm * H_del * g * H_con ) / ( realT( 1 ) + g * olX );
}

template <typename realT>
realT clGainOpt<realT>::clNTF2( int fi, realT g )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
        return 0;
#else
    if( m_f[fi] < 0 )
        return 0;
#endif

    complexT H_dm, H_del, H_con;

    complexT olX = olXfer( fi, H_dm, H_del, H_con ); // H_dm*H_wfs*H_ma*H_del*H_con;

    complexT NTF = -( H_dm * H_del * g * H_con ) / ( realT( 1 ) + g * olX );

    return norm( NTF );
}

template <typename realT>
void clGainOpt<realT>::clTF2( realT &ETF, realT &NTF, int fi, realT g )
{
#ifndef ALLOW_F_ZERO
    if( m_f[fi] <= 0 )
#else
    if( m_f[fi] < 0 )
#endif
    {
        ETF = 0;
        NTF = 0;
        return;
    }

    complexT H_dm, H_del, H_con;

    complexT olX = olXfer( fi, H_dm, H_del, H_con ); // H_dm*H_wfs*H_ma*H_del*H_con;

    if( m_f[fi] == 0 )
    {
    }

    ETF = norm( realT( 1 ) / ( realT( 1 ) + g * olX ) );
    NTF = norm( -( H_dm * H_del * g * H_con ) / ( realT( 1 ) + g * olX ) );

    /*if(m_f[fi] == 0)
    {
       std::cerr << "ETF: " << ETF << " NTF: " << NTF << "\n";
    }*/
}

template <typename realT>
realT clGainOpt<realT>::clVariance(
    realT &varErr, realT &varNoise, const std::vector<realT> &PSDerr, const std::vector<realT> &PSDnoise, realT g )
{
    if( m_f.size() != PSDerr.size() || m_f.size() != PSDnoise.size() )
    {
        std::cerr << "clVariance: Frequency grid and PSDs must be same size." << std::endl;
        return -1;
    }

    realT ETF, NTF, df;

    varErr = 0;
    varNoise = 0;

    df = m_f[1] - m_f[0];

    for( size_t i = 0; i < PSDerr.size(); ++i )
    {
        if( g == 0 )
        {
            ETF = 1;
            NTF = 0;
        }
        else
        {
            clTF2( ETF, NTF, i, g );
        }
        varErr += ETF * PSDerr[i] * df;
        varNoise += NTF * PSDnoise[i] * df;
    }

    return varErr + varNoise;
}

template <typename realT>
realT clGainOpt<realT>::clVariance( const std::vector<realT> &PSDerr, const std::vector<realT> &PSDnoise, realT g )
{
    realT varErr;
    realT varNoise;

    return clVariance( varErr, varNoise, PSDerr, PSDnoise, g );
}

template <typename realT>
realT clGainOpt<realT>::maxStableGain( realT &ll, realT &ul )
{
    static_cast<void>( ul );

    std::vector<realT> re, im;

    if( ll == 0 )
        ll = m_maxFindMin;

    nyquist( re, im, 1.0 );

    int gi_c = re.size() - 1;

    for( int gi = re.size() - 2; gi >= 0; --gi )
    {
        if( -1.0 / re[gi] < ll )
            continue;

        if( ( re[gi] < 0 ) && ( im[gi + 1] >= 0 && im[gi] < 0 ) )
        {
            // Check for loop back in Nyquist diagram
            if( re[gi] <= re[gi_c] )
                gi_c = gi;
        }
    }

    return -1.0 / re[gi_c];
}

template <typename realT>
realT maxStableGain( const realT &ll, const realT &ul )
{
    realT rll = ll;
    realT rul = ul;

    maxStableGain( rll, rul );
}

template <typename realT>
realT clGainOpt<realT>::maxStableGain()
{
    realT ul = 0;
    realT ll = m_maxFindMin;

    return maxStableGain( ll, ul );
}

// Implement the minimization, allowing pre-compiled specializations
namespace impl
{

template <typename realT>
realT optGainOpenLoop( clGainOptOptGain_OL<realT> &olgo,
                       realT &var,
                       const realT &gmax,
                       const realT &minFindMin,
                       const realT &minFindMaxFact,
                       int minFindBits,
                       uintmax_t minFindMaxIter,
                       uintmax_t &iters )
{
#ifdef MX_INCLUDE_BOOST
    realT gopt;

    try
    {
        std::pair<realT, realT> brack;
        brack = boost::math::tools::brentm_findm_minima<clGainOptOptGain_OL<realT>, realT>( olgo,
                                                                                            minFindMin,
                                                                                            minFindMaxFact * gmax,
                                                                                            minFindBits,
                                                                                            minFindMaxIter,
                                                                                            iters );
        gopt = brack.first;
        var = brack.second;
    }
    catch( ... )
    {
        std::cerr << "optGainOpenLoop: No root found\n";
        gopt = minFindMaxFact * gmax;
        var = 0;
    }

    return gopt;
#else
    staticm_assert( std::is_fundamental<realT>::value || !std::is_fundamental<realT>::value,
                    "impl::optGainOpenLoop<realT> is not specialized for type realT, and MX_INCLUDE_BOOST is not "
                    "defined, so I can't just use boost." );
    return 0;
#endif
}

template <>
float optGainOpenLoop<float>( clGainOptOptGain_OL<float> &olgo,
                              float &var,
                              const float &gmax,
                              const float &minFindMin,
                              const float &minFindMaxFact,
                              int minFindBits,
                              uintmax_t minFindMaxIter,
                              uintmax_t &iters );

template <>
double optGainOpenLoop<double>( clGainOptOptGain_OL<double> &olgo,
                                double &var,
                                const double &gmax,
                                const double &minFindMin,
                                const double &minFindMaxFact,
                                int minFindBits,
                                uintmax_t minFindMaxIter,
                                uintmax_t &iters );

template <>
long double optGainOpenLoop<long double>( clGainOptOptGain_OL<long double> &olgo,
                                          long double &var,
                                          const long double &gmax,
                                          const long double &minFindMin,
                                          const long double &minFindMaxFact,
                                          int minFindBits,
                                          uintmax_t minFindMaxIter,
                                          uintmax_t &iters );

#ifdef HASQUAD
template <>
_m_float128 optGainOpenLoop<_m_float128>( clGainOptOptGain_OL<_m_float128> &olgo,
                                          _m_float128 &var,
                                          const _m_float128 &gmax,
                                          const _m_float128 &minFindMin,
                                          const _m_float128 &minFindMaxFact,
                                          int minFindBits,
                                          uintmax_t minFindMaxIter,
                                          uintmax_t &iters );
#endif

} // namespace impl

template <typename realT>
realT clGainOpt<realT>::optGainOpenLoop( realT &var,
                                         const std::vector<realT> &PSDerr,
                                         const std::vector<realT> &PSDnoise,
                                         bool gridSearch )
{
    realT gmax = 0;
    return optGainOpenLoop( var, PSDerr, PSDnoise, gmax, gridSearch );
}

template <typename realT>
realT clGainOpt<realT>::optGainOpenLoop(
    realT &var, const std::vector<realT> &PSDerr, const std::vector<realT> &PSDnoise, realT &gmax, bool gridSearch )
{
    clGainOptOptGain_OL<realT> olgo;
    olgo.go = this;
    olgo.PSDerr = &PSDerr;
    olgo.PSDnoise = &PSDnoise;

    if( gmax <= 0 )
        gmax = maxStableGain();

    realT ming = m_minFindMin;
    realT maxg = gmax;

    if( gridSearch )
    {
        realT gstpsz = 0.05;
        realT gg = m_minFindMaxFact * gmax;
        realT var0 = clVariance( PSDerr, PSDnoise, gg );
        realT mingg = gg;

        while( gg > m_minFindMin )
        {
            gg -= gstpsz;
            realT var1 = clVariance( PSDerr, PSDnoise, gg );

            if( var1 < var0 )
            {
                var0 = var1;
                mingg = gg;
            }
        }

        ming = mingg - gstpsz;
        maxg = mingg + gstpsz;

        if( ming < m_minFindMin )
            ming = m_minFindMin;
        if( maxg > gmax )
            maxg = gmax;
    }

    uintmax_t iters;
    realT val =
        impl::optGainOpenLoop( olgo, var, maxg, ming, m_minFindMaxFact, m_minFindBits, m_minFindMaxIter, iters );

    if( iters >= m_minFindMaxIter )
    {
        // #pragma omp critical
        {
            std::cerr << "\nclGainOpt<realT>::optGainOpenLoop: minFindMaxIter (" << m_minFindMaxIter << ") reached\n";
        }
    }

    return val;
}

template <typename realT>
int clGainOpt<realT>::pseudoOpenLoop( std::vector<realT> &PSD, realT g )
{
    realT e;
    for( int f = 0; f < m_f.size(); ++f )
    {
        e = clETF2( f, g );

        if( e > 0 )
            PSD[f] = PSD[f] / e;
    }

    return 0;
}

template <typename realT>
int clGainOpt<realT>::nyquist( std::vector<realT> &re, std::vector<realT> &im, realT g )
{
    re.resize( m_f.size() );
    im.resize( m_f.size() );

    complexT etf;

    for( size_t f = 0; f < m_f.size(); ++f )
    {
        etf = g * olXfer( f ); // clETF(f, g);
        re[f] = real( etf );
        im[f] = imag( etf );
    }

    return 0;
}

//------------ Workers ---------------------

/// Bisection worker struct for finding optimum closed loop gain from open loop PSDs
template <typename realT>
struct clGainOptOptGain_OL
{
    clGainOpt<realT> *go;
    const std::vector<realT> *PSDerr;
    const std::vector<realT> *PSDnoise;

    realT operator()( const realT &g )
    {
        return go->clVariance( *PSDerr, *PSDnoise, g );
    }
};

// Explicit Instantiation
extern template class clGainOpt<float>;

extern template class clGainOpt<double>;

extern template class clGainOpt<long double>;

#ifdef HASQUAD
extern template class clGainOpt<_m_float128>;
#endif

} // namespace analysis
} // namespace AO
} // namespace mx

#endif // clGainOpt_hpp
