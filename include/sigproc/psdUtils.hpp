/** \file psdUtils.hpp
 * \brief Tools for working with PSDs
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup signal_processing_files
 *
 */

//***********************************************************************//
// Copyright 2015-2021 Jared R. Males (jaredmales@gmail.com)
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

#ifndef psdUtils_hpp
#define psdUtils_hpp

#ifndef EIGEN_NO_CUDA
#define EIGEN_NO_CUDA
#endif

#include <type_traits>

#include <Eigen/Dense>

#include <iostream>

#include "../mxError.hpp"

#include "../math/fft/fft.hpp"
#include "../math/vectorUtils.hpp"

namespace mx
{
namespace sigproc
{

/** \ingroup psds
 * @{
 */

/// Calculate the variance of a 1-D, 1-sided PSD
/** By default uses trapezoid rule integration.  This can be changed to mid-point integration.
 *
 * \returns the variance of a PSD (the integral).
 *
 * \tparam realT the real floating point type
 *
 * \test Scenario: calculating variance from a 1D PSD \ref tests_sigproc_psdUtils_psdVar_1D "[test doc]"
 */
template <typename realT>
realT psdVar1sided( realT df,         ///< [in] the frequency scale of the PSD
                    const realT *PSD, ///< [in] the PSD to integrate.
                    size_t sz,        ///< [in] the size of the PSD vector
                    realT half = 0.5  ///< [in] [optional] controls if trapezoid (0.5) or mid-point (1.0) integration is
                                      ///< used.  Do not use other values.
)
{
    realT var = 0;

    var = half * PSD[0];

    for( size_t i = 1; i < sz - 1; ++i )
        var += PSD[i];

    var += half * PSD[sz - 1];

    var *= df;

    return var;
}

/// Calculate the variance of a 1-D, 2-sided PSD
/** By default uses trapezoid rule integration.  This can be changed to mid-point integration.
 *
 * Assumes the 2-sided PSD is in standard FFT storage order, and that sz is even.
 *
 * \returns the variance of a PSD (the integral).
 *
 * \tparam realT the real floating point type
 *
 * \test Scenario: calculating variance from a 1D PSD. \ref tests_sigproc_psdUtils_psdVar_1D "[test doc]"
 */
template <typename realT>
realT psdVar2sided( realT df,         ///< [in] the frequency scale of the PSD
                    const realT *PSD, ///< [in] the PSD to integrate.
                    size_t sz,        ///< [in] the size of the PSD vector
                    realT half = 0.5  ///< [in] [optional] controls if trapezoid (0.5) or mid-point (1.0) integration is
                                      ///< used.  Do not use other values.
)
{
    realT var = 0;

    var = PSD[0];

    size_t i = 1;
    for( ; i < sz / 2; ++i )
    {
        var += PSD[i];
        var += PSD[sz - i];
    }
    var += half * PSD[i]; // The mid-point is double.  It is also the endpoint of integration from each side, so it
                          // would enter twice, hence once here.

    var *= df;

    return var;
}

/// Calculate the variance of a 1-D PSD
/** By default uses trapezoid rule integration.  This can be changed to mid-point integration.
 *
 * If f.back() < 0, then a 2-sided PSD in FFT storage order is assumed.  Otherwise, PSD is treated as 1-sided.
 *
 * \returns the variance of a PSD (the integral).
 *
 * \tparam realT the real floating point type
 *
 * \test Scenario: calculating variance from a 1D PSD. \ref tests_sigproc_psdUtils_psdVar_1D "[test doc]"
 */
template <typename realT>
realT psdVar( const std::vector<realT> &f,   ///< [in] the frequency scale of the PSD.
              const std::vector<realT> &PSD, ///< [in] the PSD to integrate.
              realT half = 0.5 ///< [in] [optional] controls if trapezoid (0.5) or mid-point (1.0) integration is used.
                               ///< Do not use other values.
)
{
    if( f.back() < 0 )
        return psdVar2sided( f[1] - f[0], PSD.data(), PSD.size(), half );
    else
        return psdVar1sided( f[1] - f[0], PSD.data(), PSD.size(), half );
}

/// Calculate the variance of a PSD
/** By default uses trapezoid rule integration.  This can be changed to mid-point integration.
 *
 * \overload
 *
 * \returns the variance of a PSD (the integral).
 *
 * \tparam realT the real floating point type
 */
template <typename eigenArrT>
typename eigenArrT::Scalar psdVarDisabled(
    eigenArrT &freq, ///< [in] the frequency scale of the PSD
    eigenArrT &PSD,  ///< [in] the PSD to integrate.
    bool trap = true ///< [in] [optional] controls if trapezoid (true) or mid-point (false) integration is used.
)
{
    typename eigenArrT::Scalar half = 0.5;
    if( !trap )
        half = 1.0;

    typename eigenArrT::Scalar var = 0;

    var = half * PSD( 0, 0 );

    for( int i = 1; i < freq.rows() - 1; ++i )
        var += PSD( i, 0 );

    var += half * PSD( freq.rows() - 1, 0 );

    var *= ( freq( 1, 0 ) - freq( 0, 0 ) );

    return var;
}

/// Calculates the frequency sampling for a grid given maximum dimension and maximum frequency.
/** The freq_sampling is
 * @f$ \Delta f = f_{max}/ (0.5*dim) @f$
 * where @f$ f_{max} = 1/(2\Delta t) @f$ is the maximum frequency and @f$ dim @f$ is the size of the grid.
 *
 * \param [in] dim is the size of the grid
 * \param [in] f_max is the maximum frequency of the grid
 *
 * \returns the sampling interval @f$ \delta f @f$
 *
 * \tparam realT is the real floating point type used for calculations.
 *
 */
template <class realT>
realT freq_sampling( size_t dim, realT f_max )
{
    return ( f_max / ( 0.5 * dim ) );
}

#if 0
///Create a 1-D frequency grid
/**
  * \param [out] vec the pre-allocated Eigen-type 1xN or Nx1 array, on return contains the frequency grid
  * \param [in] dt the temporal sampling of the time series
  * \param [in] inverse [optional] if true
  *
  * \tparam eigenArr the Eigen-like array type
  */
template<typename eigenArr>
void frequency_grid1D( eigenArr & vec,
                       typename eigenArr::Scalar dt,
                       bool inverse = false )
{
   typename eigenArr::Index dim, dim_1, dim_2;
   typename eigenArr::Scalar df;

   dim_1 = vec.rows();
   dim_2 = vec.cols();

   dim = std::max(dim_1, dim_2);

   df = freq_sampling(dim, 0.5/dt);

   if( !inverse )
   {
      for(int ii=0; ii < ceil(0.5*(dim-1) + 1); ++ii)
      {
         vec(ii) = ii*df;
      }

      for(int ii=ceil(0.5*(dim-1)+1); ii < dim_1; ++ii)
      {
         vec(ii) = (ii-dim)*df;
      }
   }
   else
   {
      for(int ii=0; ii < dim; ++ii)
      {
         vec(ii) = df * ii / dim;
      }
   }
}
#endif

/// Create a 1-D frequency grid
/**
 *
 * \tparam realT a real floating point type
 * \tparam realParamT a real floating point type, convenience to avoid double-float confusion.
 *
 * \test Verify creation of a 1D frequency grid \ref tests_sigproc_psdUtils_frequencyGrid_1D "[test doc]"
 */
template <typename realT, typename realParamT>
int frequencyGrid(
    std::vector<realT> &vec, ///< [out] vec the pre-allocated vector, on return contains the frequency grid
    realParamT dt,           ///< [in] dt the temporal sampling of the time series
    bool fftOrder = true     ///< [in] fftOrder [optional] if true the frequency grid is in FFT order
)
{
    realT dtTT = dt;

    if( fftOrder )
    {
        if( vec.size() % 2 == 1 )
        {
            mxError( "frequencyGrid", MXE_INVALIDARG, "Frequency scale can't be odd-sized for FFT order" );
            return -1;
        }

        realT df = ( 1.0 / dtTT ) / ( (realT)vec.size() );

        for( ssize_t ii = 0; ii < ceil( 0.5 * ( vec.size() - 1 ) + 1 ); ++ii )
        {
            vec[ii] = ii * df;
        }

        for( ssize_t ii = ceil( 0.5 * ( vec.size() - 1 ) + 1 ); ii < vec.size(); ++ii )
        {
            vec[ii] = ( ii - (ssize_t)vec.size() ) * df;
        }

        return 0;
    }
    else
    {
        if( vec.size() % 2 == 0 )
        {
            realT df = ( 0.5 / dtTT ) / ( (realT)vec.size() - 1 );
            for( int ii = 0; ii < vec.size(); ++ii )
            {
                vec[ii] = df * ii;
            }

            return 0;
        }
        else
        {
            realT df = ( 0.5 / dt ) / ( (realT)vec.size() );
            for( int ii = 0; ii < vec.size(); ++ii )
            {
                vec[ii] = df * ( ii + 1 );
            }

            return 0;
        }
    }
}

/// Create a 2-D frequency grid
template <typename eigenArr, typename realParamT>
void frequencyGrid( eigenArr &arr, realParamT drT, eigenArr *k_x, eigenArr *k_y )
{
    typename eigenArr::Scalar dr = drT;

    typename eigenArr::Index dim_1, dim_2;
    typename eigenArr::Scalar k_1, k_2, df;

    dim_1 = arr.rows();
    dim_2 = arr.cols();

    if( k_x )
        k_x->resize( dim_1, dim_2 );
    if( k_y )
        k_y->resize( dim_1, dim_2 );

    df = freq_sampling( std::max( dim_1, dim_2 ), 0.5 / dr );

    for( int ii = 0; ii < 0.5 * ( dim_1 - 1 ) + 1; ++ii )
    {
        k_1 = ii * df;
        for( int jj = 0; jj < 0.5 * ( dim_2 - 1 ) + 1; ++jj )
        {
            k_2 = jj * df;

            arr( ii, jj ) = sqrt( k_1 * k_1 + k_2 * k_2 );

            if( k_x )
                ( *k_x )( ii, jj ) = k_1;
            if( k_x )
                ( *k_y )( ii, jj ) = k_2;
        }

        for( int jj = 0.5 * ( dim_2 - 1 ) + 1; jj < dim_2; ++jj )
        {
            k_2 = ( jj - dim_2 ) * df;

            arr( ii, jj ) = sqrt( k_1 * k_1 + k_2 * k_2 );

            if( k_x )
                ( *k_x )( ii, jj ) = k_1;
            if( k_x )
                ( *k_y )( ii, jj ) = k_2;
        }
    }

    for( int ii = 0.5 * ( dim_1 - 1 ) + 1; ii < dim_1; ++ii )
    {
        k_1 = ( ii - dim_1 ) * df;
        for( int jj = 0; jj < 0.5 * ( dim_2 - 1 ) + 1; ++jj )
        {
            k_2 = jj * df;

            arr( ii, jj ) = sqrt( k_1 * k_1 + k_2 * k_2 );

            if( k_x )
                ( *k_x )( ii, jj ) = k_1;
            if( k_x )
                ( *k_y )( ii, jj ) = k_2;
        }

        for( int jj = 0.5 * ( dim_2 - 1 ) + 1; jj < dim_2; ++jj )
        {
            k_2 = ( jj - dim_2 ) * df;

            arr( ii, jj ) = sqrt( k_1 * k_1 + k_2 * k_2 );

            if( k_x )
                ( *k_x )( ii, jj ) = k_1;
            if( k_x )
                ( *k_y )( ii, jj ) = k_2;
        }
    }
}

/// Create a frequency grid
template <typename eigenArr>
void frequencyGrid( eigenArr &arr, typename eigenArr::Scalar dt )
{
    frequencyGrid( arr, dt, (eigenArr *)0, (eigenArr *)0 );
}

/// Create a frequency grid
template <typename eigenArr>
void frequencyGrid( eigenArr &arr, typename eigenArr::Scalar dt, eigenArr &k_x, eigenArr &k_y )
{
    frequencyGrid( arr, dt, &k_x, &k_y );
}

/// Calculate the normalization for a 1-D @f$ 1/|f|^\alpha @f$ PSD.
/**
 * \param [in] fmin is the minimum non-zero absolute value of frequency
 * \param [in] fmax is the maximum absolute value of frequencey
 * \param [in] alpha is the power-law exponent, by convention @f$ \alpha > 0 @f$.
 *
 * \returns the normalization for a 2-sided power law PSD.
 *
 * \tparam realT is the real floating point type used for calculations.
 */
template <typename realT>
realT oneoverf_norm( realT fmin, realT fmax, realT alpha )
{
    realT integ = 2 * ( pow( fmax, -1.0 * alpha + 1.0 ) - pow( fmin, -1.0 * alpha + 1.0 ) ) / ( -1.0 * alpha + 1.0 );

    return 1 / integ;
}

/// Calculate the normalization for a 2-D @f$ 1/|k|^\alpha @f$ PSD.
/**
 * \param [in] kmin is the minimum non-zero absolute value of frequency
 * \param [in] kmax is the maximum absolute value of frequencey
 * \param [in] alpha is the power-law exponent, by convention @f$ \alpha > 0 @f$.
 *
 * \returns the normalization for a 2-D, 2-sided power law PSD.
 *
 * \tparam realT is the real floating point type used for calculations.
 */
template <typename realT>
realT oneoverk_norm( realT kmin, realT kmax, realT alpha )
{
    realT integ = 2 * ( pow( kmax, -1 * alpha + 2.0 ) - pow( kmin, -1.0 * alpha + 2.0 ) ) / ( -1 * alpha + 2.0 );

    return 1 / integ;
}

/// Normalize a 1-D PSD to have a given variance
/** A frequency range can be specified to calculate the norm, otherwise f[0] to f[f.size()-1] is the range.  The entire
 * PSD is normalized regardless.
 *
 * \tparam floatT the floating point type of the PSD.
 * \tparam floatParamT a floating point type, convenience to avoid double-float cofusion.
 *
 * \test Verify scaling and normalization of augment1SidedPSD \ref tests_sigproc_psdUtils_augment1SidedPSD "[test doc]"
 */
template <typename floatT, typename floatParamT>
int normPSD( std::vector<floatT> &psd, ///< [in.out] the PSD to normalize, will be altered.
             std::vector<floatT> &f,   ///< [in] the frequency points for the PSD
             floatParamT normT,        ///< [in] the desired total variance (or integral) of the PSD.
             floatT fmin = std::numeric_limits<floatT>::min(), ///< [in] [optiona] the minimum frequency of the range
                                                               ///< over which to normalize.
             floatT fmax = std::numeric_limits<floatT>::max()  ///< [in] [optiona] the maximum frequency of the range
                                                               ///< over which to normalize.
)
{
    floatT norm = normT;

    floatT s = 0; // accumulate

    // Check if inside-the-loop branch is needed
    if( fmin != std::numeric_limits<floatT>::min() || fmax != std::numeric_limits<floatT>::max() )
    {
        for( size_t i = 0; i < psd.size(); ++i )
        {
            if( fabs( f[i] ) < fmin || fabs( f[i] ) > fmax )
                continue;
            s += psd[i];
        }
    }
    else
    {
        for( size_t i = 0; i < psd.size(); ++i )
        {
            s += psd[i];
        }
    }

    s *= ( f[1] - f[0] );

    for( size_t i = 0; i < psd.size(); ++i )
        psd[i] *= norm / s;

    return 0;
}

/// Normalize a 2-D PSD to have a given variance
/** A frequency range can be specified for calculating the norm, otherwise the entire PSD is used.  The entire PSD is
 * normalized regardless.
 *
 * \tparam floatT the floating point type of the PSD.
 * \tparam floatParamT a floating point type, convenience to avoid double-float cofusion.
 *
 * \test Verify scaling and normalization of augment1SidedPSD \ref tests_sigproc_psdUtils_augment1SidedPSD "[test doc]"
 */
template <typename floatT, typename floatParamT>
floatT
normPSD( Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> &psd, ///< [in.out] the PSD to normalize, will be altered.
         Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> &k,   ///< [in] the frequency grid for psd.
         floatParamT normT, ///< [in] the desired total variance (or integral) of the PSD.
         floatT kmin = std::numeric_limits<floatT>::min(), ///< [in] [optiona] the minimum frequency of the range over
                                                           ///< which to normalize.
         floatT kmax = std::numeric_limits<floatT>::max()  ///< [in] [optiona] the maximum frequency of the range over
                                                           ///< which to normalize.
)
{
    floatT norm = normT;

    floatT dk1, dk2;

    if( k.rows() > 1 )
        dk1 = k( 1, 0 ) - k( 0, 0 );
    else
        dk1 = 1;

    if( k.cols() > 1 )
        dk2 = k( 0, 1 ) - k( 0, 0 );
    else
        dk2 = 1;

    floatT s = 0;

    // Check if inside-the-loop branch is needed
    if( kmin != std::numeric_limits<floatT>::min() || kmax != std::numeric_limits<floatT>::max() )
    {
        for( int c = 0; c < psd.cols(); ++c )
        {
            for( int r = 0; r < psd.rows(); ++r )
            {
                if( fabs( k( r, c ) ) < kmin || fabs( k( r, c ) ) > kmax )
                    continue;
                s += psd( r, c );
            }
        }
    }
    else
    {
        for( int c = 0; c < psd.cols(); ++c )
        {
            for( int r = 0; r < psd.rows(); ++r )
            {
                s += psd( r, c );
            }
        }
    }

    s *= dk1 * dk2;

    for( int c = 0; c < psd.cols(); ++c )
    {
        for( int r = 0; r < psd.rows(); ++r )
        {
            psd( r, c ) *= norm / s;
        }
    }

    return 0;
}

/// Generates a @f$ 1/|f|^\alpha @f$ power spectrum
/**
 * Populates an Eigen array  with
 * \f[
 *  P(|f| = 0) = 0
 * \f]
 * \f[
 *  P(|f| > 0) = \frac{\beta}{|f|^{\alpha}}
 * \f]
 *
 *
 * \param [out] psd is the array to populate
 * \param [in] freq is a frequency grid, must be the same logical size as psd
 * \param [in] alpha is the power law exponent, by convention @f$ alpha > 0 @f$.
 * \param [in] beta [optional is a normalization constant to multiply the raw spectrum by.  If beta==-1 (default) then
 *                           the PSD is normalized using \ref oneoverf_norm.
 *
 * \tparam eigenArrp is the Eigen-like array type of the psd
 * \tparam eigenArrf is the Eigen-like array type of the frequency grid
 */
template <typename eigenArrp, typename eigenArrf>
void oneoverf_psd( eigenArrp &psd,
                   eigenArrf &freq,
                   typename eigenArrp::Scalar alpha,
                   typename eigenArrp::Scalar beta = -1 )
{
    typedef typename eigenArrp::Scalar Scalar;

    typename eigenArrp::Index dim_1, dim_2;
    Scalar f_x, f_y, p;

    dim_1 = psd.rows();
    dim_2 = psd.cols();

    if( beta == -1 )
    {
        Scalar fmin;
        Scalar fmax;

        fmax = freq.abs().maxCoeff();

        // Find minimum non-zero Coeff.
        fmin = ( freq.abs() > 0 ).select( freq.abs(), freq.abs() + fmax ).minCoeff();

        beta = oneoverf_norm( fmin, fmax, alpha );
    }

    for( int ii = 0; ii < dim_1; ++ii )
    {
        for( int jj = 0; jj < dim_2; ++jj )
        {
            if( freq( ii, jj ) == 0 )
            {
                p = 0;
            }
            else
            {
                p = beta / std::pow( std::abs( freq( ii, jj ) ), alpha );
            }
            psd( ii, jj ) = p;
        }
    }
}

/// Generate a 1-D von Karman power spectrum
/**
 * Populates an Eigen array  with
 *
 * \f[
 *  P(f) = \frac{\beta}{ (f^2 + (1/T_0)^2)^{\alpha/2}} e^{ - f^2 t_0^2}
 * \f]
 *
 * If you set \f$ T_0 \le 0 \f$ and \f$ t_0 = 0\f$ this reverts to a simple \f$ 1/f^\alpha \f$ law (i.e.
 * it treats this as infinite outer scale and inner scale).
 *
 *
 * \tparam floatT a floating point
 */
template <typename floatT,
          typename floatfT,
          typename alphaT,
          typename T0T = double,
          typename t0T = double,
          typename betaT = double>
int vonKarmanPSD( std::vector<floatT> &psd, ///< [out] the PSD vector, will be resized.
                  std::vector<floatfT> &f,  ///< [in] the frequency vector
                  alphaT alpha,             ///< [in] the exponent, by convention @f$ alpha > 0 @f$.
                  T0T T0 = 0,               ///< [in] the outer scale, default is 0 (not used).
                  t0T t0 = 0,               ///< [in] the inner scale, default is 0 (not used).
                  betaT beta = 1            ///< [in] the scaling constant, default is 1
)
{

#ifndef MX_VKPSD_REFACT
    static_assert( 0 * std::is_floating_point<floatT>::value,
                   "the 1D vonKarmanPSD has been refactored.  After modifying your code to match, you must define "
                   "MX_VKPSD_REFACT before including psdUtils.hpp to avoid this error." );
#endif

    floatT T02;
    if( T0 > 0 )
        T02 = 1.0 / ( T0 * T0 );
    else
        T02 = 0;

    floatT sqrt_alpha = 0.5 * alpha;

    floatT _beta;
    if( beta <= 0 )
        _beta = 1;
    else
        _beta = beta;

    psd.resize( f.size() );

    for( size_t i = 0; i < f.size(); ++i )
    {
        floatT p = _beta / pow( pow( f[i], 2 ) + T02, sqrt_alpha );
        if( t0 > 0 )
            p *= exp( -1 * pow( f[i] * static_cast<floatT>( t0 ), 2 ) );
        psd[i] = p;
    }

    return 0;
}

/// Generate a 1-D "knee" PSD
/**
 * Populates an Eigen array  with
 *
 * \f[
 *  P(f) = \frac{\beta}{ 1 + (f/f_n)^{\alpha}}
 * \f]
 *
 * If you set \f$ T_0 \le 0 \f$ and \f$ t_0 = 0\f$ this reverts to a simple \f$ 1/f^\alpha \f$ law (i.e.
 * it treats this as infinite outer scale and inner scale).
 *
 * \tparam floatT a floating point
 */
template <typename floatT>
int kneePSD( std::vector<floatT> &psd, ///< [out] the PSD vector, will be resized.
             std::vector<floatT> &f,   ///< [in] the frequency vector
             floatT beta,              ///< [in] the scaling constant
             floatT fn,                ///< [in] the knee frequency
             floatT alpha              ///< [in] the exponent, by convention @f$ alpha > 0 @f$.
)
{

    psd.resize( f.size() );

    for( int i = 0; i < f.size(); ++i )
    {
        floatT p = beta / ( 1 + pow( f[i] / fn, alpha ) );
        psd[i] = p;
    }

    return 0;
}

/// Generates a von Karman power spectrum
/**
 * Populates an Eigen array  with
 *
 * \f[
 *  P(k) = \frac{\beta}{ (k^2 + (1/L_0)^2)^{\alpha/2}} e^{ - k^2 l_0^2}
 * \f]
 *
 * If you set \f$ L_0 \le 0 \f$ and \f$ l_0 = 0\f$ this reverts to a simple \f$ 1/f^\alpha \f$ law (i.e.
 * it treats this as infinite outer scale and inner scale).
 *
 * \param [out] psd is the array to populate, allocated.
 * \param [in] freq is a frequency grid, must be the same logical size as psd
 * \param [in] alpha is the power law exponent, by convention @f$ alpha > 0 @f$.
 * \param [in] L0 [optional] is the outer scale.
 * \param [in] l0 [optional] is the inner scale.
 * \param [in] beta [optional] is a normalization constant to multiply the raw spectrum by.  If beta==-1 (default) then
 *                           the PSD is normalized using \ref oneoverf_norm.
 *
 * \tparam eigenArrp is the Eigen array type of the psd
 * \tparam eigenArrf is the Eigen array type of the frequency grid
 */
template <typename eigenArrp, typename eigenArrf, typename alphaT, typename L0T, typename l0T, typename betaT>
void vonKarmanPSD( eigenArrp &psd, eigenArrf &freq, alphaT alpha, L0T L0 = 0, l0T l0 = 0, betaT beta = -1 )
{
    typedef typename eigenArrp::Scalar Scalar;

    typename eigenArrp::Index dim_1, dim_2;
    Scalar p;

    dim_1 = psd.rows();
    dim_2 = psd.cols();

    Scalar _beta;

    if( beta == -1 )
    {
        Scalar fmin;
        Scalar fmax;

        fmax = freq.abs().maxCoeff();

        // Find minimum non-zero Coeff.
        fmin = ( freq.abs() > 0 ).select( freq.abs(), freq.abs() + fmax ).minCoeff();

        _beta = beta = oneoverf_norm( fmin, fmax, static_cast<Scalar>( alpha ) );
    }
    else
        _beta = static_cast<Scalar>( beta );

    Scalar L02;
    if( L0 > 0 )
        L02 = 1.0 / ( L0 * L0 );
    else
        L02 = 0;

    Scalar sqrt_alpha = 0.5 * alpha; // std::sqrt(alpha);

    for( int ii = 0; ii < dim_1; ++ii )
    {
        for( int jj = 0; jj < dim_2; ++jj )
        {
            if( freq( ii, jj ) == 0 && L02 == 0 )
            {
                p = 0;
            }
            else
            {
                p = _beta / pow( pow( freq( ii, jj ), 2 ) + L02, sqrt_alpha );
                if( l0 > 0 )
                    p *= exp( -1 * pow( freq( ii, jj ) * static_cast<Scalar>( l0 ), 2 ) );
            }
            psd( ii, jj ) = p;
        }
    }
}

/// Augment a 1-sided PSD to standard 2-sided FFT form.
/** Allocates psdTwoSided to hold a flipped copy of psdOneSided.
 * Default assumes that psdOneSided[0] corresponds to 0 frequency,
 * but this can be changed by setting zeroFreq to a non-zero value.
 * In this case psdTwoSided[0] is set to 0, and the augmented psd
 * is shifted by 1.
 *
 * To illustrate, the bins are re-ordered as:
 * \verbatim
 * {1,2,3,4,5} --> {0,1,2,3,4,5,-4,-3,-2,-1}
 * \endverbatim
 *
 * The output is scaled so that the total power remains the same.  The 0-freq and
 * Nyquist freq are not scaled.
 *
 *
 * Entries in psdOneSided are cast to the value_type of psdTwoSided,
 * for instance to allow for conversion to complex type.
 *
 * \test Verify scaling and normalization of augment1SidedPSD \ref tests_sigproc_psdUtils_augment1SidedPSD "[test doc]"
 */
template <typename vectorTout, typename vectorTin>
void augment1SidedPSD(
    vectorTout &psdTwoSided, ///< [out] on return contains the FFT storage order copy of psdOneSided.
    vectorTin &psdOneSided,  ///< [in] the one-sided PSD to augment
    bool addZeroFreq =
        false, ///< [in] [optional] set to true if psdOneSided does not contain a zero frequency component.
    typename vectorTin::value_type scale = 0.5 ///< [in] [optional] value to scale the input by when copying to the
                                               ///< output.  The default 0.5 re-normalizes for a 2-sided PSD.
)
{
    typedef typename vectorTout::value_type outT;

    bool needZero = 1;

    size_t N;

    if( addZeroFreq == 0 )
    {
        needZero = 0;
        N = 2 * psdOneSided.size() - 2;
    }
    else
    {
        N = 2 * psdOneSided.size();
    }

    psdTwoSided.resize( N );

    // First set the 0-freq point
    if( needZero )
    {
        psdTwoSided[0] = outT( 0.0 );
    }
    else
    {
        psdTwoSided[0] = outT( psdOneSided[0] );
    }

    // Now set all the rest.
    unsigned long i;
    for( i = 0; i < psdOneSided.size() - 1 - ( 1 - needZero ); ++i )
    {
        psdTwoSided[i + 1] = outT( psdOneSided[i + ( 1 - needZero )] * scale );
        psdTwoSided[i + psdOneSided.size() + needZero] = outT( psdOneSided[psdOneSided.size() - 2 - i] * scale );
    }
    psdTwoSided[i + 1] = outT( psdOneSided[i + ( 1 - needZero )] );
}

/// Augment a 1-sided frequency scale to standard FFT form.
/** Allocates freqTwoSided to hold a flipped copy of freqOneSided.
 * If freqOneSided[0] is not 0, freqTwoSided[0] is set to 0, and the augmented
 * frequency scale is shifted by 1.
 *
 * Example:
 *
 * {1,2,3,4,5} --> {0,1,2,3,4,5,-4,-3,-2,-1}
 *
 * \test Verify scaling and normalization of augment1SidedPSD \ref tests_sigproc_psdUtils_augment1SidedPSD "[test doc]"
 */
template <typename T>
void augment1SidedPSDFreq(
    std::vector<T> &freqTwoSided, ///< [out] on return contains the FFT storage order copy of freqOneSided.
    std::vector<T> &freqOneSided  ///< [in] the one-sided frequency scale to augment
)
{
    int needZero = 1;

    size_t N;

    if( freqOneSided[0] != 0 )
    {
        N = 2 * freqOneSided.size();
    }
    else
    {
        needZero = 0;
        N = 2 * freqOneSided.size() - 2;
    }

    freqTwoSided.resize( N );

    if( needZero )
    {
        freqTwoSided[0] = 0.0;
    }
    else
    {
        freqTwoSided[0] = freqOneSided[0]; // 0
    }

    int i;
    for( i = 0; i < freqOneSided.size() - 1 - ( 1 - needZero ); ++i )
    {
        freqTwoSided[i + 1] = freqOneSided[i + ( 1 - needZero )];
        freqTwoSided[i + freqOneSided.size() + needZero] = -freqOneSided[freqOneSided.size() - 2 - i];
    }
    freqTwoSided[i + 1] = freqOneSided[i + ( 1 - needZero )];
}

/// Rebin a PSD, including its frequency scale, to a larger frequency bin size (fewer bins)
/** The rebinning uses trapezoid integration within bins to ensure minimum signal loss.
 *
 * Maintains DFT sampling.  That is, if initial frequency grid is 0,0.1,0.2...
 * and the binSize is 1.0, the new grid will be 0,1,2 (as opposed to 0.5, 1.5, 2.5).
 *
 * This introduces a question of what to do with first half-bin, which includes 0.  It can be
 * integrated (binAtZero = true, the default).  This may cause inaccurate behavior if the value of the PSD when
 * f=0 is important (e.g. when analyzing correlated noise), so setting binAtZero=false causes the f=0 value to be
 * copied (using the nearest neighbor if no f=0 point is in the input.
 *
 * The last half bin is always integrated.
 *
 * The output is variance normalized to match the input variance.
 *
 * \tparam realT  the real floating point type
 */
template <typename realT>
int rebin1SidedPSD( std::vector<realT> &binFreq, ///< [out] the binned frequency scale, resized.
                    std::vector<realT> &binPSD,  ///< [out] the binned PSD, resized.
                    std::vector<realT> &freq,    ///< [in] the frequency scale of the PSD to bin.
                    std::vector<realT> &PSD,     ///< [in] the PSD to bin.
                    realT binSize,               ///< [in] in same units as freq
                    bool binAtZero = true ///< [in] [optional] controls whether the zero point is binned or copied.
)
{
    binFreq.clear();
    binPSD.clear();

    realT sumPSD = 0;
    realT startFreq = 0;
    realT sumFreq = 0;
    int nSum = 0;

    int i = 0;

    realT df = freq[1] - freq[0];

    // Now move to first bin
    while( freq[i] <= 0.5 * binSize + 0.5 * df )
    {
        sumPSD += PSD[i];
        ++nSum;
        ++i;
        if( i >= freq.size() )
            break;
    }

    if( !binAtZero )
    {
        binFreq.push_back( 0 );
        binPSD.push_back( PSD[0] );
    }
    else
    {
        binFreq.push_back( 0 );
        binPSD.push_back( sumPSD / nSum );
    }

    --i;
    startFreq = freq[i];
    nSum = 0;
    sumFreq = 0;
    sumPSD = 0;

    while( i < freq.size() )
    {
        realT sc = 0.5; // First one is multiplied by 1/2 for trapezoid rule.
        while( freq[i] - startFreq + 0.5 * df < binSize )
        {
            sumFreq += freq[i];
            sumPSD += sc * PSD[i];
            sc = 1.0;
            ++nSum;

            ++i;
            if( i >= freq.size() - 1 )
                break; // break 1 element early so last point is mult by 0.5
        }

        if( i < freq.size() )
        {
            sumFreq += freq[i];
            sumPSD += 0.5 * PSD[i]; // last one is multiplied by 1/2 for trapezoid rule
            ++nSum;
            ++i;
        }

        // Check if this is last point
        if( i < freq.size() )
        {
            binFreq.push_back( sumFreq / nSum );
        }
        else
        {
            // last point frequencyis not averaged.
            binFreq.push_back( freq[freq.size() - 1] );
        }

        binPSD.push_back( sumPSD / ( nSum - 1 ) );

        sumFreq = 0;
        sumPSD = 0;
        nSum = 0;
        if( i >= freq.size() )
            break;

        --i; // Step back one, so averages are edge to edge.

        startFreq = freq[i];
    }

    // Now normalize variance
    realT var = psdVar( freq, PSD );
    realT binv = psdVar( binFreq, binPSD );

    for( int i = 0; i < binFreq.size(); ++i )
        binPSD[i] *= var / binv;

    return 0;
}
///@}

} // namespace sigproc
} // namespace mx

#endif // psdUtils_hpp
