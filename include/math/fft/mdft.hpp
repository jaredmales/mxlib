/** \file mdft.hpp
 * \brief The Matrix Discrete Fourier Transform interface
 * \ingroup fft_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2024 Jared R. Males (jaredmales@gmail.com)
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

#ifndef mdft_hpp
#define mdft_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include "fft.hpp"
#include "../../math/constants.hpp"

namespace mx
{
namespace math
{
namespace fft
{

template <typename _inputT, typename _outputT, size_t _rank, int _cudaGPU = 0>
class mdftT;

/// Matrix Discrete Fourier Transforms
/** Calculates the Discrete Fourier Transform (DFT) using matrix multiplication.  This
 *  is normally much less efficient than the FFT, but for large oversampling (padding)
 *  the Matrix DFT (MDFT) will be much more space efficient and faster at the cost of
 *  "field of view".
 *
 *  This interface is modeled after the fftT interface to fftw.  Note that by default
 *  these MDFTs are normalized as in fftw, so 1/N on the forward.
 *
 *  \tparam inputT is the input type of the transform, can be either real or complex
 *  \tparam outputT is the output type of the transform, can be either real or complex
 *
 *  \ingroup fft
 */
template <typename _inputT, typename _outputT, size_t _rank>
class mdftT<_inputT, _outputT, _rank, 0>
{
    typedef _inputT inputT;
    typedef _outputT outputT;

    typedef typename fftwTypeSpec<inputT, outputT>::realT realT;

    typedef typename fftwTypeSpec<_inputT, _outputT>::complexT complexT;

    static const size_t rank = _rank;

    typedef Eigen::Array<inputT, -1, -1> eigenArrayInputT;
    typedef Eigen::Array<outputT, -1, -1> eigenArrayOutputT;

    typedef Eigen::Matrix<complexT, -1, -1> eigenMatrixT;

  protected:
    int m_dir{ MXFFT_FORWARD }; ///< Direction of this MDFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD

    int m_szX{ 0 }; ///< Size of the x dimension
    int m_szY{ 0 }; ///< Size of the y dimension
    int m_szZ{ 0 }; ///< size of the z dimension

    float m_osFac{ 1 }; ///< The oversampling factor

    realT m_xOff{ 0 }; ///< The offset in the rows direction for the center of the DFT.
    realT m_yOff{ 0 }; ///< The offset in the columns direction for the center of the DFT.

    eigenMatrixT m_dftR; ///< DFT matrix for the rows
    eigenMatrixT m_dftC; ///< DFT matrix for the columnss

  public:
    /// Default c'tor
    mdftT();

    /// Constructor for rank 1 MDFT.
    template <size_t crank = _rank>
    mdftT( int nx,                   ///< [in] the desired size of the MDFT
           int ndir = MXFFT_FORWARD, /**< [in] [optional] direction of this MDFT, either MXFFT_FORWARD (default)
                                                          or MXFFT_BACKWARD */
           realT xOff = 0,    /**< [in] [optional] the x offset of the center of the transformed array. Default 0.*/
           realT osFac = 1.0, /**< [in] [optional] the oversampling factor.  Default 1. */
           typename std::enable_if<crank == 1>::type * = 0 ) = delete;

    /// Constructor for rank 2 MDFT.
    template <size_t crank = _rank>
    mdftT( int nx,                   ///< [in] the desired x size of the MDFT
           int ny,                   ///< [in] the desired y size of the MDFT
           int ndir = MXFFT_FORWARD, /**< [in] [optional] direction of this MDFT, either MXFFT_FORWARD (default)
                                                          or MXFFT_BACKWARD */
           realT xOff = 0,    /**< [in] [optional] the x offset of the center of the transformed array. Default 0.*/
           realT yOff = 0,    /**< [in] [optional] the y offset of the center of the transformed array  Default 0.*/
           realT osFac = 1.0, /**< [in] [optional] the oversampling factor.  Default 1. */
           typename std::enable_if<crank == 2>::type * = 0 );

    /// Constructor for rank 3 MDFT.
    template <size_t crank = _rank>
    mdftT( int nx,                   ///< [in] the desired x size of the MDFT
           int ny,                   ///< [in] the desired y size of the MDFT
           int nz,                   ///< [in] the desired z size of the MDFT
           int ndir = MXFFT_FORWARD, /**< [in] [optional] direction of this MDFT, either MXFFT_FORWARD (default)
                                                          or MXFFT_BACKWARD */
           realT xOff = 0,    /**< [in] [optional] the x offset of the center of the transformed array. Default 0.*/
           realT yOff = 0,    /**< [in] [optional] the y offset of the center of the transformed array  Default 0.*/
           realT zOff = 0,    /**< [in] [optional] the z offset of the center of the transformed array  Default 0.*/
           realT osFac = 1.0, /**< [in] [optional] the oversampling factor.  Default 1. */
           typename std::enable_if<crank == 3>::type * = 0 ) = delete;

    /// Planning routine for rank 2 transforms.
    template <size_t crank = _rank>
    void plan( int nx,                   ///< [in] the desired x size of the MDFT
               int ny,                   ///< [in] the desired y size of the MDFT
               int ndir = MXFFT_FORWARD, /**< [in] [optional] direction of this MDFT, either MXFFT_FORWARD (default)
                                                              or MXFFT_BACKWARD */
               realT xOff = 0,    /**< [in] [optional] the x offset of the center of the transformed array. Default 0.*/
               realT yOff = 0,    /**< [in] [optional] the y offset of the center of the transformed array  Default 0.*/
               realT osFac = 1.0, /**< [in] [optional] the oversampling factor.  Default 1. */
               typename std::enable_if<crank == 2>::type * = 0 );

    /// Conduct the DFT
    void operator()( eigenArrayOutputT &out, /**< [out] the output of the DFT, must be pre-allocated */
                     const eigenArrayInputT &in /**< [in] the input to the DFT */ ) const;
};

template <typename inputT, typename outputT, size_t rank>
mdftT<inputT, outputT, rank, 0>::mdftT()
{
}

template <typename inputT, typename outputT, size_t rank>
template <size_t crank>
mdftT<inputT, outputT, rank, 0>::mdftT(
    int nx, int ny, int ndir, realT xoff, realT yoff, realT osFac, typename std::enable_if<crank == 2>::type * )
{
    plan( nx, ny, ndir, xoff, yoff, osFac );
}

template <typename inputT, typename outputT, size_t rank>
template <size_t crank>
void mdftT<inputT, outputT, rank, 0>::plan(
    int nx, int ny, int ndir, realT xOff, realT yOff, realT osFac, typename std::enable_if<crank == 2>::type * )
{
    m_szX = nx;
    m_szY = ny;
    m_dir = ndir;
    m_xOff = xOff;
    m_yOff = yOff;
    m_osFac = osFac;

    if( m_szX != m_szY )
    {
        throw std::invalid_argument( "MDFT of non-square size is not implemented. nx must equal ny." );
    }

    // These should depend on szX too.
    m_dftR.resize( m_szX, m_szX );
    m_dftC.resize( m_szX, m_szX );

    // There is probably an osNx and an osNy?
    int osN = ceil( m_szX * m_osFac );

    if( m_dir == MXFFT_FORWARD )
    {
        realT norm = 1.0 / ( m_szX * m_szY );
        for( int cc = 0; cc < m_szY; ++cc )
        {
            for( int rr = 0; rr < m_szX; ++rr )
            {
                realT x = static_cast<realT>( rr - m_xOff ) * static_cast<realT>( cc );

                m_dftR( rr, cc ) = norm * exp( complexT( { 0, -2 * pi<realT>() * x / ( osN ) } ) );

                x = static_cast<realT>( rr ) * static_cast<realT>( cc - m_yOff );

                m_dftC( rr, cc ) = norm * exp( complexT( { 0, -2 * pi<realT>() * x / ( osN ) } ) );
            }
        }
    }
    else
    {
        realT norm = 1.0;

        for( int rr = 0; rr < m_szX; ++rr )
        {
            for( int cc = 0; cc < m_szY; ++cc )
            {
                realT x = static_cast<realT>( rr - m_xOff ) * static_cast<realT>( cc );

                m_dftR( cc, rr ) = norm * exp( complexT( { 0, +2 * pi<realT>() * x / ( osN ) } ) );

                x = static_cast<realT>( rr ) * static_cast<realT>( cc - m_yOff );

                m_dftC( cc, rr ) = norm * exp( complexT( { 0, +2 * pi<realT>() * x / ( osN ) } ) );
            }
        }
    }
}

template <typename inputT, typename outputT, size_t rank>
void mdftT<inputT, outputT, rank, 0>::operator()( eigenArrayOutputT &out, const eigenArrayInputT &in ) const
{
    out = ( m_dftR * in.matrix() * m_dftC ).array();
}

} // namespace fft
} // namespace math
} // namespace mx

#endif // mdft_hpp
