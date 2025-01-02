/** \file mftT.hpp
 * \brief The Matrix Fourier Transform interface
 * \ingroup ft_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2024-2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef mftT_hpp
#define mftT_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include "ft.hpp"
#include "../../math/constants.hpp"

/** \addtogroup mft
 *
 * The Matrix Fourier Transform is the <a href="https://en.wikipedia.org/wiki/DFT_matrix">
 * matrix-multiplication implementation of the
 * Discrete Fourier Transform</a>.  It is slower than the FFT given the same size matrix, e.g.
 * <a href="https://pages.hmc.edu/ruye/e161/lectures/fourier/node11.html">in 2D</a>
 * \f$O(N^3)\f$ vs \f$O(N^2\ln(N))\f$.  However, compared to zero-padding,
 * it can provide advantages in both speed
 * and memory required for problems where oversampling is needed but only over a small region
 * of the output.
 */
namespace mx
{
namespace math
{
namespace ft
{

template <typename _inputT, typename _outputT, size_t _rank, int _cudaGPU = 0>
class mftT;

/// Matrix Fourier Transforms
/** Calculates the Discrete Fourier Transform (DFT) using matrix multiplication.  This
 *  is normally much less efficient than the FFT, but for large oversampling (padding)
 *  the Matrix FT (MFT) will be much more space efficient and faster at the cost of
 *  "field of view".
 *
 *  This interface is modeled after the \ref fftT<_inputT, _outputT, _rank, 0> "fftT" interface to fftw, and is interoperable with it.
 *  That means that using, e.g., mftT (unshifted, not-oversampled) for the forward transform
 *  and \ref fftT<_inputT, _outputT, _rank, 0> "fftT" for the backward
 *  transform is equivalent to using one or the other for both. Note that
 *  the MFT is normalized as in fftw, so by 1/N on the forward and by 1 on the backward.
 *
 *  Oversampling is optionally included in the transform, which is equivalent to zero padding
 *  in terms of increased resolution \cite soummer_2007.
 *  The cost is that the output domain is truncated by the oversampling factor.  I.e. if we
 *  oversample by a factor of 10, only 1/10th of the output transform is available.
 *
 *  The output can be shifted as part of the MFT calculation, which is implemented similarly
 *  to <a href="https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation">
 *  this matlab code</a>.
 *
 *  Note that when either oversampling or shifting is done on a forward (backward) transform,
 *  a subsequent backward (forward) transform will not in general be the inverse.
 *
 *  \tparam inputT is the input type of the transform, only complex types are suppored by mftT
 *  \tparam outputT is the output type of the transform, only complex types are suppored by mftT
 *  \tparam _rank is the rank of the transform. Currently only rank 2 is implemented.
 *
 *  \ingroup mft
 */
template <typename _inputT, typename _outputT, size_t _rank>
class mftT<_inputT, _outputT, _rank, 0>
{
    typedef _inputT inputT;
    typedef _outputT outputT;

    typedef typename fftwTypeSpec<inputT, outputT>::realT realT;

    typedef typename fftwTypeSpec<_inputT, _outputT>::complexT complexT;

    static const size_t rank = _rank;

    typedef Eigen::Array<inputT, -1, -1> eigenArrayInT;
    typedef Eigen::Array<outputT, -1, -1> eigenArrayOutT;

    typedef Eigen::Matrix<complexT, -1, -1> eigenMatrixT;

  protected:
    dir m_dir{ dir::forward }; /**< Direction of this MFT, either dir::forward (default)
                                     or dir::backward */

    int m_szX{ 0 }; ///< Size of the x dimension
    int m_szY{ 0 }; ///< Size of the y dimension
    int m_szZ{ 0 }; ///< size of the z dimension

    float m_osFac{ 1 }; ///< The oversampling factor

    realT m_xOff{ 0 }; ///< The offset in the rows direction for the center of the DFT.
    realT m_yOff{ 0 }; ///< The offset in the columns direction for the center of the DFT.
  public:
    eigenMatrixT m_dftR; ///< DFT matrix for the rows
    eigenMatrixT m_dftC; ///< DFT matrix for the columnss

  public:
    /// Default c'tor
    mftT();

    /// Constructor for rank 1 MFT.
    template <size_t crank = _rank>
    mftT( int nx,                   ///< [in] the desired size of the MFT
          dir ndir = dir::forward, /**< [in] [optional] direction of this MFT, either dir::forward
                                                         (default) or dir::backward */
          realT xOff = 0,           /**< [in] [optional] the x offset of the center of the
                                                         transformed array. Default 0.*/
          realT osFac = 1.0,        /**< [in] [optional] the oversampling factor.  Default 1. */
          typename std::enable_if<crank == 1>::type * = 0 ) = delete;

    /// Constructor for rank 2 MFT.
    template <size_t crank = _rank>
    mftT( int nx,                   ///< [in] the desired x size of the MFT
          int ny,                   ///< [in] the desired y size of the MFT
          dir ndir = dir::forward, /**< [in] [optional] direction of this MFT, either dir::forward
                                                         (default) or dir::backward */
          realT xOff = 0,           /**< [in] [optional] the x offset of the center of the
                                                         transformed array. Default 0.*/
          realT yOff = 0,           /**< [in] [optional] the y offset of the center of the
                                                         transformed array  Default 0.*/
          realT osFac = 1.0,        /**< [in] [optional] the oversampling factor.  Default 1. */
          typename std::enable_if<crank == 2>::type * = 0 );

    /// Constructor for rank 3 MFT.
    template <size_t crank = _rank>
    mftT( int nx,                   ///< [in] the desired x size of the MFT
          int ny,                   ///< [in] the desired y size of the MFT
          int nz,                   ///< [in] the desired z size of the MFT
          dir ndir = dir::forward, /**< [in] [optional] direction of this MFT, either dir::forward
                                                         (default) or dir::backward */
          realT xOff = 0,           /**< [in] [optional] the x offset of the center of the
                                                         transformed array. Default 0.*/
          realT yOff = 0,           /**< [in] [optional] the y offset of the center of the
                                                         transformed array  Default 0.*/
          realT zOff = 0,           /**< [in] [optional] the z offset of the center of the
                                                         transformed array  Default 0.*/
          realT osFac = 1.0,        /**< [in] [optional] the oversampling factor.  Default 1. */
          typename std::enable_if<crank == 3>::type * = 0 ) = delete;

    /// Planning routine for rank 2 transforms.
    template <size_t crank = _rank>
    void plan( int nx,                   ///< [in] the desired x size of the MFT
               int ny,                   ///< [in] the desired y size of the MFT
               dir ndir = dir::forward, /**< [in] [optional] direction of this MFT, either dir::forward
                                                              (default) or dir::backward */
               realT xOff = 0,           /**< [in] [optional] the x offset of the center of the
                                                              transformed array. Default 0.*/
               realT yOff = 0,           /**< [in] [optional] the y offset of the center of the
                                                              transformed array  Default 0.*/
               realT osFac = 1.0,        /**< [in] [optional] the oversampling factor.  Default 1. */
               typename std::enable_if<crank == 2>::type * = 0 );

    /// Conduct the MFT
    void operator()( eigenArrayOutT &out, /**< [out] the output of the DFT */
                     const eigenArrayInT &in /**< [in] the input to the DFT */ ) const;
};

template <typename inputT, typename outputT, size_t rank>
mftT<inputT, outputT, rank, 0>::mftT()
{
}

template <typename inputT, typename outputT, size_t rank>
template <size_t crank>
mftT<inputT, outputT, rank, 0>::mftT(
    int nx, int ny, dir ndir, realT xoff, realT yoff, realT osFac, typename std::enable_if<crank == 2>::type * )
{
    plan( nx, ny, ndir, xoff, yoff, osFac );
}

template <typename inputT, typename outputT, size_t rank>
template <size_t crank>
void mftT<inputT, outputT, rank, 0>::plan(
    int nx, int ny, dir ndir, realT xOff, realT yOff, realT osFac, typename std::enable_if<crank == 2>::type * )
{
    m_szX = nx;
    m_szY = ny;
    m_dir = ndir;
    m_xOff = xOff;
    m_yOff = yOff;
    m_osFac = osFac;

    if( m_szX != m_szY )
    {
        throw std::invalid_argument( "MFT of non-square size is not implemented. nx must equal ny." );
    }

    // These should depend on szX too.
    m_dftR.resize( m_szX, m_szX );
    m_dftC.resize( m_szX, m_szX );

    // There is probably an osNx and an osNy?
    realT osN = m_szX * m_osFac;

    if( m_dir == dir::forward )
    {
        realT sign = -1;
        realT norm = 1.0 / ( m_szX * m_szY );

        for( int cc = 0; cc < m_szY; ++cc )
        {
            realT ccx = cc;
            // if(ccx > m_szY/2) ccx = -1*(m_szY - ccx);

            for( int rr = 0; rr < m_szX; ++rr )
            {
                realT x = static_cast<realT>( rr - m_xOff ) * ccx / osN;

                m_dftR( rr, cc ) = norm * exp( complexT( { 0, sign * 2 * pi<realT>() * x } ) );

                realT rrx = rr;
                // if(rrx > m_szX/2) rrx = -1*(m_szX - rrx);

                x = rrx * static_cast<realT>( cc - m_yOff ) / osN;

                m_dftC( rr, cc ) = norm * exp( complexT( { 0, sign * 2 * pi<realT>() * x } ) );
            }
        }
    }
    else
    {
        realT sign = +1;
        realT norm = 1.0;

        for( int cc = 0; cc < m_szY; ++cc )
        {
            realT ccx = cc;
            if( ccx > m_szY / 2 )
                ccx = -1 * ( m_szY - ccx );

            for( int rr = 0; rr < m_szX; ++rr )
            {
                realT rrx = rr;
                if( rrx > m_szX / 2 )
                    rrx = -1 * ( m_szX - rrx );

                realT x = rrx * static_cast<realT>( cc - m_xOff ) / osN;

                m_dftR( cc, rr ) = norm * exp( complexT( { 0, sign * 2 * pi<realT>() * x } ) );

                x = static_cast<realT>( rr - m_yOff ) * ccx / osN;

                m_dftC( cc, rr ) = norm * exp( complexT( { 0, sign * 2 * pi<realT>() * x } ) );
            }
        }
    }
}

template <typename inputT, typename outputT, size_t rank>
void mftT<inputT, outputT, rank, 0>::operator()( eigenArrayOutT &out, const eigenArrayInT &in ) const
{
    out = ( m_dftR * in.matrix() * m_dftC ).array();
}

} // namespace ft
} // namespace math
} // namespace mx

#endif // mdft_hpp
