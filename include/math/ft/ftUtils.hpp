/** \file ftUtils.hpp
 * \brief Fourier Transform Utilities
 * \ingroup ft_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef ftUtils_hpp
#define ftUtils_hpp

namespace mx
{
namespace math
{
namespace ft
{

/// Add padding to the output of a real-to-complex (R2C) FFT
/**
 * The output array should be pre-sized according to:
 * \verbatim
 * aout.cols() == ain.cols() + padCols
 * aout.rows() == ain.rows() + padRows
 * \endverbatim
 * where `padRows` and `padCols` are the desired amount of padding. \p aout will be in
 * FFT order and suitable for input to the inverse transform (complext-to-real).
 *
 * \ingroup ft
 */
template <class eigenArrOutT, class eigenArrInT>
void padR2CFFTOutput( eigenArrOutT &aout, /**< [out] The padded array.  Must be pre-allocated. The pad region will
                                                    not be set to zero.  Size of the pad is determined by the
                                                    difference in size between \p aout and \p ain*/
                      eigenArrInT &ain /**< [in] The array to be paddded, in real-to-complex FFTW storage order.*/ )
{
    int ci = 0.5 * ain.cols() + 1;
    int cf = ain.cols() - ci;

    aout.block( 0, 0, ain.rows(), ci ) = ain.block( 0, 0, ain.rows(), ci );
    aout.block( 0, aout.cols() - cf, ain.rows(), cf ) = ain.block( 0, ain.cols() - cf, ain.rows(), cf );
}

/// Add padding to the output of a complex-to-complex (C2C) FFT
/** The output array should be pre-sized according to:
 * \verbatim
 * aout.row() ==  ain.rows() + padRows
 * aout.cols() == ain.cols() + padCols
 * \endverbatim
 * where `padRows` and `padCols` are the desired amount of padding. \p aout will be in
 * FFT order and suitable for input to the inverse transform (complex-to-complex).
 *
 * \ingroup ft
 */
template <class eigenArrOutT, class eigenArrInT>
void padC2CFFTOutput( eigenArrOutT &aout, /**< [out] The padded array.  Must be pre-allocated. The pad region will
                                                     not be set to zero.  Size of the pad is determined by the
                                                     difference in size between \p aout and \p ain*/
                      eigenArrInT &ain /**< [in] The array to be paddded, in complex-to-complex
                                                 FFTW storage order.*/ )
{
    int ri = 0.5 * ain.rows() + 1;
    int ci = 0.5 * ain.cols() + 1;

    int rf = ain.rows() - ri;
    int cf = ain.cols() - ci;

    aout.topLeftCorner( ri, ci ) = ain.topLeftCorner( ri, ci );
    aout.bottomLeftCorner( rf, ci ) = ain.bottomLeftCorner( rf, ci );
    aout.topRightCorner( ri, cf ) = ain.topRightCorner( ri, cf );
    aout.bottomRightCorner( rf, cf ) = ain.bottomRightCorner( rf, cf );
}

/// Augment and pad the output of a real-to-complex (R2C) forward FFT for use with a complex-to-complex (C2C) inverse
/// FFT
/** The output array should be pre-sized according to:
 *  \verbatim
 *   aout.row() == 2*(ain.rows()-1) + padRows
 *   aout.cols() == ain.cols() + padCols
 * \endverbatim
 * where padRows and padCols can be 0. \p aout will be in
 * FFT order and suitable for input to the inverse transform (complex-to-complex).
 *
 * \ingroup ft
 */
template <class eigenArrOutT, class eigenArrInT>
void augmentR2CFFTOutput( eigenArrOutT &aout, eigenArrInT &ain )
{

    int padRows = aout.rows() - ain.rows();
    int padCols = aout.cols() - ain.cols();

    // aout.resize( 2 * ( ain.rows() - 1 ) + padRows, ain.cols() + padCols );

    int r1 = ain.rows();
    int c1 = 0.5 * ain.cols() + 1;

    int r2 = aout.rows() - ( r1 + padRows );
    int c2 = aout.cols() - ( c1 + padCols );

    aout.block( 0, 0, r1, c1 ) = ain.block( 0, 0, r1, c1 );

    aout.block( 0, c1 + padCols, r1, c2 ) = ain.block( 0, c1, r1, c2 );

    // Upper Right corner
    for( int cc = 1; cc < c1 - 1; ++cc )
    {
        for( int rr = 1; rr < r1 - 1; ++rr )
        {
            aout( aout.rows() - rr, aout.cols() - cc ) = conj( ain( rr, cc ) );
        }
    }

    // Lower Right corner
    for( int rr = 1; rr < ain.rows() - 1; ++rr )
    {
        aout( aout.rows() - rr, 0 ) = conj( ain( rr, 0 ) );

        // After ain.rows()+1, for column 1 and greater, the columns are reversed and conjugate
        aout.block( aout.rows() - rr, 1, 1, aout.cols() - padCols - c1 + 1 ) =
            conj( ain.block( rr, c1 - 1, 1, ain.cols() - c1 + 1 ).reverse() );
    }
}

} // namespace ft
} // namespace math
} // namespace mx

#endif // ftUtils_hpp
