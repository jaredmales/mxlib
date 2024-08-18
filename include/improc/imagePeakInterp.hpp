/** \file imagePeakInterp.hpp
 * \brief A class to find the location of a peak using interpolation
 * \ingroup image_processing_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2023 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __imagePeakInterp_hpp__
#define __imagePeakInterp_hpp__

#include "eigenImage.hpp"
#include "imageTransforms.hpp"

// #define ICCS_OMP
namespace mx
{
namespace improc
{

/// Find the peak of an image using interpolation
/** Interpolates onto a finer grid according to m_tol
 *
 * \tparam transformT is a transformation type
 *
 * \ingroup image_reg
 */
template <typename transformT>
struct imagePeakInterp
{

    typedef typename transformT::arithT realT;

    transformT m_transform;

    realT m_tol{ 0.1 };

    eigenImage<realT> m_magIm;

    imagePeakInterp()
    {
    }

    imagePeakInterp( realT tol ) : m_tol( tol )
    {
    }

    void operator()( realT &x, realT &y, eigenImage<realT> &im )
    {
        int r = ( 1.0 * im.rows() ) / m_tol + 1;
        int c = ( 1.0 * im.cols() ) / m_tol + 1;

        m_magIm.resize( r, c );

        imageMagnify( m_magIm, im, m_transform );

        int nx, ny;
        m_magIm.maxCoeff( &nx, &ny );

        x = nx * m_tol;
        y = ny * m_tol;
    }
};

} // namespace improc
} // namespace mx

#endif //__imagePeakInterp_hpp__
