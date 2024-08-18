/** \file array2FitGaussian1D.hpp
 * \author Jared R. Males
 * \brief Wrapper for a native array to pass to \ref levmarInterface, with 1D Gaussian details.
 * \ingroup fitting_files
 *
 */

//***********************************************************************//
// Copyright 2022 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_fit_array2FitGaussian1D_hpp
#define math_fit_array2FitGaussian1D_hpp

#include "../../mxError.hpp"

namespace mx
{
namespace math
{
namespace fit
{

/// Wrapper for a native array to pass to \ref levmarInterface, with !D Gaussian details.
/** Supports fixing G0, G, x0, and sigma independently.
 * \ingroup gaussian_peak_fit
 */
template <typename realT>
struct array2FitGaussian1D
{
    realT *m_data{ nullptr };   ///< ///< Pointer to the array of y values
    realT *m_coords{ nullptr }; ///< Pointer to the array of x values (optional)
    size_t m_nx{ 0 };           ///< X dimension of the array

    realT *m_mask{ nullptr }; ///< Pointer to the (optional) mask array.  Any 0 pixels are excluded from the fit.

    realT m_G0{ 0 };
    realT m_G{ 0 };
    realT m_x0{ 0 };
    realT m_sigma{ 0 };

    int m_G0_idx{ 0 };
    int m_G_idx{ 1 };
    int m_x0_idx{ 2 };
    int m_sigma_idx{ 3 };

    int m_nparams{ 4 };
    int m_maxNparams{ 4 };

    /// Set whether each parameter is fixed.
    /** Sets the parameter indices appropriately.
     */
    void setFixed( bool G0, bool G, bool x0, bool sigma );

    realT G0( realT *p );

    void G0( realT *p, realT nG0 );

    realT G( realT *p );

    void G( realT *p, realT nG );

    realT x0( realT *p );

    void x0( realT *p, realT nx0 );

    realT sigma( realT *p );

    void sigma( realT *p, realT nsigma_x );

    int nparams();
};

template <typename realT>
void array2FitGaussian1D<realT>::setFixed( bool G0, bool G, bool x0, bool sigma )
{
    int idx = 0;

    if( G0 )
        m_G0_idx = -1;
    else
        m_G0_idx = idx++;

    if( G )
        m_G_idx = -1;
    else
        m_G_idx = idx++;

    if( x0 )
        m_x0_idx = -1;
    else
        m_x0_idx = idx++;

    if( sigma )
        m_sigma_idx = -1;
    else
        m_sigma_idx = idx++;

    m_nparams = idx;
}

template <typename realT>
realT array2FitGaussian1D<realT>::G0( realT *p )
{
    if( m_G0_idx < 0 )
    {
        return m_G0;
    }
    else
    {
        return p[m_G0_idx];
    }
}

template <typename realT>
void array2FitGaussian1D<realT>::G0( realT *p, realT nG0 )
{
    if( m_G0_idx < 0 )
    {
        m_G0 = nG0;
    }
    else
    {
        p[m_G0_idx] = nG0;
    }
}

template <typename realT>
realT array2FitGaussian1D<realT>::G( realT *p )
{
    if( m_G_idx < 0 )
    {
        return m_G;
    }
    else
    {
        return p[m_G_idx];
    }
}

template <typename realT>
void array2FitGaussian1D<realT>::G( realT *p, realT nG )
{
    if( m_G_idx < 0 )
    {
        m_G = nG;
    }
    else
    {
        p[m_G_idx] = nG;
    }
}

template <typename realT>
realT array2FitGaussian1D<realT>::x0( realT *p )
{
    if( m_x0_idx < 0 )
    {
        return m_x0;
    }
    else
    {
        return p[m_x0_idx];
    }
}

template <typename realT>
void array2FitGaussian1D<realT>::x0( realT *p, realT nx0 )
{
    if( m_x0_idx < 0 )
    {
        m_x0 = nx0;
    }
    else
    {
        p[m_x0_idx] = nx0;
    }
}

template <typename realT>
realT array2FitGaussian1D<realT>::sigma( realT *p )
{
    if( m_sigma_idx < 0 )
    {
        return m_sigma;
    }
    else
    {
        return p[m_sigma_idx];
    }
}

template <typename realT>
void array2FitGaussian1D<realT>::sigma( realT *p, realT nsigma )
{
    if( m_sigma_idx < 0 )
    {
        m_sigma = nsigma;
    }
    else
    {
        p[m_sigma_idx] = nsigma;
    }
}

template <typename realT>
int array2FitGaussian1D<realT>::nparams()
{
    return m_nparams;
}

} // namespace fit
} // namespace math

} // namespace mx

#endif // math_fit_array2FitGaussian1D_hpp
