/** \file fitEmpirical.hpp
 * \author Jared R. Males
 * \brief Tools for fitting empirical functions to data.
 * \ingroup fitting_files
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

#ifndef fitEmpirical_hpp
#define fitEmpirical_hpp

#include "levmarInterface.hpp"

#include "../../improc/eigenImage.hpp"
#include "../../improc/imageTransforms.hpp"

namespace mx
{
namespace math
{
namespace fit
{

// forward
template <typename realT>
struct array2FitEmpirical;

// forward
template <typename _realT>
struct empirical2D_sym_fitter;

/// Class to manage fitting a 2D Moffat to data via the \ref levmarInterface
/** Fits \ref gen_math_moffats to a 2-dimensional array of data.
 *
 * This class allows for treating any of the parameters as fixed.
 *
 * \tparam fitterT a type meeting the requirements specified in \ref levmarInterface.
 *
 * \ingroup moffat_peak_fit
 *
 */
template <typename fitterT>
class fitEmpirical2DGen : public levmarInterface<fitterT>
{

  public:
    typedef typename fitterT::realT realT;

  public:
    array2FitEmpirical<realT> arr;

    void initialize()
    {
        this->allocate_params( arr.nparams() );
        this->adata = &arr;
    }

  public:
    fitEmpirical2DGen()
    {
        initialize();
    }

    ~fitEmpirical2DGen()
    {
    }

    /// Set whether each parameter is fixed.
    /** Sets the parameter indices appropriately.
     */
    void setFixed( bool scale, ///< [in] if true, then scale will be not be part of the fit
                   bool dx,    ///< [in] if true, then dx will be not be part of the fit
                   bool dy     ///< [in] if true, then dy will be not be part of the fit
    )
    {
        arr.setFixed( scale, dx, dy );
        this->allocate_params( arr.nparams() );
    }

    /// Set the initial guess for the empirical fit
    void setGuess( realT scale, ///< [in] the scale factr
                   realT dx,    ///< [in] the x shift
                   realT dy     ///< [in] the y shift
    )
    {
        arr.scale( this->p, scale );
        arr.dx( this->p, dx );
        arr.dy( this->p, dy );
    }

    /// Set the data aray.
    void setArray( const eigenImage<realT> *data, ///< [in] pointer to an nx X ny array of data to be fit
                   const eigenImage<realT> *ref,   ///< [in] pointer to the empirical function to fit
                   const eigenImage<realT> *weights
    )
    {
        arr.setup( data, ref, weights );

        this->n = arr.m_nx * arr.m_ny;
    }

    /// Set the data aray.
    void setArray( const eigenImage<realT> *data, ///< [in] pointer to an nx X ny array of data to be fit
                   const eigenImage<realT> *ref   ///< [in] pointer to the empirical function to fit
    )
    {
        setArray( data, ref, nullptr);
    }

    /// Get the current value of the scale factor.
    /**
     * \returns the current value of scale
     */
    realT scale()
    {
        return arr.scale( this->p );
    }

    /// Get the current value of dx, the x shift
    /**
     * \returns the current value of dx
     */
    realT dx()
    {
        return arr.dx( this->p );
    }

    /// Get the current value of dy, the y shift
    /**
     * \returns the current value of dy
     */
    realT dy()
    {
        return arr.dy( this->p );
    }
};

/// Wrapper for a native array to pass to \ref levmarInterface, with empirical function fit details.
/** \ingroup moffat_peak_fit
 */
template <typename realT>
struct array2FitEmpirical
{
    const eigenImage<realT> *m_data{ nullptr };    ///< Pointer to the data array
    const eigenImage<realT> *m_ref{ nullptr };     ///< Pointer to the reference image to fit to the data
    const eigenImage<realT> *m_weights{ nullptr }; ///< Pointer to the weight image
    eigenImage<realT> m_refShifted;                ///< Working memory for the shifted reference image

    size_t m_nx{ 0 };                              ///< X dimension of the array
    size_t m_ny{ 0 };                              ///< Y dimension of the array

    realT m_scale{ 0 };
    realT m_dx{ 0 };
    realT m_dy{ 0 };

    int m_scale_idx{ 0 };
    int m_dx_idx{ 1 };
    int m_dy_idx{ 2 };

    int m_nparams = 3;

    void setup( const eigenImage<realT> *data, const eigenImage<realT> *ref, const eigenImage<realT> *weights )
    {
        m_nx = data->rows();
        m_ny = data->cols();

        if( ref->rows() != m_nx || ref->cols() != m_ny )
        {
            std::cerr << "ref and data not same size\n";
            exit( -1 );
        }

        if( weights )
        {
            if( weights->rows() != m_nx || weights->cols() != m_ny )
            {
                std::cerr << "ref and data not same size\n";
                exit( -1 );
            }
        }

        m_data = data;
        m_ref = ref;
        m_weights = weights;
        m_refShifted.resize( m_nx, m_ny );
    }

    void setup( const eigenImage<realT> *data, const eigenImage<realT> *ref )
    {
        setup( data, ref, nullptr );
    }

    /// Set whether each parameter is fixed.
    /** Sets the parameter indices appropriately.
     */
    void setFixed( bool scale, bool dx, bool dy )
    {
        int idx = 0;

        if( scale )
            m_scale_idx = -1;
        else
            m_scale_idx = idx++;

        if( dx )
            m_dx_idx = -1;
        else
            m_dx_idx = idx++;

        if( dy )
            m_dy_idx = -1;
        else
            m_dy_idx = idx++;

        m_nparams = idx;
    }

    realT scale( realT *p )
    {
        if( m_scale_idx < 0 )
        {
            return m_scale;
        }
        else
        {
            return p[m_scale_idx];
        }
    }

    void scale( realT *p, realT nscale )
    {
        if( m_scale_idx < 0 )
        {
            m_scale = nscale;
        }
        else
        {
            p[m_scale_idx] = nscale;
        }
    }

    realT dx( realT *p )
    {
        if( m_dx_idx < 0 )
        {
            return m_dx;
        }
        else
        {
            return p[m_dx_idx];
        }
    }

    void dx( realT *p, realT ndx )
    {
        if( m_dx_idx < 0 )
        {
            m_dx = ndx;
        }
        else
        {
            p[m_dx_idx] = ndx;
        }
    }

    realT dy( realT *p )
    {
        if( m_dy_idx < 0 )
        {
            return m_dy;
        }
        else
        {
            return p[m_dy_idx];
        }
    }

    void dy( realT *p, realT ndy )
    {
        if( m_dy_idx < 0 )
        {
            m_dy = ndy;
        }
        else
        {
            p[m_dy_idx] = ndy;
        }
    }

    int nparams()
    {
        return m_nparams;
    }
};

///\ref levmarInterface fitter structure for 2D empirical functions.
/** \ingroup moffat_peak_fit
 *
 */
template <typename _realT>
struct empirical2D_fitter
{
    typedef _realT realT;

    static const int nparams = 3;

    static void func( realT *p, realT *hx, int m, int n, void *adata )
    {
        array2FitEmpirical<realT> *arr = (array2FitEmpirical<realT> *)adata;

        size_t idx_dat;

        idx_dat = 0;

        realT scale = arr->scale( p );
        realT dx = arr->dx( p );
        realT dy = arr->dy( p );

        imageShift( arr->m_refShifted, *( arr->m_ref ), dx, dy, cubicConvolTransform<realT>() );

        arr->m_refShifted *= scale;

        if( arr->m_weights )
        {
            for( int cc = 0; cc < arr->m_ny; ++cc )
            {
                for( int rr = 0; rr < arr->m_nx; ++rr )
                {
                    hx[idx_dat] =(*arr->m_weights)(rr,cc)*( ( *arr->m_data)( rr, cc ) - arr->m_refShifted( rr, cc ) );

                    ++idx_dat;
                }
            }
        }
        else
        {
            for( int cc = 0; cc < arr->m_ny; ++cc )
            {
                for( int rr = 0; rr < arr->m_nx; ++rr )
                {
                    hx[idx_dat] = ( ( *( arr->m_data ) )( rr, cc ) - arr->m_refShifted( rr, cc ) );

                    ++idx_dat;
                }
            }
        }
    }
};

/// Alias for the fitEmpirical2D type fitting the symmetric Moffat profile.
/** \ingroup moffat_peak_fit
 */
template <typename realT>
using fitEmpirical2D = mx::math::fit::fitEmpirical2DGen<mx::math::fit::empirical2D_fitter<realT>>;

} // namespace fit
} // namespace math

} // namespace mx

#endif // fitEmpirical_hpp
