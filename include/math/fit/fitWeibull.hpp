/** \file fitWeibull.hpp
 * \author Jared R. Males
 * \brief Tools for fitting the Weibull distribution to data.
 * \ingroup fitting_files
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

#ifndef fitWeibull_hpp
#define fitWeibull_hpp

#include "levmarInterface.hpp"
#include "../func/weibull.hpp"

namespace mx
{
namespace math
{
namespace fit
{

template <typename realT>
struct array2FitWeibull;

/** \defgroup weibull_peak_fit Weibull Distribution
 * \brief Fitting the Weibull Distribution to data.
 *
 * The Weibull Distribution is fit to data.
 *
 * \ingroup peak_fit
 */

/// Class to manage fitting the Weibull Distribution to data via the \ref levmarInterface
/** In addition to the requirements on fitterT specified by \ref levmarInterface
 * this class also requires this definition in fitterT
 * \code
 * static const int nparams = 3;
 * \endcode
 * where the number 3 is replaced by the number of parameters that fitterT expects to fit.
 *
 * \tparam fitterT a type meeting the above requirements.
 *
 * \ingroup weibull_peak_fit
 *
 */
template <typename fitterT>
class fitWeibull : public levmarInterface<fitterT>
{

  public:
    typedef typename fitterT::realT realT;

    static const int nparams = fitterT::nparams;

    array2FitWeibull<realT> arr;

    void initialize()
    {
        this->allocate_params( nparams );
        this->adata = &arr;
    }

    fitWeibull()
    {
        initialize();
    }

    ~fitWeibull()
    {
    }

    /// Set the initial guess when platescale and central obscuration are fixed.
    void setGuess( realT x0,    ///< [in] the location parameter
                   realT k,     ///< [in] the shape parameter
                   realT lambda ///< [in] the scale parameter
    )
    {
        static_assert( nparams == 3, "fitWeibull: Wrong setGuess called for no location parameter." );

        this->p[2] = x0;
        this->p[0] = k;
        this->p[1] = lambda;
    }

    /// Set the initial guess when platescale and central obscuration are fixed.
    void setGuess( realT k,     ///< [in] the shape parameter
                   realT lambda ///< [in] the scale parameter
    )
    {
        static_assert( nparams == 2, "fitWeibull: Wrong setGuess called for location parameter." );

        this->p[0] = k;
        this->p[1] = lambda;
    }

    void setArray( realT *data, int n )
    {
        arr.data = data;
        arr.n = n;

        this->n = n;
    }

    void x0( realT nx0 )
    {
        arr.x0 = nx0;
        if( nparams == 3 )
        {
            this->p[2] = nx0;
        }
    }

    void k( realT nk )
    {
        arr.k = nk;
        this->p[0] = nk;
    }

    void lambda( realT nl )
    {
        arr.lambda = nl;
        this->p[1] = nl;
    }

    int fit()
    {
        return levmarInterface<fitterT>::fit();
    }

    realT x0()
    {
        if( nparams == 3 )
        {
            return this->p[2];
        }
        else
        {
            return 0;
        }
    }

    realT k()
    {
        return this->p[0];
    }

    realT lambda()
    {
        return this->p[1];
    }
};

/// Wrapper for a native array to pass to \ref levmarInterface, with Weibull details.
/** \ingroup weibull_peak_fit
 */
template <typename realT>
struct array2FitWeibull
{
    realT *data{ nullptr }; ///< Pointer to the array
    size_t n{ 0 };          ///< dimension of the array

    realT x0{ 0 };     ///< the location parameter.
    realT k{ 0 };      ///< the shape parameter
    realT lambda{ 0 }; ///< the scale parameter
};

///\ref levmarInterface fitter structure for the shifted Weibull Distribution
/**
 *
 * \ingroup weibull_peak_fit
 *
 */
template <typename _realT>
struct weibull_3param_fitter
{
    typedef _realT realT;

    static const int nparams = 3;

    static void func( realT *p, realT *hx, int m, int n, void *adata )
    {
        array2FitWeibull<realT> *arr = (array2FitWeibull<realT> *)adata;

        for( int i = 0; i < arr->n; i++ )
        {
            hx[i] = func::weibull<realT>( i, p[2], p[0], p[1] ) - arr->data[i];
        }
    }
};

///\ref levmarInterface fitter structure for the shifted Weibull Distribution
/**
 *
 * \ingroup weibull_peak_fit
 *
 */
template <typename _realT>
struct weibull_2param_fitter
{
    typedef _realT realT;

    static const int nparams = 2;

    static void func( realT *p, realT *hx, int m, int n, void *adata )
    {
        array2FitWeibull<realT> *arr = (array2FitWeibull<realT> *)adata;

        for( int i = 0; i < arr->n; i++ )
        {
            hx[i] = func::weibull<realT>( i, 0.0, p[0], p[1] ) - arr->data[i];
        }
    }
};

} // namespace fit
} // namespace math
} // namespace mx

#endif // fitWeibull_hpp
