/** \file zernike.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Working with the Zernike polynomials.
 *
 * \todo the basic zernike polys should be in math::func.
 *
 * \ingroup signal_processing_files
 *
 */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_zernike_hpp
#define math_zernike_hpp

#include <cmath>
#include <complex>
#include <vector>

#include "../mxlib.hpp"
#include "../math/func/bessel.hpp"
#include "../math/func/jinc.hpp"
#include "../math/func/factorial.hpp"
#include "../math/func/sign.hpp"
#include "../math/constants.hpp"

namespace mx
{
namespace sigproc
{

/**
 * \ingroup zernike_basis
 * @{
 */

/// Get the Zernike coefficients n,m corrresponding the Noll index j.
/** Calculates the values of (n,m) for an index j following Noll (1976) \cite noll_1976
 * See also: http://en.wikipedia.org/wiki/Zernike_polynomials
 *
 * If j is odd, this returns m <= 0.
 *
 *
 * \retval 0 on success
 * \retval -1 on error (j < 1)
 *
 */
int noll_nm( int &n, ///< [out] n the radial index of the Zernike polynomial
             int &m, ///< [out] m the azimuthal index of the Zernnike polynomial.  m < 0 if j odd.
             int j   ///< [in]  j the Noll index, j > 0.
);

/// Get the Noll index j corresponding to Zernike coefficients n,m
/** Calculates the value j for(n,m) following Noll (1976) \cite noll_1976
 * See also: http://en.wikipedia.org/wiki/Zernike_polynomials
 *
 * \retval >= 0 on success
 * \retval -1 on error (n-m odd)
 *
 */
int noll_j( unsigned n, ///< [in] n the radial index of the Zernike polynomial
            int m       ///< [in] m the azimuthal index of the Zernnike polynomial.
);

/// Get the number of Zernikes up to and including a radial order.
/** Calculates the total number of Zernike polynomials through radial order \p n.  See Noll (1976) \cite noll_1976
 * See also: http://en.wikipedia.org/wiki/Zernike_polynomials
 *
 * \retval the number of
 * \retval -1 on error (n-m odd)
 *
 */
int nZernRadOrd( unsigned n /**< [n] the radial order */ );

/// Calculate the coefficients of a Zernike radial polynomial
/**
 * \retval 0 on success
 * \retval -1 on error
 *
 * \tparam realT is a real floating type
 */
template <typename realT>
int zernikeRCoeffs(
    std::vector<realT> &c, ///< [out] allocated to length \f$ 0.5(n-m)+1\f$ and filled with the coefficients.
    int n,                 ///< [in] the radial index of the Zernike polynomial.
    int m                  ///< [in] the azimuthal index of the Zernike polynomial.
)
{
    m = abs( m );

    if( n < m )
    {
        internal::mxlib_error_report( error_t::invalidarg, "n cannot be less than m in the Zernike polynomials" );
        return -1;
    }

    // If odd, it's 0.
    if( ( n - m ) % 2 > 0 )
    {
        c.resize( 1, 0 );
        return 0;
    }

    int ul = 0.5 * ( n - m ) + 1;

    c.resize( ul );

    for( int k = 0; k < ul; ++k )
    {
        c[k] = pow( -1.0, k ) * math::func::factorial<realT>( n - k ) /
               ( math::func::factorial<realT>( k ) * math::func::factorial<realT>( 0.5 * ( n + m ) - k ) *
                 math::func::factorial<realT>( 0.5 * ( n - m ) - k ) );
    }

    return 0;
}

// Explicit instantiations:
extern template int zernikeRCoeffs<float>( std::vector<float> &c, int n, int m );

extern template int zernikeRCoeffs<double>( std::vector<double> &c, int n, int m );

extern template int zernikeRCoeffs<long double>( std::vector<long double> &c, int n, int m );

#ifdef HASQUAD
extern template int zernikeRCoeffs<__float128>( std::vector<__float128> &c, int n, int m );
#endif

/// Calculate the value of a Zernike radial polynomial at a given separation.
/**
 *
 * \retval -9999 indicates a possible error
 * \retval R the value of the Zernike radial polynomial otherwise
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calcs, should be at least double.
 */
template <typename realT, typename calcRealT>
realT zernikeR( realT rho,                ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.
                int n,                    ///< [in] the radial index of the Zernike polynomial.
                int m,                    ///< [in] the azimuthal index of the Zernike polynomial.
                std::vector<calcRealT> &c /**< [in] contains the radial polynomial coeeficients,
                                                    and must be of length \f$ 0.5(n-m)+1\f$. */
)
{
    m = abs( m );

    // If odd, it's 0.
    if( ( n - m ) % 2 > 0 )
    {
        c.resize( 1, 0 );
        return 0.0;
    }

    if( c.size() != 0.5 * ( n - m ) + 1 )
    {
        internal::mxlib_error_report( error_t::invalidarg, "c vector has incorrect length for n and m." );
        return -9999;
    }

    realT R = 0.0;
    for( size_t k = 0; k < c.size(); ++k )
    {
        R += c[k] * pow( rho, n - 2 * k );
    }

    return R;
}

extern template float zernikeR<float, double>( float rho, int n, int m, std::vector<double> &c );

extern template double zernikeR<double, double>( double rho, int n, int m, std::vector<double> &c );

extern template long double
zernikeR<long double, long double>( long double rho, int n, int m, std::vector<long double> &c );

#ifdef HASQUAD
extern template __float128 zernikeR<__float128, __float128>( __float128 rho, int n, int m, std::vector<__float128> &c );
#endif

/// Calculate the value of a Zernike radial polynomial at a given separation.
/**
 * \retval -9999 indicates a possible error
 * \retval R the value of the Zernike radial polynomial otherwise
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename realT, typename calcRealT>
realT zernikeR( realT rho, ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.
                int n,     ///< [in] the radial index of the Zernike polynomial.
                int m      ///< [in] the azimuthal index of the Zernike polynomial.
)
{
    m = abs( m );

    // If odd, it's 0.
    if( ( n - m ) % 2 > 0 )
    {
        return 0.0;
    }

    std::vector<calcRealT> c;

    if( zernikeRCoeffs( c, n, m ) < 0 )
        return -9999;

    return zernikeR<realT, calcRealT>( rho, n, m, c );
}

extern template float zernikeR<float, double>( float rho, int n, int m );

extern template double zernikeR<double, double>( double rho, int n, int m );

extern template long double zernikeR<long double, long double>( long double rho, int n, int m );

#ifdef HASQUAD
extern template __float128 zernikeR<__float128, __float128>( __float128 rho, int n, int m );
#endif

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/**
 * \retval -9999 indicates a possible error
 * \retval R the value of the Zernike radial polynomial otherwise
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename realT, typename calcRealT>
realT zernike( realT rho,                /**< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.*/
               realT phi,                /**< [in] the azimuthal angle (in radians)*/
               int n,                    /**< [in] the radial index of the Zernike polynomial.*/
               int m,                    /**< [in] the azimuthal index of the Zernike polynomial.*/
               std::vector<calcRealT> &c /**< [in] contains the radial polynomial coeeficients, and
                                                   must be of length \f$ 0.5(n-m)+1\f$.*/
)
{
    realT azt;

    if( n == 0 && m == 0 )
    {
        return 1.0;
    }

    if( m < 0 )
    {
        azt = math::root_two<realT>() * sin( -m * phi );
    }
    else if( m > 0 )
    {
        azt = math::root_two<realT>() * cos( m * phi );
    }
    else
    {
        azt = 1.0;
    }

    return sqrt( (realT)n + 1 ) * zernikeR<realT, calcRealT>( rho, n, m, c ) * azt;
}

extern template float zernike<float, double>( float rho, float phi, int n, int m, std::vector<double> &c );

extern template double zernike<double, double>( double rho, double phi, int n, int m, std::vector<double> &c );

extern template long double
zernike<long double, long double>( long double rho, long double phi, int n, int m, std::vector<long double> &c );

#ifdef HASQUAD
extern template __float128
zernike<__float128, __float128>( __float128 rho, __float128 phi, int n, int m, std::vector<__float128> &c );
#endif

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/**
 *
 * \retval -9999 indicates a possible error
 * \retval R the value of the Zernike radial polynomial otherwise
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename realT, typename calcRealT>
realT zernike( realT rho, ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.
               realT phi, ///< [in] the azimuthal angle (in radians)
               int n,     ///< [in] the radial index of the Zernike polynomial.
               int m      ///< [in] the azimuthal index of the Zernike polynomial.
)
{

    std::vector<calcRealT> c;

    if( zernikeRCoeffs<calcRealT>( c, n, m ) < 0 )
        return -9999;

    return zernike<realT, calcRealT>( rho, phi, n, m, c );
}

extern template float zernike<float, double>( float rho, float phi, int n, int m );

extern template double zernike<double, double>( double rho, double phi, int n, int m );

extern template long double zernike<long double, long double>( long double rho, long double phi, int n, int m );

#ifdef HASQUAD
extern template __float128 zernike<__float128, __float128>( __float128 rho, __float128 phi, int n, int m );
#endif

/// Calculate the value of a Zernike radial polynomial at a given radius and angle.
/**
 * \retval -9999 indicates a possible error
 * \retval R the value of the Zernike radial polynomial otherwise
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename realT, typename calcRealT>
realT zernike( realT rho, ///< [in] the radial coordinate, \f$ 0 \le \rho \le 1 \f$.
               realT phi, ///< [in] the azimuthal angle (in radians)
               int j      ///< [in] the Noll index of the Zernike polynomial.
)
{
    int n, m;

    // Get n and m from j
    if( noll_nm( n, m, j ) < 0 )
        return -9999;

    return zernike<realT, calcRealT>( rho, phi, n, m );
}

extern template float zernike<float, double>( float rho, float phi, int j );

extern template double zernike<double, double>( double rho, double phi, int j );

extern template long double zernike<long double, long double>( long double rho, long double phi, int j );

#ifdef HASQUAD
extern template __float128 zernike<__float128, __float128>( __float128 rho, __float128 phi, int j );
#endif

/// Fill in an Eigen-like array with a Zernike polynomial
/** Sets any pixel which is at rad \<= r \< rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename arrayT, typename calcRealT, int overscan = 2>
int zernike( arrayT &arr,                     /**< [out] allocated array with an Eigen-like interface. The rows()
                                                         and cols() members are used to size the polynomial. */
             int n,                           /**< [in] the radial index of the polynomial*/
             int m,                           /**< [in] the azimuthal index of the polynomial*/
             typename arrayT::Scalar xcen,    /**< [in] the x coordinate of the desired center of the polynomial,
                                                        in pixels*/
             typename arrayT::Scalar ycen,    /**< [in] the y coordinate of the desired center of the polynomial,
                                                        in pixels*/
             typename arrayT::Scalar rad = -1 /**< [in] the desired radius. If rad \<= 0, then the maximum radius
                                                        based on dimensions of m is used.*/
)
{
    typedef typename arrayT::Scalar realT;
    realT x;
    realT y;
    realT r, rho;
    realT phi;

    std::vector<calcRealT> c;

    if( zernikeRCoeffs( c, n, m ) < 0 )
        return -1;

    size_t l0 = arr.rows();
    size_t l1 = arr.cols();

    if( rad <= 0 )
        rad = 0.5 * std::min( l0 - 1, l1 - 1 );

    for( size_t i = 0; i < l0; ++i )
    {
        for( size_t j = 0; j < l1; ++j )
        {
            x = i - xcen;
            y = j - ycen;

            r = std::sqrt( x * x + y * y );

            // This is to be consistent with mx::circularPupil while still respecting the Zernike rules
            if( r > rad && r <= rad + ( 1.0 / overscan ) )
                r = rad;

            rho = r / rad;

            if( rho <= 1.0 )
            {
                phi = std::atan2( y, x );
                arr( i, j ) = zernike( rho, phi, n, m, c );
            }
            else
            {
                arr( i, j ) = 0.0;
            }
        }
    }
    return 0;
}

/// Fill in an Eigen-like array with a Zernike polynomial
/** Sets any pixel which is at rad \<= r \<= rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
 *
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename arrayT, typename calcRealT>
int zernike(
    arrayT &arr,                     /**< [out] is the allocated array with an Eigen-like interface. The rows() and
                                                cols() members are used to size the polynomial.*/
    int j,                           /**< [in] is the Noll index of the polynomial*/
    typename arrayT::Scalar xcen,    /**< [in] is the x coordinate of the desired center of the polynomial, in pixels */
    typename arrayT::Scalar ycen,    /**< [in] is the y coordinate of the desired center of the polynomial, in pixels*/
    typename arrayT::Scalar rad = -1 /**< [in] is the desired radius. If rad \<= 0, then the maximum radius
    based on dimensions of m is used.*/
)
{
    typedef typename arrayT::Scalar realT;

    int n, m;

    if( noll_nm( n, m, j ) < 0 )
        return -1;

    return zernike<arrayT, calcRealT>( arr, n, m, xcen, ycen, rad );
}

/// Fill in an Eigen-like array with a Zernike polynomial
/** The geometric center of the array, 0.5*(arr.rows()-1), 0.5*(arr.cols()-1), is used as the center.
 * Sets any pixel which is at rad \<= r \< rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
 *
 * \tparam realT is a real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename arrayT, typename calcRealT>
int zernike( arrayT &arr,                     /**< [out] allocated array with an Eigen-like interface.
                                                         The rows() and cols() members are used to size
                                                         the polynomial.*/
             int n,                           /**< [in] the radial index of the polynomial*/
             int m,                           /**< [in] the azimuthal index of the polynomial*/
             typename arrayT::Scalar rad = -1 /**< [in] [opt] the desired radius. If rad \<= 0, then the
                                                              maximum radius based on dimensions of m is used.*/
)
{
    typename arrayT::Scalar xcen = 0.5 * ( arr.rows() - 1.0 );
    typename arrayT::Scalar ycen = 0.5 * ( arr.cols() - 1.0 );

    return zernike<arrayT, calcRealT>( arr, n, m, xcen, ycen, rad );
}

/// Fill in an Eigen-like array with a Zernike polynomial
/** The geometric center of the array, 0.5*(arr.rows()-1), 0.5*(arr.cols()-1), is used as the center.
 * Sets any pixel which is at rad \<= r \< rad+0.5 pixels to rho = 1, to be consistent with mx::circularPupil
 *
 * \tparam arrayT is an Eigen-like array of real floating type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename arrayT, typename calcRealT>
int zernike( arrayT &arr,                     /**< [out] the allocated array with an Eigen-like interface.
                                                         The rows() and cols() members are used to size
                                                         the polynomial.*/
             int j,                           /**< [in] the Noll index of the polynomial*/
             typename arrayT::Scalar rad = -1 /**< [in] [opt] the desired radius. If rad \<= 0, then the maximum
                                                              radius based on dimensions of m is used.*/
)
{
    typename arrayT::Scalar xcen = 0.5 * ( arr.rows() - 1.0 );
    typename arrayT::Scalar ycen = 0.5 * ( arr.cols() - 1.0 );

    return zernike<arrayT, calcRealT>( arr, j, xcen, ycen, rad );
}

/// Fill in an Eigencube-like array with Zernike polynomials in Noll order
/** The cube is pre-allocated to set the image size and the number of modes.
 *
 * \returns 0 on success
 * \returns -1 on error
 *
 * \tparam cubeT is an Eigencube-like array with real floating point type
 * \tparam calcRealT is a real floating type used for internal calculations, should be at least double
 */
template <typename cubeT, typename calcRealT>
int zernikeBasis( cubeT &cube,                     /**< [in/out] the pre-allocated cube which will be filled
                                                                 with the Zernike basis*/
                  typename cubeT::Scalar rad = -1, /**< [in] [opt] the radius of the aperture.  If -1 then the full
                                                                   image size is used.*/
                  int minj = 2                     /**< [in] [opt] the minimum j value to include. The default is j=2,
                                                                   which skips piston (j=1).*/
)
{
    typedef typename cubeT::imageT arrayT;

    typename cubeT::imageT im;

    im.resize( cube.rows(), cube.cols() );

    int rv;
    for( int i = 0; i < cube.planes(); ++i )
    {
        rv = zernike<arrayT, calcRealT>( im, minj + i, rad );

        if( rv < 0 )
        {
            return rv;
        }
        cube.image( i ) = im;
    }

    return 0;
}

/// Calculate the square-normed Fourier transform of a Zernike polynomial at position (k,phi)
/** Implements Equation (8) of Noll (1976) \cite noll_1976.
 *
 * \todo need a more robust jinc_n function for n > 1
 *
 *
 * \returns the value of |Q(k,phi)|^2
 *
 * \tparam realT is the floating point type used for arithmetic
 */
template <typename realT>
std::complex<realT> zernikeQ( realT k,   /**< [in] the radial coordinate of normalized spatial frequency. This is in the
                                                  \cite noll_1976 convention of cycles-per-radius.*/
                              realT phi, ///< [in] the azimuthal coordinate of normalized spatial frequency
                              int n,     ///< [in] the Zernike polynomial n
                              int m      ///< [in] the Zernike polynomial m
)
{

    std::complex<realT> Q;

    // sloppy implementation of jinc_n for k ~ 0
    if( k < 1e-12 )
    {
        if( n == 0 )
            Q = 1.0;
        else
            Q = 0.0;
    }
    else
    {
        Q = math::func::bessel_j( n + 1, math::two_pi<realT>() * k ) / ( math::pi<realT>() * k );
    }

    Q = sqrt( n + 1 ) * Q;

    if( m > 0 ) // Even j (see Noll 1976)
    {
        Q = Q * pow( -1, 0.5 * ( n - m ) ) * pow( std::complex<realT>( { 0, 1 } ), m ) * sqrt( 2 ) * cos( m * phi );
    }
    else if( m < 0 ) // Odd j (see Noll 1976) , but m can't really be neg
    {
        Q = Q * pow( -1, 0.5 * ( n + m ) ) * pow( std::complex<realT>( { 0, 1 } ), -m ) * sqrt( 2 ) * sin( -m * phi );
    }
    else
    {
        Q = Q * pow( -1, 0.5 * n );
    }

    return Q;
}

/// Calculate the square-normed Fourier transform of a Zernike polynomial at position (k,phi)
/** Implements Equation (8) of Noll (1976) \cite noll_1976.
 *
 * \todo need a more robust jinc_n function for n > 1
 *
 *
 * \returns the value of |Q(k,phi)|^2
 *
 * \tparam realT is the floating point type used for arithmetic
 */
template <typename realT>
realT zernikeQNorm( realT k,   /**< [in] the radial coordinate of normalized spatial frequency. This is in the
                                         \cite noll_1976 convention of cycles-per-radius.*/
                    realT phi, ///< [in] the azimuthal coordinate of normalized spatial frequency
                    int n,     ///< [in] the Zernike polynomial n
                    int m      ///< [in] the Zernike polynomial m
)
{

    realT B;

    // sloppy implementation of jinc_n for k ~ 0
    if( k < 0.00001 )
    {
        if( n == 0 )
        {
            B = 1.0;
        }
        else
        {
            B = 0.0;
        }
    }
    else
    {
        B = math::func::bessel_j( n + 1, math::two_pi<realT>() * k ) / ( math::pi<realT>() * k );
    }

    realT Q2 = ( n + 1 ) * ( B * B );

    if( m > 0 ) // Even j (see Noll 1976)
    {
        Q2 = 2 * Q2 * pow( cos( m * phi ), 2 );
    }
    else if( m < 0 ) // Odd j (see Noll 1976)
    {
        Q2 = 2 * Q2 * pow( sin( -m * phi ), 2 );
    }

    return Q2;
}

extern template float zernikeQNorm<float>( float k, float phi, int n, int m );

extern template double zernikeQNorm<double>( double k, double phi, int n, int m );

extern template long double zernikeQNorm<long double>( long double k, long double phi, int n, int m );

#ifdef HASQUAD
extern template __float128 zernikeQNorm<__float128>( __float128 k, __float128 phi, int n, int m );
#endif

/// Calculate the square-normed Fourier transform of a Zernike polynomial at position (k,phi)
/** Implements Equation (8) of Noll (1976) \cite noll_1976.
 *
 * \returns the value of |Q(k,phi)|^2
 *
 * \tparam realT is the floating point type used for arithmetic
 *
 */
template <typename realT>
realT zernikeQNorm( realT k,   /**< [in] the radial coordinate of normalized spatial frequency. This is in the
                                         \cite noll_1976 convention of cycles-per-radius.*/
                    realT phi, ///< [in] the azimuthal coordinate of normalized spatial frequency
                    int j      ///< [in] the Zernike polynomial index j (Noll convention)
)
{
    int n, m;

    noll_nm( n, m, j );

    return zernikeQNorm( k, phi, n, m );
}

/// Fill in an Eigen-like array with the square-normed Fourier transform of a Zernike polynomial
/** The array is filled in with \f$\lvert Q(k,\phi)\rvert^2\f$ according to Equation (8) of Noll (1976).
 *
 *
 * \returns 0 on success
 * \returns -1 on error
 *
 * \tparam arrayT is the Eigen-like array type.  Arithmetic will be done in arrayT::Scalar.
 */
template <typename arrayT>
int zernikeQNorm( arrayT &arr, /**< [out] the allocated array. The rows() and cols() members are used to size
                                          the transform.*/
                  arrayT &k,   /**< [in] the normalized spatial frequency magnitude at each pixel, in the Noll (1976)
                                         convention of cycles-per-radius.*/
                  arrayT &phi, ///< [in] the spatial frequency angle at each pixel
                  int j        ///< [in] the polynomial index in the Noll convention \cite noll_1976
)
{
    if( arr.rows() != k.rows() || arr.cols() != k.cols() )
    {
        internal::mxlib_error_report( error_t::invalidarg, "output array and input k are not the same size" );
        return -1;
    }

    if( arr.rows() != phi.rows() || arr.cols() != phi.cols() )
    {
        internal::mxlib_error_report( error_t::invalidarg, "output array and input phi are not the same size" );
        return -1;
    }

    int n, m;
    if( noll_nm( n, m, j ) < 0 )
        return -1; // noll_nm will explain error

    for( size_t i = 0; i < arr.rows(); ++i )
    {
        for( size_t j = 0; j < arr.cols(); ++j )
        {
            arr( i, j ) = zernikeQNorm( k( i, j ), phi( i, j ), n, m );
        }
    }
    return 0;
}

/// Calculate the spatial power spectrum of Piston
template <typename realT>
realT zernikePPiston( const realT &kD /**< [in] Spatial frequency in diameter units, i.e. cycles per aperture.*/ )
{
    return 4 * pow( math::func::jinc( math::pi<realT>() * kD ), 2 );
}

/// Calculate the spatial power spectrum of Tip \& Tilt
template <typename realT>
realT zernikePTipTilt( const realT &kD /**< [in] Spatial frequency in diameter units, i.e. cycles per aperture.*/ )
{
    return 16 * pow( math::func::jincN( 2, math::pi<realT>() * kD ), 2 );
}

/// Calculate the spatial power spectrum of Defocus
template <typename realT>
realT zernikePDefocus( const realT &kD /**< [in] Spatial frequency in diameter units, i.e. cycles per aperture.*/ )
{
    return 12 * pow( math::func::jincN( 3, math::pi<realT>() * kD ), 2 );
}

/// Calculate the spatial power spectrum of Astigmatism
template <typename realT>
realT zernikePAstig( const realT &kD /**< [in] Spatial frequency in diameter units, i.e. cycles per aperture.*/ )
{
    return 24 * pow( math::func::jincN( 3, math::pi<realT>() * kD ), 2 );
}

/// Calculate the spatial power spectrum of Coma
template <typename realT>
realT zernikePComa( const realT &kD /**< [in] Spatial frequency in diameter units, i.e. cycles per aperture.*/ )
{
    return 32 * pow( math::func::jincN( 4, math::pi<realT>() * kD ), 2 );
}

/// Calculate the spatial power spectrum of Trefoil
template <typename realT>
realT zernikePTrefoil( const realT &kD /**< [in] Spatial frequency in diameter units, i.e. cycles per aperture.*/ )
{
    return 32 * pow( math::func::jincN( 4, math::pi<realT>() * kD ), 2 );
}

/// Get the degrees of correction coefficient for Zernike polynomials in Kolmogorov turbulence.
/** Returns the coefficient from Table IV of Noll (1976) \cite noll_1976 for the given index.
 * Given this, the total variance in radians at wavelength \f$\lambda\f$ for a diameter \f$D\f$
 * aperture after correcting \f$j\f$ modes can be calculated from
 * \f[
 * realT c = zernikeModeDOCKolmogorov( j );
 * var_lambda = c * pow(D/r_0_lambda, math::five_thirds<realT>)
 * \f]
 * \returns 0 if noll_j is 0, indicating an error
 * \returns the coefficient
 */
template <typename realT>
realT zernikeModeDOCKolmogorov( unsigned noll_j /**< [in] the mode index, must be greater than 0 */ )
{
    realT c;

    switch( noll_j )
    {
        case 1:
            c = 1.0299;
            break;
        case 2:
            c = 0.582;
            break;
        case 3:
            c = 0.134;
            break;
        case 4:
            c = 0.111;
            break;
        case 5:
            c = 0.0880;
            break;
        case 6:
            c = 0.0648;
            break;
        case 7:
            c = 0.0587;
            break;
        case 8:
            c = 0.0525;
            break;
        case 9:
            c = 0.0463;
            break;
        case 10:
            c = 0.0401;
            break;
        case 11:
            c = 0.0377;
            break;
        case 12:
            c = 0.0352;
            break;
        case 13:
            c = 0.0328;
            break;
        case 14:
            c = 0.0304;
            break;
        case 15:
            c = 0.0279;
            break;
        case 16:
            c = 0.0267;
            break;
        case 17:
            c = 0.0255;
            break;
        case 18:
            c = 0.0243;
            break;
        case 19:
            c = 0.0232;
            break;
        case 20:
            c = 0.0220;
            break;
        case 21:
            c = 0.0208;
            break;
        default:
            if( noll_j == 0 )
            {
                return 0;
            }

            c = 0.2944 * pow( static_cast<realT>( noll_j ), -1 * math::half_root_three<realT>() );
    }

    return c;
}

/// Get the degrees of correction for Zernike polynomials in Kolmogorov turbulence.
/** Returns the degree of correction from Table IV of Noll (1976) \cite noll_1976 for the given index.
 * This is the total variance in radians for Fried parameter $r_0$ for a diameter $D$.
 * Equivalent to:
 * \f[
 * realT c = zernikeModeDOCKolmogorov( j );
 * var_lambda = c * pow(D/r_0_lambda, math::five_thirds<realT>)
 * \f]
 * \returns 0 if noll_j is 0, indicating an error
 * \returns the coefficient
 */
template <typename realT>
realT zernikeModeDOCKolmogorov( unsigned noll_j, /**< [in] the mode index, must be greater than 0 */
                                realT D,         /**< [in] the aperture diameter, in same units as r_0 */
                                realT r_0        /** <[in] the Fried parameter, in same units as D */
)
{
    if( noll_j == 0 )
    {
        return 0;
    }

    return zernikeModeDOCKolmogorov<realT>( noll_j ) * pow( D / r_0, math::five_thirds<realT>() );
}

/// Get the difference in degrees of correction coefficient for Zernike polynomials in Kolmogorov turbulence.
/** Returns the difference in coefficients from Table IV of Noll (1976) \cite noll_1976 for the given index
 * from the previous mode. Given this, the variance in radians-squared at wavelength \f$\lambda\f$ for a
 * diameter \f$D\f$ aperture in the \f$j\f$-th mode can be calculated from
 * \f[
 * realT c = zernikeModeDOCDiffKolmogorov( j );
 * var_j_lambda = c * pow(D/r_0_lambda, math::five_thirds<realT>)
 * \f]
 *
 * \returns 0 if noll_j is 0 or 1, indicating an error
 * \returns the coefficient difference
 */
template <typename realT>
realT zernikeModeDOCDiffKolmogorov( unsigned noll_j /**< [in] the mode index, must be greater than 1 */ )
{
    realT c;

    switch( noll_j )
    {
        case 2:
            c = 1.0299 - 0.582;
            break;
        case 3:
            c = 0.582 - 0.134;
            break;
        case 4:
            c = 0.134 - 0.111;
            break;
        case 5:
            c = 0.111 - 0.0880;
            break;
        case 6:
            c = 0.0880 - 0.0648;
            break;
        case 7:
            c = 0.0648 - 0.0587;
            break;
        case 8:
            c = 0.0587 - 0.0525;
            break;
        case 9:
            c = 0.0525 - 0.0463;
            break;
        case 10:
            c = 0.0463 - 0.0401;
            break;
        case 11:
            c = 0.0401 - 0.0377;
            break;
        case 12:
            c = 0.0377 - 0.0352;
            break;
        case 13:
            c = 0.0352 - 0.0328;
            break;
        case 14:
            c = 0.0328 - 0.0304;
            break;
        case 15:
            c = 0.0304 - 0.0279;
            break;
        case 16:
            c = 0.0279 - 0.0267;
            break;
        case 17:
            c = 0.0267 - 0.0255;
            break;
        case 18:
            c = 0.0255 - 0.0243;
            break;
        case 19:
            c = 0.0243 - 0.0232;
            break;
        case 20:
            c = 0.0232 - 0.0220;
            break;
        case 21:
            c = 0.0220 - 0.0208;
            break;
        case 22:
            c = 0.0208 - 0.2944 * pow( static_cast<realT>( 22 ), -1 * math::half_root_three<realT>() );
            break;
        default:
            if( noll_j < 2 )
            {
                return 0;
            }

            c = 0.2944 * ( pow( static_cast<realT>( noll_j - 1 ), -1 * math::half_root_three<realT>() ) -
                           pow( static_cast<realT>( noll_j ), -1 * math::half_root_three<realT>() ) );
    }

    return c;
}

/// Get the variance for a single Zernike polynomial in Kolmogorov turbulence.
/** Returns the variance from Table IV of Noll (1976) \cite noll_1976 for the given mode index.
 * For the $j$-th mode on a diameter $D$ aperture and Fried parameter $r_0$ this is equivalent to calculating
 * \f[
 * realT c = zernikeModeDOCDiffKolmogorov( j );
 * var_j_lambda = c * pow(D/r_0_lambda, math::five_thirds<realT>)
 * \f]
 *
 * \returns 0 if noll_j is 0 or 1, indicating an error
 * \returns the variance for the \param noll_j mode in radians squared.
 */
template <typename realT>
realT zernikeModeDOCDiffKolmogorov( unsigned noll_j, /**< [in] the mode index, must be greater than 0 */
                                    realT D,         /**< [in] the aperture diameter, in same units as r_0 */
                                    realT r_0 /** <[in] the Fried parameter, in same units as D */ )
{
    if( noll_j < 2 )
    {
        return 0;
    }

    return zernikeModeDOCDiffKolmogorov<realT>( noll_j ) * pow( D / r_0, math::five_thirds<realT>() );
}

///@} signal_processing

} // namespace sigproc
} // namespace mx

#endif // math_zernike_hpp
