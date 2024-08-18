/** \file expModGaussian.hpp
 * \brief The Exponentially Modified Gaussian distribution.
 * \ingroup gen_math_files
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

#ifndef expModGaussian_hpp
#define expModGaussian_hpp

#include <boost/math/tools/minima.hpp>
#include <limits>

namespace mx
{
namespace math
{
namespace func
{

/// The Exponentially Modified Gaussian at a point.
/** Calculates the value of the Exponentially Modified Gaussian distribution at a location specified by x.
 *
 *
 * \tparam realT a real floating point type
 *
 * \returns the value of the Exponentially Modified Gaussian distribution at x.
 *
 * \ingroup gen_math_expModGaussian
 */
template <typename realT>
realT expModGaussian( realT x,     ///< [in] the location at which to calculate the distribution
                      realT mu,    ///< [in] the mean parameter
                      realT sigma, ///< [in] the standard deviation
                      realT lambda ///< [in] the rate of decay
)
{
    return ( lambda / 2 ) * exp( ( lambda / 2 ) * ( 2 * mu + lambda * sigma * sigma - 2 * x ) ) *
           std::erfc( ( mu + lambda * sigma * sigma - x ) / ( root_two<realT>() * sigma ) );
}

/// The Mean of the Exponentially Modified Gaussian.
/** Calculates the mean of the Exponentially Modified Gaussian distribution.
 *
 *
 * \tparam realT a real floating point type
 *
 * \returns the mean of the Exponentially Modified Gaussian.
 *
 * \ingroup gen_math_expModGaussian
 */
template <typename realT>
realT expModGaussianMean( realT mu,    ///< [in] the mean parameter
                          realT lambda ///< [in] the rate of decay
)
{
    return mu + 1.0 / lambda;
}

/// The Variance of the Exponentially Modified Gaussian.
/** Calculates the variance of the Exponentially Modified Gaussian distribution.
 *
 *
 * \tparam realT a real floating point type
 *
 * \returns the variance of the Exponentially Modified Gaussian.
 *
 * \ingroup gen_math_expModGaussian
 */
template <typename realT>
realT expModGaussianVariance( realT sigma, ///< [in] the standard deviation
                              realT lambda ///< [in] the rate of decay
)
{
    return sigma * sigma + 1.0 / ( lambda * lambda );
}

template <typename realT>
struct emgModeFunc
{
    realT mu;
    realT sigma;
    realT lambda;

    realT operator()( const realT &x )
    {
        return -expModGaussian( x, mu, sigma, lambda );
    }
};

/// The Mode of the Exponentially Modified Gaussian.
/** Calculates the mode of the Exponentially Modified Gaussian distribution.
 * This is done iteratively with Brent's method.
 *
 * \tparam realT a real floating point type
 *
 * \returns the mode of the Exponentially Modified Gaussian.
 *
 * \ingroup gen_math_expModGaussian
 */
template <typename realT>
realT expModGaussianMode( realT mu,    ///< [in] the mean parameter
                          realT sigma, ///< [in] the standard deviation
                          realT lambda ///< [in] the rate of decay
)
{
    realT mn = expModGaussianMean( mu, lambda );
    realT sd = sqrt( expModGaussianVariance( sigma, lambda ) );

    std::cerr << mn << " " << sd << "\n";

    emgModeFunc<realT> mf;
    mf.mu = mu;
    mf.sigma = sigma;
    mf.lambda = lambda;

    uintmax_t maxit = 1000;
    try
    {
        std::pair<realT, realT> brack;
        brack = boost::math::tools::brent_find_minima<emgModeFunc<realT>, realT>(
            mf, mn - 2 * sd, mn + 2 * sd, std::numeric_limits<realT>::digits, maxit );
        std::cerr << brack.first << " " << brack.second << " " << maxit << "\n";

        return brack.first;
    }
    catch( ... )
    {
        std::cerr << "expModGaussianMode: No mode found\n";
        return std::numeric_limits<realT>::quiet_NaN();
    }
}

} // namespace func
} // namespace math
} // namespace mx

#endif // expModGaussian_hpp
