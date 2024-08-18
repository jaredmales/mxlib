/** \file gammaDistribution.hpp
 * \brief The Gamma Distribution.
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

#ifndef gammaDistribution_hpp
#define gammaDistribution_hpp

#include <cmath>

#include "gamma.hpp"

namespace mx
{
namespace math
{
namespace func
{

/// The denominator of the Gamma Distribution
/** Can be used to avoid repeated calculations when the parameters are constant
 *
 * \tparam realT a real floating point type
 *
 * \returns the denominator of the Gamma Distribution.
 *
 * \ingroup gen_math_gammaDist
 */
template <typename realT>
realT gammaDistributionDenom( realT k, ///< [in] shape parameter
                              realT q  ///< [in] the scale parameter
)
{
    try
    {
        realT d = pow( q, k ) * tgamma<realT>( k );

        if( !std::isnormal( d ) )
            return std::numeric_limits<realT>::max();

        return d;
    }
    catch( ... )
    {
        return std::numeric_limits<realT>::max();
    }
}

/// The general shifted Gamma Distribution at a point using an arbitrary peak scaling parameter
/** Calculates the value of the Gamma Distribution at a location specified by x.
 *
 * \tparam realT a real floating point type
 *
 * \returns the value of the Gamma distribution at x.
 *
 * \ingroup gen_math_gammaDist
 */
template <typename realT>
realT gammaDistribution( realT x,    ///< [in] the location at which to calculate the distribution
                         realT x0,   ///< [in] the location parameter
                         realT k,    ///< [in] shape parameter
                         realT q,    ///< [in] the scale parameter
                         realT denom ///< [in] the denominator, or 1/peak-scale.
)
{
    if( x - x0 < 0 )
        return 0;

    realT v = pow( x - x0, k - 1 ) * exp( -( x - x0 ) / q ) / denom;

    if( !std::isnormal( v ) )
        return 0;

    return v;
}

/// The general shifted Gamma Distribution at a point.
/** Calculates the value of the Gamma Distribution at a location specified by x.
 *
 * \tparam realT a real floating point type
 *
 * \returns the value of the Gamma distribution at x.
 *
 * \ingroup gen_math_gammaDist
 */
template <typename realT>
realT gammaDistribution( realT x,  ///< [in] the location at which to calculate the distribution
                         realT x0, ///< [in] the location parameter
                         realT k,  ///< [in] shape parameter
                         realT q   ///< [in] the scale parameter
)
{

    return gammaDistribution<realT>( x, x0, k, q, gammaDistributionDenom<realT>( k, q ) );
}

/// The mean of the Gamma Distribution
/** Calculates the mean of the Gamma Distribution for the given parameters.
 *
 * \tparam realT a real floating point type
 *
 * \returns the mean of the Gamma Distribution.
 *
 * \ingroup gen_math_gammaDist
 */
template <typename realT>
realT gammaDistributionMean( realT x0,   ///< [in] the location parameter
                             realT k,    ///< [in] shape parameter
                             realT theta ///< [in] the scale parameter
)
{
    return x0 + k * theta;
}

/// The mode of the Gamma Distribution
/** Calculates the mode of the Gamma Distribution for the given parameters.
 *
 * \tparam realT a real floating point type
 *
 * \returns the mode of the Gamma Distribution.
 *
 * \ingroup gen_math_gammaDist
 */
template <typename realT>
realT gammaDistributionMode( realT x0,   ///< [in] the location parameter
                             realT k,    ///< [in] shape parameter
                             realT theta ///< [in] the scale parameter
)
{
    if( k >= 1 )
        return x0 + ( k - 1 ) * theta;

    return 0;
}

/// The variance of the Gamma Distribution
/** Calculates the variance of the Gamma Distribution for the given parameters.
 *
 * \tparam realT a real floating point type
 *
 * \returns the variance of the Gamma Distribution.
 *
 * \ingroup gen_math_gammaDist
 */
template <typename realT>
realT gammaDistributionVariance( realT k,    ///< [in] shape parameter
                                 realT theta ///< [in] the scale parameter
)
{
    return k * theta * theta;
}

} // namespace func
} // namespace math
} // namespace mx

#endif // gammaDistribution_hpp
