/** \file clGainOpt.cpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief clGainOpt specializations for pre-compiled boost brent minimization
 * \ingroup mxAO_files
 *
 */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "ao/analysis/clGainOpt.hpp"

#include <boost/math/tools/minima.hpp>

#define aaaa

namespace mx
{
namespace AO
{
namespace analysis
{
namespace impl
{

template <typename realT>
mx::error_t _optGainOpenLoop( realT &gain,
                              realT &var,
                              clGainOptOptGain_OL<realT> &olgo,
                              const realT &minimumGain,
                              const realT &maximumGain,
                              int minFindBits,
                              uintmax_t minFindMaxIter,
                              uintmax_t &iters )
{
    gain = std::numeric_limits<realT>::quiet_NaN();
    var = std::numeric_limits<realT>::quiet_NaN();
    iters = minFindMaxIter;

    try
    {
        std::pair<realT, realT> brack;
        brack = boost::math::tools::brent_find_minima<clGainOptOptGain_OL<realT>, realT>( olgo,
                                                                                          minimumGain,
                                                                                          maximumGain,
                                                                                          minFindBits,
                                                                                          iters );
        gain = brack.first;
        var = brack.second;
    }
    catch( ... )
    {
        return error_t::exception;
    }

    if( iters >= minFindMaxIter )
    {
        return error_t::timeout;
    }

    return error_t::noerror;
}

template <>
mx::error_t optGainOpenLoop<float>( float &gain,
                                    float &var,
                                    clGainOptOptGain_OL<float> &olgo,
                                    const float &minimumGain,
                                    const float &maximumGain,
                                    int minFindBits,
                                    uintmax_t minFindMaxIter,
                                    uintmax_t &iters )
{
    return _optGainOpenLoop<float>( gain, var, olgo, minimumGain, maximumGain, minFindBits, minFindMaxIter, iters );
}

template <>
mx::error_t optGainOpenLoop<double>( double &gain,
                                     double &var,
                                     clGainOptOptGain_OL<double> &olgo,
                                     const double &minimumGain,
                                     const double &maximumGain,
                                     int minFindBits,
                                     uintmax_t minFindMaxIter,
                                     uintmax_t &iters )
{
    return _optGainOpenLoop<double>( gain, var, olgo, minimumGain, maximumGain, minFindBits, minFindMaxIter, iters );
}

template <>
mx::error_t optGainOpenLoop<long double>( long double &gain,
                                          long double &var,
                                          clGainOptOptGain_OL<long double> &olgo,
                                          const long double &minimumGain,
                                          const long double &maximumGain,
                                          int minFindBits,
                                          uintmax_t minFindMaxIter,
                                          uintmax_t &iters )
{
    return _optGainOpenLoop<long double>( gain,
                                          var,
                                          olgo,
                                          minimumGain,
                                          maximumGain,
                                          minFindBits,
                                          minFindMaxIter,
                                          iters );
}

#ifdef HASQUAD
template <>
mx::error_t optGainOpenLoop<__float128>( __float128 &gain,
                                         __float128 &var,
                                         clGainOptOptGain_OL<__float128> &olgo,
                                         const __float128 &minimumGain,
                                         const __float128 &maximumGain,
                                         int minFindBits,
                                         uintmax_t minFindMaxIter,
                                         uintmax_t &iters )
{
    return _optGainOpenLoop<__float128>( gain,
                                         var,
                                         olgo,
                                         minimumGain,
                                         maximumGain,
                                         minFindBits,
                                         minFindMaxIter,
                                         iters );
}
#endif

} // namespace impl

// Explicit Instantiation
template class clGainOpt<float>;

template class clGainOpt<double>;

template class clGainOpt<long double>;

#ifdef HASQUAD
template class clGainOpt<__float128>;
#endif

} // namespace analysis
} // namespace AO
} // namespace mx
