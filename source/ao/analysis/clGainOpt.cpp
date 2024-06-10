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
   
template<typename realT>
realT _optGainOpenLoop( clGainOptOptGain_OL<realT> & olgo,
                        realT & var,
                        const realT & gmax,
                        const realT & minFindMin,
                        const realT & minFindMaxFact,
                        int minFindBits,
                        uintmax_t minFindMaxIter,
                        uintmax_t & iters
                      )
{
   realT gopt;

   iters = minFindMaxIter;
   
   try
   {
      std::pair<realT,realT> brack;
      brack = boost::math::tools::brent_find_minima<clGainOptOptGain_OL<realT>, realT>(olgo, minFindMin, minFindMaxFact*gmax, minFindBits, iters);
      gopt = brack.first;
      var = brack.second;
   }
   catch(...)
   {
      std::cerr << "optGainOpenLoop: No root found\n";
      gopt = minFindMaxFact*gmax;
      var = 0;
   }
   
   return gopt;
}

template<>
float optGainOpenLoop<float>( clGainOptOptGain_OL<float> & olgo,
                              float & var,
                              const float & gmax,
                              const float & minFindMin,
                              const float & minFindMaxFact,
                              int minFindBits,
                              uintmax_t minFindMaxIter,
                              uintmax_t & iters
                            )
{
   return _optGainOpenLoop<float>(olgo, var, gmax, minFindMin, minFindMaxFact, minFindBits, minFindMaxIter, iters);
}

template<>
double optGainOpenLoop<double>( clGainOptOptGain_OL<double> & olgo,
                                double & var,
                                const double & gmax,
                                const double & minFindMin,
                                const double & minFindMaxFact,
                                int minFindBits,
                                uintmax_t minFindMaxIter,
                                uintmax_t & iters
                              )
{
   return _optGainOpenLoop<double>(olgo, var, gmax, minFindMin, minFindMaxFact, minFindBits, minFindMaxIter, iters);
}
 
template<>
long double optGainOpenLoop<long double>( clGainOptOptGain_OL<long double> & olgo,
                                          long double & var,
                                          const long double & gmax,
                                          const long double & minFindMin,
                                          const long double & minFindMaxFact,
                                          int minFindBits,
                                          uintmax_t minFindMaxIter,
                                          uintmax_t & iters
                                        )
{
   return _optGainOpenLoop<long double>(olgo, var, gmax, minFindMin, minFindMaxFact, minFindBits, minFindMaxIter, iters);
}

#ifdef HASQUAD
template<>
__float128 optGainOpenLoop<__float128>( clGainOptOptGain_OL<__float128> & olgo,
                                          __float128 & var,
                                          __float128 & gmax,
                                          __float128 & minFindMin,
                                          __float128 & minFindMaxFact,
                                          int minFindBits,
                                          uintmax_t minFindMaxIter,
                                          uintmax_t iters
                                        )
{
   return _optGainOpenLoop<__float128>(olgo, var, gmax, minFindMin, minFindMaxFact, minFindBits, minFindMaxIter, iters);
}
#endif

} //namespace impl

   
//Explicit Instantiation
template
class clGainOpt<float>;

template
class clGainOpt<double>;

template
class clGainOpt<long double>;

#ifdef HASQUAD
template
class clGainOpt<__float128>;
#endif
   
   
   
} //namespace analysis
} //namespace ao
} //namespace mx

