/** \file aoAtmosphere.cpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Implementation of the AO Atmosphere.
 * \ingroup mxAO_files
 *
 */

#include "ao/analysis/aoAtmosphere.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

template class aoAtmosphere<float>;

template class aoAtmosphere<double>;

#ifdef MXLIB_USE_LONGDOUBLE
template class aoAtmosphere<long double>;
#endif

#ifdef MXLIB_HAS_QUAD
template class aoAtmosphere<__float128>;
#endif

} // namespace analysis
} // namespace AO
} // namespace mx
