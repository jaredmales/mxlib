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

template class aoAtmosphere<long double>;

#ifdef HASQUAD
template class aoAtmosphere<__float128>;
#endif

} // namespace analysis
} // namespace AO
} // namespace mx
