
#include "math/func/precision.hpp"

#include <boost/math/tools/precision.hpp>

namespace mx
{
namespace math
{
namespace func
{

template <>
float root_epsilon<float>()
{
    return boost::math::tools::root_epsilon<float>();
}

template <>
double root_epsilon<double>()
{
    return boost::math::tools::root_epsilon<double>();
}

} // namespace func
} // namespace math
} // namespace mx
