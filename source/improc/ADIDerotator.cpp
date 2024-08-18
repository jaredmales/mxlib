/** \file ADIDerotator.cpp
 * \author Jared R. Males
 * \brief Implements a generic ADI derotator class.
 * \ingroup hc_imaging_files
 *
 */

#include "improc/ADIDerotator.hpp"

namespace mx
{
namespace improc
{

template struct ADIDerotator<float>;
template struct ADIDerotator<double>;

} // namespace improc
} // namespace mx
