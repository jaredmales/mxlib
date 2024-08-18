/** \file ADIobservation.cpp
 * \author Jared R. Males
 * \brief Instantiates the ADI high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#include "improc/ADIDerotator.hpp"
#include "improc/ADIobservation.hpp"

namespace mx
{
namespace improc
{
namespace HCI
{
std::string fakeMethodsStr( int method )
{
    if( method == single )
    {
        return "single";
    }
    else if( method == list )
    {
        return "list";
    }
    else
    {
        return "unknown";
    }
}

int fakeMethodFmStr( const std::string &method )
{
    if( method == "single" )
    {
        return single;
    }
    else if( method == "list" )
    {
        return list;
    }
    else
    {
        return -1;
    }
}
} // namespace HCI

template class ADIobservation<float, ADIDerotator<float>>;
template class ADIobservation<double, ADIDerotator<double>>;

} // namespace improc
} // namespace mx
