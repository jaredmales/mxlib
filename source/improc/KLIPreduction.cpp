/** \file KLIPreduction.cpp
 * \author Jared R. Males
 * \brief Instantiations of an implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#include "improc/ADIDerotator.hpp"
#include "improc/KLIPreduction.hpp"

namespace mx
{
namespace improc
{
namespace HCI
{


std::string excludeMethodStr( int method )
{
    if( method == excludeNone )
        return "excludeNone";
    else if( method == excludePixel )
        return "excludePixel";
    else if( method == excludeAngle )
        return "excludeAngle";
    else if( method == excludeImno )
        return "excludeImno";
    else
        return "UNKNOWN";
}

int excludeMethodFmStr( const std::string &method )
{
    if( method == "excludeNone" )
        return excludeNone;
    else if( method == "excludePixel" )
        return excludePixel;
    else if( method == "excludeAngle" )
        return excludeAngle;
    else if( method == "excludeImno" )
        return excludeImno;
    else
        return -1;
}

std::string includeMethodStr( int method )
{
    if( method == includeAll )
        return "includeAll";
    else if( method == includeCorr )
        return "includeCorr";
    else if( method == includeTime )
        return "includeTime";
    else if( method == includeAngle )
        return "includeAngle";
    else if( method == includeImno )
        return "includeImno";
    else
        return "UNKNOWN";
}

int includeMethodFmStr( const std::string &method )
{
    if( method == "includeAll" )
        return includeAll;
    else if( method == "includeCorr" )
        return includeCorr;
    else if( method == "includeTime" )
        return includeTime;
    else if( method == "includeAngle" )
        return includeAngle;
    else if( method == "includeImno" )
        return includeImno;
    else
        return -1;
}
} // namespace HCI

template <typename realT>
class ADIDerotator;

template struct KLIPreduction<float, ADIDerotator<float>, float>;
template struct KLIPreduction<float, ADIDerotator<float>, double>;
template struct KLIPreduction<double, ADIDerotator<double>, double>;

} // namespace improc
} // namespace mx
