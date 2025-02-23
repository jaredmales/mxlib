/** \file HCIobservation.cpp
 * \author Jared R. Males
 * \brief Instantiation of the basic high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#include "improc/HCIobservation.hpp"

namespace mx
{

namespace improc
{

namespace HCI
{

std::string meanSubMethodStr( meanSubMethod method )
{
    if( method == meanSubMethod::none )
    {
        return "none";
    }
    else if( method == meanSubMethod::meanImage )
    {
        return "meanImage";
    }
    else if( method == meanSubMethod::medianImage )
    {
        return "medianImage";
    }
    else if( method == meanSubMethod::imageMean )
    {
        return "imageMean";
    }
    else if( method == meanSubMethod::imageMedian )
    {
        return "imageMedian";
    }
    else if( method == meanSubMethod::imageMode )
    {
        return "imageMode";
    }
    else
    {
        return "UNKNOWN";
    }
}

meanSubMethod meanSubMethodFmStr( const std::string &method )
{
    if( method == "none" )
    {
        return meanSubMethod::none;
    }
    else if( method == "meanImage" )
    {
        return meanSubMethod::meanImage;
    }
    else if( method == "medianImage" )
    {
        return meanSubMethod::medianImage;
    }
    else if( method == "imageMean" )
    {
        return meanSubMethod::imageMean;
    }
    else if( method == "imageMedian" )
    {
        return meanSubMethod::imageMedian;
    }
    else if( method == "imageMode" )
    {
        return meanSubMethod::imageMode;
    }
    else
    {
        return meanSubMethod::unknown;
    }
}

std::string pixelTSNormMethodStr( pixelTSNormMethod method )
{
    if( method == pixelTSNormMethod::none )
    {
        return "none";
    }
    else if( method == pixelTSNormMethod::rms )
    {
        return "rms";
    }
    else if( method == pixelTSNormMethod::rmsSigmaClipped )
    {
        return "rmsSigmaClipped";
    }
    else
    {
        return "UNKNOWN";
    }
}

pixelTSNormMethod pixelTSNormMethodFmStr( const std::string &method )
{
    if( method == "none" )
    {
        return pixelTSNormMethod::none;
    }
    else if( method == "rms" )
    {
        return pixelTSNormMethod::rms;
    }
    else if( method == "rmsSigmaClipped" )
    {
        return pixelTSNormMethod::rmsSigmaClipped;
    }
    else
    {
        return pixelTSNormMethod::unknown;
    }
}

std::string combineMethodStr( int method )
{
    if( method == noCombine )
        return "noCombine";
    else if( method == medianCombine )
        return "medianCombine";
    else if( method == meanCombine )
        return "meanCombine";
    else if( method == sigmaMeanCombine )
        return "sigmaMeanCombine";
    else
        return "UNKNOWN";
}

int combineMethodFmStr( const std::string &method )
{
    if( method == "noCombine" )
        return noCombine;
    else if( method == "medianCombine" )
        return medianCombine;
    else if( method == "meanCombine" )
        return meanCombine;
    else if( method == "sigmaMeanCombine" )
        return sigmaMeanCombine;
    else
        return -1;
}
} // namespace HCI

template class HCIobservation<float>;
template class HCIobservation<double>;

} // namespace improc
} // namespace mx
