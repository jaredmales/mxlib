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
   std::string combineMethodStr( int method )
   {
      if(method == noCombine) return "noCombine";
      else if (method == medianCombine) return "medianCombine";
      else if (method == meanCombine) return "meanCombine";
      else if (method == sigmaMeanCombine) return "sigmaMeanCombine";
      else return "UNKNOWN";
   }
   
   int combineMethodFmStr( const std::string & method )
   {
      if(method == "noCombine") return noCombine;
      else if (method == "medianCombine") return medianCombine;
      else if (method == "meanCombine") return meanCombine;
      else if (method == "sigmaMeanCombine") return sigmaMeanCombine;
      else return -1;
   }
}

template class HCIobservation<float>;
template class HCIobservation<double>;

} //namespace improc
} //namespace mx

