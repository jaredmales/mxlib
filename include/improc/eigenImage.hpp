

#ifndef __eigenImage__
#define __eigenImage__

#pragma GCC system_header
#include <Eigen/Dense>

#include "eigenCube.hpp"

namespace mx
{
namespace improc 
{
   
///Definition of the eigenImage type, which is an alias for Eigen::Array.
/** \ingroup eigen_image_processing
  */
template<typename scalarT>
using eigenImage = Eigen::Array<scalarT, -1, -1>;


#if 0
///Function object to get the number of planes in an eigenImage or eigenCube
/** The Eigen types don't have a planes() member.  This uses the \ref is_eigenCube compile time test.
  * \ingroup eigen_image_processing 
  */
template<typename eigenT, bool iscube = is_eigenCube<eigenT>::value>
struct imPlanes
{
   typename eigenT::Index operator()(const eigenT & im)
   {
      return im.planes();
   }
};

///Specialization of \ref imPlanes for a 2D image. 
/** If this is used then the \ref is_eigenCube test is false.
  * * \ingroup eigen_image_processing 
  */
template<typename eigenT>
struct imPlanes<eigenT, false>
{
   typename eigenT::Index operator()(const eigenT &im)
   {
      return 1;
   }
};
#endif

} //namespace improc
} //namespace mx

#endif //__eigenImage__

