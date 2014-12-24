

#ifndef __eigenImage__
#define __eigenImage__


#include <Eigen/Dense>
#include "fitsUtils.hpp"
#include "eigenCube.hpp"

#include "ds9_interface.h"


namespace mx
{

/** \addtogroup eigen_image_processing
  * @{
  */

typedef Eigen::Array<short, Eigen::Dynamic, Eigen::Dynamic> eigenImages;

typedef Eigen::Array<unsigned short, Eigen::Dynamic, Eigen::Dynamic> eigenImageus;

typedef Eigen::ArrayXXi eigenImagei;

typedef Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic> eigenImageui;

typedef Eigen::Array<long, Eigen::Dynamic, Eigen::Dynamic> eigenImagel;

typedef Eigen::Array<unsigned long, Eigen::Dynamic, Eigen::Dynamic> eigenImageul;

typedef Eigen::ArrayXXf eigenImagef;

typedef Eigen::ArrayXXd eigenImaged;

typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> eigenImageld;

typedef Eigen::ArrayXXd eigenImagecf;

typedef Eigen::ArrayXXd eigenImagecd;

typedef Eigen::Array<std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic> eigenImagecld;


///Function object to get the number of planes in an eigenImage or eigenCube
/** The Eigen types don't have a planes() member.  This uses the \ref is_eigenCube compile time test.
 */
template<typename eigenT, bool iscube = is_eigenCube<eigenT>::value>
struct imPlanes
{
   typename eigenT::Index operator()(const eigenT & im)
   {
      return im.planes();
   }
};

///Specialization of \ref imPlanes for a 1D image. 
/** If this is used then the \ref is_eigenCube test is false.
 */
template<typename eigenT>
struct imPlanes<eigenT, false>
{
   typename eigenT::Index operator()(const eigenT &im)
   {
      return 1;
   }
};


/** \ingroup plotting
  * @{
  */

///Display an eigenImage or eigenCube in ds9
/** Is a wrapper for the \ref ds9_interface facility
 */
template <typename eigenT>
void ds9( const eigenT & im, int frame =1)
{   
   ::ds9_display(frame, im.data(), im.rows(), im.cols(), imPlanes<eigenT>()(im), getFitsBITPIX<typename eigenT::Scalar>());
}

///@}

///@}

} //namespace mx

#endif //__eigenImage__

