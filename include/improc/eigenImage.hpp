

#ifndef __eigenImage__
#define __eigenImage__

#pragma GCC system_header
#include <Eigen/Dense>

#include "../vectorUtils.hpp"
//#include "eigenCube.hpp"

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

///Test whether a type is an eigenCube by testing whether it has a typedef of "is_eigenCube"
/** Used for compile-time determination of type
  * Example usage:
  * \code
  * bool is_eC = is_eigenCube<eigenCube<float> >; //Evaluates to true
  * bool is_not_eC = is_eigenCube<eigenImagef>; //Evaluates to false
  * \endcode
  * 
  * This was taken directly from the example at http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
  * 
  * \ingroup eigen_image_processing
  */

template <typename T>
struct is_eigenCube 
{
   // Types "yes" and "no" are guaranteed to have different sizes,
   // specifically sizeof(yes) == 1 and sizeof(no) == 2.
   typedef char yes[1];
   typedef char no[2];
 
   template <typename imageT>
   static yes& test(typename imageT::is_eigenCube*);
 
   template <typename>
   static no& test(...);
 
   // If the "sizeof" of the result of calling test<T>(0) would be equal to sizeof(yes),
   // the first overload worked and T has a nested type named "is_mmatrix".
   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};


///Function object to retun the number of planes for any Eigen like object, whether 2D or a 3D cube.
/** Uses SFINAE to check for 3D eigenCube.
  *
  * \ingroup eigen_image_processing 
  */
template<typename arrT, bool isCube=is_eigenCube<arrT>::value>
struct eigenArrPlanes
{
   //If it's an eigenCube, call planes planes()
   int operator()(const arrT & arr)
   {
      return arr.planes(); 
   }
};

template<typename arrT>
struct eigenArrPlanes<arrT, false>
{
   //If it's not an eigenCube, never call planes()
   int operator()(const arrT & arr)
   {
      return 1;
   }
};

/// Calculate the median of an Eigen-like array.
/** Calculates the median of the entire array, allowing for some pixels to be ignored using a mask.
  * Working memory can be retained between calls.
  * 
  * \tparam imageT is an Eigen-like type
  * \tparam maskT is an Eigen-like type
  * 
  * \returns the median of the unmasked pixels of mat, using \ref vectorMedianInPlace().
  * 
  * \ingroup eigen_image_processing
  * 
  */ 
template<typename imageT, typename maskT=imageT>
typename imageT::Scalar imageMedian( const imageT & mat,  ///< [in] the image to take the median of
                                     const maskT * mask,  ///< [in] if non-0, a 1/0 mask where 0 pixels are ignored.
                                     std::vector<typename imageT::Scalar> * work =0 ///< [in] [optional] working memory can be retained and re-passed.
                                   )
{
   typename imageT::Scalar med;
   
   bool localWork = false;
   if(work == 0) 
   {
      work = new std::vector<typename imageT::Scalar>;
      localWork = true;
   }
   
   int sz = mat.size();
   
   if(mask)
   {
      sz = mask->sum();
   }
   
   work->resize(sz);
   
   int ii = 0;
   for(int i=0;i<mat.rows();++i)
   {
      for(int j=0; j<mat.cols();++j)
      {
         if(mask)
         {
            if( (*mask)(i,j) == 0) continue;
         }
         
         (*work)[ii] = mat(i,j);
         ++ii;
      }
   }

   
   med = vectorMedianInPlace(*work);
   
   if(localWork) delete work;
   
   return med;
} 

/// Calculate the median of an Eigen-like array.
/** Calculates the median of the entire array.
  * Working memory can be retained between calls.
  * 
  * \tparam imageT is an Eigen-like type
  * 
  * \returns the median of the unmasked pixels of mat, using \ref vectorMedianInPlace().
  * 
  * \ingroup eigen_image_processing
  * 
  */ 
template<typename imageT>
typename imageT::Scalar imageMedian( const imageT & mat, ///< [in] the image to take the median of
                                     std::vector<typename imageT::Scalar> * work =0 ///< [in] [optional] working memory can be retained and re-passed.
                                   )
{
   return imageMedian( mat, (Eigen::Array<typename imageT::Scalar,-1,-1> *) 0, work);
}


} //namespace improc
} //namespace mx

#endif //__eigenImage__
