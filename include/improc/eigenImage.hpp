/** \file eigenImage.hpp
  * \brief Tools for using the eigen library for image processing
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef improc_eigenImage_hpp
#define improc_eigenImage_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include "../math/vectorUtils.hpp"

namespace mx
{
namespace improc 
{
   
///Definition of the eigenImage type, which is an alias for Eigen::Array.
/** \ingroup eigen_image_processing
  */
template<typename scalarT>
using eigenImage = Eigen::Array<scalarT, -1, -1>;


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
   // the first overload worked and T has a nested type named "is_eigenCube".
   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};


///Function object to return the number of planes for any Eigen like object, whether 2D or a 3D cube.
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
typename imageT::Scalar imageMedian( const imageT & mat,                            ///< [in] the image to take the median of
                                     const maskT * mask,                            ///< [in] if non-0, a 1/0 mask where 0 pixels are ignored.
                                     std::vector<typename imageT::Scalar> * work =0 ///< [in] [optional] working memory can be retained and re-passed.  Is resized.
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

   
   med = math::vectorMedianInPlace(*work);
   
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
typename imageT::Scalar imageMedian( const imageT & mat,                            ///< [in] the image to take the median of
                                     std::vector<typename imageT::Scalar> * work =0 ///< [in] [optional] working memory can be retained and re-passed.
                                   )
{
   return imageMedian( mat, (Eigen::Array<typename imageT::Scalar,-1,-1> *) 0, work);
}


} //namespace improc
} //namespace mx

#endif //improc_eigenImage_hpp

