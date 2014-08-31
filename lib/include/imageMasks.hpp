/** \file imageMasks.hpp
  * \brief Declares and defines functions to work with image masks
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __imageMasks_hpp__
#define __imageMasks_hpp__

namespace mx
{

template<class eigenT> 
void setMask(eigenT & mask, std::vector<size_t> & idx, typename eigenT::Scalar maskval = 1)
{
   for(size_t i=0; i< idx.size(); ++i)
   {
      mask(idx[i]) = maskval;
   }
}


} //namespace mx

#endif // __imageMasks_hpp__
