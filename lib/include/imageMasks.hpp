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
void applyMask(eigenT & mask, std::vector<size_t> & idx, typename eigenT::Scalar maskval = 1)
{
   for(size_t i=0; i< idx.size(); ++i)
   {
      mask(idx[i]) = maskval;
   }
}

inline void genRectMask(std::vector<size_t> & idx, size_t rows, size_t cols, size_t xmin, size_t xmax, size_t ymin, size_t ymax)
{
      
   if(xmin < 0) xmin = 0;
   if(xmax > rows-1) xmax = rows-1;
   
   if(ymin < 0) ymin = 0;
   if(ymax > cols-1) ymax = cols-1;
   
   idx.reserve( (xmax-xmin+1)*(ymax-ymin + 1) );
   
   for(int i=xmin; i<=xmax; ++i)
   {
      for(int j=ymin;j<=ymax; ++j)
      {
         idx.push_back( j*rows + i);
      }
   }
}

template<class eigenT>
void genRectMask(std::vector<size_t> & idx, eigenT & mask, size_t xmin, size_t xmax, size_t ymin, size_t ymax)
{
   genRectMask(idx, (size_t) mask.rows(), (size_t) mask.cols(), xmin, xmax, ymin, ymax);  
}


} //namespace mx

#endif // __imageMasks_hpp__
