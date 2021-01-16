/** \file imageUtils.cpp
  * \author Jared R. Males
  * \brief  Implementation of image processing utilities
  * \ingroup image_processing_files
  *
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#include "improc/imageUtils.hpp"
#include <cstring>

namespace mx
{
namespace improc
{

void * imcpy( void * dest,   
              void * src,    
              size_t width,  
              size_t height, 
              size_t szof    
            )
{
   return memcpy(dest, src, width*height*szof);
}

void * imcpy_flipUD( void * dest,   
                     void * src,    
                     size_t width,  
                     size_t height, 
                     size_t szof    
                   )
{
   for(size_t rr=0; rr< height; ++rr)
   {
      memcpy((char *)dest + rr*width*szof, (char *)src + (height-1-rr)*width*szof, width*szof);
   }
   
   return dest;
}

void * imcpy_flipLR( void * dest,   
                     void * src,    
                     size_t width,  
                     size_t height, 
                     size_t szof    
                   )
{
   if(szof == 1)
   {
      uint8_t * d = (uint8_t *) dest;
      uint8_t * s = (uint8_t *)src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            d[rr*width + cc] =  s[rr*width + (width-1-cc)];
         }   
      }
   }
   else if(szof == 2)
   {
      uint16_t * d = (uint16_t *) dest;
      uint16_t * s = (uint16_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            d[rr*width + cc] =  s[rr*width + (width-1-cc)];
         }   
      }
   }
   else if(szof == 4)
   {
      uint32_t * d = (uint32_t *) dest;
      uint32_t * s = (uint32_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            d[rr*width + cc] =  s[rr*width + (width-1-cc)];
         }   
      }
   }
   else if(szof == 8)
   {
      uint64_t * d = (uint64_t *) dest;
      uint64_t * s = (uint64_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            d[rr*width + cc] =  s[rr*width + (width-1-cc)];
         }   
      }
   }
   else
   {
      uint8_t * d = (uint8_t *) dest;
      uint8_t * s = (uint8_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width*szof; cc += szof)
         {
            for(size_t pp=0; pp<szof;++pp)
            {
               d[rr*width*szof + cc + pp] =  s[rr*width*szof + (width*szof-szof-cc) +pp];
            }  
         }
      }  
   }  
   return dest;
}

void * imcpy_flipUDLR( void * dest,   
                       void * src,    
                       size_t width,  
                       size_t height, 
                       size_t szof    
                     )
{
   if(szof == 1)
   {
      uint8_t * d = (uint8_t *) dest;
      uint8_t * s = (uint8_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            {
               d[rr*width + cc] =  s[(height-1-rr)*width + (width-1-cc)];
            }
         }
      }
   }
   else if(szof == 2)
   {
      uint16_t * d = (uint16_t *) dest;
      uint16_t * s = (uint16_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            {
               d[rr*width + cc] =  s[(height-1-rr)*width + (width-1-cc)];
            }
         }
      }
   }
   else if(szof == 4)
   {
      uint32_t * d = (uint32_t *) dest;
      uint32_t * s = (uint32_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            {
               d[rr*width + cc] =  s[(height-1-rr)*width + (width-1-cc)];
            }
         }
      }
   }
   else if(szof == 8)
   {
      uint64_t * d = (uint64_t *) dest;
      uint64_t * s = (uint64_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width; ++cc)
         {
            {
               d[rr*width + cc] =  s[(height-1-rr)*width + (width-1-cc)];
            }
         }
      }
   }
   else
   {
      uint8_t * d = (uint8_t *) dest;
      uint8_t * s = (uint8_t *) src;
      for(size_t rr=0; rr< height; ++rr)
      {
         for(size_t cc=0; cc<width*szof; cc += szof)
         {
            for(size_t pp=0; pp<szof;++pp)
            {
               d[rr*width*szof + cc + pp] =  s[(height-1-rr)*width*szof + (width*szof-szof-cc) +pp];
            }
         }
      }
   }
   
   return dest;
}

} //namespace math
} //namespace mx

