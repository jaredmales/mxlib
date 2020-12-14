/** \file randomSeed.hpp
  * \author Jared R. Males
  * \brief Defines a random number seed generator
  * \ingroup gen_math_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef mx_math_randomSeed_hpp
#define mx_math_randomSeed_hpp


#include <unistd.h>
#include <fcntl.h>


#include "../mxError.hpp"



namespace mx
{
namespace math 
{

///Get a value to use as a random seed 
/** On Linux systems, uses /dev/urandom to populate the value with sizeof(intT) bytes.  
  * Otherwise, uses time(0) to get time since the epoch.
  * 
  * \returns 0 on success.
  * \returns -1 on error.
  * 
  * \tparam intT is the integer type of seeval. 
  * 
  * \ingroup random
  * 
  */ 
template<typename intT>
int randomSeed(intT & seedval /**< [out] will be populated with the seed.*/ )
{   
   #ifdef __linux__
   
   int fd;
   
   errno = 0;
   fd = open("/dev/urandom", O_RDONLY);

   if(fd < 0)
   {
      mxPError("randomSeed", errno, "error opening /dev/urandom");
         
      return -1;
   }
   
   seedval = 0;

   errno = 0;
   int rv = ::read(fd, &seedval, sizeof(intT));
      
      
   if(rv < 0)
   {
      mxPError("randomSeed", errno, "Error on read from /dev/urandom.");
      close(fd);
      return -1;
   }
   
   close(fd);

   int sz = sizeof(intT);
   
   if(rv < sz)
   {
      mxError("randomSeed", MXE_FILERERR, "Read from /dev/urandom did not return enough bytes");
         
      return -1;
   }
   
   return 0;

   
   #endif //__linux__
   
   
   seedval = time(0);
   
   return 0;

}

} //namespace math
} //namespace mx


#endif //mx_math_randomSeed_hpp
