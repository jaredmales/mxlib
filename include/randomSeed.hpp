/** \file randomSeed.hpp
  * \author Jared R. Males
  * \brief Defines a random number seed generator
  * \ingroup random
  *
  */

#ifndef __randomSeed_hpp__
#define __randomSeed_hpp__


#include <unistd.h>
#include <fcntl.h>


#include "mxError.hpp"



namespace mx
{

///Get a value to use as a random seed 
/** On Linux systems, uses /dev/urandom to populate the value with sizeof(intT) bytes.  
  * Otherwise, uses time(0) to get time since the epoch.
  *
  * \param [out] seedval will be populated with the seed.
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
int randomSeed(intT & seedval)
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

   if(rv < sizeof(intT))
   {
      mxError("randomSeed", MXE_FILERERR, "Read from /dev/urandom did not return enough bytes");
         
      return -1;
   }
   
   return 0;

   
   #endif //__linux__
   
   
   seedval = time(0);
   
   return 0;

}

} //namespace mx


#endif //__randomSeed_hpp__
