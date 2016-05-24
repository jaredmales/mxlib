/** \file randomSeed.hpp
  * \author Jared R. Males
  * \brief Defines a random number seed generator
  * \ingroup random
  *
  */

#ifndef __randomSeed_hpp__
#define __randomSeed_hpp__

#include <fcntl.h>


namespace mx
{

///Get a value to use as a random seed 
/** On Linux systems, uses /dev/urandom to populate the value with sizeof(intT) bytes.  
  * Otherwise, uses time(0) to get time since the epoch.
  *
  * \param seedval will be populated with the seed.
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
      
   fd = open("/dev/urandom", O_RDONLY);

   seedval = 0;

   int rv =read(fd, &seedval, sizeof(intT));
      
   close(fd);
      
    
   if(rv < 0)
   {
      std::cerr << "mx::randomSeed: read from /dev/urandom returned error\n";
         
      return rv;
   }
     
   return 0;
#endif //__linux__
   
   
   seedval = time(0);
   
   return 0;

}

} //namespace mx


#endif //__randomSeed_hpp__
