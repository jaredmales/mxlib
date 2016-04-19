/** \file environment.hpp
  * \author Jared R. Males
  * \brief Utilities for working with the environment
  *
  */

#ifndef __environment_hpp__
#define __environment_hpp__


#include <locale>
#include <string>

namespace mx
{

   
/// Return the value of an environment variable
/** Call the standard getenv function, but handles the null pointer as an empty string.
  *
  * \param estr is the name of the environment variable to query
  *
  * \returns the value of the environment varialbe, or empty string if it doesn't exist
  */   
inline
std::string getEnv(const std::string & estr)
{
   char * e = getenv(estr.c_str());
   
   if(e) return std::string(e);
   else return std::string("");
}


}

#endif //__environment_hpp__
