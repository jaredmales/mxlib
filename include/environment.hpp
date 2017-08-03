/** \file environment.hpp
  * \author Jared R. Males
  * \brief Utilities for working with the environment
  * \ingroup utils_files
  */

#ifndef environment_hpp
#define environment_hpp


#include <locale>
#include <string>

namespace mx
{

   
/// Return the value of an environment variable
/** Call the standard getenv function, but handles the null pointer as an empty string.
  *
  * \returns the value of the environment varialbe, or empty string if it doesn't exist
  * 
  * \ingroup system
  */   
inline
std::string getEnv(const std::string & estr /**< [in] is the name of the environment variable to query */)
{
   char * e = getenv(estr.c_str());
   
   if(e) return std::string(e);
   else return std::string("");
}


}

#endif //environment_hpp
