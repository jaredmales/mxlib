/** \file mxlib.hpp
 * \author Jared R. Males
 * \brief Declarations of some libarary wide utilities
 *
 */

#ifndef mxlib_hpp
#define mxlib_hpp


#include "mxlib_uncomp_version.h"

#ifndef UNUSED
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif
#endif

namespace mx
{
   
///Dump the current git status of the library to a stream
/** Prints the current SHA1 hash and whether or not the library
  * has been modified since that commit.
  *
  * \tparam iosT a std::ostream-like type 
  * \tparam comment a character to print at the beginning of each line.  Default is '#'.
  */ 
template<typename iosT, char comment='#'>
iosT & dumpGitStatus( iosT & ios /**< [in] a std::ostream-like stream. */)
{
   //This causes the stream to not output a '\0' if comment = '\0'
   char c[] = {comment, '\0'};
      
   ios << c << "--------------------------------------------\n";
   ios << c << "mxlib git SHA1: " << MXLIB_UNCOMP_CURRENT_SHA1 << "\n";
   ios << c << "mxlib git modified flag: " << std::boolalpha << (bool) MXLIB_UNCOMP_REPO_MODIFIED << "\n";
   ios << c << "--------------------------------------------\n";
   
}

}//namespace mx
   
#endif // mxlib_hpp
