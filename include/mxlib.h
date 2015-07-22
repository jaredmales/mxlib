/** \file mxlib.h
 * \author Jared R. Males
 * \brief Declarations of some libarary wide utilities
 *
 */

#ifndef __mxlib_h__
#define __mxlib_h__


#include "mxlib_uncomp_version.h"

#ifdef __cplusplus
extern "C"
{
#endif


   
/// Get the compile-time git repository SHA-1
/** Returns a pointer to a static string containing the SHA-1 hash of the
  * mxlib git repository at compile time.  This pointer should not be modified.
  * 
  * \retval char* which points to a string containing the SHA-1 hash
  * 
  */
char * mxlib_compiled_git_sha1();

/// Get the compile-time git repository modification status
/** Returns a pointer to a static int, if 1 then at compile 
  * time there were modifications not committed.
  * 
  * \retval int denoting whether the repo contained uncommitted modifications at compile time.
  */
int mxlib_compiled_git_repo_modified();

#ifdef __cplusplus
} //extern "C"
#endif


#endif // __mxlib_h__
