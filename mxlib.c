/** \file mxlib.c
 * \author Jared R. Males
 * \brief Definitions of some libarary wide utilities
 *
 */

#include "mxlib.h"


//Include the build-time git version header
#include "mxlib_comp_version.h"

char * mxlib_compiled_git_sha1()
{
   static char sha1[] = MXLIB_COMP_CURRENT_SHA1;

   return sha1;
}


int mxlib_compiled_git_repo_modified()
{
   return MXLIB_COMP_REPO_MODIFIED;
}



