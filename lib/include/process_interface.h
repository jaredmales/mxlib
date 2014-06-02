/** \file process_interface.h
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib c process interface facilities
  * \ingroup IPC
  * 
*/

#ifndef __process_interface_h__
#define __process_interface_h__

#include <string.h>
#include <stdio.h>
#include "IPC.h"

#ifdef __cplus_plus
extern "C"
{
#endif

///Run a process and copy the output to a string
/** Uses popen, so the attendant precautions about privileges apply.
  *
  * \param cmd is the command to run
  * \param resp is the string to copy the response to.
  * \param respsz is the available size of the string.
  *
  * \retval 0 on success 
  * \retval -1 on error
  * 
  */
int command_response(const char * cmd, char * resp, size_t respsz);


#ifdef __cplus_plus
} //extern "C"
#endif
   
#endif //__process_interface_h__
