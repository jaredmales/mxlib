/** \file IPC.h
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib interprocess communication (IPC) tools
  * 
*/


#ifndef __mx_IPC_h__
#define __mx_IPC_h__

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>
#include <string.h>

#include <unistd.h>
   
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>

/** \addtogroup IPC
  * @{
  */
///The maximum length of the IPC key string
#define MX_IPC_KEYLEN 1024

///The process interface buffer size
#define MX_IPC_PI_BUFSZ 128


//Include rest of IPC facility here.
#include "msgq.h"
#include "sharedmem_segment.h"
#include "process_interface.h"

#ifdef __cplusplus
} //extern "C"
#endif

///@}

#endif //__mx_IPC_h__


