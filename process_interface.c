/** \file process_interface.c
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions for the mxlib c process interface facilities
  * \ingroup IPC
  * 
*/

#include "process_interface.h"


int command_response(const char * cmd, char * resp, size_t respsz)
{
   size_t written;
   char buf[MX_IPC_PI_BUFSZ];
   FILE *cmd_ptr;

   written = 0;

   if ( (cmd_ptr = popen(cmd, "r")) != NULL) 
   {
      while (fgets(buf, MX_IPC_PI_BUFSZ, cmd_ptr) != NULL && written+1 < respsz)
      {
         strcat(resp, buf);
         written = strlen(resp);                  
      }
      pclose(cmd_ptr);
   }
   else
   {
      return -1;
   }

   return 0;
}


