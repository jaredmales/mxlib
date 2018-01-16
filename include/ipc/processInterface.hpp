/** \file processInterface.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Process interface facilities
  * \ingroup IPC
  * 
*/

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef ipc_processInterface_hpp
#define ipc_processInterface_hpp

#include <string.h>
#include <stdio.h>

#include "ipc.hpp"

namespace mx
{
namespace ipc 
{
   
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
inline
int command_response(const char * cmd, char * resp, size_t respsz)
{
   FILE *cmd_ptr;

   if ( (cmd_ptr = popen(cmd, "r")) != NULL) 
   {
      size_t written = 0;
         
      char buf[MX_IPC_PI_BUFSZ];
      
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

}//namespace ipc
} //namespace mx

#endif //ipc_processInterface_hpp
