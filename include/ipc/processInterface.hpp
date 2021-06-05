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

#include <vector>
#include <string>

#include "../mxlib.hpp"
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
  * \ingroup IPC 
  */
int command_response(const char * cmd, char * resp, size_t respsz);

/// Runs a command (with parameters) passed in using fork/exec
/** A new process is fork()-ed, and the child runs execvp with command provided.  The output of the process
  * is captured in \p commandOutput, and error messages in \p commandStderr.
  * 
  * The process return value in \p retVal is only meaningful if this function returns 0.
  * 
  * If this function returns -1, the last entry in \p commandStderr will contain a message,
  * and errno may be useful.
  * 
  * \returns 0 on success
  * \returns -1 on error
  * 
  * \ingroup IPC
  */
int runCommand( int & retVal,                                ///< [out] the return value of the process. Only meaningful if this returns 0.
                std::vector<std::string> & commandOutput,    ///< [out] the output, line by line.  If an error, first entry contains the message.
                std::vector<std::string> & commandStderr,    ///< [out] the output of stderr.
                const std::vector<std::string> & commandList ///< [in] command to be run, with one entry per command line word
              );

}//namespace ipc
} //namespace mx

#endif //ipc_processInterface_hpp
