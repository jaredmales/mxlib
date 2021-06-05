/** \file processInterface.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Process interface facilities
  * \ingroup IPC
  * 
*/

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#include "ipc/processInterface.hpp"

#include <cstring>
#include <sstream>
#include <cstdio>
#include <iostream>

#include <unistd.h>
#include <sys/wait.h>

namespace mx
{
namespace ipc 
{
   
int command_response( const char * cmd, 
                      char * resp, 
                      size_t respsz
                    )
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

int runCommand( std::vector<std::string> & commandOutput,    // [out] the output, line by line.  If an error, first entry contains the message.
                std::vector<std::string> & commandStderr,    // [out] the output of stderr.
                const std::vector<std::string> & commandList // [in] command to be run, with one entry per command line word
              )
{
   int link[2];
   int errlink[2];
   
   pid_t pid;
   
   if (pipe(link)==-1) 
   {
      commandOutput.push_back(std::string("Pipe error stdout: ") + strerror(errno));
      return -1;
   }

   if (pipe(errlink)==-1) 
   {
      commandOutput.push_back(std::string("Pipe error stderr: ") + strerror(errno));
      return -1;
   }
   
   if ((pid = fork()) == -1) 
   {
      commandOutput.push_back(std::string("Fork error: ") + strerror(errno));
      return -1;
   }

   if(pid == 0) 
   {
      dup2 (link[1], STDOUT_FILENO);
      close(link[0]);
      close(link[1]);
      
      dup2 (errlink[1], STDERR_FILENO);
      close(errlink[0]);
      close(errlink[1]);
      
      std::vector<const char *>charCommandList( commandList.size()+1, NULL);
      for(int index = 0; index < (int) commandList.size(); ++index)
      {
         charCommandList[index]=commandList[index].c_str();
      }
      execvp( charCommandList[0], const_cast<char**>(charCommandList.data()));
      commandOutput.push_back(std::string("execvp returned: ") + strerror(errno));
      return -1;
   }
   else 
   {
      char commandOutput_c[4096];
      
      wait(NULL);
      
      close(link[1]);
      close(errlink[1]);
      
      int rd;
      if ( (rd = read(link[0], commandOutput_c, sizeof(commandOutput_c))) < 0) 
      {
         commandOutput.push_back(std::string("Read error: ") + strerror(errno));  
         close(link[0]);
         return -1;
      }
      close(link[0]);
      
      std::string line;
      
      commandOutput_c[rd] = '\0';
      std::string commandOutputString(commandOutput_c);
      
      std::istringstream iss(commandOutputString);
      
      while (getline(iss, line)) 
      {
         commandOutput.push_back(line);
      }
      
      //----stderr
      if ( (rd = read(errlink[0], commandOutput_c, sizeof(commandOutput_c))) < 0) 
      {
         commandStderr.push_back(std::string("Read error on stderr: ") + strerror(errno));  
         close(errlink[0]);
         return -1;
      }
      close(errlink[0]);
      
      commandOutput_c[rd] = '\0';
      commandOutputString = commandOutput_c;
      
      std::istringstream iss2(commandOutputString);
      
      while (getline(iss2, line)) 
      {
         commandStderr.push_back(line);
      }
      
      wait(NULL);
      return 0;
   }
}

}//namespace ipc
} //namespace mx

