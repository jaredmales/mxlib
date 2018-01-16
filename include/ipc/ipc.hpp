/** \file ipc.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib interprocess communication (IPC) tools
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

#ifndef ipc_ipc_hpp
#define ipc_ipc_hpp



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


///@}

#endif //mx_IPC_hpp


