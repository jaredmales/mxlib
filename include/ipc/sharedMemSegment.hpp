/** \file sharedMemSegment.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib shared memory facility
  * \ingroup IPC_sharedmem
  * \ingroup IPC
  * 
*/

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018, 2021 Jared R. Males (jaredmales@gmail.com)
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

#ifndef ipc_sharedMemSegment_hpp
#define ipc_sharedMemSegment_hpp

#include <sys/shm.h>
#include <stdint.h>
#include <stdio.h>

#include "ipc.hpp"

namespace mx
{
namespace ipc 
{

/** \addtogroup IPC_sharedmem
  * @{
  */

/// A c++ class to manage a System V shared memory segment with memory mapping
/** The associated functions create and/or attach to a shared memory segment.  If a segment is created, a block of
  * size sizeof(uintptr_t) is reserved where the address is stored.  This address is used when subsequently
  * attaching to the segment so that pointers stored in the block are valid, etc.
  * 
  */ 
class sharedMemSegment
{
public:
   ///The path to use for key creation
   char key_path[MX_IPC_KEYLEN];
   
   ///The id to use for key creation
   int key_id;
   
   ///The shared memory key 
   key_t key;

   ///The shared memory id associated with the key
   int shmemid;

   ///The base address of the segment
   void * addr;
   
   ///The size of the segment
   size_t size;

   ///Flag indicating whether or not the segment is attached
   int attached;

public:
   
   ///Initialize the class
   void initialize();

   ///Set the key
   /** If path == 0, then the key is set directly to id, without using ftok.
     * 
     * \param path is the full path to use in a call to ftok
     * \param id is the id number to use in a call to ftok
     * 
     * \returns the key value, which is also set in the msgq
     */
   key_t setKey( const char * path, 
                 const int id
               );
   
   ///Create and attach to the segment
   /** A segment of size = sz + sizeof(uintptr_t) is actually created.
     * 
     * \note key must be set before calling this function.
     * 
     * \param sz the size of the segment to create
     * 
     */
   int create(size_t sz);

   ///Attach to a segment without creating it.
   /** If donot_set_addr == false, then the address is mapped to match that stored in the first uintptr_t of the segment.
     * 
     * \param donot_set_addr if > 0 (true) then the address is not mapped to match that stored in the segment
     */
   int attach( bool donot_set_addr = false);
   
   ///Detach from the segment
   /**
     *
     */ 
   int detach();

};

/// @}

}//namespace ipc 
}//namespace mx

#endif //ipc_sharedMemSegment_hpp

