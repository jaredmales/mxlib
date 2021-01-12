/** \file sharedMemSegment.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions for the mxlib shared memory facility
  * \ingroup IPC_sharedmem
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

#include "ipc/sharedMemSegment.hpp"

namespace mx
{
namespace ipc 
{

void sharedMemSegment::initialize()
{
   key_id = 0;
   key = 0;
   shmemid = 0;
   addr = 0;
   size = 0;
   
   attached = 0;

}

key_t sharedMemSegment::setKey( const char * path, 
                                const int id
                              )
{
   if(path != 0)
   {
      strncpy( key_path, path, MX_IPC_KEYLEN);
      key_id = id;
   
      key = ftok( key_path, key_id);
   }
   else
   {
      key_path[0] = 0;
      key_id = id;
      
      key = id;
   }
   
   return key;
}

int sharedMemSegment::create( size_t sz )
{
   if( ( shmemid = shmget( key, sz + 1*sizeof(uintptr_t), IPC_CREAT | 0666))<0)
   {
      //If it failed, try to remove the shmem block and recreate it.
      shmemid  = shmget( key, 1, 0666);
      if(shmctl( shmemid, IPC_RMID, 0) < 0)
      {
         fprintf(stderr, "Could not remove shared memory with key %i\n", key);
         return -1;
      }
          
      //removal successful, now try to create again
      if(( shmemid = shmget( key, sz + 1*sizeof(uintptr_t), IPC_CREAT | 0666))<0)
      {
         fprintf(stderr, "Could not create shared memory with key %i\n", key);
         return -1;
      }
   }
    
   attach( true );
      
   //Since we created this segment, we set the address field.
   *((uintptr_t *) addr) = (uintptr_t) addr;
   
   return 0;
}

int sharedMemSegment::attach( bool donot_set_addr)
{
   struct shmid_ds shmstats;
   void * new_addr;
   
   if( shmemid == 0 )
   {
      if(( shmemid = shmget(key, 0, 0666))<0)
      {
         fprintf(stderr, "Could not remove shared memory with key %i\n", key);
         return -1;
      }
   }
   
   if ((new_addr = shmat( shmemid, 0, 0)) == (char *) -1) 
   {
      fprintf(stderr, "Could not attach to shared memory with key %i\n", key);
      return -1;
   }
   
   attached = 1;
   
   if (shmctl( shmemid, IPC_STAT, &shmstats) < 0)
   {
      fprintf(stderr, "Could not get shared memory stats with key %i\n", key);
      return -1;
   }

   size = shmstats.shm_segsz;
   

   //Here we first read in the address from the first unitptr_t size block
   //then detach, then re-attach specifying the address.   
   if(!donot_set_addr)
   {
      addr =  (void *) *((uintptr_t *) new_addr); //read the address from the segment itself
      
      //now detach
      if(shmdt(new_addr) != 0)
      {
         fprintf(stderr, "Unable to detach from shared memory\n");
         return -1;
      }
      
      attached = 0;
      
      //and then re-attach, but now specifying an address
      if ((new_addr = shmat(shmemid, addr, 0)) == (char *) -1) 
      {
         fprintf(stderr, "Could not re-attach shared memory with key %i\n",key);
         return -1;
      }
   }
   else
   {
      addr = new_addr;
   }
   
   attached = 1;
   
   return 0;
}

int sharedMemSegment::detach()
{
   
   if(attached) return 0;
   
   if(addr == 0) return 0;

   //now detach
   if(shmdt(addr) != 0)
   {
      fprintf(stderr, "Unable to detach from shared memory\n");
      return -1;
   }      
   
   return 0;
}
   
}//namespace ipc 
}//namespace mx
