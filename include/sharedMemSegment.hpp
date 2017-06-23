/** \file sharedmem_segment.h
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib shared memory facility
  * \ingroup IPC_sharedmem
  * \ingroup IPC
  * 
*/

#ifndef __sharedMemSegment_hpp__
#define __sharedMemSegment_hpp__

#include <sys/shm.h>
#include <stdint.h>
#include <stdio.h>

//#include "IPC.h"


/** \addtogroup IPC_sharedmem
  * @{
  */

/// A c++ class to manage a shared memory segment with memory mapping
/** The associated functions create and/or attach to a shared memory segment.  If a segment is created, a block of
  * size sizeof(uintptr_t) is reserved where the address is stored.  This address is used when subsequently
  * attaching to the segment so that pointers stored in the block are valid, etc.
  * 
  */ 
class sharedMemSegment
{
protected:
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
     * \param seg is the \ref sharedmem_segment to use
     *
     */ 
   int detach();

};

inline
void sharedMemSegment::initialize()
{
   key_id = 0;
   key = 0;
   shmemid = 0;
   addr = 0;
   size = 0;
   
   attached = 0;

}

inline
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

inline
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

inline
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

int sharedmem_segment_detach(sharedmem_segment *seg)
{
   
   if(!seg->attached) return 0;
   
   if(seg->addr == 0) return 0;

   //now detach
   if(shmdt(seg->addr) != 0)
   {
      fprintf(stderr, "Unable to detach from shared memory\n");
      return -1;
   }      
   
   return 0;
}
   
/// @}


#endif //__sharedMemSegment_hpp__

