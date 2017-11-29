/** \file sharedmem_segment.c
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions for managing a shared memory segment.
  * 
*/

#include "sharedmem_segment.h"

void sharedmem_segment_initialize(sharedmem_segment * seg)
{
   seg->key_id = 0;
   seg->key = 0;
   seg->shmemid = 0;
   seg->addr = 0;
   seg->size = 0;
   
   seg->attached = 0;

}

key_t sharedmem_segment_set_key(sharedmem_segment * seg, const char * path, const int id)
{
   if(path != 0)
   {
      strncpy(seg->key_path, path, MX_IPC_KEYLEN);
      seg->key_id = id;
   
      seg->key = ftok(seg->key_path, seg->key_id);
   }
   else
   {
      seg->key_path[0] = 0;
      seg->key_id = id;
      
      seg->key = id;
   }
   return seg->key;
}


int sharedmem_segment_create(sharedmem_segment * seg, size_t sz)
{
   if( (seg->shmemid = shmget(seg->key, sz + 1*sizeof(uintptr_t), IPC_CREAT | 0666))<0)
   {
      //If it failed, try to remove the shmem block and recreate it.
      seg->shmemid  = shmget(seg->key, 1, 0666);
      if(shmctl(seg->shmemid, IPC_RMID, 0) < 0)
      {
         fprintf(stderr, "Could not remove shared memory with key %i\n", seg->key);
         return -1;
      }
          
      //removal successful, now try to create again
      if((seg->shmemid = shmget(seg->key, sz + 1*sizeof(uintptr_t), IPC_CREAT | 0666))<0)
      {
         fprintf(stderr, "Could not create shared memory with key %i\n", seg->key);
         return -1;
      }
   }
    
   sharedmem_segment_attach(seg, 1);
      
   //Since we created this segment, we set the address field.
   *((uintptr_t *) seg->addr) = (uintptr_t) seg->addr;
   
   return 0;
}

int sharedmem_segment_attach(sharedmem_segment *seg, int donot_set_addr)
{
   struct shmid_ds shmstats;
   void * new_addr;
   
   if(seg->shmemid == 0)
   {
      if((seg->shmemid = shmget(seg->key, 0, 0666))<0)
      {
         fprintf(stderr, "Could not remove shared memory with key %i\n", seg->key);
         return -1;
      }
   }
   
   if ((new_addr = shmat(seg->shmemid, 0, 0)) == (char *) -1) 
   {
      fprintf(stderr, "Could not attach to shared memory with key %i\n", seg->key);
      return -1;
   }
   
   seg->attached = 1;
   
   if (shmctl(seg->shmemid, IPC_STAT, &shmstats) < 0)
   {
      fprintf(stderr, "Could not get shared memory stats with key %i\n", seg->key);
      return -1;
   }

   seg->size = shmstats.shm_segsz;
   

   //Here we first read in the address from the first unitptr_t size block
   //then detach, then re-attach specifying the address.   
   if(!donot_set_addr)
   {
      seg->addr =  (void *) *((uintptr_t *) new_addr); //read the address from the segment itself
      
      //now detach
      if(shmdt(new_addr) != 0)
      {
         fprintf(stderr, "Unable to detach from shared memory\n");
         return -1;
      }
      
      seg->attached = 0;
      
      //and then re-attach, but now specifying an address
      if ((new_addr = shmat(seg->shmemid, seg->addr, 0)) == (char *) -1) 
      {
         fprintf(stderr, "Could not re-attach shared memory with key %i\n",seg->key);
         return -1;
      }
   }
   else
   {
      seg->addr = new_addr;
   }
   
   seg->attached = 1;
   
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

