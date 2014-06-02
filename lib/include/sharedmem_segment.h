/** \file sharedmem_segment.h
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib c shared memory facilities
  * \ingroup IPC_sharedmem
  * \ingroup IPC
  * 
*/

#ifndef __sharedmem_segment_h__
#define __sharedmem_segment_h__

#include <sys/shm.h>
#include <stdint.h>
#include <stdio.h>

#include "IPC.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** \addtogroup IPC_sharedmem
  * @{
  */

/// A c structure to manage a shared memory segment with memory mapping
/** The associated functions create and/or attach to a shared memory segment.  If a segment is created, a block of
  * size sizeof(uintptr_t) is reserved where the address is stored.  This address is used when subsequently
  * attaching to the segment so that pointers stored in the block are valid, etc.
  * 
  * \sa Functions for working with sharedmem_segment include: \ref sharedmem_segment_initialize, 
  * \ref sharedmem_segment_set_key, \ref sharedmem_segment_create, \ref sharedmem_segment_attach, 
  * and \ref sharedmem_segment_detach
  * 
  */ 
typedef struct
{
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

} sharedmem_segment;

///Initialize a \ref sharedmem_segment structure
/**
  * \param seg is the \ref sharedmem_segment to initialize
  * 
  */
void sharedmem_segment_initialize(sharedmem_segment * seg);

///Set the key for a \ref sharedmem_segment structure
/** If path == 0, then the key is set directly to id, without using ftok.
  * 
  * \param seg is a pointer to a \ref sharedmem_segment structure
  * \param path is the full path to use in a call to ftok
  * \param id is the id number to use in a call to ftok
  * 
  * \returns the key value, which is also set in the msgq
  */
key_t sharedmem_segment_set_key(sharedmem_segment * seg, const char * path, const int id);
   
///Create and attach to the segment described by a \ref sharedmem_segment
/** A segment of size = sz + sizeof(uintptr_t) is actually created.
  * 
  * \note seg->key must be set before calling this function (see \ref sharedmem_segment_set_key).
  * 
  * \param seg is the \ref sharedmem_segment meta to use
  * \param sz the size of the segment to create
  * 
  */
int sharedmem_segment_create(sharedmem_segment * seg, size_t sz);

///Attach to a segment without creating it.
/** If donot_set_addr == 0, then the address is mapped to match that stored in the first uintptr_t of the segment.
  * 
  * \param seg is the \ref sharedmem_segment to use
  * \param donot_set_addr if > 0 (true) then the address is not mapped to match that stored in the segment
  */
int sharedmem_segment_attach(sharedmem_segment *seg, int donot_set_addr);
   
///Detach from the segment
/**
  * \param seg is the \ref sharedmem_segment to use
  *
  */ 
int sharedmem_segment_detach(sharedmem_segment *seg);
   
/// @}

#ifdef __cplusplus
} //extern "C"
#endif

#endif //__sharedmem_segment_h__

