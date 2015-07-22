/** \file sharedMemSegment.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions for a class to manage a shared memory segment.
  * 
*/

#include "sharedMemSegment"

namespace mx
{
   
sharedMemSegment::sharedMemSegment()
{
   sharedmem_segment_initialize(this);
}

sharedMemSegment::sharedMemSegment(const std::string &path, const int & id)
{
   sharedmem_segment_initialize(this);
   setKey(path, id);
}

key_t sharedMemSegment::setKey(const std::string &path, const int & id)
{
   return sharedmem_segment_set_key(this, path.c_str(), id);
}


int sharedMemSegment::create(size_t sz)
{
   return sharedmem_segment_create(this, sz);
}


int sharedMemSegment::attach(int doNotSetAddr)
{
   return sharedmem_segment_attach(this, doNotSetAddr);
}

int sharedMemSegment::detach()
{
   return sharedmem_segment_detach(this);
}

void * sharedMemSegment::getRawAddr()
{
   return addr;
}

void * sharedMemSegment::getAddr()
{
   return (void *) ((uintptr_t) addr + sizeof(uintptr_t));
}

size_t sharedMemSegment::getRawSize()
{
   return size;
}

size_t sharedMemSegment::getSize()
{
   return size - sizeof(uintptr_t);
}



} //namespace mx



// int main()
// {
//    sharedMemSegment seg(5000);
//    
//    //seg.createSegment(100);
//    
//    //seg.detachSegment();
//    
//    seg.attachSegment();
//    
//    std::cout << seg.getRawAddr() << "\n";
//    std::cout << seg.getAddr() << "\n";
//    std::cout << (uintptr_t) seg.getAddr() - (uintptr_t) seg.getRawAddr() << "\n";
//    return 0;
// }
