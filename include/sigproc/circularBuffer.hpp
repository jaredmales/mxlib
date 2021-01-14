/** \file circularBuffer.hpp
  * \author Jared R. Males
  * \brief A circular buffer class
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef sigproc_circularBuffer
#define sigproc_circularBuffer

#include <vector>

namespace mx
{
namespace sigproc 
{
   
/** \defgroup circular_buffer Circular Buffer
  * \brief A circular buffer class
  * 
  * Three options for the circular buffer are provided, each inheriting the same underlying interface.  The choices
  * vary the way in which the wrapping is handled at the end of the storage.  These are:
  * - circularBufferMod: uses the mod operator, this is slowest in all situations and is provided for comparison
  * - circularBufferBranch: uses an if statement (branch) to test for location in memory.  Slower in sequential access, possibly slightly faster in random access.
  * - circularBufferIndex: uses a pre-populated index array, which adds some memory overhead.  Fastest in sequential access, possibly slightly slower in random access.
  * 
  * Benchmarks are clear that circularBufferIndex is fastest for sequential acceess, e.g. one element after the other in sequential order, by something like 30%.  
  * For random access circularBufferBranch is possibly slightly faster, but not enough tests were performed to be conclusive.  
  *  
  * circularBufferMod is always much slower due to use of the `%` operator.
  * 
  * The memory overhead of circularBufferIndex is `2*maxEntries*sizeof(indexT)`, where `maxEntries` is the maximum length of the buffer, and indexT is the 
  * type used for indexing and sizes.
  * 
  * \todo perform circular buffer testing on an isolated core, with one test per executable
  * \todo implement a circular buffer with fixed power-of-two size to test `&` modulo
  * 
  * \ingroup signal_processing 
  */

/// CRTP base class for all circular buffers, providing the underlying memory management and accessors.
/** The CRTP derived classes implement a standard interface based on how they handle wrapping from the end 
  * to the beginning of the buffer.
  *
  * \tparam _derivedT the child class type
  * \tparam _storedT the type stored in teh circular buffer
  * \tparam _indexT the index type, also used for sizes (can be unsigned) 
  * 
  * \ingroup circular_buffer
  */ 
template<typename _derivedT, typename _storedT, typename _indexT>
class circularBufferBase
{
public:
   
   typedef _derivedT derivedT; ///< The child class
   
   typedef _storedT storedT; ///< The type stored in the circular buffer
   typedef _indexT indexT; ///< The index type, also used for sizes
   
protected:
      
   std::vector<storedT> m_buffer;  ///< The circular buffer storage
   
   indexT m_maxEntries {0}; ///< The maximum number of entries to allow in the buffer before wrapping
   indexT m_nextEntry {0}; ///< Index into m_buffer of the next entry.  This is the oldest entry in the buffer.
   
public:
   
   /// Default c'tor
   circularBufferBase();

   /// Sizing constructor
   /** Sets the maximum size of the buffer.  Note that this will not be the size until 
     * a full set of entries have been added to the buffer.
     */ 
   explicit circularBufferBase(indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/);
   
   /// Set the maximum size of the buffer.
   /** Note that this will not be the size until 
     * a full set of entries have been added to the buffer.
     */ 
   void maxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/);
   
   /// Get the maximum size of the buffer.
   /** Note that this is not the current size of the buffer until
     * at least maxEntries have been added.  Use size() to check current size.
     *
     * \returns the maximum size of the buffer, m_maxEntries
     */ 
   indexT maxEntries();

   /// Get the number of entries.
   /** Note that this is the current size.  It will not 
     * be the same as maxEntries() until at least maxEntries
     * have been added.
     */ 
   indexT size();
   
   /// Add the next entry to the circular buffer
   /** Adds the entry (will incur a deep copy) and updates
     * the wrapping system.
     */ 
   void nextEntry( storedT & newEnt /**< [in] the new entry to add to the buffer*/);
   
   /// Returns the index of the next entry
   /** This is particularly useful for accessing entries with the at() function
     * over an unchanging sequence if asynchronous updates are possible.
     */ 
   indexT nextEntry();
   
   /// Get the entry at a given index
   /** `idx=0` is the earliest entry in the buffer. `idx=1` is the one after that.
     * I.e., this counts forward from the oldest data
     * 
     * \returns a reference to the indicated entry in the circular buffer.
     */ 
   storedT & operator[](indexT idx /**< [in] the index of the entry to access*/);
   
   /// Get the entry at a given index relative a fixed start entry 
   /** `idx=0` is the entry at startEntry. `idx=1` is the one after that.
     * I.e., this counts forward from the oldest data
     * 
     * \returns a reference to the indicated entry in the circular buffer.
     */
   storedT & at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
               );
private:
   derivedT & derived()
   {
      return *static_cast<derivedT *>(this);
   }
};

template<class derivedT, typename storedT, typename indexT>
circularBufferBase<derivedT, storedT, indexT>::circularBufferBase()
{
}

template<class derivedT, typename storedT, typename indexT>
circularBufferBase<derivedT, storedT, indexT>::circularBufferBase(indexT maxEnt)
{
   maxEntries(maxEnt);
}

template<class derivedT, typename storedT, typename indexT>
void circularBufferBase<derivedT, storedT, indexT>::maxEntries(indexT maxEnt)
{
   m_buffer.clear();
   m_nextEntry = 0;
   
   m_maxEntries = maxEnt;
   m_buffer.reserve(m_maxEntries);
   
   derived().setMaxEntries(maxEnt);
   
   
}

template<class derivedT, typename storedT, typename indexT>
typename circularBufferBase<derivedT, storedT, indexT>::indexT circularBufferBase<derivedT, storedT, indexT>::maxEntries()
{
   return m_maxEntries;
}

template<class derivedT, typename storedT, typename indexT>
typename circularBufferBase<derivedT, storedT, indexT>::indexT circularBufferBase<derivedT, storedT, indexT>::size()
{
   return m_buffer.size();
}

template<class derivedT, typename storedT, typename indexT>
void circularBufferBase<derivedT, storedT, indexT>::nextEntry(storedT & newEnt)
{
   if( m_buffer.size() < m_maxEntries )
   {
      m_buffer.push_back(newEnt);
      m_nextEntry = 0;
      derived().setWrapStartup();
   }
   else
   {
      m_buffer[m_nextEntry] = newEnt;
      ++m_nextEntry;
      if(m_nextEntry >= m_buffer.size()) m_nextEntry = 0;
      derived().setWrap();
   }
}

template<class derivedT, typename storedT, typename indexT>
indexT circularBufferBase<derivedT, storedT, indexT>::nextEntry()
{
   return m_nextEntry;
}

template<class derivedT, typename storedT, typename indexT>
typename circularBufferBase<derivedT, storedT, indexT>::storedT & circularBufferBase<derivedT, storedT, indexT>::operator[](indexT idx)
{
   return derived().at(m_nextEntry, idx);
}

template<class derivedT, typename storedT, typename indexT>
storedT & circularBufferBase<derivedT, storedT, indexT>::at( indexT refEntry,
                                                             indexT idx
                                                           )
{
   return derived().at(refEntry, idx);
}

/// Circular buffer which wraps with an if statement (branching) [slower, less memory]
/**
  * \ingroup circular_buffer
  */
template<typename _storedT, typename _indexT>
class circularBufferBranch : public circularBufferBase<circularBufferBranch<_storedT,_indexT>, _storedT, _indexT>
{
public:
   
   typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
   typedef _indexT indexT; ///< The index type, also used for sizes
   
protected:
   indexT m_wrapEntry{0};
   
public:
   
   /// Default c'tor
   circularBufferBranch() : circularBufferBase<circularBufferBranch,storedT,indexT>()
   {
   }
   
   /// Sizing constructor
   /** Sets the maximum size of the buffer.  Note that this will not be the size until 
     * a full set of entries have been added to the buffer.
     */
   explicit circularBufferBranch(indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/) : circularBufferBase<circularBufferBranch,storedT,indexT>(maxEnt)
   {
   }
   
   /// Interface implementation for maxEntries.
   /** Resets the wrap entry to 0.
     */
   void setMaxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/)
   {
      m_wrapEntry = 0;
   }

   /// Interface implementation for wrapping setup during the startup phase
   /** This is called before maxEntries have been added.
     */
   void setWrapStartup()
   {
      m_wrapEntry = this->m_buffer.size();
   }
   
   /// Interface implementation for wrapping setup after the startup phase
   /** This is called after maxEntries have been added.
     */
   void setWrap()
   {
      m_wrapEntry = this->m_maxEntries - this->m_nextEntry - 1; 
   }
   
   /// Interface implementation for entry access
   /** Accesses the idx-th element relative to refEntry, using a branch (if-statement) to wrap
     *
     * \returns a reference to the idx-th element
     */
   storedT & at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
               )
   {
      //Option 1: branch
      if(idx > this->m_wrapEntry) return this->m_buffer[refEntry + idx - this->m_buffer.size()];
      else return this->m_buffer[refEntry + idx];
   }
   

};

/// Circular buffer which wraps with the mod opoerator [very slow]
/**
  * \ingroup circular_buffer
  */
template<typename _storedT, typename _indexT>
class circularBufferMod : public circularBufferBase<circularBufferMod<_storedT,_indexT>, _storedT, _indexT>
{
public:
   typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
   typedef _indexT indexT; ///< The index type, also used for sizes
   
   /// Default c'tor
   circularBufferMod() : circularBufferBase<circularBufferMod,storedT,indexT>()
   {
   }
   
   /// Sizing constructor
   /** Sets the maximum size of the buffer.  Note that this will not be the size until 
     * a full set of entries have been added to the buffer.
     */
   explicit circularBufferMod(indexT maxEnt  /**< [in] the maximum number of entries this buffer will hold*/) : circularBufferBase<circularBufferMod,storedT,indexT>(maxEnt)
   {
   }
   
   /// Interface implementation for maxEntries.
   /** A no-op
     */
   void setMaxEntries( indexT maxEnt  /**< [in] the maximum number of entries this buffer will hold*/)
   {
   }

   /// Interface implementation for wrapping setup during the startup phase
   /** This is called before maxEntries have been added.
     */
   void setWrapStartup()
   {
   }
   
   /// Interface implementation for wrapping setup after the startup phase
   /** This is called after maxEntries have been added.
     */
   void setWrap()
   {
   }
   
   /// Interface implementation for entry access
   /** Accesses the idx-th element relative to refEntry, using the mod operator to wrap
     *
     * \returns a reference to the idx-th element
     */
   storedT & at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
               )
   {
      return this->m_buffer[(refEntry + idx) % this->m_buffer.size()];
   }
   

};

/// Circular buffer which wraps with a pre-populated indices array [generally fastest]
/**
  * \ingroup circular_buffer
  */
template<typename _storedT, typename _indexT>
class circularBufferIndex : public circularBufferBase<circularBufferIndex<_storedT,_indexT>, _storedT, _indexT>
{
public:
   typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
   typedef _indexT indexT; ///< The index type, also used for sizes
   
protected:

   std::vector<size_t> m_indices; ///< Vector of indices for fast indexing into parent's m_buffer
     
public:
   
   /// Default c'tor
   circularBufferIndex() : circularBufferBase<circularBufferIndex,storedT,indexT>()
   {
   }
   
   /// Sizing constructor
   /** Sets the maximum size of the buffer.  Note that this will not be the size until 
     * a full set of entries have been added to the buffer.
     */
   explicit circularBufferIndex(indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/) : circularBufferBase<circularBufferIndex,storedT,indexT>(maxEnt)
   {
   }
   
   /// Interface implementation for maxEntries.
   /** Resizes and populates the indices array.
     */
   void setMaxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/)
   {
      m_indices.resize(2*this->m_maxEntries);
      for(size_t i=0; i< this->m_maxEntries; ++i)
      {
         m_indices[i] = i;
         m_indices[this->m_maxEntries + i] = i;
      }
   }

   /// Interface implementation for wrapping setup during the startup phase
   /** This is called before maxEntries have been added.
     */
   void setWrapStartup()
   {
   }
   
   /// Interface implementation for wrapping setup after the startup phase
   /** This is called after maxEntries have been added.
     */
   void setWrap()
   {
   }
   
   /// Interface implementation for entry access
   /** Accesses the idx-th element relative to refEntry, using the pre-populated indices to wrap
     *
     * \returns a reference to the idx-th element
     */
   storedT & at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
               )
   {
      return this->m_buffer[m_indices[refEntry + idx]];
   }
   

};

} //namespace sigproc
} //namespace circularBuffer

#endif //sigproc_circularBuffer
