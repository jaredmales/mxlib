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

#include <atomic>
#include <vector>

namespace mx
{
namespace sigproc
{

/** \defgroup circular_buffer Circular Buffer
 * \brief A circular buffer class
 *
 * Three options for the circular buffer are provided, each inheriting the same underlying interface.  The choices
 * vary the way in which the wrapping is handled at the beginning and end of the storage.  These are:
 * - circularBufferMod: uses the mod operator, this is slowest in all situations and is provided for comparison
 * - circularBufferBranch: uses an if statement (branch) to test for location in memory.  Slower in sequential access,
 * possibly slightly faster in random access.
 * - circularBufferIndex: uses a pre-populated index array, which adds some memory overhead.  Fastest in sequential
 * access, possibly slightly slower in random access.
 *
 * Benchmarks are clear that circularBufferIndex is fastest for sequential acceess, e.g. one element after the other in
 * sequential order, by something like 30%. For random access circularBufferBranch is possibly slightly faster, but not
 * enough tests were performed to be conclusive.
 *
 * circularBufferMod is always much slower due to use of the `%` operator.  Don't use this for real work.
 *
 * The memory overhead of circularBufferIndex is `3*maxEntries*sizeof(indexT)`, where `maxEntries` is the maximum length
 * of the buffer, and indexT is the type used for indexing and sizes.  The factor of 3 enables negative offsets, i.e.
 * for traversal in reverse.
 *
 * \note circularBufferIndex will not properly wrap until it is full.  Attempts to access wrapped entries prior
 *       to inserting maxEntries entries will segfault.
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
 * \tparam _indexT the index type, also used for sizes (must be signed)
 *
 * \ingroup circular_buffer
 */
template <typename _derivedT, typename _storedT, typename _indexT>
class circularBufferBase
{
  public:
    typedef _derivedT derivedT; ///< The child class

    typedef _storedT storedT;   ///< The type stored in the circular buffer
    typedef _indexT indexT;     ///< The index type, also used for sizes

    static_assert( std::is_signed_v<_indexT> == true, "circularBuffer indexT must be signed" );

  protected:
    std::vector<storedT> m_buffer; ///< The circular buffer storage

    indexT m_maxEntries{ 0 };      ///< The maximum number of entries to allow in the buffer before wrapping
    indexT m_nextEntry{ 0 };       ///< Index into m_buffer of the next entry.  This is the oldest entry in the buffer.
    indexT m_latest{ 0 };          ///< Index into m_buff of the latest entry.  This is the newest entry in the buffer.

    uint64_t m_mono{ 0 };          /**< A monotonic counter, which is incremented for each new entry,
                                        and reset on call to maxEntries.*/

  public:
    /// Default c'tor
    circularBufferBase();

    /// Sizing constructor
    /** Sets the maximum size of the buffer.  Note that this will not be the size until
     * a full set of entries have been added to the buffer.
     */
    explicit circularBufferBase( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ );

    /// Set the maximum size of the buffer.
    /** Note that this will not be the size until
     * a full set of entries have been added to the buffer.
     */
    void maxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ );

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
    indexT size() const;

    /// Add the next entry to the circular buffer
    /** Adds the entry (will incur a deep copy) and updates
     * the wrapping system.
     */
    void nextEntry( const storedT &newEnt /**< [in] the new entry to add to the buffer*/ );

    /// Move to the next entry in the circular buffer
    /** If not yet at maxEntries entries, this emplaces a new default constructed entry.
     *  Otherwise it increments counters so that the old latest entry is now the newest entry.
     */
    void nextEntry();

    /// Returns the index of the earliest entry
    /** This is particularly useful for accessing entries with the at() function
     * over an unchanging sequence if asynchronous updates are possible.
     */
    indexT earliest();

    /// Returns the index of the latest entry
    /**
     */
    indexT latest();

    uint64_t mono();

    /// Get the entry at a given index
    /** `idx=0` is the earliest entry in the buffer. `idx=1` is the one after that.
     * I.e., this counts forward from the oldest data
     *
     * \returns a reference to the indicated entry in the circular buffer.
     */
    storedT &operator[]( indexT idx /**< [in] the index of the entry to access*/ );

    const storedT &operator[]( indexT idx /**< [in] the index of the entry to access*/ ) const;

    /// Get the entry at a given index relative a fixed reference entry
    /** `idx=0` is the entry at refEntry. `idx=1` is the one after that.
     * I.e., this counts forward from the oldest data
     *
     * \returns a reference to the indicated entry in the circular buffer.
     */
    storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
    );

    /// Get the entry at a given index relative a fixed start entry, const version
    /** `idx=0` is the entry at refEntry. `idx=1` is the one after that.
     * I.e., this counts forward from the oldest data
     *
     * \overload
     *
     * \returns a const reference to the indicated entry in the circular buffer.
     */
    const storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                       indexT idx       ///< [in] the index of the entry to access
    ) const;

  private:
    derivedT &derived()
    {
        return *static_cast<derivedT *>( this );
    }

    const derivedT &derived() const
    {
        return *static_cast<const derivedT *>( this );
    }
};

template <class derivedT, typename storedT, typename indexT>
circularBufferBase<derivedT, storedT, indexT>::circularBufferBase()
{
}

template <class derivedT, typename storedT, typename indexT>
circularBufferBase<derivedT, storedT, indexT>::circularBufferBase( indexT maxEnt )
{
    maxEntries( maxEnt );
}

template <class derivedT, typename storedT, typename indexT>
void circularBufferBase<derivedT, storedT, indexT>::maxEntries( indexT maxEnt )
{
    m_buffer.clear();
    m_nextEntry = 0;
    m_latest = 0;

    m_maxEntries = maxEnt;
    m_buffer.reserve( m_maxEntries );

    derived().setMaxEntries( maxEnt );

    m_mono = 0;
}

template <class derivedT, typename storedT, typename indexT>
indexT circularBufferBase<derivedT, storedT, indexT>::maxEntries()
{
    return m_maxEntries;
}

template <class derivedT, typename storedT, typename indexT>
indexT circularBufferBase<derivedT, storedT, indexT>::size() const
{
    return m_buffer.size();
}

template <class derivedT, typename storedT, typename indexT>
void circularBufferBase<derivedT, storedT, indexT>::nextEntry( const storedT &newEnt )
{
    if( m_buffer.size() < m_maxEntries )
    {
        m_buffer.push_back( newEnt );
        m_latest = m_buffer.size() - 1;
        m_nextEntry = 0;
    }
    else
    {
        m_buffer[m_nextEntry] = newEnt;
        m_latest = m_nextEntry;
        ++m_nextEntry;
        if( m_nextEntry > m_buffer.size() - 1 )
        {
            m_nextEntry = 0;
        }
    }

    ++m_mono;
}

template <class derivedT, typename storedT, typename indexT>
void circularBufferBase<derivedT, storedT, indexT>::nextEntry()
{
    if( m_buffer.size() < m_maxEntries )
    {
        m_buffer.emplace_back();
        m_latest = m_buffer.size() - 1;
        m_nextEntry = 0;
    }
    else
    {
        m_latest = m_nextEntry;
        ++m_nextEntry;
        if( m_nextEntry > m_buffer.size() - 1 )
        {
            m_nextEntry = 0;
        }
    }

    ++m_mono;
}

template <class derivedT, typename storedT, typename indexT>
indexT circularBufferBase<derivedT, storedT, indexT>::earliest()
{
    return m_nextEntry;
}

template <class derivedT, typename storedT, typename indexT>
indexT circularBufferBase<derivedT, storedT, indexT>::latest()
{
    return m_latest;
}

template <class derivedT, typename storedT, typename indexT>
uint64_t circularBufferBase<derivedT, storedT, indexT>::mono()
{
    return m_mono;
}

template <class derivedT, typename storedT, typename indexT>
storedT &circularBufferBase<derivedT, storedT, indexT>::operator[]( indexT idx )
{
    return derived().at( m_nextEntry, idx );
}

template <class derivedT, typename storedT, typename indexT>
const storedT &circularBufferBase<derivedT, storedT, indexT>::operator[]( indexT idx ) const
{
    return derived().at( m_nextEntry, idx );
}

template <class derivedT, typename storedT, typename indexT>
storedT &circularBufferBase<derivedT, storedT, indexT>::at( indexT refEntry, indexT idx )
{
    return derived().at( refEntry, idx );
}

template <class derivedT, typename storedT, typename indexT>
const storedT &circularBufferBase<derivedT, storedT, indexT>::at( indexT refEntry, indexT idx ) const
{
    return derived().at( refEntry, idx );
}

/// Circular buffer which wraps with an if statement (branching) [faster than mod, less memory than index]
/**
 * \ingroup circular_buffer
 */
template <typename _storedT, typename _indexT>
class circularBufferBranch : public circularBufferBase<circularBufferBranch<_storedT, _indexT>, _storedT, _indexT>
{
  public:
    typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
    typedef _indexT indexT;   ///< The index type, also used for sizes

  public:
    /// Default c'tor
    circularBufferBranch() : circularBufferBase<circularBufferBranch, storedT, indexT>()
    {
    }

    /// Sizing constructor
    /** Sets the maximum size of the buffer.  Note that this will not be the size until
     * a full set of entries have been added to the buffer.
     */
    explicit circularBufferBranch( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
        : circularBufferBase<circularBufferBranch, storedT, indexT>( maxEnt )
    {
    }

    /// Interface implementation for maxEntries.
    /** Resets the wrap entry to 0.
     */
    void setMaxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
    {
    }

    /// Interface implementation for entry access
    /** Accesses the idx-th element relative to refEntry, using a branch (if-statement) to wrap
     *
     * \returns a reference to the idx-th element
     */
    storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
    )
    {
        indexT _idx = refEntry + idx;
        if( _idx < 0 )
        {
            return this->m_buffer[_idx + this->m_buffer.size()];
        }
        else if( _idx > this->m_buffer.size() - 1 )
        {
            return this->m_buffer[_idx - this->m_buffer.size()];
        }
        else
        {
            return this->m_buffer[_idx];
        }
    }

    /// Interface implementation for entry access, const version
    /** Accesses the idx-th element relative to refEntry, using a branch (if-statement) to wrap
     *
     * \overload
     *
     * \returns a const reference to the idx-th element
     */
    const storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                       indexT idx       ///< [in] the index of the entry to access
    ) const
    {
        indexT _idx = refEntry + idx;
        if( _idx < 0 )
        {
            return this->m_buffer[_idx + this->m_buffer.size()];
        }
        else if( _idx > this->m_buffer.size() - 1 )
        {
            return this->m_buffer[_idx - this->m_buffer.size()];
        }
        else
        {
            return this->m_buffer[_idx];
        }
    }
};

/// Circular buffer which wraps with the mod operator [very slow]
/**
 * \ingroup circular_buffer
 */
template <typename _storedT, typename _indexT>
class circularBufferMod : public circularBufferBase<circularBufferMod<_storedT, _indexT>, _storedT, _indexT>
{
  public:
    typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
    typedef _indexT indexT;   ///< The index type, also used for sizes

    /// Default c'tor
    circularBufferMod() : circularBufferBase<circularBufferMod, storedT, indexT>()
    {
    }

    /// Sizing constructor
    /** Sets the maximum size of the buffer.  Note that this will not be the size until
     * a full set of entries have been added to the buffer.
     */
    explicit circularBufferMod( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
        : circularBufferBase<circularBufferMod, storedT, indexT>( maxEnt )
    {
    }

    /// Interface implementation for maxEntries.
    /** A no-op
     */
    void setMaxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
    {
    }

    /// Interface implementation for entry access
    /** Accesses the idx-th element relative to refEntry, using the mod operator to wrap
     *
     * \returns a reference to the idx-th element
     */
    storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
    )
    {
        return this->m_buffer[( refEntry + idx ) % this->m_buffer.size()];
    }

    /// Interface implementation for entry access, const version
    /** Accesses the idx-th element relative to refEntry, using the mod operator to wrap
     *
     * \overload
     *
     * \returns a const reference to the idx-th element
     */
    const storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                       indexT idx       ///< [in] the index of the entry to access
    ) const
    {
        return this->m_buffer[( refEntry + idx ) % this->m_buffer.size()];
    }
};

/// Circular buffer which wraps with a pre-populated indices array [generally fastest]
/**
 * \ingroup circular_buffer
 */
template <typename _storedT, typename _indexT, bool fixed_size = false>
class circularBufferIndex;

template <typename _storedT, typename _indexT>
class circularBufferIndex<_storedT, _indexT, false>
    : public circularBufferBase<circularBufferIndex<_storedT, _indexT, false>, _storedT, _indexT>
{
  public:
    typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
    typedef _indexT indexT;   ///< The index type, also used for sizes

  protected:
    std::vector<size_t> m_indices; ///< Vector of indices for fast indexing into parent's m_buffer

  public:
    /// Default c'tor
    circularBufferIndex() : circularBufferBase<circularBufferIndex, storedT, indexT>()
    {
    }

    /// Sizing constructor
    /** Sets the maximum size of the buffer.  Note that this will not be the size until
     * a full set of entries have been added to the buffer.
     */
    explicit circularBufferIndex( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
        : circularBufferBase<circularBufferIndex, storedT, indexT>( maxEnt )
    {
    }

    /// Interface implementation for maxEntries.
    /** Resizes and populates the indices array.
     */
    void setMaxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
    {
        m_indices.resize( 3 * this->m_maxEntries );
        for( size_t i = 0; i < this->m_maxEntries; ++i )
        {
            m_indices[i] = i;
            m_indices[this->m_maxEntries + i] = i;
            m_indices[2 * this->m_maxEntries + i] = i;
        }
    }

    /// Interface implementation for entry access
    /** Accesses the idx-th element relative to refEntry, using the pre-populated indices to wrap
     *
     * \returns a reference to the idx-th element
     */
    storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                 indexT idx       ///< [in] the index of the entry to access
    )
    {
        return this->m_buffer[m_indices[this->m_maxEntries + refEntry + idx]];
    }

    /// Interface implementation for entry access, const version
    /** Accesses the idx-th element relative to refEntry, using the pre-populated indices to wrap
     *
     * \overload
     *
     * \returns a const reference to the idx-th element
     */
    const storedT &at( indexT refEntry, ///< [in] the entry to start counting from
                       indexT idx       ///< [in] the index of the entry to access
    ) const
    {
        return this->m_buffer[this->m_maxEntries + m_indices[refEntry + idx]];
    }
};

/// Fixed-size circular buffer specialization for single-producer/single-consumer use.
/**
 * \ingroup circular_buffer
 */
template <typename _storedT, typename _indexT>
class circularBufferIndex<_storedT, _indexT, true>
{
  public:
    typedef _storedT storedT; ///< The maximum number of entries to allow in the buffer before wrapping
    typedef _indexT indexT;   ///< The index type, also used for sizes

    static_assert( std::is_signed_v<_indexT> == true, "circularBuffer indexT must be signed" );
    static_assert( std::is_trivially_copyable_v<_storedT> == true,
                   "fixed_size circularBufferIndex requires trivially copyable storedT" );

    /// Snapshot of the readable state of the buffer.
    struct snapshotT
    {
        indexT earliest{ 0 };     ///< Earliest readable slot in the current snapshot.
        indexT latest{ 0 };       ///< Latest published slot in the current snapshot.
        indexT validEntries{ 0 }; ///< Number of currently valid entries.
        indexT maxEntries{ 0 };   ///< Fixed capacity.
        uint64_t mono{ 0 };       ///< Publication generation for retry validation.
        bool full{ false };       ///< Whether the buffer has reached fixed capacity.
    };

  protected:
    std::vector<storedT> m_buffer;           ///< The circular buffer storage.
    std::vector<size_t> m_indices;           ///< Vector of indices for fast indexing into m_buffer.

    std::atomic<indexT> m_nextEntry{ 0 };    ///< The slot that will be written next by the producer.
    std::atomic<indexT> m_latest{ 0 };       ///< The latest slot published by the producer.
    std::atomic<indexT> m_validEntries{ 0 }; ///< The number of valid published entries.
    std::atomic<uint64_t> m_mono{ 0 };       ///< A monotonic publication counter.

    indexT m_maxEntries{ 0 };                ///< Fixed capacity of the buffer.

  public:
    /// Default c'tor.
    circularBufferIndex() = default;

    /// Sizing constructor.
    explicit circularBufferIndex( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
    {
        maxEntries( maxEnt );
    }

    /// Set the maximum size of the buffer.
    void maxEntries( indexT maxEnt /**< [in] the maximum number of entries this buffer will hold*/ )
    {
        m_buffer.clear();
        m_indices.clear();

        m_maxEntries = maxEnt;

        if( m_maxEntries > 0 )
        {
            m_buffer.resize( m_maxEntries );
            m_indices.resize( 3 * m_maxEntries );

            for( size_t i = 0; i < static_cast<size_t>( m_maxEntries ); ++i )
            {
                m_indices[i] = i;
                m_indices[m_maxEntries + i] = i;
                m_indices[2 * m_maxEntries + i] = i;
            }
        }

        m_nextEntry.store( 0, std::memory_order_release );
        m_latest.store( 0, std::memory_order_release );
        m_validEntries.store( 0, std::memory_order_release );
        m_mono.store( 0, std::memory_order_release );
    }

    /// Get the fixed capacity of the buffer.
    indexT maxEntries()
    {
        return m_maxEntries;
    }

    /// Get the number of valid entries.
    indexT size() const
    {
        return m_validEntries.load( std::memory_order_acquire );
    }

    /// Add the next entry to the circular buffer.
    void nextEntry( const storedT &newEnt /**< [in] the new entry to add to the buffer*/ )
    {
        if( m_maxEntries <= 0 )
        {
            return;
        }

        indexT next = m_nextEntry.load( std::memory_order_relaxed );

        m_buffer[next] = newEnt;

        indexT latest = next;
        ++next;
        if( next >= m_maxEntries )
        {
            next = 0;
        }

        indexT validEntries = m_validEntries.load( std::memory_order_relaxed );
        if( validEntries < m_maxEntries )
        {
            ++validEntries;
        }

        m_latest.store( latest, std::memory_order_relaxed );
        m_validEntries.store( validEntries, std::memory_order_relaxed );
        m_nextEntry.store( next, std::memory_order_relaxed );
        m_mono.fetch_add( 1, std::memory_order_release );
    }

    /// Move to the next entry using a default constructed value.
    void nextEntry()
    {
        nextEntry( storedT{} );
    }

    /// Returns the index of the earliest entry.
    indexT earliest() const
    {
        snapshotT sn = snapshot();
        return sn.earliest;
    }

    /// Returns the index of the latest entry.
    indexT latest() const
    {
        return m_latest.load( std::memory_order_acquire );
    }

    /// Returns the publication counter.
    uint64_t mono() const
    {
        return m_mono.load( std::memory_order_acquire );
    }

    /// Get the entry at a given index.
    storedT operator[]( indexT idx /**< [in] the index of the entry to access*/ ) const
    {
        snapshotT sn = snapshot();

        if( sn.validEntries <= 0 )
        {
            return storedT{};
        }

        return at( sn.earliest, idx );
    }

    /// Get the entry at a given index relative a fixed reference entry.
    storedT at( indexT refEntry, ///< [in] the entry to start counting from
                indexT idx       ///< [in] the index of the entry to access
    ) const
    {
        return m_buffer[physicalIndex( refEntry, idx )];
    }

    /// Return a snapshot of the current readable region.
    snapshotT snapshot() const
    {
        snapshotT sn;

        sn.mono = m_mono.load( std::memory_order_acquire );
        sn.latest = m_latest.load( std::memory_order_relaxed );
        sn.validEntries = m_validEntries.load( std::memory_order_relaxed );
        sn.maxEntries = m_maxEntries;
        sn.full = ( sn.validEntries == sn.maxEntries && sn.maxEntries > 0 );

        if( sn.full )
        {
            sn.earliest = m_nextEntry.load( std::memory_order_relaxed );
        }
        else
        {
            sn.earliest = 0;
        }

        return sn;
    }

    /// Determine whether a logical window is readable from the provided snapshot.
    bool windowReadable( const snapshotT &sn, ///< [in] the snapshot to validate against
                         indexT refEntry,     ///< [in] the entry to start counting from
                         indexT count         ///< [in] the number of entries requested
    ) const
    {
        if( count <= 0 || refEntry < 0 || count > sn.validEntries )
        {
            return false;
        }

        if( !sn.full )
        {
            return refEntry + count <= sn.validEntries;
        }

        for( indexT n = 0; n < count; ++n )
        {
            if( physicalIndex( refEntry, n ) == sn.earliest )
            {
                return false;
            }
        }

        return true;
    }

    /// Copy a logical sequence into caller-owned storage, retrying if the producer advances.
    bool loadSequence( indexT refEntry,   ///< [in] the entry to start counting from
                       indexT count,      ///< [in] the number of entries requested
                       storedT *dest,     ///< [out] destination storage of size count
                       snapshotT &sn,     ///< [out] the validated snapshot used for the copy
                       int maxRetries = 3 /**< [in] the maximum number of retries*/
    ) const
    {
        if( dest == nullptr )
        {
            return false;
        }

        for( int retry = 0; retry < maxRetries; ++retry )
        {
            sn = snapshot();

            if( !windowReadable( sn, refEntry, count ) )
            {
                return false;
            }

            for( indexT n = 0; n < count; ++n )
            {
                dest[n] = m_buffer[physicalIndex( refEntry, n )];
            }

            if( m_mono.load( std::memory_order_acquire ) == sn.mono )
            {
                return true;
            }
        }

        return false;
    }

    /// Copy a logical sequence from a caller-supplied snapshot, validating that the producer has not advanced.
    bool loadSequence( const snapshotT &sn, ///< [in] the snapshot to validate against
                       indexT refEntry,     ///< [in] the entry to start counting from
                       indexT count,        ///< [in] the number of entries requested
                       storedT *dest        ///< [out] destination storage of size count
    ) const
    {
        if( dest == nullptr || !windowReadable( sn, refEntry, count ) )
        {
            return false;
        }

        for( indexT n = 0; n < count; ++n )
        {
            dest[n] = m_buffer[physicalIndex( refEntry, n )];
        }

        return m_mono.load( std::memory_order_acquire ) == sn.mono;
    }

    /// Copy the latest logical sequence into caller-owned storage, retrying if the producer advances.
    bool loadLatestSequence( indexT count,      ///< [in] the number of entries requested
                             storedT *dest,     ///< [out] destination storage of size count
                             snapshotT &sn,     ///< [out] the validated snapshot used for the copy
                             int maxRetries = 3 /**< [in] the maximum number of retries*/
    ) const
    {
        for( int retry = 0; retry < maxRetries; ++retry )
        {
            sn = snapshot();

            if( count <= 0 || count > sn.validEntries )
            {
                return false;
            }

            indexT refEntry;
            if( sn.latest >= count - 1 )
            {
                refEntry = sn.latest - ( count - 1 );
            }
            else
            {
                refEntry = sn.maxEntries + sn.latest - ( count - 1 );
            }

            if( !windowReadable( sn, refEntry, count ) )
            {
                return false;
            }

            for( indexT n = 0; n < count; ++n )
            {
                dest[n] = m_buffer[physicalIndex( refEntry, n )];
            }

            if( m_mono.load( std::memory_order_acquire ) == sn.mono )
            {
                return true;
            }
        }

        return false;
    }

  protected:
    /// Convert a logical index into a physical slot index.
    indexT physicalIndex( indexT refEntry, ///< [in] the entry to start counting from
                          indexT idx       ///< [in] the index of the entry to access
    ) const
    {
        return static_cast<indexT>( m_indices[m_maxEntries + refEntry + idx] );
    }
};

} // namespace sigproc
} // namespace mx

#endif // sigproc_circularBuffer
