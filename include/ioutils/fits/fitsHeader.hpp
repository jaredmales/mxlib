/** \file fitsHeader.hpp
 * \brief Declares and defines a class to work with a FITS header
 * \ingroup fits_processing_files
 *
 */

//***********************************************************************//
// Copyright 2015-2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef ioutils_fits__fitsHeader_hpp
#define ioutils_fits__fitsHeader_hpp

#include <list>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <optional>

#include "fitsHeaderCard.hpp"

namespace mx
{
namespace fits
{

namespace fitsHeaderDetail
{

/** \cond */
/// Mutation stages exposed to deterministic failure tests.
enum class mutation
{
    appendList,
    appendMap,
    insertBeforeList,
    insertBeforeMap,
    insertAfterList,
    insertAfterMap
};

/// Signature of the resettable mutation hook used by fitsHeader.
using mutationHookT = error_t ( * )( mutation );

/// Access the process-wide fitsHeader mutation hook.
mutationHookT &mutationHook();

/// Restore the fitsHeader mutation hook to its production no-op.
void resetMutationHook();
/** \endcond */

} // namespace fitsHeaderDetail

/// Class to manage a FITS file metadata header and provide fast access to the cards by keyword
/** Manages tasks such as insertion (avoiding duplicates), and keyword lookup.
 *
 * The \ref fitsHeaderCard "cards" are stored in a `std::list` to preserve order, but
 * a `std::unordered_multimap` is used to provide fast keyword lookup into the list.
 *
 * \tparam verboseT sets the error reporting \ref mx::verbose "verbosity"
 *
 * \ingroup fits_processing
 */
template <class verboseT = verbose::d>
class fitsHeader
{
  public:
    /// The list type
    /** We use a list,  rather than forward_list, so that append (insert at end) is constant time.
     *
     */
    typedef std::list<fitsHeaderCard<verboseT>> cardListT;

    /// The iterator type for the cards list
    typedef cardListT::iterator headerIteratorT;

    /// The map type
    /** Use unordered_multimap to allow multiple HISTORY and COMMENT properly, but be as efficient as possible.
     */
    typedef std::unordered_multimap<std::string, headerIteratorT> cardMapT;

    /// The value type for the card map
    typedef cardMapT::value_type cardMapValueT;

    /// The iterator type for the card map
    typedef cardMapT::iterator mapIteratorT;

  protected:
    /// The storage for the FITS header cards
    cardListT m_cardList;

    /// This multimap allows for fast lookup by keyword.
    cardMapT m_cardMap;

    /// Card with empty keyword used as an error sentinel.
    mutable fitsHeaderCard<verboseT> m_emptyCard;

  public:
    /// Default c'tor
    fitsHeader();

    /// Copy constructor
    fitsHeader( const fitsHeader &head /**< The fitsHeader to copy */ );

    /// Destructor
    ~fitsHeader();

    /// Assignment operator
    fitsHeader &operator=( const fitsHeader &head /**< The fitsHeader to copy */ );

    /// Get iterator to the beginning of the cards list
    headerIteratorT begin();

    /// Get iterator to the end of the cards list
    headerIteratorT end();

    /// Get iterator pointing to a specific element
    headerIteratorT iterator( const std::string &keyword /**< The keyword to look up*/ );

    /// Test whether the header is empty.
    bool empty();

    /// Get number of cards currently stored in the header.
    size_t size();

    /// Clear all cards from the header
    error_t clear();

    /// Get number of cards with a given keyword
    /** Returns the result of the count() method of the header map.
     *
     * \retval the number of cards with  keyword.
     */
    size_t count( const std::string &keyword /**< [in] the keyword to look up*/ );

    /// Erase card by keyword
    /** This can not be used to erase COMMENT or HISTORY cards.
     *
     */
    error_t erase( const std::string &keyword /**< [in] the keyword of the card to delete*/ );

    /// Erase card by iterator
    /** This handles COMMENT and HISTORY cards, deleting only the one pointed to by it
     *  using all the contents of the card (not just the keyword)
     */
    error_t erase( headerIteratorT it /**< [in] iterator pointing to the card to delete. */ );

    /// Erase the standard entries at the top of the header
    /** Erases each entry down to BSCALE.  This is useful for appending
     * a header from a previous image to a newly created file. Also erases boilerplate comments,
     * such as for long string.
     */
    error_t eraseStandardTop();

    /// Append a fitsHeaderCard to the end of the header
    /**
     */
    error_t append( const fitsHeaderCard<verboseT> &card /**< [in] the card to append*/ );

    /// Append a string card to the end of the header, from the three components of a card.
    /**
     */
    error_t append( const std::string &k, /**< [in] the keyword of the new card*/
                    const char *v,        /**< [in] the value of the new card*/
                    const std::string &c  /**< [in] the comment of the new card*/
    );

    /// Append a card to the end of the header, from the three components of a card.
    /**
     * \tparam typeT is the data type of the value
     *
     */
    template <typename typeT>
    error_t append( const std::string &k, /**< [in] the keyword of the new card*/
                    const typeT &v,       /**< [in] the value of the new card*/
                    const std::string &c  /**< [in] the comment of the new card*/
    );

    /// Append a card to the end of the header, from the components of a card with no comment.
    /**
     * \tparam typeT is the data type of the value
     *
     */
    template <typename typeT>
    error_t append( const std::string &k, /**< [in] the keyword of the new card*/
                    const typeT &v        /**< [in] the value of the new card*/
    );

    /// Append a card to the end of the header, with just a keyword.
    /** Appends a headerCard with unknownType
     *
     */
    error_t append( const std::string &k /**< [in] the keyword of the new card*/ );

    /// Append a fitsHeader to the end of the header
    /**
     */
    error_t append( fitsHeader &head /**< [in]  the fitsHeader to append*/ );

    /// Insert a card before another card.
    /**
     */
    error_t insert_before( headerIteratorT it,           /**< [in] iterator pointing to the
                                                                   element before which to insert*/
                           fitsHeaderCard<verboseT> card /**< [in] the card to insert*/
    );

    /// Insert a card before another card, specifying the card by its components.
    /**
     * \tparam typeT is the type of the value, which is converted to string for insertion
     *
     */
    template <typename typeT>
    error_t insert_before( headerIteratorT it,   /**< [in] iterator pointing to the element before which to insert*/
                           const std::string &k, /**< [in] the keyword of the new card*/
                           typeT v,              /**< [in] the value of the new card*/
                           const std::string &c  /**< [in] the comment of the new card*/
    );

    /// Insert a card before another card, specifying the card by its components.
    /**
     * \tparam typeT is the type of the value, which is converted to string for insertion
     *
     */
    template <typename typeT>
    error_t insert_before( headerIteratorT it,   /**< [in] iterator pointing to the element before which to insert*/
                           const std::string &k, /**< [in] the keyword of the new card*/
                           typeT v               /**< [in] the value of the new card*/
    );

    /// Insert a card after another card.
    /**
     */
    error_t insert_after( headerIteratorT it,           /**< [in] iterator pointing to the element
                                                                  after which to insert*/
                          fitsHeaderCard<verboseT> card /**< [in] the card to insert*/
    );

    /// Insert a card after another card, specifying the card by its components.
    /**
     * \tparam typeT is the type of the value, which is converted to string for insertion
     *
     */
    template <typename typeT>
    error_t insert_after( headerIteratorT it,   /**< [in] iterator pointing to the element after which to insert*/
                          const std::string &k, /**< [in] the keyword of the new card*/
                          typeT v,              /**< [in] the value of the new card*/
                          const std::string &c  /**< [in] the comment of the new card*/
    );

    /// Insert a card after another card, specifying the card by its components.
    /**
     * \tparam typeT is the type of the value, which is converted to string for insertion
     *
     */
    template <typename typeT>
    error_t insert_after( headerIteratorT it,   /**< [in] interator pointing to the element after which to insert*/
                          const std::string &k, /**< [in] the keyword of the new card*/
                          typeT v               /**< [in] the value of the new card*/
    );

    /// Card access by keyword operator
    /** Looks up the card by its keyword, and returns a reference to it.
     *
     * \returns on success: fitsHeaderCard& reference to the \ref fitsHeaderCard
     * \returns on error: a fitsHeaderCard& to a card with an empty keyword
     */
    fitsHeaderCard<verboseT> &operator[]( const std::string &keyword /**< [in] the header keyword to look up*/ );

    /// Card access by keyword operator (const version)
    /** Looks up the card by its keyword, and returns a reference to it.
     *
     * \returns on success: fitsHeaderCard& reference to the \ref fitsHeaderCard
     * \returns on error: a fitsHeaderCard& to a card with an empty keyword
     */
    const fitsHeaderCard<verboseT> &
    operator[]( const std::string &keyword /**< [in] the header keyword to look up*/ ) const;

}; // fitsHeader

template <class verboseT>
fitsHeader<verboseT>::fitsHeader()
{
}

template <class verboseT>
fitsHeader<verboseT>::fitsHeader( const fitsHeader &head )
{
    operator=( head );
}

template <class verboseT>
fitsHeader<verboseT>::~fitsHeader()
{
    clear();
}

template <class verboseT>
fitsHeader<verboseT> &fitsHeader<verboseT>::operator=( const fitsHeader &head )
{
    m_cardList = head.m_cardList;

    headerIteratorT it = m_cardList.begin();

    m_cardMap.clear();
    while( it != m_cardList.end() )
    {
        m_cardMap.insert( cardMapValueT( it->keyword(), it ) );
        ++it;
    }

    return *this;
}

template <class verboseT>
fitsHeader<verboseT>::headerIteratorT fitsHeader<verboseT>::begin()
{
    return m_cardList.begin();
}

template <class verboseT>
fitsHeader<verboseT>::headerIteratorT fitsHeader<verboseT>::end()
{
    return m_cardList.end();
}

template <class verboseT>
fitsHeader<verboseT>::headerIteratorT fitsHeader<verboseT>::iterator( const std::string &keyword )
{
    auto mit = m_cardMap.find( keyword );
    if( mit == m_cardMap.end() )
    {
        return m_cardList.end();
    }

    return mit->second;
}

template <class verboseT>
bool fitsHeader<verboseT>::empty()
{
    return m_cardList.empty();
}

template <class verboseT>
size_t fitsHeader<verboseT>::size()
{
    return m_cardList.size();
}

template <class verboseT>
error_t fitsHeader<verboseT>::clear()
{
    m_cardList.clear();
    m_cardMap.clear();

    return error_t::noerror;
}

template <class verboseT>
size_t fitsHeader<verboseT>::count( const std::string &keyword )
{
    return m_cardMap.count( keyword );
}

template <class verboseT>
error_t fitsHeader<verboseT>::erase( const std::string &keyword )
{
    if( keyword == "COMMENT" )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "can't erase COMMENT by keyword" );
    }

    if( keyword == "HISTORY" )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "can't erase HISTORY by keyword" );
    }

    auto mit = m_cardMap.find( keyword );
    if( mit == m_cardMap.end() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::notfound, "keyword not found in map:" + keyword );
    }

    headerIteratorT it = mit->second;

    m_cardMap.erase( mit );
    m_cardList.erase( it );

    return error_t::noerror;
}

template <class verboseT>
error_t fitsHeader<verboseT>::erase( headerIteratorT it )
{
    if( it == m_cardList.end() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "invalid list iterator" );
    }

    std::string keyword = it->keyword();

    auto range = m_cardMap.equal_range( keyword );
    mapIteratorT mit;
    for( mit = range.first; mit != range.second; ++mit )
    {
        if( mit->second == it )
        {
            break;
        }
    }

    m_cardMap.erase( mit ); // does not throw
    m_cardList.erase( it );

    return error_t::noerror;
}

template <class verboseT>
error_t fitsHeader<verboseT>::eraseStandardTop()
{

    headerIteratorT it = begin(), nit;

    int n = 0;
    while( it != end() )
    {
        nit = it;
        ++nit;
        if( it->keyword() == "SIMPLE" || it->keyword() == "BITPIX" || it->keyword() == "NAXIS" ||
            it->keyword() == "NAXIS1" || it->keyword() == "NAXIS2" || it->keyword() == "NAXIS3" ||
            it->keyword() == "EXTEND" || it->keyword() == "BZERO" || it->keyword() == "BSCALE" ||
            it->keyword() == "LONGSTRN" )
        {
            mxlib_error_check( erase( it ) );
        }

        if( it->keyword() == "COMMENT" )
        {
            if( it->comment().find( "FITS (Flexible Image" ) != std::string::npos )
            {
                mxlib_error_check( erase( it ) );
            }
            else if( it->comment().find( "and Astrophysics'" ) != std::string::npos )
            {
                mxlib_error_check( erase( it ) );
            }
        }

        if( nit == end() )
        {
            break;
        }
        it = nit;
        ++n;
    }

    return error_t::noerror;
}

template <class verboseT>
error_t fitsHeader<verboseT>::append( const fitsHeaderCard<verboseT> &card )
{
    if( card.keyword() == "CONTINUE" )
    {
        if( m_cardList.empty() )
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "CONTINUE card requires a preceding card" );
        }

        headerIteratorT backIt = m_cardList.end();
        --backIt;
        return backIt->appendContinue( card );
    }

    // First check if duplicate key
    if( m_cardMap.count( card.keyword() ) > 0 )
    {
        if( card.type() != fitsType<fitsCommentType>() && card.type() != fitsType<fitsHistoryType>() )
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "attempt to duplicate keyword " + card.keyword() );
        }
    }

    headerIteratorT insertedIt;

    // Now insert in list
    try
    {
        error_t errc = fitsHeaderDetail::mutationHook()( fitsHeaderDetail::mutation::appendList );
        if( errc != error_t::noerror )
        {
            return internal::mxlib_error_report<verboseT>( errc, "inserting " + card.keyword() );
        }

        m_cardList.push_back( card );
        insertedIt = m_cardList.end();
        --insertedIt;
    }
    catch( const std::bad_alloc &e )
    {
        internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return error_t::std_bad_alloc;
        #else
            throw;
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return error_t::std_exception;
        #else
            throw;
        #endif
        // clang-format on
    }

    // Then add to the Map.
    try
    {
        error_t errc = fitsHeaderDetail::mutationHook()( fitsHeaderDetail::mutation::appendMap );
        if( errc != error_t::noerror )
        {
            m_cardList.erase( insertedIt );
            return internal::mxlib_error_report<verboseT>( errc, "inserting " + card.keyword() );
        }

        m_cardMap.insert( cardMapValueT( card.keyword(), insertedIt ) );
    }
    catch( const std::bad_alloc &e )
    {
        m_cardList.erase( insertedIt );
        internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return error_t::std_bad_alloc;
        #else
            throw;
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        m_cardList.erase( insertedIt );
        internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return error_t::std_exception;
        #else
            throw;
        #endif
        // clang-format on
    }

    return error_t::noerror;
}

template <class verboseT>
error_t fitsHeader<verboseT>::append( fitsHeader &head )
{
    if( this == &head )
    {
        fitsHeader copy( head );
        return append( copy );
    }

    headerIteratorT it;

    for( it = head.begin(); it != head.end(); ++it )
    {
        mxlib_error_check( append( *it ) );
    }

    return error_t::noerror;
}

template <class verboseT>
error_t fitsHeader<verboseT>::append( const std::string &k )
{
    mxlib_error_return( append( fitsHeaderCard<verboseT>( k ) ) );
}

template <class verboseT>
error_t fitsHeader<verboseT>::append( const std::string &k, const char *v, const std::string &c )
{
    mxlib_error_return( append( fitsHeaderCard<verboseT>( k, v, c ) ) );
}

template <class verboseT>
template <typename typeT>
error_t fitsHeader<verboseT>::append( const std::string &k, const typeT &v, const std::string &c )
{
    mxlib_error_return( append( fitsHeaderCard<verboseT>( k, v, c ) ) );
}

template <class verboseT>
template <typename typeT>
error_t fitsHeader<verboseT>::append( const std::string &k, const typeT &v )
{
    mxlib_error_return( append( fitsHeaderCard<verboseT>( k, v, "" ) ) );
}

template <class verboseT>
error_t fitsHeader<verboseT>::insert_before( headerIteratorT it, fitsHeaderCard<verboseT> card )
{
    // First check if duplicate key
    if( m_cardMap.count( card.keyword() ) > 0 )
    {
        if( card.type() != fitsType<fitsCommentType>() && card.type() != fitsType<fitsHistoryType>() )
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "duplicate keyword: " + card.keyword() );
        }
    }

    // Now insert in list
    headerIteratorT insertedIt;

    try
    {
        error_t errc = fitsHeaderDetail::mutationHook()( fitsHeaderDetail::mutation::insertBeforeList );
        if( errc != error_t::noerror )
        {
            return internal::mxlib_error_report<verboseT>( errc, "inserting " + card.keyword() );
        }

        insertedIt = m_cardList.insert( it, card );
    }
    catch( const std::bad_alloc &e )
    {
        internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return error_t::std_bad_alloc;
        #else
            throw;
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return error_t::std_exception;
        #else
            throw;
        #endif
        // clang-format on
    }

    // Then add to the Map.
    try
    {
        error_t errc = fitsHeaderDetail::mutationHook()( fitsHeaderDetail::mutation::insertBeforeMap );
        if( errc != error_t::noerror )
        {
            m_cardList.erase( insertedIt );
            return internal::mxlib_error_report<verboseT>( errc, "inserting " + card.keyword() );
        }

        m_cardMap.insert( cardMapValueT( card.keyword(), insertedIt ) );
    }
    catch( const std::bad_alloc &e )
    {
        m_cardList.erase( insertedIt );
        internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return error_t::std_bad_alloc;
        #else
            throw;
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        m_cardList.erase( insertedIt );
        internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return error_t::std_exception;
        #else
            throw;
        #endif
        // clang-format on
    }

    return error_t::noerror;
}

template <class verboseT>
template <typename typeT>
error_t fitsHeader<verboseT>::insert_before( headerIteratorT it, const std::string &k, typeT v, const std::string &c )
{
    mxlib_error_return( insert_before( it, fitsHeaderCard<verboseT>( k, v, c ) ) );
}

template <class verboseT>
template <typename typeT>
error_t fitsHeader<verboseT>::insert_before( headerIteratorT it, const std::string &k, typeT v )
{
    mxlib_error_return( insert_before( it, fitsHeaderCard<verboseT>( k, v, "" ) ) );
}

template <class verboseT>
template <typename typeT>
error_t fitsHeader<verboseT>::insert_after( headerIteratorT it, const std::string &k, typeT v, const std::string &c )
{
    mxlib_error_return( insert_after( it, fitsHeaderCard<verboseT>( k, v, c ) ) );
}

template <class verboseT>
error_t fitsHeader<verboseT>::insert_after( headerIteratorT it, fitsHeaderCard<verboseT> card )
{
    if( it == m_cardList.end() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "invalid list iterator" );
    }

    // First check if duplicate key
    if( m_cardMap.count( card.keyword() ) > 0 )
    {
        if( card.type() != fitsType<fitsCommentType>() && card.type() != fitsType<fitsHistoryType>() )
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "duplicate keyword: " + card.keyword() );
        }
    }

    // Now insert in list
    headerIteratorT insertedIt;

    try
    {
        error_t errc = fitsHeaderDetail::mutationHook()( fitsHeaderDetail::mutation::insertAfterList );
        if( errc != error_t::noerror )
        {
            return internal::mxlib_error_report<verboseT>( errc, "inserting " + card.keyword() );
        }

        insertedIt = m_cardList.insert( ++it, card );
    }
    catch( const std::bad_alloc &e )
    {
        internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return error_t::std_bad_alloc;
        #else
            throw;
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return error_t::std_exception;
        #else
            throw;
        #endif
        // clang-format on
    }

    // Then add to the Map.
    try
    {
        error_t errc = fitsHeaderDetail::mutationHook()( fitsHeaderDetail::mutation::insertAfterMap );
        if( errc != error_t::noerror )
        {
            m_cardList.erase( insertedIt );
            return internal::mxlib_error_report<verboseT>( errc, "inserting " + card.keyword() );
        }

        m_cardMap.insert( cardMapValueT( card.keyword(), insertedIt ) );
    }
    catch( const std::bad_alloc &e )
    {
        m_cardList.erase( insertedIt );
        internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS )
            return error_t::std_bad_alloc;
        #else
            throw;
        #endif
        // clang-format on
    }
    catch( const std::exception &e )
    {
        m_cardList.erase( insertedIt );
        internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                "inserting " + card.keyword() + " :" + e.what() );
        // clang-format off
        #if defined( MXLIB_CATCH_ALL_EXCEPTIONS ) || defined( MXLIB_CATCH_NONALLOC_EXCEPTIONS )
            return error_t::std_exception;
        #else
            throw;
        #endif
        // clang-format on
    }

    return error_t::noerror;
}

template <class verboseT>
template <typename typeT>
error_t fitsHeader<verboseT>::insert_after( headerIteratorT it, const std::string &k, typeT v )
{
    mxlib_error_return( insert_after( it, fitsHeaderCard<verboseT>( k, v, "" ) ) );
}

template <class verboseT>
fitsHeaderCard<verboseT> &fitsHeader<verboseT>::operator[]( const std::string &keyword )
{
    headerIteratorT it;

    auto mit = m_cardMap.find( keyword );

    // If not found, append it.
    if( mit == m_cardMap.end() )
    {
        error_t errc = append( keyword );
        if( errc != error_t::noerror )
        {
            internal::mxlib_error_report<verboseT>( errc, "appending " + keyword );
            m_emptyCard.keyword( "" ); // we have to reset this every time b/c someone could have changed it
            return m_emptyCard;
        }

        // have to do new search
        mit = m_cardMap.find( keyword );
    }

    it = mit->second;

    return *it;
}

template <class verboseT>
const fitsHeaderCard<verboseT> &fitsHeader<verboseT>::operator[]( const std::string &keyword ) const
{
    auto mit = m_cardMap.find( keyword );
    if( mit == m_cardMap.end() )
    {
        internal::mxlib_error_report<verboseT>( error_t::notfound, keyword + " not found" );
        m_emptyCard.keyword( "" ); // we have to reset this every time b/c someone could have changed it
        return m_emptyCard;
    }

    headerIteratorT it = mit->second;

    return *it;
}

/** \addtogroup fits_utils
 * @{
 */

/// Convert the values in a std::vector of \ref fitsHeader "fits headers" into a std::vector of values.
/** Resizes the vector of the appropriate type.
 *
 * \returns error_t::noerror on success. \p bad will be empty.
 * \returns error_t::error if the keyword fails to convert for any header.  \p bad will contain the indices of
 *          \p heads for which the extraction of a value for \p keyw failed
 *
 * \tparam dataT is the type of the header value
 * \tparam fitsHeaderT is the fitsHeader type
 *
 */
template <typename dataT, class fitsHeaderT>
error_t headersToValues( std::vector<dataT> &v,           /**< [out] will contain the converted values*/
                         std::vector<size_t> &bad,        /**<[out] will contain the indices of any
                                                                    headers that failed conversion */
                         std::vector<fitsHeaderT> &heads, /**< [in] contains the headers */
                         const std::string &keyw          /**< [in] contains the keyword designating which
                                                                    value to convert*/
)
{
    v.resize( heads.size() );
    bad.clear();

    for( size_t i = 0; i < heads.size(); ++i )
    {
        error_t errc;
        v[i] = heads[i][keyw].template value<dataT>( &errc ); // convertFromString<dataT>(heads[i][keyw].value);

        if( errc != error_t::noerror )
        {
            bad.push_back( i );
            v[i] = std::numeric_limits<dataT>::max();
        }
    }

    if( bad.size() > 0 )
    {
        return error_t::error;
    }
    else
    {
        return error_t::noerror;
    }
}

/// Write the status of a Git repository to HISTORY in a FITS header.
/**
 * \param [in,out] head the HISTORY cards will be appended to this header
 * \param [in] repoName the name of the repository
 * \param [in] sha1 is the SHA-1 hash string of the repository
 * \param [in] modified whether or not the repository has been modified after the
 *                      commit referred to by sha1
 */
template <class fitsHeaderT>
void fitsHeaderGitStatus( fitsHeaderT &head, const std::string &repoName, const char *sha1, int modified )
{
    std::string hist = "Git status for " + repoName + ":";
    head.append( "", fitsHistoryType(), hist );

    hist = "   sha1=";
    hist += sha1;
    if( modified )
        hist += ", modified";

    head.append( "", fitsHistoryType(), hist );
}

///@}

extern template class fitsHeader<verbose::d>;

} // namespace fits
} // namespace mx

#endif // ioutils_fits__fitsHeader_hpp
