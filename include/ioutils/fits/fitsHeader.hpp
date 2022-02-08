/** \file fitsHeader.hpp
  * \brief Declares and defines a class to work with a FITS header
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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

#include "fitsHeaderCard.hpp"

namespace mx
{
namespace fits 
{
   

/// Class to manage a complete fits header
/** Manages tasks such as insertion (avoiding duplicates), and keyword lookup.
  * The \ref fitsHeaderCard "cards" are stored in a std::list to preserve order, but
  * a std::unordered_multimap is used to provide fast keyword lookup.
  * 
  * \ingroup fits_processing
  */ 
class fitsHeader
{

public:

   /// The iterator type for the cards list
   typedef std::list<fitsHeaderCard>::iterator headerIterator ;

   /// The iterator type for the card map
   typedef std::unordered_multimap<std::string, headerIterator>::iterator mapIterator;
   
protected:
   
   /// The storage for the FITS header cards
   /** We use a list,  rather than forward_list, so that append (insert at end) is constant time.
     * 
     */
   std::list<fitsHeaderCard> cards;
   
   /// This multimap allows for fast lookup by keyword.
   /** Use unordered_multimap to handle HISTORY and COMMENT properly, and be as efficient as possible.
     */
   std::unordered_multimap<std::string, headerIterator> cardMap;
   
public:

   ///Default c'tor
   fitsHeader();
   
   ///Copy constructor
   /** Must be explicitly defined to handle creation of new iterators in the cardMap
    */
   fitsHeader(const fitsHeader & head /**< The fitsHeader to copy */);
   
   ///Destructor
   ~fitsHeader();
   
   /// Assignment
   /** Must be explicitly defined to handle creation of new iterators in the cardMap
     */
   fitsHeader & operator=(const fitsHeader & head /**< The fitsHeader to copy */);
   
   /// Get iterator to the beginning of the cards list
   headerIterator begin();
   
   /// Get iterator to the end of the cards list
   headerIterator end();
   
   /// Get iterator pointing to a specific element
   headerIterator iterator(const std::string & keyword /**< The keyword to look up*/);
   
   /// Test whether the header is empty.
   bool empty();
   
   /// Get number of cards currently stored in the header.
   size_t size();
   
   /// Clear all cards from the header
   void clear();
   
   /// Get number of cards with a given keyword 
   /** Reeturns the result of the count() method of the header map.
     * 
     * \retval the number of cards with  keyword.
     */
   size_t count( const std::string & keyword /**< [in] the keyword to loop up*/);
   
   /// Erase card by keyword
   /** This can not be used to erase COMMENT or HISTORY cards.
     *
     * \param keyword the keyword of the card to delete 
     * 
     */
   void erase(const std::string & keyword);
   
   /// Erase card by iterator
   /** This handles COMMENT and HISTORY cards, deleting only the one pointed to by it
     *
     * \param it is a headerIterator pointing to the card to delete.
     * 
     */
   void erase(headerIterator it);
   
   /// Erase the standard entries at the top of the header
   /** Erases each entry down to BSCALE.  This is useful for appending
     * a header to a newly created file.
     */  
   void eraseStandardTop();
   
   /// Append a fitsHeaderCard to the end of the header
   /**
     * \param card is a fitsHeaderCard already populated
     */ 
   void append(fitsHeaderCard card);
   
   /// Append a card to the end of the header, from the three components of a card.
   /**
     * \tparam typeT is the data type of the value
     * 
     * \param k is the keyword string
     * \param v is the value of typeT
     * \param c is the comment string
     */ 
   template<typename typeT> 
   void append(const std::string &k, const typeT &v, const std::string &c);
   
   /// Append a card to the end of the header, from the components of a card with no comment.
   /**
     * \tparam typeT is the data type of the value
     * 
     * \param k is the keyword string
     * \param v is the value of typeT
     */
   template<typename typeT> 
   void append(const std::string &k, const typeT &v);
          
   /// Append a card to the end of the header, with just a keyword.
   /** Appens a headerCard with unknownType
     * 
     * \param k is the keyword string
     */
   void append(const std::string &k);
   
   /// Append a fitsHeader to the end of the header
   /**
     * \param head is a populated fitsHeader
     */ 
   void append(fitsHeader & head);
   
   /// Insert a card before another card.
   /** 
     * \param it points to the element before which to insert
     * \param card contains the card to insert
     */ 
   void insert_before(headerIterator it, fitsHeaderCard card);
   
   /// Insert a card before another card, specifying the card by its components.
   /** 
     * \tparam typeT is the type of the value, which is converted to string for insertion
     * 
     * \param it points to the element before which to insert
     * \param k is the keyword
     * \param v is the value
     * \param c is the comment
     * 
     */ 
   template<typename typeT> 
   void insert_before(headerIterator it, const std::string &k, typeT v, const std::string &c);
   
   /// Insert a card before another card, specifying the card by its components.
   /** 
     * \tparam typeT is the type of the value, which is converted to string for insertion
     *
     * \param it points to the element before which to insert
     * \param k is the keyword
     * \param v is the value
     * 
     */
   template<typename typeT> 
   void insert_before(headerIterator it, const std::string &k, typeT v);
   
   /// Insert a card after another card.
   /** 
     * \param it points to the element after which to insert
     * \param card contains the card to insert
     */ 
   void insert_after(headerIterator it, fitsHeaderCard card);
   
   /// Insert a card after another card, specifying the card by its components.
   /** 
     * \tparam typeT is the type of the value, which is converted to string for insertion
     * 
     * \param it points to the element after which to insert
     * \param k is the keyword
     * \param v is the value
     * \param c is the comment
     * 

     */ 
   template<typename typeT> 
   void insert_after(headerIterator it, const std::string &k, typeT v, const std::string &c);
   
   /// Insert a card after another card, specifying the card by its components.
   /** 
     * \tparam typeT is the type of the value, which is converted to string for insertion
     * 
     * \param it points to the element after which to insert
     * \param k is the keyword
     * \param v is the value
     * 
     */
   template<typename typeT> 
   void insert_after(headerIterator it, const std::string &k, typeT v);
   
   /// Card access by keyword operator
   /** Looks up the card by its keyword, and returns a reference to it.
     *
     * \param keyword is the header keyword to look up
     *
     * \retval fitsHeaderCard& reference to the \ref fitsHeaderCard 
     */
   fitsHeaderCard & operator[](const std::string & keyword);
   
   /// Card access by keyword operator (const version)
   /** Looks up the card by its keyword, and returns a reference to it.
     *
     * \param keyword is the header keyword to look up
     *
     * \retval fitsHeaderCard& const reference to the \ref fitsHeaderCard 
     */
   const fitsHeaderCard & operator[](const std::string & keyword) const;
   
   
   
};  // fitsHeader



template<typename typeT> 
void fitsHeader::append( const std::string &k, 
                         const typeT &v, 
                         const std::string &c
                       )
{
   append(fitsHeaderCard(k,v,c));
}


template<typename typeT> 
void fitsHeader::append( const std::string &k, 
                         const typeT &v
                       )
{
   append(fitsHeaderCard(k,v));
}

template<typename typeT> 
void fitsHeader::insert_before( headerIterator it, 
                                const std::string &k, 
                                typeT v, 
                                const std::string &c
                              )
{
   insert_before(it, fitsHeaderCard(k,v,c));
}


template<typename typeT> 
void fitsHeader::insert_before( headerIterator it, 
                                const std::string &k, 
                                typeT v
                              )
{
   insert_before(it, fitsHeaderCard(k,v));
}

template<typename typeT> 
void fitsHeader::insert_after( headerIterator it, 
                               const std::string &k, 
                               typeT v, 
                               const std::string &c
                             )
{
   insert_after(it, fitsHeaderCard(k,v,c));
}


template<typename typeT> 
void fitsHeader::insert_after( headerIterator it, 
                               const std::string &k, 
                               typeT v
                             )
{
   insert_after(it, fitsHeaderCard(k,v));
}


/** \addtogroup fits_utils
  * @{
  */

///Convert the values in a std::vector of \ref fitsHeader "fits headers" into a std::vector of values.
/** Resizes the vector of the appropriate type.
  *
  * \tparam dataT is the type of the header value
  * 
  * \param[out] v will contain the converted values
  * \param[in] heads contains the headers 
  * \param[in] keyw contains the keyword designating which value to convert
  * 
  */
template<typename dataT>
void headersToValues( std::vector<dataT> & v, 
                      std::vector<fitsHeader> & heads, 
                      const std::string &keyw
                    )
{
   v.resize(heads.size());
   
   for(size_t i=0;i<heads.size(); ++i)
   {
      v[i] = heads[i][keyw].value<dataT>();//convertFromString<dataT>(heads[i][keyw].value);
   }

}

///Convert the values in a std::vector of \ref fitsHeader "fits headers" into a std::vector of values.
/** Creates a vector of the appropriate type and size.
  *
  * \tparam dataT is the type of the header value
  * 
  * \param[in] heads contains the headers 
  * \param[in] keyw contains the keyword designating which value to convert
  * 
  * \retval std::vector<dataT> containing the converted values
  */
template<typename dataT>
std::vector<dataT> headersToValues( std::vector<fitsHeader> & heads, 
                                    const std::string &keyw
                                  )
{
   std::vector<dataT> v(heads.size());
   
   headersToValues(v,heads, keyw);
   
   return v;
}

   
///Write the status of a Git repository to HISTORY in a FITS header.
/**
  * \param [in,out] head the HISTORY cards will be appended to this header
  * \param [in] repoName the name of the repository
  * \param [in] sha1 is the SHA-1 hash string of the repository 
  * \param [in] modified whether or not the repository has been modified after the 
  *                      commit referred to by sha1
  */
void fitsHeaderGitStatus( fitsHeader & head, 
                          const std::string & repoName,
                          const char * sha1,
                          int modified
                        );

///@}


} //namespace fits   
} //namespace mx



#endif //ioutils_fits__fitsHeader_hpp
