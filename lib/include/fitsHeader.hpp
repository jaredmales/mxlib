/** \file fitsHeader.hpp
  * \brief Declares and defines a class to work with a FITS header
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __fitsHeader_hpp__
#define __fitsHeader_hpp__


#include <list>
#include <unordered_map>
#include <iostream>
#include <vector>

#include "fitsHeaderCard.hpp"

namespace mx
{
   
/** \addtogroup fits_processing
  * @{
  */

/// Class to manage a complete fits header
/** Manages tasks such as insertion (avoiding duplicates), and keyword lookup.
  * The \ref fitsHeaderCard "cards" are stored in a std::list to preserve order, but
  * a std::unordered_multimap is used to provide fast keyword lookup.
  * 
  */ 
class fitsHeader
{

public:

   /// The iterator type for the cards list
   typedef std::list<fitsHeaderCard>::iterator headerIterator ;

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
   fitsHeader()
   {
   }
   
   ///Copy constructor
   /** Must be defined to handle creation of new iterators in the cardMap
    */
   fitsHeader(const fitsHeader & head);
   
   ///Destructor
   ~fitsHeader()
   {
      clear();
   }
   
   /// Get iterator to the beginning of the cards list
   headerIterator begin();
   
   /// Get iterator to the end of the cards list
   headerIterator end();
   
   /// Get iterator pointing to a specific element
   headerIterator iterator(const std::string & keyword)
   {
      return cardMap.find(keyword)->second;
      
   }
   
   /// Test whether the header is empty.
   bool empty();
   
   /// Get number of cards currently stored in the header.
   size_t size();
   
   /// Clear all cards from the header
   void clear();
   
   void erase(const std::string & keyword)
   {
      headerIterator it = cardMap.find(keyword)->second;
      cardMap.erase(keyword);
      cards.erase(it);
   }
      
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

//@}


template<typename typeT> 
void fitsHeader::append(const std::string &k, const typeT &v, const std::string &c)
{
   append(fitsHeaderCard(k,v,c));
}


template<typename typeT> 
void fitsHeader::append(const std::string &k, const typeT &v)
{
   append(fitsHeaderCard(k,v));
}



template<typename typeT> 
void fitsHeader::insert_before(headerIterator it, const std::string &k, typeT v, const std::string &c)
{
   insert_before(it, fitsHeaderCard(k,v,c));
}


template<typename typeT> 
void fitsHeader::insert_before(headerIterator it, const std::string &k, typeT v)
{
   insert_before(it, fitsHeaderCard(k,v));
}




template<typename typeT> 
void fitsHeader::insert_after(headerIterator it, const std::string &k, typeT v, const std::string &c)
{
   insert_after(it, fitsHeaderCard(k,v,c));
}


template<typename typeT> void fitsHeader::insert_after(headerIterator it, const std::string &k, typeT v)
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
void headersToValues(std::vector<dataT> & v, const std::vector<fitsHeader> & heads, const std::string &keyw)
{
   v.resize(heads.size());
   
   for(int i=0;i<heads.size(); ++i)
   {
      v[i] = convertFromString<dataT>(heads[i][keyw].value);
   }

   return v;
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
std::vector<dataT> headersToValues(const std::vector<fitsHeader> & heads, const std::string &keyw)
{
   std::vector<dataT> v(heads.size());
   
   headersToValues(v,heads, keyw);
   
   return v;
}

///@}

} //namespace mx



#endif //__fitsHeader_hpp__
