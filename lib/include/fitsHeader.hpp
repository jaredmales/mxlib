/** \file fitsHeader.hpp
  * \brief Declares and defines a class to work with a FITS header
  * \ingroup image_processing
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
   
/** \addtogroup image_processing
  * @{
  */

/// Class to manage a complete fits header
/** Manages tasks such as insertion (avoiding duplicates), and keyword lookup.
  * The \ref fitsHeaderCard "cards" are stored in a std::list to preserver order, but
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

   fitsHeader()
   {
   }
   ///Copy constructor
   /** Must be defined to handle creation of new iterators in the cardMap
    */
   fitsHeader(const fitsHeader & head);
   
   ~fitsHeader()
   {
      clear();
   }
   
   /// Get iterator to the beginning of the cards list
   headerIterator begin();
   
   /// Get iterator to the end of the cards list
   headerIterator end();
   
   /// Test whether the header is empty.
   bool empty();
   
   /// Get number of cards currently stored in the header.
   size_t size();
   
   /// Clear al cards from the header
   void clear();
   
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
   template<typename typeT> void append(const std::string &k, typeT v, const std::string &c);
   
   /// Append a card to the end of the header, from the components of a card with no comment.
   /**
     * \tparam typeT is the data type of the value
     * 
     * \param k is the keyword string
     * \param v is the value of typeT
     */
   template<typename typeT> void append(const std::string &k, typeT v);
      
   /// Insert a card before another card.
   void insert_before(headerIterator it, fitsHeaderCard card);
   template<typename typeT> void insert_before(headerIterator it, const std::string &k, typeT v, const std::string &c);
   
   template<typename typeT> void insert_before(headerIterator it, const std::string &k, typeT v);
   
   /// Insert a card after another card.
   void insert_after(headerIterator it, fitsHeaderCard card);
   
   template<typename typeT> void insert_after(headerIterator it, const std::string &k, typeT v, const std::string &c);
   
   template<typename typeT> void insert_after(headerIterator it, const std::string &k, typeT v);
   
   /// Card access by keyword operator
   /** Looks up the card by its keyword, and returns a reference to it.
    */
   fitsHeaderCard & operator[](const std::string & keyword);
   
   const fitsHeaderCard & operator[](const std::string & keyword) const;
   
};  // fitsHeader




template<typename typeT> void fitsHeader::append(const std::string &k, typeT v, const std::string &c)
{
   append(fitsHeaderCard(k,v,c));
}


template<typename typeT> void fitsHeader::append(const std::string &k, typeT v)
{
   append(fitsHeaderCard(k,v));
}



template<typename typeT> void fitsHeader::insert_before(headerIterator it, const std::string &k, typeT v, const std::string &c)
{
   insert_before(it, fitsHeaderCard(k,v,c));
}


template<typename typeT> void fitsHeader::insert_before(headerIterator it, const std::string &k, typeT v)
{
   insert_before(it, fitsHeaderCard(k,v));
}




template<typename typeT> void fitsHeader::insert_after(headerIterator it, const std::string &k, typeT v, const std::string &c)
{
   insert_after(it, fitsHeaderCard(k,v,c));
}


template<typename typeT> void fitsHeader::insert_after(headerIterator it, const std::string &k, typeT v)
{
   insert_after(it, fitsHeaderCard(k,v));
}

///Convert the values in a vector of \ref fitsHeader "fits headers" into a vector of values.
template<typename dataT>
std::vector<dataT> headersToValues(const std::vector<fitsHeader> & heads, const std::string &keyw)
{
   std::vector<dataT> v(heads.size());
   
   for(int i=0;i<heads.size(); ++i)
   {
      v[i] = convertFromString<dataT>(heads[i][keyw].value);
   }

   return v;
}

///@}

} //namespace mx



#endif //__fitsHeader_hpp__
