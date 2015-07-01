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
   fitsHeader(const fitsHeader & head);
   
   ///Destructor
   ~fitsHeader();
   
   /// Assignment
   /** Must be explicitly defined to handle creation of new iterators in the cardMap
     */
   fitsHeader & operator=(const fitsHeader & head);
   
   /// Get iterator to the beginning of the cards list
   headerIterator begin();
   
   /// Get iterator to the end of the cards list
   headerIterator end();
   
   /// Get iterator pointing to a specific element
   headerIterator iterator(const std::string & keyword);
   
   /// Test whether the header is empty.
   bool empty();
   
   /// Get number of cards currently stored in the header.
   size_t size();
   
   /// Clear all cards from the header
   void clear();
   
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


inline 
fitsHeader::fitsHeader()
{
}
   
inline 
fitsHeader::fitsHeader(const fitsHeader & head)
{
   operator=(head);
   
}

inline 
fitsHeader::~fitsHeader()
{
   clear();
}

inline
fitsHeader & fitsHeader::operator=(const fitsHeader & head)
{
   cards = head.cards;
   
   headerIterator it = cards.begin();
   
   cardMap.clear();
   while(it != cards.end())
   {
      cardMap.insert(std::pair<std::string, headerIterator>(it->keyword, it));         
      ++it;
   }
   
   return *this;
}


inline 
fitsHeader::headerIterator fitsHeader::begin() 
{
   return cards.begin();
}

inline 
fitsHeader::headerIterator fitsHeader::end() 
{
   return cards.end();
}

inline 
fitsHeader::headerIterator fitsHeader::iterator(const std::string & keyword)
{
   return cardMap.find(keyword)->second;   
}


inline 
bool fitsHeader::empty()
{
   return cards.empty();
}
   
inline 
size_t fitsHeader::size()
{
   return cards.size();
}
   
inline 
void fitsHeader::clear()
{
   cards.clear();
   cardMap.clear();
}

inline
void fitsHeader::erase(const std::string & keyword)
{
   if(keyword == "COMMENT" || keyword == "HISTORY")
   {
      return;
   }
   headerIterator it = cardMap.find(keyword)->second;
   cardMap.erase(keyword);
   cards.erase(it);
}

inline
void fitsHeader::erase(headerIterator it)
{
   mapIterator mit = cardMap.find(it->keyword);
              
   if(it->keyword == "COMMENT" || it->keyword == "HISTORY")
   {
      while(mit->second->keyword == it->keyword && it->comment != mit->second->comment) ++mit;
   }
   cardMap.erase(mit);
   cards.erase(it);
}

inline 
void fitsHeader::eraseStandardTop()
{
   
   headerIterator it = begin();
   
   while(it->keyword != "BSCALE" && it != end())
   {
      erase(it);
      
      it = begin();
   }
   
   erase("BSCALE");   
}

            
inline
void fitsHeader::append(fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword) > 0)
   {
      if(card.type !=  fitsTCOMMENT && card.type != fitsTHISTORY)
      {
         std::cerr << "attempt to duplicate keyword: " << card.keyword << "\n";
         return;
      }
   }
   
   //Now insert in list
   cards.push_back(card);
   
   //Then add to the Map.
   headerIterator insertedIt = cards.end();
   --insertedIt;
   cardMap.insert( std::pair<std::string, headerIterator>(card.keyword, insertedIt) );
   
}

inline
void fitsHeader::append(fitsHeader & head)
{
   headerIterator it;
   
   for(it = head.begin(); it != head.end(); ++it)
   {
      append(*it);
   }
}


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

inline
void fitsHeader::insert_before(headerIterator it, fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword) > 0)
   {
      if(card.type !=  fitsTCOMMENT && card.type != fitsTHISTORY)
      {
         std::cerr << "attempt to duplicate keyword: " << card.keyword << "\n";
         return;
      }
   }
   
   //Now insert in list
   headerIterator insertedIt = cards.insert(it, card);
   
   //Then add to the Map.
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword, insertedIt));
      
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

inline
void fitsHeader::insert_after(headerIterator it, fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword) > 0)
   {
      if(card.type !=  fitsTCOMMENT && card.type != fitsTHISTORY)
      {
         std::cerr << "attempt to duplicate keyword:" << card.keyword << "\n";
         return;
      }
   }
   
   //Now insert in list
   headerIterator insertedIt = cards.insert(++it, card);
   
   //Then add to the Map.
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword, insertedIt));
   
   
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

inline
fitsHeaderCard & fitsHeader::operator[](const std::string & keyword)
{
   headerIterator it = cardMap.find(keyword)->second;
   
   return *it;
}

inline
const fitsHeaderCard & fitsHeader::operator[](const std::string & keyword) const
{
   headerIterator it = cardMap.find(keyword)->second;
   
   return *it;
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

   
inline void fitsHeaderGitStatus(fitsHeader & head, 
                                const std::string & repoName,
                                const char * sha1,
                                int modified)
{
   std::string hist = "Git status for " + repoName + ":";
   head.append("", fitsHistoryType(), hist);
   
   hist = "   sha1=";
   hist += sha1;
   if(modified) hist += ", modified.";
   else hist += ".";
   head.append("", fitsHistoryType(), hist);
}

///@}

} //namespace mx



#endif //__fitsHeader_hpp__
