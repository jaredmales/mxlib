/** \file fitsHeader.cpp
  * \brief Definitions for a class to work with a FITS header
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  

#include "fitsHeader.hpp"

namespace mx
{

fitsHeader::fitsHeader(const fitsHeader & head)
{
   cards = head.cards;
   
   headerIterator it = cards.begin();
   
   cardMap.clear();
   while(it != cards.end())
   {
      cardMap.insert(std::pair<std::string, headerIterator>(it->keyword, it));         
      ++it;
   }
}

fitsHeader::headerIterator fitsHeader::begin() 
{
   return cards.begin();
}

fitsHeader::headerIterator fitsHeader::end() 
{
   return cards.end();
}

bool fitsHeader::empty()
{
   return cards.empty();
}
   
size_t fitsHeader::size()
{
   return cards.size();
}
   
void fitsHeader::clear()
{
   cards.clear();
   cardMap.clear();
}

   
void fitsHeader::append(fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword) > 0)
   {
      if(card.keyword != "HISTORY" && card.keyword != "COMMENT")
      {
         std::cerr << "attempt to duplicate keyword\n";
         return;
      }
   }
   
   //Now insert in list
   cards.push_back(card);
   
   //Then add to the Map.
   headerIterator insertedIt = cards.end();
   insertedIt--;
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword, insertedIt));
   
}



void fitsHeader::insert_before(headerIterator it, fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword) > 0)
   {
      if(card.keyword != "HISTORY" && card.keyword != "COMMENT")
      {
         std::cerr << "attempt to duplicate keyword\n";
         return;
      }
   }
   
   //Now insert in list
   headerIterator insertedIt = cards.insert(it, card);
   
   //Then add to the Map.
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword, insertedIt));
   
   
}


void fitsHeader::insert_after(headerIterator it, fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword) > 0)
   {
      if(card.keyword != "HISTORY" && card.keyword != "COMMENT")
      {
         std::cerr << "attempt to duplicate keyword\n";
         return;
      }
   }
   
   //Now insert in list
   headerIterator insertedIt = cards.insert(++it, card);
   
   //Then add to the Map.
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword, insertedIt));
   
   
}

fitsHeaderCard & fitsHeader::operator[](const std::string & keyword)
{
   headerIterator it = cardMap.find(keyword)->second;
   
   return *it;
}

const fitsHeaderCard & fitsHeader::operator[](const std::string & keyword) const
{
   headerIterator it = cardMap.find(keyword)->second;
   
   return *it;
}

} //namespace mx


