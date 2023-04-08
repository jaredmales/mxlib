/** \file fitsHeader.cpp
  * \brief Implementation of a class to work with a FITS header
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

#include "ioutils/fits/fitsHeader.hpp"

namespace mx
{
namespace fits 
{
   
fitsHeader::fitsHeader()
{
}
   
fitsHeader::fitsHeader(const fitsHeader & head)
{
   operator=(head);
   
}

fitsHeader::~fitsHeader()
{
   clear();
}

fitsHeader & fitsHeader::operator=(const fitsHeader & head)
{
   cards = head.cards;
   
   headerIterator it = cards.begin();
   
   cardMap.clear();
   while(it != cards.end())
   {
      cardMap.insert(std::pair<std::string, headerIterator>(it->keyword(), it));         
      ++it;
   }
   
   return *this;
}


fitsHeader::headerIterator fitsHeader::begin() 
{
   return cards.begin();
}

fitsHeader::headerIterator fitsHeader::end() 
{
   return cards.end();
}

fitsHeader::headerIterator fitsHeader::iterator(const std::string & keyword)
{
   return cardMap.find(keyword)->second;   
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

size_t fitsHeader::count( const std::string & keyword)
{
   return cardMap.count(keyword);
}

void fitsHeader::erase(const std::string & keyword)
{
   if(keyword == "COMMENT" || keyword == "HISTORY")
   {
      return;
   }
   headerIterator it;
   it = (cardMap.find(keyword))->second; //trying to fool codacy
   cardMap.erase(keyword);
   cards.erase(it);
}

void fitsHeader::erase(headerIterator it)
{
   mapIterator mit = cardMap.find(it->keyword());
              
   if(it->keyword() == "COMMENT" || it->keyword() == "HISTORY")
   {
      while(mit->second->keyword() == it->keyword() && it->comment() != mit->second->comment()) ++mit;
   }
   cardMap.erase(mit);
   cards.erase(it);
}

void fitsHeader::eraseStandardTop()
{
  
   headerIterator it = begin(), nit;
  
   int n =0;
   while(it != end())
   {
      nit = it;
      ++nit;
      if(it->keyword() == "SIMPLE" || it->keyword() == "BITPIX" || it->keyword() == "NAXIS" 
           || it->keyword() == "NAXIS1" || it->keyword() == "NAXIS2" || it->keyword() == "NAXIS3" || it->keyword() == "EXTEND"
              || it->keyword() == "BZERO" || it->keyword() == "BSCALE" || it->keyword() == "LONGSTRN")
      {
         erase(it);
      }
         
      if(it->keyword() == "COMMENT")
      {
         if(it->comment().find("FITS (Flexible Image") != std::string::npos) erase(it);
         else if(it->comment().find("and Astrophysics'") != std::string::npos) erase(it);
      }
      
      if(nit == end()) break;
      it = nit;
      ++n;
   }
   
}

void fitsHeader::append(fitsHeaderCard card)
{
   //First check if duplicate key
   if(cardMap.count(card.keyword()) > 0)
   {
      if(card.type() !=  fitsType<fitsCommentType>() && card.type() != fitsType<fitsHistoryType>())
      {
         std::cerr << "attempt to duplicate keyword: " << card.keyword() << "\n";
         return;
      }
   }
   
   //Now insert in list
   cards.push_back(card);
   
   //Then add to the Map.
   headerIterator insertedIt = cards.end();
   --insertedIt;
   cardMap.insert( std::pair<std::string, headerIterator>(card.keyword(), insertedIt) );
   
}

void fitsHeader::append(fitsHeader & head)
{
   headerIterator it;
   
   for(it = head.begin(); it != head.end(); ++it)
   {
      append(*it);
   }
}

void fitsHeader::append(const std::string &k)
{
   append(fitsHeaderCard(k));
}

void fitsHeader::insert_before( headerIterator it, 
                                fitsHeaderCard card
                              )
{
   //First check if duplicate key
   if(cardMap.count(card.keyword()) > 0)
   {
      if(card.type() !=  fitsType<fitsCommentType>() && card.type() != fitsType<fitsHistoryType>())
      {
         std::cerr << "attempt to duplicate keyword: " << card.keyword() << "\n";
         return;
      }
   }
   
   //Now insert in list
   headerIterator insertedIt = cards.insert(it, card);
   
   //Then add to the Map.
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword(), insertedIt));
      
}

void fitsHeader::insert_after( headerIterator it, 
                               fitsHeaderCard card
                             )
{
   //First check if duplicate key
   if(cardMap.count(card.keyword()) > 0)
   {
      if(card.type() !=  fitsType<fitsCommentType>() && card.type() != fitsType<fitsHistoryType>())
      {
         std::cerr << "attempt to duplicate keyword:" << card.keyword() << "\n";
         return;
      }
   }
   
   //Now insert in list
   headerIterator insertedIt = cards.insert(++it, card);
   
   //Then add to the Map.
   cardMap.insert(std::pair<std::string, headerIterator>(card.keyword(), insertedIt));
   
   
}

fitsHeaderCard & fitsHeader::operator[](const std::string & keyword)
{
   headerIterator it;

   //If not found, append it.
   if(cardMap.find(keyword) == cardMap.end())
   {
      append(keyword);
   }

   it = cardMap.find(keyword)->second;

   return *it;
   
}

const fitsHeaderCard & fitsHeader::operator[](const std::string & keyword) const
{
   if(cardMap.find(keyword) == cardMap.end())
   {
      std::cerr << "Fits header card with keyword: " << keyword << " not found.\n";
   }

   
   headerIterator it = cardMap.find(keyword)->second;
   
   return *it;
}

void fitsHeaderGitStatus( fitsHeader & head, 
                          const std::string & repoName,
                          const char * sha1,
                          int modified
                        )
{
   std::string hist = "Git status for " + repoName + ":";
   head.append("", fitsHistoryType(), hist);
   
   hist = "   sha1=";
   hist += sha1;
   if(modified) hist += ", modified";
   
   head.append("", fitsHistoryType(), hist);
}

} //namespace fits
} //namespace mx
