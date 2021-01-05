/** \file stringUtils.cpp
  * \brief Implementation of utilities for working with strings
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup stringutils
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

#include "ioutils/stringUtils.hpp"

namespace mx
{
namespace ioutils
{


// Specialization of convertToString to avoid converting a string to a string
template<> 
std::string convertToString<std::string>( const std::string & value,
                                          int precision
                                        )
{
   static_cast<void>(precision);
   return value;
}

// Template specialization of convertFromString for char
template<>  
char convertFromString<char>(const std::string & str )
{
   return (char) atoi(str.c_str());
}

// Template specialization of convertFromString for char16_t
template<>  
char16_t convertFromString<char16_t>(const std::string & str )
{
   return (char16_t) atoi(str.c_str());
}

// Template specialization of convertFromString for char32_t
template<>  
char32_t convertFromString<char32_t>(const std::string & str )
{
   return (char32_t) atoi(str.c_str());
}

// Template specialization of convertFromString for char32_t
template<>  
wchar_t convertFromString<wchar_t>(const std::string & str )
{
   return (wchar_t) atoi(str.c_str());
}

// Template specialization of convertFromString for unsigned char
template<>
signed char convertFromString<signed char>(const std::string & str )
{
   return (signed char) atoi(str.c_str());
}

// Template specialization of convertFromString for unsigned char
template<>
unsigned char convertFromString<unsigned char>(const std::string & str )
{
   return (unsigned char) atoi(str.c_str());
}


// Template specialization of convertFromString for short
template<>
short convertFromString<short>(const std::string & str )
{
   return (short) atoi(str.c_str());
}

// Template specialization of convertFromString for unsigned short
template<>
unsigned short convertFromString<unsigned short>(const std::string & str )
{
   return (unsigned short) atoi(str.c_str());
}

// Template specialization of convertFromString for int
template<> 
int convertFromString<int>(const std::string & str )
{
   return atoi(str.c_str());
}

// Template specialization of convertFromString for unsigned int
template<>
unsigned int convertFromString<unsigned int>(const std::string & str )
{
   return (unsigned int) strtoul(str.c_str(),0,0);
}

// Template specialization of convertFromString for long
template<> 
long convertFromString<long>(const std::string & str )
{
   return strtol(str.c_str(), 0, 0);
}

// Template specialization of convertFromString for unsigned long
template<>
unsigned long convertFromString<unsigned long>(const std::string & str)
{
   return strtoul(str.c_str(), 0, 0);
}

// Template specialization of convertFromString for long long
template<>
long long convertFromString<long long>(const std::string & str )
{
   return strtoll(str.c_str(), 0, 0);
}

// Template specialization of convertFromString for unsigned long long
template<>
unsigned long long convertFromString<unsigned long long>(const std::string & str )
{
   return strtoull(str.c_str(), 0, 0);
}

// Template specialization of convertFromString for float
template<>
float convertFromString<float>(const std::string & str )
{
   return strtof(str.c_str(), 0);
}

// Template specialization of convertFromString for double
template<>
double convertFromString<double>(const std::string & str )
{
   return strtod(str.c_str(), 0);
}

// Template specialization of convertFromString for long double
template<>
long double convertFromString<long double>(const std::string & str )
{
   return strtold(str.c_str(), 0);
}

// Template specialization of convertFromString for bool
template<>
bool convertFromString<bool>(const std::string & str )
{
   char c = str[0];
   size_t i=0;
   while(isspace(c) && i < str.length()) c = str[i++];

   if(c == '0' || c == 'f' || c == 'F') return false;
   if(c == '1' || c == 't' || c == 'T') return true;

   return (bool) convertFromString<int>(str);
}

// Convert a string to all lower case.
void toLower( std::string &outstr, 
              const std::string & instr 
            )
{
   outstr.resize(instr.size());

   for(size_t i=0; i < instr.size(); ++i) outstr[i] = tolower(instr[i]);

}

// Convert a string to all lower case.
std::string toLower(const std::string & instr )
{
   std::string outstr;

   toLower(outstr, instr);

   return outstr;
}

//vConvert a string to all upper case.
void toUpper( std::string &outstr, 
              const std::string & instr 
            )
{
   outstr.resize(instr.size());

   for(size_t i=0; i < instr.size(); ++i) outstr[i] = toupper(instr[i]);

}

// Convert a string to all upper case.
std::string toUpper(const std::string & instr )
{
   std::string outstr;

   toUpper(outstr, instr);

   return outstr;
}

// Remove all white space from a string.
void removeWhiteSpace( std::string & outstr,  
                       const std::string & instr 
                     )
{
   outstr = instr;

   outstr.erase(std::remove_if(outstr.begin(), outstr.end(), ::isspace), outstr.end());
}

// Remove all white space from a string.
std::string removeWhiteSpace( const std::string & instr )
{
   std::string outstr;

   removeWhiteSpace(outstr, instr);

   return outstr;
}

// Wrap a string by breaking it into smaller sized portions of a desired width
int stringWrap( std::vector<std::string> & lines, 
                const std::string & str,
                int width 
              )
{
   int L = str.length();

   if(L == 0) return 0;
   int startPos, tmpPos, endPos;

   bool done = false;

   startPos = 0;

   while( !done )
   {
      if(startPos == L) --startPos; //Just to make sure

      endPos = startPos + width;

      if(endPos >= L)
      {
         lines.push_back( str.substr( startPos, L-startPos ));
         done = true;
      }
      else
      {
         //Backup to nearest space
         tmpPos = endPos;
         while( !isspace(str[tmpPos]) && tmpPos >= startPos ) --tmpPos;

         //If we aren't at the beginning (i.e. splitting consecutive characters) we use new end position
         if(tmpPos > startPos) endPos = tmpPos;


         lines.push_back( str.substr( startPos, endPos-startPos) );


         startPos = endPos;

         //Clear 1 space
         if( str[startPos] == ' ') ++startPos;
      }
   }

   return 0;
}


} //namespace ioutils
} //namespace mx
