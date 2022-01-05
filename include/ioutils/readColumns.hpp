/** \file readColumns.hpp
  * \author Jared R. Males
  * \brief A utility to read in columns from a text file.
  * \ingroup asciiutils
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __readColumns_hpp__
#define __readColumns_hpp__

#include <fstream>
#include <string>
#include <cstring>
#include <iostream>

#include "../mxlib.hpp"
#include "../mxError.hpp"

#include "stringUtils.hpp"

#define MX_READCOL_MISSINGVALSTR "-99"

namespace mx
{
namespace ioutils
{

template<char delim=' ', char eol='\n'>
void readcol(char * sin, int sz)
{
   static_cast<void>(sin);
   static_cast<void>(sz);
   
   return;
}

template<char delim=' ', char eol='\n', typename arrT, typename... arrTs>
void readcol(char * sin, int sz, arrT & array, arrTs &... arrays)
{
   //static const unsigned short int nargs = sizeof...(arrTs);
   std::string str;

   int i=0;
   int l = strlen(sin);

   if(l < 1) return;

   //Eat white space
   while( isspace(sin[i]) && sin[i] != eol && i < l) ++i;
   sin = sin + i;
   sz = sz -i;

   //If there's nothing here, we still need to populate the vector
   if(sz <= 1)
   {
      array.push_back(convertFromString<typename arrT::value_type>(""));
      return;
   }

   std::stringstream sinstr(sin);

   std::getline(sinstr, str, delim);

   //Last entry in line might contain eol
   if( str[str.size()-1] == eol)
   {
      str.erase(str.size()-1);
   }

   if( str.size() == 0 )
   {
      array.push_back(convertFromString<typename arrT::value_type>(MX_READCOL_MISSINGVALSTR));
   }
   else
   {
      array.push_back(convertFromString<typename arrT::value_type>(str));
   }

   sin += ( str.size()+1)*sizeof(char);
   sz -= ( str.size()+1)*sizeof(char);

   readcol<delim,eol>(sin, sz, arrays...);

}



///Read in columns from a text file
/** This function opens a file containing data formatted in columns and reads in the data row by row.
  * The data are stored in std::vectors, which should not be pre-allocated (though they could be reserve()-ed).
  *
  * Example:
  * \code
  * std::vector<int> i1;
  * std::vector<float> f1;
  * std::vector<double> d1;
  *
  * readColumns("data_file.txt", i1, f1, d1);
  * \endcode
  *
  * Note that the types of the vectors do not need to be specified as template arguments.
  *
  * The format of the file can be specified with template arguments like
  * \code
  * readColumns<',', ';', '\r'>("data_file.csv", i1, f1, d1);
  * \endcode
  * which sets the delimmiter to comma, the comment character to ;, and the end-of-line to \\r.
  *
  * Columns can be skipped using mx::ioutils::skipCol.
  *
  * \tparam delim is the character separating columns,  by default this is space.
  * \tparam comment is the character starting a comment.  by default this is #
  * \tparam eol is the end of line character.  by default this is \n
  * \tparam arrTs a variadic list of array types. this is not specified by the user.
  *
  * \todo lineSize should be configurable
  * 
  * \ingroup asciiutils
  */
template<char delim=' ', char comment='#', char eol='\n', typename... arrTs>
int readColumns( const std::string & fname, ///< [in] is the file name to read from
                 arrTs &... arrays ///< [out] a variadic list of std::vectors. Any number with mixed value_type can be specified. Neither allocated nor cleared, so repeated calls will append data.
               )
{
   //open file
   errno = 0;
   std::ifstream fin;
   fin.open(fname);

   if(!fin.good())
   {
      if(errno != 0)
      {
         mxPError("readColumns", errno, "Occurred while opening " + fname + " for reading.");
      }
      else
      {
         mxError("readColumns", MXE_FILEOERR, "Occurred while opening " + fname + " for reading.");
      }
      return -1;
   }

   int lineSize = 4096;
   char * line = new char[lineSize];

   while(fin.good())
   {
      //Save one space for adding eol
      fin.getline(line, lineSize-1, eol);

      int i=0;
      int l = strlen(line);

      if(l <= 0) break;

      //std::cerr << line << "\n";

      //Find start of comment and end line at that point.
      while(line[i] != comment )
      {
         ++i;
         if( i == l ) break;
      }

      if(i <= l-1)
      {
         line[i] = '\0';
      }

      l = strlen(line);

      if(l == 0) continue;

      //Make sure line ends with eol
      line[l] = eol;
      ++l;
      line[l] = '\0';

      readcol<delim,eol>(line, strlen(line), arrays...);
   }

   delete[] line;

   //getline will have set fail if there was no new line on the last line.
   if(fin.bad() && !fin.fail())
   {
      if(errno != 0)
      {
         mxPError("readColumns", errno, "Occurred while reading from " + fname + ".");
      }
      else
      {
         mxError("readColumns", MXE_FILERERR, "Occurred while reading from " + fname + ".");
      }
      return -1;
   }

   fin.clear(); //Clear the fail bit which may have been set by getline
   fin.close();

   if(fin.fail())
   {
      if(errno != 0)
      {
         mxPError("readColumns", errno, "Occurred while closing " + fname + ".");
      }
      else
      {
         mxError("readColumns", MXE_FILECERR, "Occurred while closing " + fname + ".");
      }
      return -1;
   }



   return 0;
}

///A dummy class to allow mx::readColumns to skip a column(s) in a file without requiring memory allocation.
/** The alternative is to use dummy vectors, which result in excess memory allocations and deallocations.
  * Usage:
  \code
  std::vector<T> col1, col5;
  skipCol sk;
  readColumns("filename.txt", col1, sk, sk, sk, col5); //This results in only columns 1 and 5 being stored.
  \endcode
  *
  * \ingroup asciiutils
  */
struct skipCol
{
   typedef std::string value_type; ///< value_type is defined as std::string so that no conversions take place.

   template<typename T>
   void push_back( const T & arg )
   {
      return;
   }
};


} //namespace ioutils
} //namespace mx

#endif //__readColumns_hpp__
