/** \file textTable.hpp
  * \brief declares and defines a simple text table manager
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup stringutils
  *
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

#ifndef textTable_hpp
#define textTable_hpp

#include <string>
#include <sstream>


#include "stringUtils.hpp"

/// The mxlib c++ namespace
namespace mx
{

/** \addtogroup stringutils
  * @{
  */

///An ascii table formatter
/** Manages a table of string content, including word-wrapping within cells.
  * Contents are output to an ostream line by line with a table structure.
  * 
  * \todo have ability to add whole row, with a variadic template to handle different types
  * \todo allow non-wrapping (i.e. unlimited width) columns
  * \todo headers and footers
  * \todo add centering (h and v) and right and bottom justification
  * \todo have a latex() setup function
  */ 
struct textTable
{
   std::vector<int> colWidths; ///< The widths of each column, not including the separator.
   
   std::vector<std::vector<std::vector<std::string>>> rows; ///<The table cells.  A cell could be multiple lines.
 
   std::string lineStart; ///< Text to print at the beginning of the line.
   
   std::string lineEnd; ///< Text to print at the end of the line.
   
   std::string colSep; ///< Text to print between each column.
   
   std::string rowSep; ///< Text to print between each row.

   ///Add one cell to the table, overwriting if it already exists.
   void addCell( int row, ///< [in] the row of the cell.
                 int col, ///< [in] the column of the cell.
                 const std::string & cell  ///< [in] the new contents of the cell, will be wrapped.
               );
   
   void addCell( int row, ///< [in] the row of the cell.
                 int col, ///< [in] the column of the cell.
                 const char * cell  ///< [in] the new contents of the cell, will be wrapped.
               );
   
   ///Add one cell to the table, overwriting if it already exists.
   template<typename typeT>
   void addCell( int row, ///< [in] the row of the cell.
                 int col, ///< [in] the column of the cell.
                 const typeT & cell,  ///< [in] the new contents of the cell, will be wrapped.
                 int precision=0
               );
   
   ///Output the table to a stream
   /** Prints the table to the stream, with separators, etc.
     *
     * \tparam iosT is a std::ostream-like type.
     */ 
   template<typename iosT>
   void outPut( iosT & ios /**< [in] a std::ostream-like stream.*/ );

};

inline
void textTable::addCell( int row,
                         int col,
                         const std::string & cell 
                       )
{
   //Increase size if needed.
   if( row >= rows.size() )
   {
      int N = rows.size();
      for(int i=0; i < row - N +1; ++i)
      {
         rows.push_back( std::vector<std::vector<std::string>>( colWidths.size()));
      }
   }
   
   rows[row][col].clear();
   stringWrap(rows[row][col], cell, colWidths[col]); 
}

void textTable::addCell(int row,
                        int col,
                        const char * cell
                       )
{
   addCell(row, col, std::string(cell));
}

template<typename typeT>
void textTable::addCell(int row,
                        int col,
                        const typeT & cell,
                        int precision
                       )
{
   addCell(row, col, convertToString<typeT>(cell, precision));
}

template<typename iosT>
void textTable::outPut( iosT & ios )
{
   std::string line;
   
   int width = 0;
   for(int i=0;i<colWidths.size();++i) width += colWidths[i];
   width += colWidths.size()*(colSep.length()-1); //+100;
   
   line.resize(width, ' ');
   
   for(int i=0; i< rows.size(); ++i)
   {
      int rowL = 0;
      for(int j=0; j< colWidths.size(); ++j)
      {
         if( rows[i][j].size() > rowL ) rowL = rows[i][j].size();
      }
      
      for(int k=0; k< rowL; ++k)
      {
         //line.insert( 0, line.length(), ' ');
         line.clear();
         line.resize(width, ' ');
         //line.replace(0, width, width, ' ');
         
         int startPos = 0;
         for(int j=0; j< colWidths.size(); ++j)
         {
            if(rows[i][j].size() > k) line.replace(startPos, rows[i][j][k].length(), rows[i][j][k]);
            startPos += colWidths[j];
           
            if( j < colWidths.size()-1) line.replace(startPos, colSep.length(), colSep);
            startPos += colSep.length();

         }
         
         ios << lineStart << line << lineEnd << "\n";
      }
      
     if(rowSep.length() > 0) ios << rowSep << "\n";

   }
}

/// @}

} //namespace mx

#endif //textTable_hpp
