/** \file textTable.hpp
  * \brief declares and defines a simple text table manager
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup asciiutils
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

#include "../mxlib.hpp"

#include "stringUtils.hpp"

/// The mxlib c++ namespace
namespace mx
{
namespace ioutils
{

/** \addtogroup asciiutils
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
  * \todo add sort capability, with vector of column choices to sort first by column N, then M, etc.
  */
struct textTable
{
   std::vector<int> m_colWidths; ///< The widths of each column, not including the separator.

   ///The table cells.
   /** Organized as rows.cells.lines, where a cell could be multiple lines.
     */
   std::vector<std::vector<std::vector<std::string>>> m_rows;

   std::string m_lineStart; ///< Text to print at the beginning of the line.

   std::string m_lineEnd; ///< Text to print at the end of the line.

   std::string m_colSep; ///< Text to print between each column.

   std::string m_rowSep; ///< Text to print between each row.

   ///Add one cell to the table, overwriting if it already exists.
   void addCell( size_t row, ///< [in] the row of the cell.
                 size_t col, ///< [in] the column of the cell.
                 const std::string & cell  ///< [in] the new contents of the cell, will be wrapped.
               );

   void addCell( size_t row, ///< [in] the row of the cell.
                 size_t col, ///< [in] the column of the cell.
                 const char * cell  ///< [in] the new contents of the cell, will be wrapped.
               );

   ///Add one cell to the table, overwriting if it already exists.
   template<typename typeT>
   void addCell( size_t row, ///< [in] the row of the cell.
                 size_t col, ///< [in] the column of the cell.
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

template<typename typeT>
void textTable::addCell(size_t row,
                        size_t col,
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

   size_t width = 0;
   for(size_t i=0;i<m_colWidths.size();++i) width += m_colWidths[i];
   width += m_colWidths.size()*(m_colSep.length()-1); //+100;

   line.resize(width, ' ');

   for(size_t i=0; i< m_rows.size(); ++i)
   {
      size_t rowL = 0;
      for(size_t j=0; j< m_colWidths.size(); ++j)
      {
         if( m_rows[i][j].size() > rowL ) rowL = m_rows[i][j].size();
      }

      for(size_t k=0; k< rowL; ++k)
      {
         //line.insert( 0, line.length(), ' ');
         line.clear();
         line.resize(width, ' ');
         //line.replace(0, width, width, ' ');

         size_t startPos = 0;
         for(size_t j=0; j< m_colWidths.size(); ++j)
         {
            if(m_rows[i][j].size() > k) line.replace(startPos, m_rows[i][j][k].length(), m_rows[i][j][k]);
            startPos += m_colWidths[j];

            if( j < m_colWidths.size()-1) line.replace(startPos, m_colSep.length(), m_colSep);
            startPos += m_colSep.length();

         }

         ios << m_lineStart << line << m_lineEnd << "\n";
      }

     if(m_rowSep.length() > 0) ios << m_rowSep << "\n";

   }
}

extern template
void textTable::outPut<std::ostream>(std::ostream & ios);

/// @}

} //namespace ioutils
} //namespace mx

#endif //textTable_hpp
