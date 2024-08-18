/** \file textTable.cpp
 * \brief implementation for a simple text table manager
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 * \ingroup asciiutils
 *
 */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#include "ioutils/textTable.hpp"

/// The mxlib c++ namespace
namespace mx
{
namespace ioutils
{
void textTable::addCell( size_t row, size_t col, const std::string &cell )
{
    // Increase size if needed.
    if( row >= m_rows.size() )
    {
        size_t N = m_rows.size();
        for( size_t i = 0; i < row - N + 1; ++i )
        {
            m_rows.push_back( std::vector<std::vector<std::string>>( m_colWidths.size() ) );
        }
    }

    m_rows[row][col].clear();
    stringWrap( m_rows[row][col], cell, m_colWidths[col] );
}

void textTable::addCell( size_t row, size_t col, const char *cell )
{
    addCell( row, col, std::string( cell ) );
}

template void textTable::outPut<std::ostream>( std::ostream &ios );

} // namespace ioutils
} // namespace mx
