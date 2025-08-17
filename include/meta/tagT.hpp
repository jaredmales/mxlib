/** \file tagT.hpp
 * \brief Declares and defines an empty struct for tag dispatching
 * \ingroup utils_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef meta_tagT_hpp
#define meta_tagT_hpp

namespace mx
{

namespace meta
{

/// Empty type for tag dispatching
/**
 * \tparam T can be any type
 *
 * \ingroup meta
 */
template <typename T>
struct tagT
{
};


} // namespace meta
} // namespace mx

#endif // meta_tagT_hpp
