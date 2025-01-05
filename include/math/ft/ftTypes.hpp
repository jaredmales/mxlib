/** \file ftTypes.hpp
 * \brief Fourier Transform Types
 * \ingroup ft_files
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

#ifndef ftTypes_hpp
#define ftTypes_hpp

namespace mx
{
namespace math
{
namespace ft
{

/// Directions of the Fourier Transform
/** \ingroup ft
 *
 */
enum class dir{ forward,  ///< Specifies the forward transform
                backward, ///< Specifies the backward transform
              };


} // namespace ft
} // namespace math
} // namespace mx


#endif // ftTypes_hpp
