/** \file trueFalseT.hpp
  * \brief Declares and defines a true-false virtual type used for boolean template overrides
  * \ingroup utils_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef trueFalseT_hpp
#define trueFalseT_hpp

namespace mx
{

namespace meta
{
   
///Template declaration of a trueFalseT type.
/** The specialized types can be used for tag dispatching, that is achieving SFINAE-like behavior but with function overloading when
  * identical function signatures are desired.  This works because trueFalseT\<true\> and trueFalseT\<false\> are
  * dfferent types.
  * 
  * This is a slightly simpler construct than std::true_type and std::false_type, which are derived from std::integral_constant.
  * 
  * This template type is not defined, only the specializations are. The specializations have a typedef of True or False, accordingly,
  * and a static const member value which evaluates to true or false accordingly.
  * 
  * \tparam trueFalse bool to choose which of the specializations to invoke.
  * 
  * \ingroup meta
  */   
template<bool trueFalse>
struct trueFalseT;

///The true specialization of trueFalseT.
/** \ingroup meta
  */
template<>
struct trueFalseT<true>
{
   typedef bool True; ///< Typedef which can be used for SFINAE
   static const bool value = true; ///< bool member value = true
};

///The false specialization of trueFalseT.
/** \ingroup meta
  */
template<>
struct trueFalseT<false>
{
   typedef bool False; ///< Typedef which can be used for SFINAE
   static const bool value = false; ///< bool member value = false
};

} //namespace meta
} //namespace mx


#endif //trueFalseT_hpp
