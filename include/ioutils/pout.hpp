/** \file pout.hpp
 * \author Jared R. Males
 * \brief Declaration and definition of a simple formatted output function.
 * \ingroup utils_files
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

#ifndef ioutils_pout_hpp
#define ioutils_pout_hpp

#include <iostream>

#include "../mxlib.hpp"

namespace mx
{
namespace ioutils
{

template<char space, bool flush, char eol>
void pout()
{
   if(eol) std::cout << eol;
   if(flush) std::cout.flush();

   return;
}

///A simple formatted output function.
/** This function writes its arguments, of any type and of any number, to stdout.  By default, the
  * arguments are separated by a space, the new line '\\n' is written at the end, and the std::cout
  * stream is flushed.  These behaviors can be altered via template parameters.
  *
  * Example:
  * \code
  * std::string s="output:";
  * double d=2.567;
  * int i=3;
  *
  * pout(s,d,i);
  * \endcode
  * Note that the types of the values do not need to be specified as templated arguments.
  * When run, this code results in
  * \verbatim
  * $ output: 2.567 3
  * $
  * \endverbatim
  *
  * The behavior can be changed with template arguments like
  * \code
  * pout<'\t', false, 0>(s,d,i);
  * \endcode
  * which would use the tab instead of space, and neither flush nor print a newline at the end.
  *
  * \tparam space is the character to place between the values,  by default this is space.
  * \tparam flush controls whether std::cout.flush() is called, by default it is true.
  * \tparam eol is the character to print at end of line, by default it is '\\n'.
  * \tparam valT a type which can be output by std::cout
  * \tparam valTs a variadic list of additional types which can be output by std::cout
  *
  * \ingroup ioutils
  */
template<char space=' ', bool flush=true, char eol='\n', typename valT, typename... valTs>
void pout( valT value,            ///<  [in] a value to print.
           const valTs&... values ///< [in] a variadic list of additional values. Any number of values can be specified, with any type handled by std::cout.
         )
{
   static const unsigned short int nargs = sizeof...(valTs);

   std::cout << value;

   if(nargs > 0) std::cout << space;

   pout<space,flush,eol>(values...);

}

} //namespace ioutils 
} //namespace mx

#endif //ioutils_pout_hpp
