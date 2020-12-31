/** \file fft.cpp
  * \brief The fast Fourier transform interface
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#include "math/fft/fft.hpp"

namespace mx
{
namespace math
{
namespace fft 
{
   
template<>
std::vector<int> fftwDimVec<1>( int szX, 
                                int szY, 
                                int szZ
                              )
{
   static_cast<void>(szY);
   static_cast<void>(szZ);
   
   std::vector<int> v({szX});
   return v;
}

template<>
std::vector<int> fftwDimVec<2>( int szX, 
                                int szY, 
                                int szZ
                              )
{
   static_cast<void>(szZ);
   
   std::vector<int> v({szX, szY});
   return v;
}

template<>
std::vector<int> fftwDimVec<3>( int szX, 
                                int szY, 
                                int szZ
                              )
{
   std::vector<int> v({szX, szY, szZ});
   return v;
}




}//namespace ffit 
}//namespace math
}//namespace mx


