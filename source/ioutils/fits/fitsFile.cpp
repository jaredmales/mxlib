/** \file fitsFile.cpp
 * \brief Definitions for a class to work with FITS files
 * \ingroup fits_processing_files
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

#include "ioutils/fits/fitsFile.hpp"

namespace mx
{
namespace fits
{

template class fitsFile<short,verbose::o>;
template class fitsFile<short,verbose::v>;
template class fitsFile<short,verbose::vv>;
template class fitsFile<short,verbose::vvv>;

template class fitsFile<unsigned short,verbose::o>;
template class fitsFile<unsigned short,verbose::v>;
template class fitsFile<unsigned short,verbose::vv>;
template class fitsFile<unsigned short,verbose::vvv>;

template class fitsFile<int,verbose::o>;
template class fitsFile<int,verbose::v>;
template class fitsFile<int,verbose::vv>;
template class fitsFile<int,verbose::vvv>;

template class fitsFile<unsigned int,verbose::o>;
template class fitsFile<unsigned int,verbose::v>;
template class fitsFile<unsigned int,verbose::vv>;
template class fitsFile<unsigned int,verbose::vvv>;

template class fitsFile<float,verbose::o>;
template class fitsFile<float,verbose::v>;
template class fitsFile<float,verbose::vv>;
template class fitsFile<float,verbose::vvv>;

template class fitsFile<double,verbose::o>;
template class fitsFile<double,verbose::v>;
template class fitsFile<double,verbose::vv>;
template class fitsFile<double,verbose::vvv>;

} // namespace fits
} // namespace mx
