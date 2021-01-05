/** \file rawBinary.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides functions for working with raw binary files.
  * \ingroup utils_files
  *
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef rawBinary_hpp
#define rawBinary_hpp

#include <cstdio>
#include <string>

#include "../mxlib.hpp"
#include "../mxError.hpp"

namespace mx
{
namespace ioutils
{

/// Write an array of data to file as raw binary.
/** Here raw binary means no formatting or metadata,
  * just the bytes pointed to by the array are written to
  * disk.
  *
  * \ingroup ioutils
  *
  * \returns 0 on success
  * \returns -1 on error
  */
template<typename T>
int writeRawBinary( const std::string & fileName, ///< [in] the file to write to
                    T * data,                     ///< [in] the data pointer
                    size_t szData)                ///< [in] number of elements of sizeof(T) to write
{
   FILE * fout;

   fout = fopen(fileName.c_str(), "wb");
   if(fout == 0)
   {
      mxPError("writeRawBinary", errno, "Error from fopen [" + fileName + "]");
      return -1;
   }


   int nwr = fwrite( data, sizeof(T), szData, fout);

   if(nwr != szData)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("writeRawBinary", errno, "Error writing to file [" + fileName + "]");
      }
      else
      {
         mxError("writeRawBinary", MXE_FILERERR, "Error writing to file, did not write all elements. [" + fileName+ "]");
      }
      fclose(fout);
      return -1;
   }

   int res = fclose(fout);

   if(res != 0)
   {
      mxPError("writeRawBinary", errno, "Error closing file [" + fileName + "]");
      return -1;
   }

   return 0;
}

} //namespace ioutils
} //namespace mx

#endif //rawBinary_hpp
