/** \file binVector.hpp
  * \author Jared R. Males
  * \brief A utility to read/write vectors of data from/to a binary file.
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

#ifndef binVector_hpp
#define binVector_hpp


#include <vector>
#include <complex>

#include "../mxlib.hpp"
#include "../mxError.hpp"

namespace mx
{
namespace ioutils
{

/** \defgroup binvector binVector Binary File Format
  * \ingroup ioutils
  * \brief A simple binary file format for storing vectors of data on disk.
  *
  * The binVector file format is very simple: the first 8 bytes contain a unit64_t integer which specifies
  * the type of the data. The second 8 bytes contain a uint64_t integer which specifies the length L of the data vector.
  * The remaining L*sizeof(dataT) bytes contain the data.
  *
  * The suggested extension for BinVector files is ".binv".
  *
  *
  */

/// The type used for binVector sizes.
/** \ingroup binVector
  */
typedef uint64_t binVSizeT;

/// The type of binVector type codes.
/** \ingroup binvector
  */
typedef uint64_t binVTypeT;

///Namespace fo the bin-vector type-codes.
/** \ingroup binvector
  */
namespace binVTypes
{

/// The pre-defined type codes for binVector.
/** \ingroup binvector
  */
enum : binVTypeT { Bool = 0,
                   SChar = 1,
                   UChar = 2,
                   Char = 3,
                   WChar = 4,
                   Char16 = 5,
                   Char32 = 6,
                   Int = 7,
                   UInt = 8,
                   SInt = 9,
                   SUInt = 10,
                   LInt = 11,
                   LUInt = 12,
                   LLInt = 13,
                   LLUInt = 14,
                   Float = 15,
                   Double = 16,
                   LDouble = 17,
                   Quad = 18,
                   CFloat = 19,
                   CDouble = 20,
                 };

} //namespace binVTyes

///Get the integer type code corresponding to the type.
/**
  * \returns an integer which uniquely identifies the type.
  *
  * \tparam dataT is the type
  *
  * \ingroup binvector
  */
template<typename dataT>
binVTypeT binVectorTypeCode();

template<>
binVTypeT binVectorTypeCode<bool>()
{
   return binVTypes::Bool;
}

template<>
binVTypeT binVectorTypeCode<signed char>()
{
   return binVTypes::SChar;
}

template<>
binVTypeT binVectorTypeCode<unsigned char>()
{
   return binVTypes::UChar;
}

template<>
binVTypeT binVectorTypeCode<char>()
{
   return binVTypes::Char;
}

template<>
binVTypeT binVectorTypeCode<wchar_t>()
{
   return binVTypes::WChar;
}

template<>
binVTypeT binVectorTypeCode<char16_t>()
{
   return binVTypes::Char16;
}

template<>
binVTypeT binVectorTypeCode<char32_t>()
{
   return binVTypes::Char32;
}

template<>
binVTypeT binVectorTypeCode<int>()
{
   return binVTypes::Int;
}

template<>
binVTypeT binVectorTypeCode<unsigned int>()
{
   return binVTypes::UInt;
}

template<>
binVTypeT binVectorTypeCode<short int>()
{
   return binVTypes::SInt;
}

template<>
binVTypeT binVectorTypeCode<short unsigned int>()
{
   return binVTypes::SUInt;
}

template<>
binVTypeT binVectorTypeCode<long int>()
{
   return binVTypes::LInt;
}

template<>
binVTypeT binVectorTypeCode<long unsigned int>()
{
   return binVTypes::LUInt;
}

template<>
binVTypeT binVectorTypeCode<long long int>()
{
   return binVTypes::LLInt;
}

template<>
binVTypeT binVectorTypeCode<long long unsigned int>()
{
   return binVTypes::LLUInt;
}


template<>
binVTypeT binVectorTypeCode<float>()
{
   return binVTypes::Float;
}

template<>
binVTypeT binVectorTypeCode<double>()
{
   return binVTypes::Double;
}

template<>
binVTypeT binVectorTypeCode<long double>()
{
   return binVTypes::LDouble;
}

template<>
binVTypeT binVectorTypeCode<std::complex<float>>()
{
   return binVTypes::CFloat;
}

template<>
binVTypeT binVectorTypeCode<std::complex<double>>()
{
   return binVTypes::CDouble;
}

/// Read a BinVector file from disk.
/**
  * \note dataT must match what was stored in the file.
  *
  * \returns 0 on success.
  * \returns -1 if an error occurs.
  *
  * \ingroup binvector
  */
template<typename dataT>
int readBinVector( std::vector<dataT> & vec, ///< [out] vec is a vector which will be resized and populated.
                   const std::string & fname ///<  [in] fname is the name (full-path) of the file.
                 )
{
   FILE *fin;
   binVTypeT typecode;
   binVSizeT sz;
   size_t nrd;

   fin = fopen(fname.c_str(), "r");
   if(fin == 0)
   {
      mxPError("readBinVector", errno, "Error from fopen [" + fname + "]");
      return -1;
   }

   errno = 0;
   nrd = fread(&typecode, sizeof(binVTypeT), 1, fin);

   if(nrd != 1)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("readBinVector", errno, "Error reading data size [" + fname + "]");
      }
      else
      {
         mxError("readBinVector", MXE_FILERERR, "Error reading data size, did not read enough bytes. [" + fname + "]");
      }
      fclose(fin);
      return -1;
   }

   if( typecode != binVectorTypeCode<dataT>() )
   {
      mxError("readBinVector", MXE_SIZEERR, "Mismatch between type dataT and type in file [" + fname + "]");
      fclose(fin);
      return -1;
   }


   errno = 0;
   nrd = fread(&sz, sizeof(binVSizeT), 1, fin);

   if(nrd != 1)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("readBinVector", errno, "Error reading vector size [" + fname + "]");
      }
      else
      {
         mxError("readBinVector", MXE_FILERERR, "Error reading vector size, did not read enough bytes [" + fname + "]");
      }
      fclose(fin);
      return -1;
   }

   vec.resize(sz);

   errno = 0;
   nrd = fread(vec.data(), sizeof(dataT), sz, fin);

   if(nrd != sz)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("readBinVector", errno, "Error reading data [" + fname + "]");
      }
      else
      {
         mxError("readBinVector", MXE_FILERERR, "Did not read enough data [" + fname + "]");
      }
      fclose(fin);
      return -1;
   }

   fclose(fin);

   return 0;
}

/// Write a BinVector file to disk.
/**
  *
  * \returns 0 on success.
  * \returns -1 if an error occurs.
  *
  * \ingroup binvector
  */
template<typename dataT>
int writeBinVector( const std::string & fname, ///< [in] fname is the name (full-path) of the file.
                    std::vector<dataT> & vec ///< [in] vec is the vector which will be written to disk.
                  )
{

   FILE *fout;
   size_t nwr;
   binVTypeT typecode = binVectorTypeCode<dataT>();
   binVSizeT sz = vec.size();

   fout = fopen(fname.c_str(), "wb");
   if(fout == 0)
   {
      mxPError("writeBinVector", errno, "Error from fopen [" + fname + "]");
      return -1;
   }

   nwr = fwrite( &typecode, sizeof(binVTypeT), 1, fout);
   if(nwr != 1)
   {
      mxPError("writeBinVector", errno, "Error writing typecode [" + fname + "]");
      fclose(fout);
      return -1;
   }

   nwr = fwrite( &sz, sizeof(binVSizeT), 1, fout);
   if(nwr != 1)
   {
      mxPError("writeBinVector", errno, "Error writing vector size [" + fname + "]");
      fclose(fout);
      return -1;
   }

   nwr = fwrite( vec.data(), sizeof(dataT), vec.size(), fout);
   if(nwr != sz)
   {
      mxPError("writeBinVector", errno, "Error writing data [" + fname + "]");
      fclose(fout);
      return -1;
   }

   fclose(fout);

   return 0;
}


} //namespace ioutils
} //namespace mx

#endif //binVector_hpp
