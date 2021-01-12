/** \file fitsHeaderCard.cpp
  * \brief Definitiions for a class to work with a FITS header card
  * \ingroup fits_processing_files
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

#include "ioutils/fits/fitsHeaderCard.hpp"

namespace mx
{
namespace fits
{
   
fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                void * v,
                                const std::string & c
                              )
{
   keyword = k;
   value = (char *) v;
   type = getFitsType<fitsUnknownType>();
   comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                const void * v,
                                const std::string & c
                              )
{
   keyword = k;
   value = (char *) v;
   type = getFitsType<fitsUnknownType>();
   comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                fitsCommentType v,
                                const std::string & c
                              )
{
   static_cast<void>(v); //be unused
   
   keyword = k;
   value = "";
   type = getFitsType<fitsCommentType>();
   comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                fitsHistoryType v,
                                const std::string & c)
{
   static_cast<void>(v); //be unused
   
   keyword = k;
   value = "";
   type = getFitsType<fitsHistoryType>();
   comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                void * v
                              )
{
   keyword = k;
   value = (char *) v;
   type = getFitsType<fitsUnknownType>();
}

fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                const void * v
                              )
{
   keyword = k;
   value = (char *) v;
   type = getFitsType<fitsUnknownType>();
}

fitsHeaderCard::fitsHeaderCard( const std::string & k)
{
   keyword = k;
   value = "";
   type = getFitsType<fitsUnknownType>();
}

//Specialization for strings, removing the apostrophe.
template<>
std::string fitsHeaderCard::Value<std::string>() const
{
   std::string s = value;
   fitsStripApost(s);
   return s;
}

std::string fitsHeaderCard::String()
{
   std::string s = value;
   fitsStripApost(s);
   return s;
}

int fitsHeaderCard::Int()
{
   return ioutils::convertFromString<int>(value);
}

unsigned int fitsHeaderCard::unsignedInt()
{
   return ioutils::convertFromString<unsigned int>(value);
}

long fitsHeaderCard::Long()
{
   return ioutils::convertFromString<long>(value);
}

unsigned long fitsHeaderCard::unsignedLong()
{
   return ioutils::convertFromString<unsigned long>(value);
}

long long fitsHeaderCard::longLong()
{
   return ioutils::convertFromString<long long>(value);
}

unsigned long long fitsHeaderCard::unsignedLongLong()
{
   return ioutils::convertFromString<unsigned long long>(value);
}

float fitsHeaderCard::Float()
{
   return ioutils::convertFromString<float>(value);
}

double fitsHeaderCard::Double()
{
   return ioutils::convertFromString<double>(value);
}

long double fitsHeaderCard::longDouble()
{
   return ioutils::convertFromString<long double>(value);
}

int fitsHeaderCard::write(fitsfile *fptr)
{
   switch(type)
   {
      case TSTRING:
      {
         return fits_write_key<char *>(fptr, (char *) keyword.c_str(), (void *)value.c_str(), (char *) comment.c_str());
      }
      case TBYTE:
      {
         unsigned char b = Int();
         return fits_write_key<unsigned char>(fptr, (char *) keyword.c_str(), &b, (char *) comment.c_str());
      }
      case TSBYTE:
      {
         signed char sc = Int();
         return fits_write_key<signed char>(fptr, (char *) keyword.c_str(), &sc, (char *) comment.c_str());
      }
      case TSHORT:
      {
         short s = Int();
         return fits_write_key<short>(fptr, (char *) keyword.c_str(), &s, (char *) comment.c_str());
      }
      case TUSHORT:
      {
         unsigned short us = Int();
         return fits_write_key<unsigned short>(fptr, (char *) keyword.c_str(), &us, (char *) comment.c_str());
      }
      case TINT:
      {
         int i = Int();
         return fits_write_key<int>(fptr, (char *) keyword.c_str(), &i, (char *) comment.c_str());
      }
      case TLONG:
      {
         long l = Long();
         return fits_write_key<long>(fptr, (char *) keyword.c_str(), &l, (char *) comment.c_str());
      }
      case TLONGLONG:
      {
         long long ll = longLong();
         return fits_write_key<long long>(fptr, (char *) keyword.c_str(), &ll, (char *) comment.c_str());
      }
      case TULONG:
      {
         unsigned long ul = unsignedLong();
         return fits_write_key<unsigned long>(fptr, (char *) keyword.c_str(), &ul, (char *) comment.c_str());
      }
      case TFLOAT:
      {
         float f = Float();
         return fits_write_key<float>(fptr, (char *) keyword.c_str(), &f, (char *) comment.c_str());
      }
      case TDOUBLE:
      {
         double d = Double();
         return fits_write_key<double>(fptr, (char *) keyword.c_str(), &d, (char *) comment.c_str());
      }
      case fitsTCOMMENT:
      {
         return fits_write_comment(fptr, (char *) comment.c_str());
      }
      case fitsTHISTORY:
      {
         return fits_write_history(fptr, (char *) comment.c_str());
      }
      default:
      {
         return fits_write_key<fitsUnknownType>(fptr, (char *) keyword.c_str(), (void *)value.c_str(), (char *) comment.c_str());
      }
   }
}

} //namespace fits
} //namespace mx


