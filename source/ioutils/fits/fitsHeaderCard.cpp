/** \file fitsHeaderCard.cpp
  * \brief Definitiions for a class to work with a FITS header card
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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
                                const std::string & v,     
                                const std::string & c
                              )
{
   m_keyword = k;
   m_valueStr.str(v);
   m_valueGood = false;
   m_valueStrGood = true;
   m_type = fitsType<std::string>();
   m_comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                const std::string & v,
                                const int & type,      
                                const std::string & c 
                              )
{
   m_keyword = k;
   m_valueStr.str(v);
   m_valueGood = false;
   m_valueStrGood = true;
   m_type = type;
   m_comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k, 
                                fitsCommentType v,     
                                const std::string & c  
                              )
{
   m_keyword = k;
   m_valueGood = false;
   m_valueStrGood = false;
   m_type = fitsType<fitsCommentType>();
   m_comment = c;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k, 
                                fitsHistoryType v,     
                                const std::string & c  
                              )
{
   m_keyword = k;
   m_valueGood = false;
   m_valueStrGood = false;
   m_type = fitsType<fitsHistoryType>();
   m_comment = c;
}

fitsHeaderCard::fitsHeaderCard(const std::string & k)
{
   m_keyword = k;
}

fitsHeaderCard::fitsHeaderCard( const std::string & k, 
                                const int type         
                              )
{
   m_keyword = k;
   m_type = type;
}


fitsHeaderCard::fitsHeaderCard( const fitsHeaderCard & card )
{
   m_keyword = card.m_keyword;
   m_type = card.m_type;
   memcpy( &m_value, &card.m_value, sizeof(values));
   m_valueStr.str(card.m_valueStr.str());
   m_valueGood = card.m_valueGood;
   m_valueStrGood = card.m_valueStrGood;
   m_comment = card.m_comment;
}

fitsHeaderCard & fitsHeaderCard::operator=(const fitsHeaderCard & card )
{
   m_keyword = card.m_keyword;
   m_type = card.m_type;
   memcpy( &m_value, &card.m_value, sizeof(values));
   m_valueStr.str(card.m_valueStr.str());
   m_valueGood = card.m_valueGood;
   m_valueStrGood = card.m_valueStrGood;
   m_comment = card.m_comment;

   return *this;
}

void fitsHeaderCard::convertToString()
{
   if(!m_valueGood)
   {
      mxThrowException(err::paramnotset, "fitsHeaderCard::convertToString()", "no value to convert for " + m_keyword);
      return;
   }

   if(m_type == fitsType<char *>()|| m_type == fitsType<std::string>())
   {
      m_valueStrGood = true; //It should be hard to get here, but just in case.
      return;
   }

   m_valueStr.str("");
   m_valueStr.precision(10);

   switch(m_type)
   {
      case fitsType<char>():
         m_valueStr << static_cast<int>(m_value.Char);
         break;
      case fitsType<unsigned char>():
         m_valueStr << static_cast<int>(m_value.UChar);
         break;
      case fitsType<short>():
         m_valueStr << m_value.Short;
         break;
      case fitsType<unsigned short>():
         m_valueStr << m_value.UShort;
         break;
      case fitsType<int>():
         m_valueStr << m_value.Int;
         break;
      case fitsType<unsigned int>():
         m_valueStr << m_value.UInt;
         break;
      case fitsType<long>():
         m_valueStr << m_value.Long;
         break;
      case fitsType<unsigned long>():
         m_valueStr << m_value.ULong;
         break;
      case fitsType<long long>():
         m_valueStr << m_value.LongLong;
         break;
      case fitsType<unsigned long long>():
         m_valueStr << m_value.ULongLong;
         break;
      case fitsType<float>():
         m_valueStr << m_value.Float;
         break;
      case fitsType<std::complex<float>>():
         m_valueStr << m_value.complexFloat;
         break;
      case fitsType<double>():
         m_valueStr << m_value.Double;
         break;
      case fitsType<std::complex<double>>():
         m_valueStr << m_value.complexDouble;
         break;
      case fitsType<fitsCommentType>():
         return;
      case fitsType<fitsHistoryType>():
         return;   
      default:
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertToString()", "Unknown FITS type");
      
   }

   m_valueStrGood = true;
   return;
}

template<>
void fitsHeaderCard::convertFromString<char>()
{
   m_type = fitsType<char>();
   m_value.Char = std::stoi(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<unsigned char>()
{
   m_type = fitsType<unsigned char>();
   m_value.UChar = std::stoi(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<short>()
{
   m_type = fitsType<short>();
   m_value.Short = std::stoi(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<unsigned short>()
{
   m_type = fitsType<unsigned short>();
   m_value.UShort = std::stoi(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<int>()
{
   m_type = fitsType<int>();
   m_value.Int = std::stoi(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<unsigned int>()
{
   m_type = fitsType<unsigned int>();
   m_value.UInt = std::stol(m_valueStr.str());
   m_valueGood = true;
}


template<>
void fitsHeaderCard::convertFromString<long>()
{
   m_type = fitsType<long>();
   m_value.Long = std::stol(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<unsigned long>()
{
   m_type = fitsType<unsigned long>();
   m_value.ULong = std::stoul(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<long long>()
{
   m_type = fitsType<long long>();
   m_value.LongLong = std::stoll(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<unsigned long long>()
{
   m_type = fitsType<unsigned long long>();
   m_value.ULongLong = std::stoull(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<float>()
{
   m_type = fitsType<float>();
   m_value.Float = std::stof(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<std::complex<float>>()
{
   mxThrowException(err::notimpl, "fitsHeaderCard::convertFromString<std::complex<float>>",  "no conversion from string to std::complex<float>");
   m_type = fitsType<std::complex<float>>();
   m_value.complexFloat = std::stof(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<double>()
{
   m_type = fitsType<double>();
   m_value.Double = std::stof(m_valueStr.str());
   m_valueGood = true;
}

template<>
void fitsHeaderCard::convertFromString<std::complex<double>>()
{
   mxThrowException(err::notimpl, "fitsHeaderCard::convertFromString<std::complex<double>>",  "no conversion from string to std::complex<double>");
   m_type = fitsType<std::complex<double>>();
   m_value.complexDouble = std::stof(m_valueStr.str());
   m_valueGood = true;
}

template<typename typeT>
typeT fitsHeaderCard::convertedValue()
{
   switch(m_type)
   {
      case fitsType<unsigned char>():
      {
         return m_value.UChar;
      }
      case fitsType<char>():
      {
         return m_value.Char;
      }
      case fitsType<short>():
      {
         return m_value.Short;
      }
      case fitsType<unsigned short>():
      {
         return m_value.UShort;
      }
      case fitsType<int>():
      {
         return m_value.Int;
      }
      case fitsType<unsigned int>():
      {
         return m_value.UInt;
      }
      case fitsType<long>():
      {
         return m_value.Long;
      }
      case fitsType<unsigned long>():
      {
         return m_value.ULong;
      }
      case fitsType<long long>():
      {
         return m_value.LongLong;
      }
      case fitsType<unsigned long long>():
      {
         return m_value.ULongLong;
      }
      case fitsType<float>():
      {
         return m_value.Float;
      }
      case fitsType<std::complex<float>>():
      {
         mxThrowException(err::notimpl, "fitsHeaderCard::convertedValue<typeT>", std::string("conversions no supported for complex types in ") + m_keyword);
      }
      case fitsType<double>():
      {
         return m_value.Double;
      }
      case fitsType<std::complex<double>>():
      {
         mxThrowException(err::notimpl, "fitsHeaderCard::convertedValue<typeT>", std::string("conversions no supported for complex types in ") + m_keyword);
      }
      case fitsType<fitsCommentType>():
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertedValue<typeT>", "cannot convert comment to numeric type");
      }
      case fitsType<fitsHistoryType>():
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertedValue<typeT>", "cannot convert history to numeric type");
      }
      case TSTRING:
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertedValue<typeT>", std::string("cannot convert string to numeric type in ") + m_keyword);
      }
      default:
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertedValue<typeT>", std::string("invalid FITS type for conversion in ") + m_keyword);
      }
   }
}

void fitsHeaderCard::convertValue(int newtype)
{
   if(!m_valueGood)
   {
      m_type = newtype;
      return;
   }

   switch(newtype)
   {
      case fitsType<unsigned char>():
      {
         m_value.UChar = convertedValue<unsigned char>();
         break;
      }
      case fitsType<char>():
      {
         m_value.Char = convertedValue<char>();
         break;
      }
      case fitsType<short>():
      {
         m_value.Short = convertedValue<short>();
         break;
      }
      case fitsType<unsigned short>():
      {
         m_value.UShort = convertedValue<unsigned short>();
         break;
      }
      case fitsType<int>():
      {
         m_value.Int = convertedValue<int>();
         break;
      }
      case fitsType<unsigned int>():
      {
         m_value.UInt = convertedValue<unsigned int>();
         break;
      }
      case fitsType<long>():
      {
         m_value.Long = convertedValue<long>();
         break;
      }
      case fitsType<unsigned long>():
      {
         m_value.ULong = convertedValue<unsigned long>();
         break;
      }
      case fitsType<long long>():
      {
         m_value.LongLong = convertedValue<long long>();
         break;
      }
      case fitsType<unsigned long long>():
      {
         m_value.ULongLong = convertedValue<unsigned long long>();
         break;
      }
      case fitsType<float>():
      {
         m_value.Float = convertedValue<float>();
         break;
      }
      case fitsType<std::complex<float>>():
      {
         mxThrowException(err::notimpl, "fitsHeaderCard::convertValue", std::string("conversions not supported for complex types in ") + m_keyword);
      }
      case fitsType<double>():
      {
         m_value.Double = convertedValue<double>();
         break;
      }
      case fitsType<std::complex<double>>():
      {
         mxThrowException(err::notimpl, "fitsHeaderCard::convertValue", std::string("conversions not supported for complex types in ") + m_keyword);
      }
      case fitsType<fitsCommentType>():
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertValue", "cannot convert comment to numeric type");
      }
      case fitsType<fitsHistoryType>():
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertValue", "cannot convert history to numeric type");
      }
      case TSTRING:
      {
         convertToString();
         m_type = newtype;
         m_valueGood = false;
         return;
      }
      default:
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::convertValue", std::string("invalid FITS type for conversion in ") + m_keyword);
      }
   }

   m_type = newtype;
   m_valueGood = true;
}

const std::string & fitsHeaderCard::keyword()
{
   return m_keyword;
}

void fitsHeaderCard::keyword(const std::string & kw)
{
   m_keyword = kw;
}

template<>
std::string fitsHeaderCard::value<std::string>()
{
   if(m_valueStrGood == false)
   {
      convertToString();
   }

   //Strip ' from beginning and end if present
   std::string str = m_valueStr.str();

   if(str[0] == '\'')
   {
      str.erase(0,1);
   }
   
   if(str[str.size()-1] == '\'')
   {
      str.erase(str.size()-1,1);
   }

   return str;
}

template<>
char fitsHeaderCard::value<char>()
{
   if(m_valueGood == false)
   {
      convertFromString<char>();
   }

   if(m_type != fitsType<char>())
   {
      return convertedValue<char>();
   }

   return m_value.Char;
}

template<>
unsigned char fitsHeaderCard::value<unsigned char>()
{
   if(m_valueGood == false)
   {
      convertFromString<unsigned char>();
   }

   if(m_type != fitsType<unsigned char>())
   {
      return convertedValue<unsigned char>();
   }

   return m_value.UChar;
}

template<>
short fitsHeaderCard::value<short>()
{
   if(m_valueGood == false)
   {
      convertFromString<short>();
   }

   if(m_type != fitsType<short>())
   {
      return convertedValue<short>();
   }

   return m_value.Short;
}

template<>
unsigned short fitsHeaderCard::value<unsigned short>()
{
   if(m_valueGood == false)
   {
      convertFromString<unsigned short>();
   }

   if(m_type != fitsType<unsigned short>())
   {
      return convertedValue<unsigned short>();
   }

   return m_value.UShort;
}

template<>
int fitsHeaderCard::value<int>()
{
   if(m_valueGood == false)
   {
      convertFromString<int>();
   }

   if(m_type != fitsType<int>())
   {
      return convertedValue<int>();
   }

   return m_value.Int;
}

template<>
unsigned int fitsHeaderCard::value<unsigned int>()
{
   if(m_valueGood == false)
   {
      convertFromString<unsigned int>();
   }

   if(m_type != fitsType<unsigned int>())
   {
      return convertedValue<unsigned int>();
   }

   return m_value.UInt;
}

template<>
long fitsHeaderCard::value<long>()
{
   if(m_valueGood == false)
   {
      convertFromString<long>();
   }

   if(m_type != fitsType<long>())
   {
      return convertedValue<long>();
   }

   return m_value.Long;
}

template<>
unsigned long fitsHeaderCard::value<unsigned long>()
{
   if(m_valueGood == false)
   {
      convertFromString<unsigned long>();
   }

   if(m_type != fitsType<unsigned long>())
   {
      return convertedValue<unsigned long>();
   }

   return m_value.ULong;
}

template<>
long long fitsHeaderCard::value<long long>()
{
   if(m_valueGood == false)
   {
      convertFromString<long long>();
   }

   if(m_type != fitsType<long long>())
   {
      return convertedValue<long long>();
   }

   return m_value.LongLong;
}

template<>
unsigned long long fitsHeaderCard::value<unsigned long long>()
{
   if(m_valueGood == false)
   {
      convertFromString<unsigned long long>();
   }

   if(m_type != fitsType<unsigned long long>())
   {
      return convertedValue<unsigned long long>();
   }

   return m_value.ULongLong;
}

template<>
float fitsHeaderCard::value<float>()
{
   if(m_valueGood == false)
   {
      convertFromString<float>();
   }

   if(m_type != fitsType<float>())
   {
      return convertedValue<float>();
   }

   return m_value.Float;
}

template<>
std::complex<float> fitsHeaderCard::value<std::complex<float>>()
{
   if(m_valueGood == false)
   {
      convertFromString<std::complex<float>>();
   }

   if(m_type != fitsType<std::complex<float>>())
   {
      return convertedValue<std::complex<float>>();
   }

   return m_value.complexFloat;
}

template<>
double fitsHeaderCard::value<double>()
{
   if(m_valueGood == false)
   {
      convertFromString<double>();
   }

   if(m_type != fitsType<double>())
   {
      return convertedValue<double>();
   }

   return m_value.Double;
}

template<>
std::complex<double> fitsHeaderCard::value<std::complex<double>>()
{
   if(m_valueGood == false)
   {
      convertFromString<std::complex<double>>();
   }

   if(m_type != fitsType<std::complex<double>>())
   {
      return convertedValue<std::complex<double>>();
   }

   return m_value.complexDouble;
}

std::string fitsHeaderCard::String()
{
   return value<std::string>();
}

char fitsHeaderCard::Char()
{
   return value<char>();
}

unsigned char fitsHeaderCard::UChar()
{
   return value<unsigned char>();
}

short fitsHeaderCard::Short()
{
   return value<short>();
}

unsigned short fitsHeaderCard::UShort()
{
   return value<unsigned short>();
}

int fitsHeaderCard::Int()
{
   return value<int>();
}

unsigned int fitsHeaderCard::UInt()
{
   return value<unsigned int>();
}

long fitsHeaderCard::Long()
{
   return value<long>();
}

unsigned long fitsHeaderCard::ULong()
{
   return value<unsigned long>();
}

long long fitsHeaderCard::LongLong()
{
   return value<long long>();
}

unsigned long long fitsHeaderCard::ULongLong()
{
   return value<unsigned long long>();
}

float fitsHeaderCard::Float()
{
   return value<float>();
}

std::complex<float> fitsHeaderCard::complexFloat()
{
   return value<std::complex<float>>();
}

double fitsHeaderCard::Double()
{
   return value<double>();
}

std::complex<double> fitsHeaderCard::complexDouble()
{
   return value<std::complex<double>>();
}

void fitsHeaderCard::value(const char * v)
{
   //Strip ' from beginning and end if present
   std::string str = v;

   if(str[0] == '\'')
   {
      str.erase(0,1);
   }
   
   if(str[str.size()-1] == '\'')
   {
      str.erase(str.size()-1,1);
   }

   m_valueStr.str(str);

   m_valueGood = false;
   m_valueStrGood = true;
   m_type = fitsType<char *>();
}

void fitsHeaderCard::value(const std::string & v)
{
   //Strip ' from beginning and end if present
   std::string str = v;

   if(str[0] == '\'')
   {
      str.erase(0,1);
   }
   
   if(str[str.size()-1] == '\'')
   {
      str.erase(str.size()-1,1);
   }

   m_valueStr.str(v);
   m_valueGood = false;
   m_valueStrGood = true;
   m_type = fitsType<std::string>();
}

void fitsHeaderCard::value(const char & v)
{
   m_value.Char = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<char>();
}

void fitsHeaderCard::value(const unsigned char & v)
{
   m_value.UChar = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<unsigned char>();
}

void fitsHeaderCard::value(const short int & v)
{
   m_value.Short = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<short>();
}

void fitsHeaderCard::value(const unsigned short & v)
{
   m_value.UShort = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<unsigned short>();
}

void fitsHeaderCard::value(const int & v)
{
   m_value.Int = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<int>();
}

void fitsHeaderCard::value(const unsigned int & v)
{
   m_value.UInt = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<unsigned int>();
}

void fitsHeaderCard::value(const long & v)
{
   m_value.Long = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<long>();
}

void fitsHeaderCard::value(const unsigned long int & v)
{
   m_value.ULong = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<unsigned long>();
}

void fitsHeaderCard::value(const long long & v)
{
   m_value.LongLong = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<long long>();
}

void fitsHeaderCard::value(const unsigned long long int & v)
{
   m_value.ULongLong = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<unsigned long long>();
}

void fitsHeaderCard::value(const float & v)
{
   m_value.Float = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<float>();
}


void fitsHeaderCard::value(const std::complex<float> & v)
{
   m_value.complexFloat = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<std::complex<float>>();
}

void fitsHeaderCard::value(const double & v)
{
   m_value.Double= v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<double>();
}

void fitsHeaderCard::value(const std::complex<double> & v)
{
   m_value.complexDouble = v;
   m_valueGood = true;
   m_valueStrGood = false;
   m_type = fitsType<std::complex<double>>();
}

int fitsHeaderCard::type()
{
   return m_type;
}

void fitsHeaderCard::type( const int & t)
{
   if(t == m_type) return;

   if(m_valueGood)
   {
      convertValue(t);
   }
   else m_type = t;

   //Need to reconvert, always favor the actual value.
   if(m_valueGood && m_valueStrGood)
   {
      m_valueStrGood = false;
   }
}

std::string fitsHeaderCard::valueStr()
{
   if(!m_valueStrGood)
   {
      convertToString();
   }

   std::string s = m_valueStr.str();
   return s;
}

bool fitsHeaderCard::valueGood()
{
   return m_valueGood;
}

bool fitsHeaderCard::valueStrGood()
{   
   return m_valueStrGood;
}

const std::string & fitsHeaderCard::comment()
{
   return m_comment;
}

void fitsHeaderCard::comment( const std::string & c)
{
   m_comment = c;
}

int fitsHeaderCard::write(fitsfile *fptr)
{
   if(m_type == fitsType<char *>() || m_type == fitsType<std::string>())
   {
      return fits_write_key<char *>(fptr, (char *) m_keyword.c_str(), (void *)m_valueStr.str().c_str(), (char *) m_comment.c_str());
   }

   //If the string is good, meaning already converted.
   if(m_valueStrGood == true)
   {
      //This populates the card directly.
      return fits_write_key<fitsUnknownType>(fptr, (char *) m_keyword.c_str(), (void *)m_valueStr.str().c_str(), (char *) m_comment.c_str());
   }

   //Ok, now we write the type directly using fitsio routines because it hasn't been converted.
   switch(m_type)
   {
      case fitsType<unsigned char>():
      {
         return fits_write_key<unsigned char>(fptr, (char *) m_keyword.c_str(), &m_value.UChar, (char *) m_comment.c_str());
      }
      case fitsType<char>():
      {
         return fits_write_key<char>(fptr, (char *) m_keyword.c_str(), &m_value.Char, (char *) m_comment.c_str());
      }
      case fitsType<short>():
      {
         return fits_write_key<short>(fptr, (char *) m_keyword.c_str(), &m_value.Short, (char *) m_comment.c_str());
      }
      case fitsType<unsigned short>():
      {
         return fits_write_key<unsigned short>(fptr, (char *) m_keyword.c_str(), &m_value.UShort, (char *) m_comment.c_str());
      }
      case fitsType<int>():
      {
         return fits_write_key<int>(fptr, (char *) m_keyword.c_str(), &m_value.Int, (char *) m_comment.c_str());
      }
      case fitsType<unsigned int>():
      {
         return fits_write_key<unsigned int>(fptr, (char *) m_keyword.c_str(), &m_value.UInt, (char *) m_comment.c_str());
      }
      case fitsType<long>():
      {
         return fits_write_key<long>(fptr, (char *) m_keyword.c_str(), &m_value.Long, (char *) m_comment.c_str());
      }
      case fitsType<unsigned long>():
      {
         return fits_write_key<unsigned long>(fptr, (char *) m_keyword.c_str(), &m_value.ULong, (char *) m_comment.c_str());
      }
      case fitsType<long long>():
      {
         return fits_write_key<long long>(fptr, (char *) m_keyword.c_str(), &m_value.LongLong, (char *) m_comment.c_str());
      }
      case fitsType<unsigned long long>():
      {
         return fits_write_key<unsigned long long>(fptr, (char *) m_keyword.c_str(), &m_value.ULongLong, (char *) m_comment.c_str());
      }
      case fitsType<float>():
      {
         return fits_write_key<float>(fptr, (char *) m_keyword.c_str(), &m_value.Float, (char *) m_comment.c_str());
      }
      case fitsType<double>():
      {
         return fits_write_key<double>(fptr, (char *) m_keyword.c_str(), &m_value.Double, (char *) m_comment.c_str());
      }
      case fitsType<fitsCommentType>():
      {
         return fits_write_comment(fptr, (char *) m_comment.c_str());
      }
      case fitsType<fitsHistoryType>():
      {
         return fits_write_history(fptr, (char *) m_comment.c_str());
      }
      default:
      {
         mxThrowException(err::invalidarg, "fitsHeaderCard::write", std::string("invalid FITS type for ") + m_keyword);
      }
   }
}

} //namespace fits
} //namespace mx


