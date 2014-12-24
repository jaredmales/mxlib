/** \file fitsHeaderCard.cpp
  * \brief Definitions for a class to work with FITS header cards
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#include "fitsHeaderCard.hpp"

#include <iostream>

namespace mx
{

int fitsHeaderCard::Int()
{
   return convertFromString<int>(value);
}

unsigned int fitsHeaderCard::unsignedInt()
{
   return convertFromString<unsigned int>(value);
}

long fitsHeaderCard::Long()
{
   return convertFromString<long>(value);
}

unsigned long fitsHeaderCard::unsignedLong()
{
   return convertFromString<unsigned long>(value);
}

long long fitsHeaderCard::longLong()
{
   return convertFromString<long long>(value);
}
   
unsigned long long fitsHeaderCard::unsignedLongLong()
{
   return convertFromString<unsigned long long>(value);
}

float fitsHeaderCard::Float()
{
   return convertFromString<float>(value);
}

double fitsHeaderCard::Double()
{
   return convertFromString<double>(value);
}   

long double fitsHeaderCard::longDouble()
{
   return convertFromString<long double>(value);
}



int fitsHeaderCard::write(fitsfile *fptr)
{
   int fstatus;
   switch(type)
   {
      case TSTRING:
      {
         return fits_update_key<char *>(fptr, (char *) keyword.c_str(), (void *)value.c_str(), (char *) comment.c_str());
      }
      case TBYTE:
      {
         unsigned char b = Int();
         return fits_update_key<unsigned char>(fptr, (char *) keyword.c_str(), &b, (char *) comment.c_str());
      }  
      case TSBYTE:
      {
         signed char sc = Int();
         return fits_update_key<signed char>(fptr, (char *) keyword.c_str(), &sc, (char *) comment.c_str());   
      }
      case TSHORT:
      {
         short s = Int();
         return fits_update_key<short>(fptr, (char *) keyword.c_str(), &s, (char *) comment.c_str());   
      }   
      case TUSHORT:
      {
         unsigned short us = Int();
         return fits_update_key<unsigned short>(fptr, (char *) keyword.c_str(), &us, (char *) comment.c_str());
      }   
      case TINT:
      {
         int i = Int();
         return fits_update_key<int>(fptr, (char *) keyword.c_str(), &i, (char *) comment.c_str());
      }   
      case TLONG:
      {
         long l = Long();
         return fits_update_key<long>(fptr, (char *) keyword.c_str(), &l, (char *) comment.c_str());
      }
      case TLONGLONG:
      {
         long long ll = longLong();
         return fits_update_key<long long>(fptr, (char *) keyword.c_str(), &ll, (char *) comment.c_str());
      }
      case TULONG:
      {
         unsigned long ul = unsignedLong();
         return fits_update_key<unsigned long>(fptr, (char *) keyword.c_str(), &ul, (char *) comment.c_str());
      }   
      case TFLOAT:
      {
         float f = Float();
         return fits_update_key<float>(fptr, (char *) keyword.c_str(), &f, (char *) comment.c_str());   
      }   
      case TDOUBLE:
      {
         double d = Double();
         return fits_update_key<double>(fptr, (char *) keyword.c_str(), &d, (char *) comment.c_str());   
      }
      default:
      {
         return fits_update_key<char *>(fptr, (char *) keyword.c_str(), (void *)value.c_str(), (char *) comment.c_str());
      }
   }
}
         
} //namespace mx
