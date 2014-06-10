/** \file fitsHeaderCart.cpp
  * \brief Definitoins for a class to work with FITS header cards
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#include "fitsHeaderCard.hpp"

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

} //namespace mx
