/** \file fitsUtils.cpp
  * \brief Defines utilities to work with FITS files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#include "fitsUtils.hpp"

namespace mx
{
   
template<> int getFitsType<unsigned char>()
{
   return TBYTE;
}

template<> int getFitsType<signed char>()
{
   return TSBYTE;
}

template<> int getFitsType<short>()
{
   return TSHORT;
}

template<> int getFitsType<unsigned short>()
{
   return TUSHORT;
}

template<> int getFitsType<int>()
{
   return TINT;
}

template<> int getFitsType<long>()
{
   return TLONG;
}

template<> int getFitsType<long long>()
{
   return TLONGLONG;
}

template<> int getFitsType<unsigned long>()
{
   return TULONG;
}

template<> int getFitsType<float>()
{
   return TFLOAT;
}

template<> int getFitsType<double>()
{
   return TDOUBLE;
}

template<> int getFitsBITPIX<char>()
{
   return SBYTE_IMG;
}

template<> int getFitsBITPIX<unsigned char>()
{
   return BYTE_IMG;
}

template<> int getFitsBITPIX<short>()
{
   return SHORT_IMG;
}

template<> int getFitsBITPIX<unsigned short>()
{
   return USHORT_IMG;
}

template<> int getFitsBITPIX<int>()
{
   return LONG_IMG; //Yes, this is right.  This returns 32
}

template<> int getFitsBITPIX<unsigned int>()
{
   return ULONG_IMG; //Yes, this is right, this returns 40
}

template<> int getFitsBITPIX<long>()
{
   return LONGLONG_IMG; //Yes, this is right, this returns 64
}

template<> int getFitsBITPIX<float>()
{
   return FLOAT_IMG;
}

template<> int getFitsBITPIX<double>()
{
   return DOUBLE_IMG;
}

}

