/** \file fitsUtils
  * \brief Declares and defines utilities to work with FITS files
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __fitsUtils__
#define __fitsUtils__


#include "fitsio.h"

namespace mx
{

/** \addtogroup image_processing
  * @{
  */

/** Return the cfitsio constant for a given data type.
  *
  * \tparam scalarT is the type 
  * \returns a constant define in cfitsio corresponding to the native type
  * \retval -1 if not a define type in cfitsio
  */
template<typename scalarT> int getFitsType()
{
   return -1;
}

template<> int getFitsType<unsigned char>();

template<> int getFitsType<signed char>();

template<> int getFitsType<short>();

template<> int getFitsType<unsigned short>();

template<> int getFitsType<int>();

template<> int getFitsType<long>();

template<> int getFitsType<long long>();

template<> int getFitsType<unsigned long>();

template<> int getFitsType<float>();

template<> int getFitsType<double>();

/** Return the cfitsio BITPIX value for a given data type.
  *
  * \tparam scalarT is the type 
  * \returns a constant defined in cfitsio corresponding to the native type
  * \retval -1 if not a defined type in cfitsio
  */
template<typename scalarT> int getFitsBITPIX()
{
   return -1;
}

template<> int getFitsBITPIX<char>();

template<> int getFitsBITPIX<unsigned char>();

template<> int getFitsBITPIX<short>();

template<> int getFitsBITPIX<unsigned short>();

template<> int getFitsBITPIX<int>();

template<> int getFitsBITPIX<unsigned int>();

template<> int getFitsBITPIX<long>();

template<> int getFitsBITPIX<float>();

template<> int getFitsBITPIX<double>();

///@}

} //namespace mx

#endif //__fitsUtils__

