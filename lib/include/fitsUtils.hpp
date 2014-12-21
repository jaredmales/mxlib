/** \file fitsUtils
  * \brief Declares and defines utilities to work with FITS files
  * \ingroup image_processing
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __fitsUtils__
#define __fitsUtils__


#include "fitsio.h"
#include <iostream>

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

template<> int getFitsType<char *>();

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


/** Update a header keyword. 
  * Templatized wrapper for the cfitsio routine.
  *
  * \tparam typeT is the type of the value
  *
  * \param fptr is a pointer to an open fits file
  * \param keyword is a c-string containing the keyword
  * \param value is a pointer to the memory location of the value
  * \param comment is a c-string, possibly NULL, containing a comment string
  *  
  * \returns the status returned by the cfitsio routine.
  */
template<typename typeT> int fits_update_key(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   int fstatus = 0;
 
   //std::cout << "writing " << keyword << " " << comment << "\n";
   
   fits_update_key(fptr, getFitsType<typeT>(), keyword, value, comment,  &fstatus);
   
   //std::cout << fstatus << "\n";
   return fstatus;
}

///@}

} //namespace mx

#endif //__fitsUtils__

