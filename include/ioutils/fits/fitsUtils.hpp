/** \file fitsUtils.hpp
  * \brief Declares and defines utilities to work with FITS files
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

#ifndef ioutils_fits_fitsUtils_hpp
#define ioutils_fits_fitsUtils_hpp


#include <iostream>
#include <cstdlib>
#include <cstring>

#include <fitsio.h>

#include "../stringUtils.hpp"

namespace mx
{
namespace fits
{


///The standard width of the value entry in a header card
#define stdValWidth (20)

#define fitsTUNKNOWN (-5000)
struct fitsUnknownType
{
};

#define fitsTCOMMENT (-5001)
struct fitsCommentType
{
   fitsCommentType()
   {
   }
   
   explicit fitsCommentType(char * v)
   {
      static_cast<void>(v);
   }
   
   explicit fitsCommentType(const char * v)
   {
      static_cast<void>(v);
   }
};

#define fitsTHISTORY (-5002)

struct fitsHistoryType
{
   fitsHistoryType()
   {
   }

   explicit fitsHistoryType(char * v)
   {
      static_cast<void>(v);
   }
   
   explicit fitsHistoryType(const char * v)
   {
      static_cast<void>(v);
   }
   
};


/** \ingroup fits_utils
  * @{
  */

/// Return the cfitsio constant for a given data type.
/** 
  *
  * \tparam scalarT is the type 
  * 
  * \returns a constant defined in cfitsio corresponding to the native type
  * \returns -1 if not a define type in cfitsio
  * 
  * \ingroup fits_utils
  */
template<typename scalarT> 
int getFitsType()
{
   return -1;
}

template<> 
int getFitsType<char *>();

template<> 
int getFitsType<std::string>();

template<> 
int getFitsType<unsigned char>();

template<> 
int getFitsType<signed char>();

template<> 
int getFitsType<char>();

template<> 
int getFitsType<short>();

template<> 
int getFitsType<unsigned short>();

template<> 
int getFitsType<int>();

template<> 
int getFitsType<unsigned int>();

template<> 
int getFitsType<long>();

template<> 
int getFitsType<long long>();

template<> 
int getFitsType<unsigned long>();

template<> 
int getFitsType<float>();

template<> 
int getFitsType<double>();

template<> 
int getFitsType<fitsUnknownType>();

template<> 
int getFitsType<fitsCommentType>();

template<> 
int getFitsType<fitsHistoryType>();

/** Return the cfitsio BITPIX value for a given data type.
  *
  * \tparam scalarT is the type 
  * \retval int > 0 if a constant is defined in cfitsio corresponding to the native type
  * \retval -1 if not a defined type in cfitsio
  */
template<typename scalarT> 
int getFitsBITPIX();

template<> 
int getFitsBITPIX<char>();

template<> 
int getFitsBITPIX<signed char>();

template<> 
int getFitsBITPIX<unsigned char>();

template<> 
int getFitsBITPIX<short>();

template<> 
int getFitsBITPIX<unsigned short>();

template<> 
int getFitsBITPIX<int>();

template<> 
int getFitsBITPIX<unsigned int>();

template<> 
int getFitsBITPIX<long>();

template<> 
int getFitsBITPIX<float>();

template<> 
int getFitsBITPIX<double>();

/// Strip the apostrophes from a FITS value string
/** The argument is modified if the first and/or last non-whitespace character is '
  *
  * \param s is the string from which to strip apostrophes
  * 
  * \retval int containing the number of stripped apostrophes
  */ 
int fitsStripApost(std::string & s);

///Populate a fits header card with the value string copied verbatim
/** 
  * \param headStr is a c-string which must be 81 characters in length, including the '\n'
  * \param keyword is the keyword name 
  * \param value is the value string
  * \param comment is the comment string
  */ 
void fitsPopulateCard(char headStr[81], char *keyword, char *value, char *comment);

/// Write a header card to a file
/** This is a templatized wrapper for the cfitsio routine fits_write_key.
  *
  * \tparam typeT is the type of the value
  *
  * \param fptr is a pointer to an open fits file
  * \param keyword is a c-string containing the keyword
  * \param value is a pointer to the memory location of the value
  * \param comment is a c-string, possibly NULL, containing a comment string
  *  
  * \retval int containing the status returned by the cfitsio routine.
  * 
  * \ingroup fits_utils
  */
template<typename typeT> 
int fits_write_key( fitsfile * fptr, 
                    char * keyword, 
                    void * value, 
                    char * comment
                  )
{
   int fstatus = 0;
   
   fits_write_key(fptr, getFitsType<typeT>(), keyword, value, comment,  &fstatus);
   
   return fstatus;
}

template<>
int fits_write_key<char *>( fitsfile * fptr, 
                            char * keyword, 
                            void * value, 
                            char * comment
                          );

template<>
int fits_write_key<std::string>( fitsfile * fptr, 
                                 char * keyword, 
                                 void * value, 
                                 char * comment
                               );

template<> 
int fits_write_key<fitsUnknownType>( fitsfile * fptr, 
                                     char * keyword, 
                                     void * value, 
                                     char * comment
                                   );

int fits_write_comment( fitsfile *fptr, 
                        char * comment
                      );

int fits_write_history( fitsfile *fptr, 
                        char * history
                      );

/// Generate a rich error meesage from a FITS status code.
void fitsErrText( std::string & explan,         ///< [out] the explanatory message
                  const std::string & filename, ///< [in] the FITS file's name which generated the problem
                  int fstatus                   ///< [in] the cfitstio status code
                );
///@}

} //namespace fits
} //namespace mx

#endif //ioutils_fits_fitsUtils_hpp

