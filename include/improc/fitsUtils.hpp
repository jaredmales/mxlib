/** \file fitsUtils.hpp
  * \brief Declares and defines utilities to work with FITS files
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __fitsUtils__
#define __fitsUtils__


#include <iostream>
#include <cstdlib>
#include <cstring>

#include "fitsio.h"


namespace mx
{
namespace improc
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
   
   explicit fitsCommentType(char *v)
   {
   }
   
   explicit fitsCommentType(const char *v)
   {
   }
};

#define fitsTHISTORY (-5002)
struct fitsHistoryType
{
   fitsHistoryType()
   {
   }

   explicit fitsHistoryType(char *v)
   {
   }
   
   explicit fitsHistoryType(const char *v)
   {
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
  * \retval int contains a constant defined in cfitsio corresponding to the native type
  * \retval -1 if not a define type in cfitsio
  * 
  * \ingroup fits_utils
  */
template<typename scalarT> 
int getFitsType()
{
   return -1;
}

template<> 
inline int getFitsType<char *>()
{
   return TSTRING;
}

template<> 
inline int getFitsType<std::string>()
{
   return TSTRING;
}


template<> 
inline int getFitsType<unsigned char>()
{
   return TBYTE;
}

template<> 
inline int getFitsType<signed char>()
{
   return TSBYTE;
}

template<> 
inline int getFitsType<short>()
{
   return TSHORT;
}

template<> 
inline int getFitsType<unsigned short>()
{
   return TUSHORT;
}

template<> 
inline int getFitsType<int>()
{
   return TINT;
}

template<> 
inline int getFitsType<unsigned int>()
{
   return TUINT;
}

template<> 
inline int getFitsType<long>()
{
   return TLONG;
}

template<> 
inline int getFitsType<long long>()
{
   return TLONGLONG;
}

template<> 
inline int getFitsType<unsigned long>()
{
   return TULONG;
}

template<> 
inline int getFitsType<float>()
{
   return TFLOAT;
}

template<> 
inline int getFitsType<double>()
{
   return TDOUBLE;
}

template<> 
inline int getFitsType<std::complex<double>>()
{
   return TDBLCOMPLEX;
}

template<> 
inline int getFitsType<fitsUnknownType>()
{
   return fitsTUNKNOWN;
}

template<> 
inline int getFitsType<fitsCommentType>()
{
   return fitsTCOMMENT;
}

template<> 
inline int getFitsType<fitsHistoryType>()
{
   return fitsTHISTORY;
}


/** Return the cfitsio BITPIX value for a given data type.
  *
  * \tparam scalarT is the type 
  * \retval int > 0 if a constant is defined in cfitsio corresponding to the native type
  * \retval -1 if not a defined type in cfitsio
  */
template<typename scalarT> int getFitsBITPIX()
{
   return -1;
}

template<> 
inline int getFitsBITPIX<char>()
{
   return SBYTE_IMG;
}

template<> 
inline int getFitsBITPIX<unsigned char>()
{
   return BYTE_IMG;
}

template<> 
inline int getFitsBITPIX<short>()
{
   return SHORT_IMG;
}

template<> 
inline int getFitsBITPIX<unsigned short>()
{
   return USHORT_IMG;
}

template<> 
inline int getFitsBITPIX<int>()
{
   return LONG_IMG; //Yes, this is right.  This returns 32
}

template<> 
inline int getFitsBITPIX<unsigned int>()
{
   return ULONG_IMG; //Yes, this is right, this returns 40
}

template<> 
inline int getFitsBITPIX<long>()
{
   return LONGLONG_IMG; //Yes, this is right, this returns 64
}

template<> 
inline int getFitsBITPIX<float>()
{
   return FLOAT_IMG;
}

template<> 
inline int getFitsBITPIX<double>()
{
   return DOUBLE_IMG;
}

/// Strip the apostrophes from a FITS value string
/** The argument is modified if the first and/or last non-whitespace character is '
  *
  * \param s is the string from which to strip apostrophes
  * 
  * \retval int containing the number of stripped apostrophes
  */ 
inline int fitsStripApost(std::string & s)
{  
   int stripped = 0;
   int p = s.find_first_not_of(" \t\n\t");
   if(s[p] == '\'')  
   {
      s.erase(0,p+1);
      ++stripped;
   }
   
   p = s.find_last_not_of(" \t\n\t");
   if(s[p] == '\'')
   {
      s.erase(p);
      ++stripped;
   }
   
   --p;

   while(s[p] == ' ' && p >=0)
   {
      s.erase(p);
      --p;
      ++stripped;
   }
   
   return stripped;
}

///Populate a fits header card with the value string copied verbatim
/** 
  * \param headStr is a c-string which must be 81 characters in length, including the '\n'
  * \param keyword is the keyword name 
  * \param value is the value string
  * \param comment is the comment string
  */ 
inline void fitsPopulateCard(char headStr[81], char *keyword, char *value, char *comment)
{
   memset(headStr, ' ', 80);
   headStr[80] = '\0';
   
   int rpos = 0;
   
   snprintf(headStr, 9, "%-8s", keyword);
   headStr[8] = '=';
   rpos = 10;
   
   if(strlen(value) < stdValWidth)
   {
      char fstr[10];
      snprintf(fstr, 10, "%%+%ds", stdValWidth);
      snprintf(headStr + rpos, 81-rpos, fstr, value);
      
   }
   else
   {
      snprintf(headStr + rpos, 81-rpos, "%s", value);
   }
   
   rpos = strlen(headStr);
   
   headStr[rpos] = ' ';
   ++rpos;
   
   headStr[rpos] = '/';
   ++rpos;
   
   headStr[rpos] = ' ';
   ++rpos;
   
   snprintf(headStr + rpos, 81-rpos, "%s", comment);
}

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
int fits_write_key(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   int fstatus = 0;
   
   fits_write_key(fptr, getFitsType<typeT>(), keyword, value, comment,  &fstatus);
   
   return fstatus;
}

template<>
inline int fits_write_key<char *>(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   int fstatus = 0;
  
   fits_write_key_longwarn(fptr, &fstatus);
   
      
   fits_write_key_longstr(fptr, keyword, (const char *)  value, comment,  &fstatus);
   
   return fstatus;
   
}

template<>
inline int fits_write_key<std::string>(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   return fits_write_key<char *>(fptr, keyword, value, comment);
}


template<> 
inline int fits_write_key<fitsUnknownType>(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   
   
   int fstatus = 0;

   char *vstr = (char *) value;
   
   //If type is unkown, do it verbatim
   char headStr[81];
   
   fitsPopulateCard(headStr, keyword, vstr, comment);
      
   fits_write_record(fptr, headStr, &fstatus);
   return fstatus;
}





inline int fits_write_comment(fitsfile *fptr, char * comment)
{
   int fstatus = 0;
   
   fits_write_comment(fptr, comment, &fstatus);
   
   return fstatus;
}

inline int fits_write_history(fitsfile *fptr, char * history)
{
   int fstatus = 0;
   
   fits_write_history(fptr, history, &fstatus);
   
   return fstatus;
}




///@}

} //namespace improc
} //namespace mx

#endif //__fitsUtils__

