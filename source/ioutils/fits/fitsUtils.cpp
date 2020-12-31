/** \file fitsUtils.cpp
  * \brief Implementation of utilities to work with FITS files
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

#include "ioutils/fits/fitsUtils.hpp"

namespace mx
{
namespace fits 
{
   
template<> 
int getFitsType<char *>()
{
   return TSTRING;
}

template<> 
int getFitsType<std::string>()
{
   return TSTRING;
}


template<> 
int getFitsType<unsigned char>()
{
   return TBYTE;
}

template<> 
int getFitsType<signed char>()
{
   return TSBYTE;
}

template<> 
int getFitsType<char>()
{
   return TSBYTE;
}

template<> 
int getFitsType<short>()
{
   return TSHORT;
}

template<> 
int getFitsType<unsigned short>()
{
   return TUSHORT;
}

template<> 
int getFitsType<int>()
{
   return TINT;
}

template<> 
int getFitsType<unsigned int>()
{
   return TUINT;
}

template<> 
int getFitsType<long>()
{
   return TLONG;
}

template<> 
int getFitsType<long long>()
{
   return TLONGLONG;
}

template<> 
int getFitsType<unsigned long>()
{
   return TULONG;
}

template<> 
int getFitsType<float>()
{
   return TFLOAT;
}

template<> 
int getFitsType<double>()
{
   return TDOUBLE;
}


template<> 
int getFitsType<fitsUnknownType>()
{
   return fitsTUNKNOWN;
}

template<> 
int getFitsType<fitsCommentType>()
{
   return fitsTCOMMENT;
}

template<> 
int getFitsType<fitsHistoryType>()
{
   return fitsTHISTORY;
}

template<> 
int getFitsBITPIX<char>()
{
   return SBYTE_IMG;
}

template<> 
int getFitsBITPIX<signed char>()
{
   return SBYTE_IMG;
}

template<> 
int getFitsBITPIX<unsigned char>()
{
   return BYTE_IMG;
}

template<> 
int getFitsBITPIX<short>()
{
   return SHORT_IMG;
}

template<> 
int getFitsBITPIX<unsigned short>()
{
   return USHORT_IMG;
}

template<> 
int getFitsBITPIX<int>()
{
   return LONG_IMG; //Yes, this is right.  This returns 32
}

template<> 
int getFitsBITPIX<unsigned int>()
{
   return ULONG_IMG; //Yes, this is right, this returns 40
}

template<> 
int getFitsBITPIX<long>()
{
   return LONGLONG_IMG; //Yes, this is right, this returns 64
}

template<> 
int getFitsBITPIX<float>()
{
   return FLOAT_IMG;
}

template<> 
int getFitsBITPIX<double>()
{
   return DOUBLE_IMG;
}

int fitsStripApost(std::string & s)
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

void fitsPopulateCard( char headStr[81], 
                       char *keyword, 
                       char *value, 
                       char *comment
                     )
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

template<>
int fits_write_key<char *>(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   int fstatus = 0;
  
   fits_write_key_longwarn(fptr, &fstatus);
   
      
   fits_write_key_longstr(fptr, keyword, (const char *)  value, comment,  &fstatus);
   
   return fstatus;
   
}

template<>
int fits_write_key<std::string>(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   return fits_write_key<char *>(fptr, keyword, value, comment);
}


template<> 
int fits_write_key<fitsUnknownType>(fitsfile * fptr, char * keyword, void * value, char * comment)
{
   
   
   int fstatus = 0;

   char *vstr = (char *) value;
   
   //If type is unkown, do it verbatim
   char headStr[81];
   
   fitsPopulateCard(headStr, keyword, vstr, comment);
      
   fits_write_record(fptr, headStr, &fstatus);
   return fstatus;
}

int fits_write_comment(fitsfile *fptr, char * comment)
{
   int fstatus = 0;
   
   fits_write_comment(fptr, comment, &fstatus);
   
   return fstatus;
}

int fits_write_history(fitsfile *fptr, char * history)
{
   int fstatus = 0;
   
   fits_write_history(fptr, history, &fstatus);
   
   return fstatus;
}

void fitsErrText( std::string & explan, 
                  const std::string & filename, 
                  int fstatus 
                )
{
   char emnem[31];

   fits_get_errstatus(fstatus, emnem);

   explan += ": ";
   explan += filename;
   explan += ". CFITSIO: ";

   explan += emnem;
   explan += " (";
   explan += ioutils::convertToString(fstatus);
   explan += ")";

}

} //namespace fits
} //namespace mx

