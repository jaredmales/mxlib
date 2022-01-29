/** \file mxException.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines the mxlib exception class.
  * \ingroup error_handling_files
  * 
*/

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __mxException__
#define __mxException__

#include <exception>
#include <sstream>

#include "mxError.hpp"


namespace mx
{
namespace err
{

/// The mxlib exception class
/** Provides a rich error report via the standard what().
  * \ingroup exceptions
  */
class mxException : public std::exception
{

protected:
   ///Contains the what() string
   char m_whatstr[4096];

   ///The source of the exception, such as stdlib or cfitisio or the function name
   std::string m_source {""};
   
   ///The mxlib error code
   int  m_code {0};
   
   ///The name of the error code
   std::string m_codeName {""};
   
   ///The errno error code (only used if non-zero)
   int m_errno {0};

   ///The source file where the exception originated
   std::string m_file {""};
   
   ///The line number of the file where the exception originated
   int  m_line {0};
   
   ///The long explanation of the error
   std::string m_explanation {""};

public:

   ///Default constructor
   mxException () noexcept
   {      
      build_what();
   }

   ///Copy constructor
   mxException (const mxException & e) noexcept : m_source(e.m_source), m_code(e.m_code), m_codeName(e.m_codeName), m_file(e.m_file), m_line(e.m_line), m_explanation(e.m_explanation)
   {
      build_what();
   }

   ///Construct and fill in each of the values, except errno
   mxException( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
                const int & ec,            ///< [in] the error code
                const std::string & emnem, ///< [in] the name of the error code
                const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
                const int & line,          ///< [in] the line number where the exception was thrown
                const std::string & expl   ///< [in] the explanation for why the exception was thrown
              ) : m_source(esrc), m_code(ec), m_codeName(emnem), m_errno(0), m_file(efile), m_line(line), m_explanation(expl)
   {
      build_what();
   }

   ///Construct and fill in each of the values, including errno
   mxException( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
                const int & ec,            ///< [in] the error code
                const std::string & emnem, ///< [in] the name of the error code
                const int & en,            ///< [in] the value of errno
                const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
                const int & line,          ///< [in] the line number where the exception was thrown
                const std::string & expl   ///< [in] the explanation for why the exception was thrown
              ) : m_source(esrc), m_code(ec), m_codeName(emnem), m_errno(en), m_file(efile), m_line(line), m_explanation(expl)
   {
      build_what();
   }

   ///Assignment operator
   mxException & operator=(const mxException & e) noexcept
   {
      m_source = e.m_source;
      m_code = e.m_code;
      m_codeName = e.m_codeName;
      m_errno = e.m_errno;
      m_file = e.m_file;
      m_line = e.m_line;
      m_explanation = e.m_explanation;
      
      build_what();
      
      return *this;
   }

   ///Destructor
   virtual ~mxException() throw()
   {
   }

   ///Build the what string.
   /** Must be called after updating any values, since the what() method is const const.
     */
   virtual void build_what()
   {
      std::ostringstream s;
      s.str("");
      
      s << "An exception has been thrown in an mxlib component.\n";
      s << "      source: " << m_source << "\n";
      if(m_code != 0)
      {
         s << "        code: " << m_code;
         if(m_codeName != "") s << " ("<< m_codeName << ")";
         s << "\n";
      }
      if(m_errno != 0)
      {
         s << "       errno: " << m_errno << " (" << errno_CodeToName(m_errno) << ")\n";
      }
      s << "     in file: " << m_file << "\n";
      s << "     at line: " << m_line << "\n";
      if(m_explanation != "") s << " explanation: " << m_explanation << "\n";
      if(m_errno != 0)
      {
         s << "              " << strerror(m_errno) << "\n";
      }
      snprintf(m_whatstr, sizeof(m_whatstr), "%s", s.str().c_str());
   }
   
   ///Return the details of the exception as a single string.
   virtual const char * what() const noexcept
   {   
      return m_whatstr;
   }
   
};

/// mxException for invalid arguments
/** 
  * \ingroup exceptions
  */ 
class invalidarg : public mxException 
{
public:
   invalidarg( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
               const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
               const int & line,          ///< [in] the line number where the exception was thrown
               const std::string & expl   ///< [in] the explanation for why the exception was thrown
             ) : mxException(esrc, MXE_INVALIDARG,MXE_INVALIDARG_NAME, efile, line, expl)
   {
   }
};

/// mxException for invalid config settings
/** 
  * \ingroup exceptions
  */ 
class invalidconfig : public mxException 
{
public:
   invalidconfig( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
                  const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
                  const int & line,          ///< [in] the line number where the exception was thrown
                  const std::string & expl   ///< [in] the explanation for why the exception was thrown
                ) : mxException(esrc, MXE_INVALIDARG,MXE_INVALIDARG_NAME, efile, line, expl)
   {
   }
};

/// mxException for not implemented features
/** 
  * \ingroup exceptions
  */
class notimpl : public mxException 
{
public:
   notimpl( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
            const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
            const int & line,          ///< [in] the line number where the exception was thrown 
            const std::string & expl   ///< [in] the explanation for why the exception was thrown
          ) : mxException(esrc, MXE_NOTIMPL,MXE_NOTIMPL_NAME, efile, line, expl)
   {
   }
};

/// mxException for parameters which aren't set
/** 
  * \ingroup exceptions
  */
class paramnotset : public mxException 
{
public:
   paramnotset( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
                const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
                const int & line,          ///< [in] the line number where the exception was thrown 
                const std::string & expl   ///< [in] the explanation for why the exception was thrown
              ) : mxException(esrc, MXE_PARAMNOTSET,MXE_PARAMNOTSET_NAME, efile, line, expl)
   {
   }
};

/// mxException for a size error
/** 
  * \ingroup exceptions
  */
class sizeerr : public mxException 
{
public:
   sizeerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
            const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
            const int & line,          ///< [in] the line number where the exception was thrown 
            const std::string & expl   ///< [in] the explanation for why the exception was thrown
          ) : mxException(esrc, MXE_SIZEERR,MXE_SIZEERR_NAME, efile, line, expl)
   {
   }
};

/// mxException for an allocation error
/** 
  * \ingroup exceptions
  */
class allocerr : public mxException 
{
public:
   allocerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown 
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_ALLOCERR,MXE_ALLOCERR_NAME, efile, line, expl)
   {
   }
};

/// mxException for errors on opening a file
/** 
  * \ingroup exceptions
  */
class fileoerr : public mxException 
{
public:
   fileoerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILEOERR,MXE_FILEOERR_NAME, efile, line, expl)
   {
   }
   fileoerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const int & en,            ///< [in] the value of errno
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILEOERR,MXE_FILEOERR_NAME, en, efile, line, expl)
   {
   }
};

/// mxException for errors reading from a file
/** 
  * \ingroup exceptions
  */
class filererr : public mxException 
{
public:
   filererr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILERERR,MXE_FILERERR_NAME, efile, line, expl)
   {
   }
   filererr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const int & en,            ///< [in] the value of errno
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILERERR,MXE_FILERERR_NAME, en, efile, line, expl)
   {
   }
};

/// mxException for errors writing to a file
/** 
  * \ingroup exceptions
  */
class filewerr : public mxException 
{
public:
   filewerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILEWERR,MXE_FILEWERR_NAME, efile, line, expl)
   {
   }
   filewerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const int & en,            ///< [in] the value of errno
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILEWERR,MXE_FILEWERR_NAME, en, efile, line, expl)
   {
   }
};

/// mxException for errors closing a file
/** 
  * \ingroup exceptions
  */
class filecerr : public mxException 
{
public:
   filecerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILECERR,MXE_FILECERR_NAME, efile, line, expl)
   {
   }
   filecerr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
             const int & en,            ///< [in] the value of errno
             const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
             const int & line,          ///< [in] the line number where the exception was thrown
             const std::string & expl   ///< [in] the explanation for why the exception was thrown
           ) : mxException(esrc, MXE_FILECERR,MXE_FILECERR_NAME, en, efile, line, expl)
   {
   }
};

/// mxException for errors returned by a library call
/** 
  * \ingroup exceptions
  */
class liberr : public mxException 
{
public:
   liberr( const std::string & esrc,  ///< [in] the source of the exception, typically the class and function
           const std::string & efile, ///< [in] the source file in which the exception occurred, normally __FILE__
           const int & line,          ///< [in] the line number where the exception was thrown
           const std::string & expl   ///< [in] the explanation for why the exception was thrown
         ) : mxException(esrc, MXE_LIBERR,MXE_LIBERR_NAME, efile, line, expl)
   {
   }
};

/// Throw an exception.  This macro takes care of the file and line.
/** \ingroup exceptions
  */
#define mxThrowException( extype, src, expl ) throw extype(src,  __FILE__, __LINE__, expl);

#define mxThrowExceptionErrno( extype, en, src, expl ) throw extype(src,  en, __FILE__, __LINE__, expl);

} //err
} //namespace mx

#endif //mxException_hpp
