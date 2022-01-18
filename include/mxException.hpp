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
  * \ingroup error_handling
  */
class mxException : public std::exception
{

protected:
   ///Contains the what() string
   char m_whatstr[4096];

   ///The source of the exception, such as stdlib or cfitisio or the function name
   std::string m_source {""};
   
   ///The source specific error code
   int  m_code {0};
   
   ///The mnemonic associated with the error code
   std::string m_codeName {""};
   
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

   ///Construct and fill in each of the values.
   mxException( const std::string & esrc, 
                const int & ec, 
                const std::string & emnem, 
                const std::string & efile, 
                const int & line, 
                const std::string & expl
              ) : m_source(esrc), m_code(ec), m_codeName(emnem), m_file(efile), m_line(line), m_explanation(expl)
   {
      build_what();
   }

   ///Assignment operator
   mxException & operator=(const mxException & e) noexcept
   {
      m_source = e.m_source;
      m_code = e.m_code;
      m_codeName = e.m_codeName;
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
      s << "        code: " << m_code;
      if(m_codeName != "") s << " ("<< m_codeName << ")";
      s << "\n";
      s << "     in file: " << m_file << "\n";
      s << "     at line: " << m_line << "\n";
      if(m_explanation != "") s << " explanation: " << m_explanation << "\n";
      
      snprintf(m_whatstr, 4096, "%s", s.str().c_str());
   }
   
   ///Return the details of the exception as a single string.
   virtual const char * what() const noexcept
   {   
      return m_whatstr;
   }
   
};

class paramnotset : public mxException 
{
public:
   paramnotset( const std::string & esrc, 
                const std::string & efile, 
                const int & line, 
                const std::string & expl
              ) : mxException(esrc, MXE_PARAMNOTSET,MXE_PARAMNOTSET_NAME, efile, line, expl)
   {
   }
};

class sizeerr : public mxException 
{
public:
   sizeerr( const std::string & esrc, 
            const std::string & efile, 
            const int & line, 
            const std::string & expl
          ) : mxException(esrc, MXE_SIZEERR,MXE_SIZEERR_NAME, efile, line, expl)
   {
   }
};

class allocerr : public mxException 
{
public:
   allocerr( const std::string & esrc, 
             const std::string & efile, 
             const int & line, 
             const std::string & expl
           ) : mxException(esrc, MXE_ALLOCERR,MXE_ALLOCERR_NAME, efile, line, expl)
   {
   }
};

class liberr : public mxException 
{
public:
   liberr( const std::string & esrc, 
           const std::string & efile, 
           const int & line, 
           const std::string & expl
         ) : mxException(esrc, MXE_LIBERR,MXE_LIBERR_NAME, efile, line, expl)
   {
   }
};


#define mxThrowException( extype, src, expl ) throw extype(src,  __FILE__, __LINE__, expl);

} //err
} //namespace mx

#endif //mxException_hpp
