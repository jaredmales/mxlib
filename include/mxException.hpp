/** \file mxException.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines the mxlib exception class.
  * \ingroup error_handling_files
  * 
*/

#ifndef __mxException__
#define __mxException__

#include <exception>
#include <sstream>

   
/// The mxlib exception class
/** Provides a rich error report via the standard what().
  * \ingroup error_handling
  */
class mxException : public std::exception
{

protected:
   ///Contains the what() string
   char whatstr[4096];
   
public:

   ///The source of the exception, such as stdlib or cfitisio
   std::string except_code_source;
   
   ///The source specific error code
   int  except_code;
   
   ///The mnemonic associated with the error code
   std::string except_code_mnem;
   
   ///The source file where the exception originated
   std::string except_file;
   
   ///The line number of the file where the exception originated
   int  except_line;
   
   ///The long explanation of the error
   std::string except_explanation;

   ///Default constructor
   mxException () noexcept
   {
      except_code = 0;
      except_line = 0;
      
      build_what();
   }

   ///Copy constructor
   mxException (const mxException & e) noexcept
   {
      except_code_source = e.except_code_source;
      except_code = e.except_code;
      except_code_mnem = e.except_code_mnem;
      except_file = e.except_file;
      except_line = e.except_line;
      except_explanation = e.except_explanation;
      
      build_what();
   }

   ///Construct and fill in each of the values.
   mxException(const std::string & esrc, const int & ec, const std::string & emnem, 
               const std::string & efile, const int & line, const std::string & expl)
   {
      except_code_source = esrc;
      except_code = ec;
      except_code_mnem = emnem;
      except_file = efile;
      except_line = line;
      except_explanation = expl;
      
      build_what();
   }

   ///Assignment operator
   mxException & operator=(const mxException & e) noexcept
   {
      except_code_source = e.except_code_source;
      except_code = e.except_code;
      except_code_mnem = e.except_code_mnem;
      except_file = e.except_file;
      except_line = e.except_line;
      except_explanation = e.except_explanation;
      
      build_what();
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
      s << "      source: " << except_code_source << "\n";
      s << "        code: " << except_code << "\n";
      s << "   code mnem: " << except_code_mnem << "\n";
      s << "     in file: " << except_file << "\n";
      s << "     at line: " << except_line << "\n";
      s << " explanation: " << except_explanation << "\n";
      
      snprintf(whatstr, 4096, "%s", s.str().c_str());
   }
   
   ///Return the details of the exception as a single string.
   virtual const char * what() const noexcept
   {   
      return whatstr;//s.str().c_str();
   }
   
};

#endif //__mxException__
