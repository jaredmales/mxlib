/** \file readColumns.hpp
  * \author Jared R. Males
  * \brief A utility to read in columns from a text file.
  * \ingroup ioutils
  */

#ifndef __readColumns_hpp__
#define __readColumns_hpp__

#include <fstream>
#include <string>
#include <cstring>

#include "stringUtils.hpp"


namespace mx
{
   
template<char delim=' ', char eol='\n'> 
void readcol(char * sin, int sz)
{
   return;
}

template<char delim=' ', char eol='\n', typename arrT, typename... arrTs> 
void readcol(char * sin, int sz, arrT & array, arrTs &... arrays)
{
   static const unsigned short int nargs = sizeof...(arrTs);
   std::string str;
   
   int i=0;
   int l = strlen(sin);
   
   if(l < 1) return;
   
   //Eat white space
   while( isspace(sin[i]) && sin[i] != eol && i < l) ++i;
   sin = sin + i;
   sz = sz -i;

   
   if(nargs >= 0) 
   {
      std::stringstream sinstr(sin);
      
      std::getline(sinstr, str, delim);
      
      array.push_back(convertFromString<typename arrT::value_type>(str));
   }
      
   sin += ( str.size()+1)*sizeof(char);
   sz -= ( str.size()+1)*sizeof(char);
   
   readcol<delim,eol>(sin, sz, arrays...);
 
}

/** \addtogroup ioutils
  * @{
  */

///Read in columns from a text file
/** This function opens a file containing data formatted in columns and reads in the data row by row.  
  * The data are stored in std::vectors, which should not be pre-allocated (though they could be reserve()-ed). 
  * 
  * Example:
  * \code
  * std::vector<int> i1;
  * std::vector<float> f1;
  * std::vector<double> d1;
  * 
  * readcol("data_file.txt", i1, f1, d1);
  * \endcode
  * 
  * Note that the types of the vectors do not need to be specified as template arguments.  
  * 
  * The format of the file can be specified with template arguments like
  * \code
  * pout<',', ';', '\r'>(readcol("data_file.csv", i1, f1, d1););
  * \endcode
  * which sets the delimmiter to comma, the comment character to ;, and the end-of-line to \\r.
  * 
  * \tparam delim is the character separating columns,  by default this is space.
  * \tparam comment is the character starting a comment.  by default this is #
  * \tparam eol is the end of line character.  by default this is \n
  * \tparam arrTs a variadic list of array types. this is not specified by the user.
  * 
  * \param fname is the file name to read from
  * \param arrays a variadic list of std::vectors. Any number with mixed value_type can be specified. Neither allocated nor
  *               cleared, so repeated calls will append data.
  * 
  */
template<char delim=' ', char comment='#', char eol='\n', typename... arrTs> 
void readColumns(const std::string & fname, arrTs &... arrays)
{
   //open file
   std::ifstream fin;
   fin.open(fname);
   
   int lineSize = 1024;
   char * line = new char[lineSize];
  
   while(fin.good())
   {
      fin.getline(line, lineSize, eol);
            
      int i=0;
      int l = strlen(line);
      
      //Find start of comment.
      while(line[i] != comment && i < l ) ++i;
      line[i] = '\0';
      
      readcol<delim,eol>(line, strlen(line), arrays...);      
   }
   
   fin.close();
   
   delete line;
}

/// @}

} //namespace mx

#endif //__readColumns_hpp__