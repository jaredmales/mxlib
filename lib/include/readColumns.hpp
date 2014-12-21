/** \file readColumns.hpp
  * \author Jared R. Males
  * \brief Declaration and definition of a utility to read in columns from a text file.
  * \ingroup ioutils
  */

#ifndef __readColumns_hpp__
#define __readColumns_hpp__

#include <fstream>
#include <string>

namespace mx
{
   
template<char delim=' ', char eol='\n'> 
void readcol(std::ifstream & fin, int i)
{
   return;
}

template<char delim=' ', char eol='\n', typename arrT, typename... arrTs> 
void readcol(std::ifstream & fin, int i, arrT* & array, arrTs*&... arrays)
{
   static const unsigned short int nargs = sizeof...(arrTs);
   
   if(nargs >= 0) 
   {
      std::string str;
      fin >> str;
      
      array[i] = convertFromString<arrT>(str);
   }
      
   readcol<delim,eol>(fin, i, arrays...);
 
}

/** \addtogroup ioutils
  * @{
  */

///Read in columns from a text file
/** This function opens a file containing data formatted in columns and reads in the data row by raw.  
  * The data are stored in previously allocated arrays. 
  * 
  * Example:
  * \code
  * int i1 = new int(5);
  * float f1 = new float(5);
  * double d1 = new double(5);
  * 
  * readcol("data_file.txt", 5, i1, f1, d1);
  * \endcode
  * 
  * Note that the types of the arrays do not need to be specified as template arguments.  
  * 
  * The format of the file can be specified with template arguments like
  * \code
  * pout<',', '\r'>(readcol("data_file.csv", 5, i1, f1, d1););
  * \endcode
  * which sets the delimmiter to comma, and the end-of-line to \r.
  * 
  * \tparam delim is the character separating columns,  by default this is space.
  * \tparam eol is the end of line character
  * \tparam arrTs... a variadic list of array types
  * 
  * \param fname is the file name to read from
  * \param n_lines is the number of lines in the file.
  * \param arrays... a variadic list of allocate arrays of the appropriate type. Any number of values with mixed type can be specified.  These must be pre-allocated with at least n_lines length.
  * 
  */
template<char delim=' ', char eol='\n', typename... arrTs> 
void readcol(const std::string & fname, int n_lines,  arrTs*&... arrays)
{
   //open file
   std::ifstream fin;
   fin.open(fname);
   
   for(int j = 0; j< n_lines; j++)
   {
      readcol<delim,eol>(fin, j, arrays...);
      
   }
   
   fin.close();
}

/// @}

} //namespace mx

#endif //__readColumns_hpp__