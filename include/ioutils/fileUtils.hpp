
/** \file fileUtils.hpp
  * \brief Declarations of utilities for working with files
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup fileutils
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

#ifndef filtUtils_hpp
#define filtUtils_hpp

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <libgen.h>
#include <cmath>

#include <boost/filesystem.hpp>


using namespace boost::filesystem;

namespace mx
{
namespace ioutils
{

/** \addtogroup fileutils
  * @{
  */

///Get a list of file names from the specified directory, specifying a prefix, a substring to match, and an extension
/**
  * \returns a std::vector\<std::string\> which contains the matching file names.
  */
inline
std::vector<std::string> getFileNames( const std::string & directory, ///< [in] the path to the directory to search. Can not be empty.
                                       const std::string & prefix,  ///< [in] the file name prefix (the beginning characters of the file name) to search for, if "" then not used.
                                       const std::string & substr,   ///< [in] a substring of the filename to earch for, if "" then not used.
                                       const std::string & extension  ///< [in] the file name extension to search for, if "" then not used.  Note that this must include the ".", as in".ext".
                                     )
{
   typedef std::vector<path> vec;             // store paths,

   std::vector<std::string> vect;
   if( exists(directory) )
   {
      if(is_directory(directory) )
      {
         vec v;                                // so we can sort them later

         copy(directory_iterator(directory), directory_iterator(), back_inserter(v));

         sort(v.begin(), v.end());             // sort, since directory iteration
                                              // is not ordered on some file systems

         auto it = v.begin();
         auto it_end = v.end();

         while(it != it_end)
         {
            bool inc = true;

            if(extension != "")
            {
               if(it->extension() != extension)
               {
                  inc = false;
               }
            }

            if(prefix != "" && inc)
            {
               std::string p = it->filename().generic_string();

               if( p.size() < prefix.size() )
               {
                  inc = false;
               }
               else
               {
                  if(p.compare(0, prefix.size(), prefix) != 0)
                  {
                     inc = false;
                  }
               }
            }

            if(substr != "" && inc)
            {
               std::string p = it->filename().generic_string();
               if(p.find(substr) == std::string::npos)
               {
                  inc = false;
               }
            }

            if(inc)
            {
               vect.push_back(it->native());
            }

            ++it;
         }
      }
      else
      {
         std::cerr << "is not a directory\n";
      }

   }
   else
   {
      std::cerr << "directory does not exist\n";
   }

   return vect;
}

///Get a list of file names from the specified directory, specifying the extension
/** \overload
  *
  * \returns a std::vector\<std::string\> which contains the matching file names.
  */
inline
std::vector<std::string> getFileNames( const std::string & directory, ///< [in] the path to the directory to search. Can not be empty.
                                       const std::string & extension ///< [in] the file name extension to search for, if "" then not used.  Note that this must include the ".", as in ".ext".
                                     )
{
   return getFileNames(directory, "", "", extension);
}

///Get a list of file names from the specified directory
/** \overload
  *
  * \returns a std::vector\<std::string\> which contains the matching file names.
  */
inline
std::vector<std::string> getFileNames( const std::string & directory /**< [in] the path to the directory to search.  Can not be empty. */)
{
   return getFileNames(directory, "", "", "");
}


///Prepend and/or append strings to a file name, leaving the directory and extension unaltered.
/**
  * \returns the new file name
  */
inline
std::string  fileNamePrependAppend( const std::string & fname,  ///< [in] the original file name, possibly including a directory and extension
                                    const std::string & prepend, ///< [in] is the string to insert at the beginning of the file name after the path
                                    const std::string & append ///< [in] is the string to insert at the end of the file name, before the extension
                                  )
{
   std::string dir, base, ext;

   path p = fname;
   dir = p.parent_path().string();
   base = p.stem().string();
   ext = p.extension().string();


   return dir +'/' + prepend + base + append + ext;


}

///Append a string to a file name, leaving the directory and extension unaltered.
/**
  * \param fname the original file name, possibly including a directory and extension
  * \param append is the string to insert at the end of the file name, before the extension
  *
  * \returns the new file name
  */
inline
std::string  fileNameAppend(const std::string & fname, const std::string & append)
{
   return fileNamePrependAppend(fname, "", append);
}

///Prepend strings to a file name, leaving the directory and extension unaltered.
/**
  * \param fname the original file name, possibly including a directory and extension
  * \param prepend is the string to insert at the beginning of the file name after the path
  *
  * \returns the new file name
  */
inline
std::string  fileNamePrepend(const std::string & fname, const std::string & prepend)
{
   return fileNamePrependAppend(fname, prepend, "");
}

///Get the next file in a numbered sequence
/** Searches for files in the path designated by basename of the form basenameXXXXextension
  * where the number of digits in XXXX is set by the \a ndigit template parameter.
  *
  * \param[in] basename  path and initial name of the file
  * \param[in] extension [optional] extension to append after the number. Default is empty.
  * \param[in] startat [optional] number to start the search from.  Default is 0.
  *
  * \retval std::string containing the next filename.
  *
  * \tparam ndigit [optional] number of digits in string representation of the number.  Default is 4.
  */
template<int ndigit = 4>
std::string getSequentialFilename( const std::string & basename,
                                   const std::string & extension = "",
                                   const int startat = 0)
{
   //int maxdig = 1;
   //for(int j=0;j<ndigit;++j) maxdig *= 10;
   int maxdig = pow(10, ndigit);

   char digstr[ndigit+1];
   int i = startat;

   std::stringstream outn;

   snprintf(digstr,ndigit+1,"%04d", i);

   outn << basename;
   outn << digstr;
   outn << extension;

   while(boost::filesystem::exists(outn.str()) && i < maxdig)
   {
      ++i;
      outn.str("");

      snprintf(digstr,ndigit+1,"%04d", i);

      outn << basename;
      outn << digstr;

      outn << extension;
   }

   return outn.str();
}

///Get the next file in a numbered sequence
/** \overload
  *
  * Searches for files in the path designated by basename of the form basenameXXXX
  * where the number of digits in XXXX is set by the \a ndigit template parameter.
  *
  * \param[in] basename  path and initial name of the file
  * \param[in] startat number to start the search from.
  *
  * \retval std::string containing the next filename.
  *
  * \tparam ndigit [optional] number of digits in string representation of the number.  Default is 4.
  */
template<int ndigit = 4>
std::string getSequentialFilename( const std::string & basename,
                                   const int startat )
{
   return getSequentialFilename<ndigit>(basename, "", startat);
}

///Get the size in bytes of a file
/** Uses fstat.
  *
  * \returns the file size if fd is valid and no errors occur
  * \returns -1 on an error 
  */ 
inline
off_t fileSize( int fd /**< [in] an open file descriptor */)
{
   if (fd == -1) 
   {
      return -1;
   }
 
   struct stat stbuf;
 
   if ((fstat(fd, &stbuf) != 0) || (!S_ISREG(stbuf.st_mode))) 
   {
      return -1; 
   }
  
   return stbuf.st_size;
  
}

///@} -fileutils

} //namespace ioutils
} //namespace mx

#endif //fileUtils_hpp
