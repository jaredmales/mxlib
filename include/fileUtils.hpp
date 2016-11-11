
/** \file fileUtils.hpp
  * \brief Declarations of utilities for working with files
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup fileutils
  *
  */

#ifndef __fileUtils_hpp__
#define __fileUtils_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <libgen.h>
 
#include <boost/filesystem.hpp>


using namespace boost::filesystem;

namespace mx
{

/** \addtogroup fileutils
  * @{
  */

///Get a list of file names from the specified directory, specifying a prefix, a substring to match, and an extension
/** 
  * \param[in] directory the path to the directory to search
  * \param[in] prefix the file name prefix (the beginning characters of the file name) to search for, if "" then not used.
  * \param[in] extension the file name extension to search for, if "" then not used 
  *
  * \retval std::vector<std::string> which contains the matching file names.
  */ 
inline
std::vector<std::string> getFileNames(const std::string & directory, const std::string & prefix, const std::string & substr, const std::string & extension)
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
               if(p.find(prefix) != 0)
               {
                  inc = false;
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
/** 
  * \param[in] directory the path to the directory to search
  * \param[in] extension the file name extension to search for, if "" then not used 
  *
  * \retval std::vector<std::string> which contains the matching file names.
  */ 
inline
std::vector<std::string> getFileNames(const std::string & directory, const std::string & extension)
{
   return getFileNames(directory, "", "", extension);
}

///Get a list of file names from the specified directory
/** \overload
  * \param directory the path to the directory to search
  *
  * \retval std::vector<std::string> which contains the matching file names.
  */ 
inline
std::vector<std::string> getFileNames(const std::string & directory)
{
   return getFileNames(directory, "", "", "");
}


///Prepend and/or append strings to a file name, leaving the directory and extension unaltered.
/**
  * \param fname the original file name, possibly including a directory and extension
  * \param prepend is the string to insert at the beginning of the file name after the path
  * \param append is the string to insert at the end of the file name, before the extension
  * 
  * \returns the new file name
  */
inline
std::string  fileNamePrependAppend(const std::string & fname, const std::string & prepend, const std::string & append)
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
   int maxdig = 1;
   for(int j=0;j<ndigit;++j) maxdig *= 10;
   
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

///@} -fileutils

} //namespace mx

#endif //__fileUtils_hpp__
