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
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

namespace mx
{

/** \addtogroup fileutils
  * @{
  */

///Get a list of file names from the specified directory, specifying a prefix and an extension
/** 
  * \param directory the path to the directory to search
  * \param extension the file name extension to search for, if "" then not used 
  * \param prefix the file name prefix (the beginning characters of the file name) to search for, if "" then not used.
  *
  * \returns a std::vector<std::string> containing the matching file names.
  */ 
std::vector<std::string> getFileNames(const std::string & directory, const std::string & prefix, const std::string & extension);


///Get a list of file names from the specified directory, specifying the extension
/** 
  * \param directory the path to the directory to search
  * \param extension the file name extension to search for, if "" then not used 
  *
  * \returns a std::vector<std::string> containing the matching file names.
  */ 
std::vector<std::string> getFileNames(const std::string & directory, const std::string & extension);

///Get a list of file names from the specified directory
/** 
  * \param directory the path to the directory to search
  *
  * \returns a std::vector<std::string> containing the matching file names.
  */ 
std::vector<std::string> getFileNames(const std::string & directory);


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
                                   const std::string & extension ="", 
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

///@} -fileutils

} //namespace mx

#endif //__fileUtils_hpp__
