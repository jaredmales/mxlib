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

#include <string>
#include <vector>
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

///@} -fileutils

} //namespace mx

#endif //__fileUtils_hpp__
