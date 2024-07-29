/** \file fileUtils.cpp
  * \brief Definitions of utilities for working with files
  *
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  * \ingroup fileutils
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "ioutils/fileUtils.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <libgen.h>
#include <cmath>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <boost/filesystem.hpp>

using namespace boost::filesystem;

namespace mx
{
namespace ioutils
{

int createDirectories( const std::string & path )
{
   //Use the non throwing version and silently ignore EEXIST errors 
   boost::system::error_code ec;
   boost::filesystem::create_directories(path, ec);
   if(ec.value() != boost::system::errc::success && ec.value() != boost::system::errc::file_exists)
   {
      return -1;
   }

   return 0;
}

std::string pathStem(const std::string & fname)
{
   boost::filesystem::path p(fname);
   return p.stem().string();
}

std::string pathFilename( const std::string & fname)
{
   boost::filesystem::path p(fname);
   return p.filename().string();
}

std::string parentPath(const std::string & fname)
{
   boost::filesystem::path p(fname);
   return p.parent_path().string();
}



std::vector<std::string> getFileNamesOld( const std::string & directory, 
                                       const std::string & prefix,    
                                       const std::string & substr,    
                                       const std::string & extension  
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

         std::sort(v.begin(), v.end());             // sort, since directory iteration
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
         std::cerr << directory << " is not a directory\n";
      }

   }
   else
   {
      std::cerr << "directory " << directory << " does not exist\n";
   }

   return vect;
}

std::vector<std::string> getFileNames( const std::string & directory, 
                                       const std::string & prefix,    
                                       const std::string & substr,    
                                       const std::string & extension  
                                     )
{
   //typedef std::vector<path> vec;             // store paths,

   std::vector<std::string> vect;
   if( exists(directory) )
   {
      if(is_directory(directory) )
      {
         directory_iterator it{directory};
         auto it_end = directory_iterator{};
         for(it; it != it_end; ++it)
         {
            if(extension != "")
            {
               if(it->path().extension() != extension)
               {
                  continue;
               }
            }

            std::string p = it->path().filename().generic_string();

            if(prefix != "")
            {
               if( p.size() < prefix.size() )
               {
                  continue;
               }
               else
               {
                  if(p.compare(0, prefix.size(), prefix) != 0)
                  {
                     continue;
                  }
               }
            }

            if(substr != "")
            {
               if(p.find(substr) == std::string::npos)
               {
                  continue;
               }
            }

            //If here then it passed all checks
            vect.push_back(it->path().native());

         }

         sort(vect.begin(), vect.end());
      }
      else
      {
         std::cerr << directory << " is not a directory\n";
      }
   }
   else
   {
      std::cerr << "directory " << directory << " does not exist\n";
   }

   return vect;
}


std::vector<std::string> getFileNames( const std::string & directory, 
                                       const std::string & extension  
                                     )
{
   return getFileNames(directory, "", "", extension);
}

std::vector<std::string> getFileNames( const std::string & directory )
{
   return getFileNames(directory, "", "", "");
}


std::string  fileNamePrependAppend( const std::string & fname,   
                                    const std::string & prepend, 
                                    const std::string & append   
                                  )
{
   std::string dir, base, ext;

   path p = fname;
   dir = p.parent_path().string();
   base = p.stem().string();
   ext = p.extension().string();


   return dir +'/' + prepend + base + append + ext;


}

std::string  fileNameAppend( const std::string & fname, 
                             const std::string & append
                           )
{
   return fileNamePrependAppend(fname, "", append);
}

std::string  fileNamePrepend( const std::string & fname, 
                              const std::string & prepend
                            )
{
   return fileNamePrependAppend(fname, prepend, "");
}

std::string getSequentialFilename( const std::string & basename,
                                   const std::string & extension,
                                   const int startat,
                                   const int ndigit
                                 )
{
   //int maxdig = 1;
   //for(int j=0;j<ndigit;++j) maxdig *= 10;
   int maxdig = pow(10, ndigit);

   char formstr[64];
   snprintf(formstr, sizeof(formstr), "%%0%dd", ndigit);
   
   char digstr[64];
   int i = startat;

   std::stringstream outn;

   snprintf(digstr,sizeof(digstr),formstr, i);

   outn << basename;
   outn << digstr;
   outn << extension;

   while(boost::filesystem::exists(outn.str()) && i < maxdig)
   {
      ++i;
      outn.str("");

      snprintf(digstr,sizeof(digstr), formstr, i);

      outn << basename;
      outn << digstr;

      outn << extension;
   }

   return outn.str();
}

off_t fileSize( int fd )
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

off_t fileSize( FILE * f )
{
   return fileSize(fileno(f));  
}


} //namespace ioutils
} //namespace mx

