/** \file fileUtils.cpp
  * \brief Definitions of utilities for working with files
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#include "fileUtils.hpp"

namespace mx
{

std::vector<std::string> getFileNames(const std::string & directory, const std::string & prefix, const std::string & extension)
{

   std::vector<std::string> vect;


   if( exists(directory) )
   {
      if(is_directory(directory) )
      {

//         auto it = directory_iterator(directory);
//         auto it_end = directory_iterator();
         
         typedef std::vector<path> vec;             // store paths,
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

std::vector<std::string> getFileNames(const std::string & directory, const std::string & extension)
{
   return getFileNames(directory, "", extension);
}

std::vector<std::string> getFileNames(const std::string & directory)
{
   return getFileNames(directory, "", "");
}


   
   
} //namespace mx

