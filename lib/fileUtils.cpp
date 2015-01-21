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

         auto it = directory_iterator(directory);
         auto it_end = directory_iterator();
         
         while(it != it_end)
         {
            bool inc = true;
            
            if(extension != "")
            {
               if(it->path().extension() != extension)
               {
                  inc = false;
               }
            }
            
            if(prefix != "" && inc)
            {
               std::string p = it->path().filename().generic_string();
               if(p.find(prefix) != 0)
               {
                  inc = false;
               }
            }
                  
            if(inc)
            {
               vect.push_back(it->path().native());
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

