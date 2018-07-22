

#ifndef iniFile_hpp
#define iniFile_hpp

#include <unordered_map>

#include "ini.hpp"

namespace mx
{

///A wrapper for the ini functions.
struct iniFile
{
   typedef std::unordered_map<std::string, std::string> nameMapT;

   nameMapT names;

   std::string makeKey(const std::string & section, const std::string & name)
   {
      return section + "=" + name;
   }

   int parse(const std::string & fname)
   {
      return ini_parse(fname.c_str(), handler, this);
   }

   int insert(const char *section, const char *name, const char * value)
   {
      std::string nkey = makeKey(section, name);

      names[nkey];
//       if (names[nkey].size() > 0) //This is where the insertion actually occurs.
//       {
//         names[nkey] += ""; //We actually do nothing.
//       }

      names[nkey] += value; //This is just an update to the existing value.

      return 0;
   }

   static int handler( void* user,
                       const char* section,
                       const char* name,
                       const char* value
                     )
   {
      iniFile * iF = static_cast<iniFile *>(user);
      return iF->insert(section, name, value);
   }

   int count(const std::string &section, const std::string & name)
   {
      return names.count(makeKey(section, name));
   }

   std::string operator()(const std::string &section, const std::string & name)
   {
      std::string key = makeKey(section, name);
      if(names.count(key) > 0)
      {
         return names[key];
      }
      else
      {
         return std::string("");
      }
   }

   std::string operator()(const std::string & name)
   {
      return names[makeKey("", name)];
   }


};

} //namespace mx

#endif
