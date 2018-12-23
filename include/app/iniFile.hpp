/** \file iniFile.hpp
 * \author Jared R. Males
 * \brief Declares and defines an ini-style (toml) file parser
 *
 * \ingroup mxApp_files
 *
 */

#ifndef iniFile_hpp
#define iniFile_hpp

#include <unordered_map>

#include "ini.hpp"

namespace mx
{

/// A wrapper for the ini functions.
/** Places results of the ini parser in an unordered map, with keys of "section=name", and value
  * strings containting the value of the config item.
  * 
  * \ingroup mxApp
  */
struct iniFile
{
   typedef std::unordered_map<std::string, std::string> nameMapT; ///< the unordered map type used for storing values.

   nameMapT names; ///< The map of section=name keys to values.

   /// Return a key generated from the section and name.
   /** Constructs the key as "section=name".
    * 
    * \returns the created key.
    */
   std::string makeKey( const std::string & section, ///< [in] The section for the key
                        const std::string & name     ///< [in] the name for the key
                      )
   {
      return section + "=" + name;
   }

   /// Calls the inih parse function with this->handler.
   /** This returns the result of the ini_parse function.
     * 
     * \returns 0 on success
     * \returns >0 is the line number of the first parsing error
     * \returns -1 on file open error
     * \return -2 on memory allocation error 
     */
   int parse( const std::string & fname /**< [in] the full path of the file*/ )
   { 
      return ini_parse(fname.c_str(), handler, this);
   }

   /// Insert a config value.  Appends if the key already exists in the map.
   /** \todo need error checking.
     * 
     * \returns 0 on success
     * \returns -1 on error.
     */ 
   int insert( const char *section, ///< [in] The section for the key
               const char *name,    ///< [in] the name for the key
               const char * value   ///< [in] the value to insert
             )
   {
      std::string nkey = makeKey(section, name);

      names[nkey]; //Inserts with empty value if doesn't exist, otherwise does nothing.
      
      names[nkey] += value; //This is just an update to the existing value.

      return 0;
   }

   /// Config entry handler for the parse.
   /** Calls insert, and returns its result.  Any non-zero return will cause ini_parse to report the current
     * line number as an error. 
     * 
     * \returns 1 on success
     * \returns 0 on error  
     */
   static int handler( void* user,          ///< [in] a pointer to this.
                       const char* section, ///< [in] the section of the config entry
                       const char* name,    ///< [in] the name of the config entry
                       const char* value    ///< [in] the value of the config entry
                     )
   {
      iniFile * iF = static_cast<iniFile *>(user);
      return (iF->insert(section, name, value) == 0);
   }

   /// Get the number of entries for the given section and name.
   /** Is either 1 or 0, depending on if this cconfig key exists.
     * 
     * \returns the number of entries for the section=name key.
     */
   size_t count( const std::string &section, ///< [in] The section for the key
                 const std::string & name    ///< [in] the name for the key
               )
   {
      return names.count(makeKey(section, name));
   }

   /// Erase the entry for the given section and name.
   /** 
     * \returns 0 on sucess
     * \returns -1 if no section=name entry in the map.
     */ 
   int erase( const std::string & section, ///< [in] The section for the key
              const std::string & name     ///< [in] the name for the key
            )
   {
      if ( names.erase(makeKey(section,name)) > 0) return 0;
      else return -1;
   }
   
   /// Get the value associated with this section=name pair.
   /**
     * \returns the value if the section=name key exists
     * \returns and empty string if the section=name key does not exist.
     */ 
   std::string operator()( const std::string &section, ///< [in] The section for the key
                           const std::string & name    ///< [in] the name for the key
                         )
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

   /// Get the value associated with this name with an empty section.
   /**
     * \returns the value if the section=name key exists, with section empty
     * \returns and empty string if the section=name key does not exist, with section empty.
     */
   std::string operator()( const std::string & name /**< [in] the name for the key*/ )
   {
      return operator()("", name);
   }


};

} //namespace mx

#endif
