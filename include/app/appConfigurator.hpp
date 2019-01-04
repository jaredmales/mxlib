/** \file appConfigurator.hpp
 * \author Jared R. Males
 * \brief An application configuration manager
 *
 * \ingroup mxApp_files
 *
 */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef appConfigurator_hpp
#define appConfigurator_hpp

#include <list>
#include <fstream>

#include "../mxlib.hpp"
#include "../mxError.hpp"


#include "clOptions.hpp"
#include "iniFile.hpp"
#include "configTarget.hpp"


namespace mx
{
namespace app 
{

/// Class to manage a set of configurable values, and read their values from config/ini files and the command line.
/** 
  * The configuration files are ini-style, with sections.  That is
  \verbatim
  key1=value1
  key2=value2

  [section1]
  key3=value3
  key4=value4,value4.1, value4.2, value4.3

  [section2]
  key3=value5
  key3=value5.1
  key4=value6_over_
       multiple_lines

  \endverbatim
  * such that section1.key3 is distinct from section2.key3  (they must have different config-target names though).
  * 
  * Additional syntax rules:
  * - Leading whitespace is stripped from the value, so `key=val` and `key= val` are equivalent.
  * - Additional entries within one file with the same section and key are appended to the previous entry.
  *   So the value of section2.key3 is "value5value5.1".   
  * - Multi-line values are handled such that in the above example the result is key4=value6_over_multiple_lines.  
  * - Vectors are input as comma separated lists, as in section1.key4 above.  Leading whitespace is stripped from each 
  *   component of the vector.
  *
  * \todo add handling of += in subsequent files.
  *
  * The command line parser handles both short-opt ("-h -vArg -n Arg") and long-opt ("--help --value=Arg --number=Arg") styles.
  * 
  * 
  * \bug a config=value pair listed in a conf file twice seems to cause a failure, even if they are the same value.
  *
  * \ingroup mxApp
  *
  */
struct appConfigurator
{
   ///Iterator for the targets unordered_map
   typedef std::unordered_map<std::string, configTarget>::iterator targetIterator;

   ///Iterator for the clOnlyTargets list.
   typedef std::list<configTarget>::iterator clOnlyTargetIterator;

   ///The targets are stored in an unordered_map for fast access by key.
   std::unordered_map<std::string, configTarget> m_targets;

   /// Config file entries present in the file(s), but not corresponding to a target when parsed.   Set aside for possible analysis.
   std::unordered_map<std::string, configTarget> m_unusedConfigs;

   
   ///Targets which are only for the command line are stored separately in a list.
   std::list<configTarget> clOnlyTargets;

   /// Non-option arguments from the command line.
   std::vector<std::string> nonOptions;

   
   /// Running count of options added, used to track order.
   int nAdded {0};

   /// Flag controlling whether or not to record config sources
   bool m_sources {false};
   
   /// Clear the containers and free up the associated memory.
   void clear();

   /// Add a configTarget
   /** Note that if name is a duplicate but the section and keyword are empty, it is handled as command-line only.
     */
   void add( const configTarget & tgt /**< [in] The configuration target to add */);

   /// Add a configTarget
   /**
     * \overload
     */
   void add( const std::string &n,        ///< [in] The name of the target
             const std::string &so,       ///< [in] The command-line short option (e.g. "f" for -f)
             const std::string &lo,       ///< [in] The command-line long option (e.g. "file" for --file)
             int clt,                     ///< [in] The command-line option type, argType::false, argType::true, argType::optional, argType::required
             const std::string & s,       ///< [in] The config file section name, can be empty ""
             const std::string & kw,      ///< [in] The config file keyword, read in a "keyword=value" pair
             bool isReq = false,          ///< [in] Whether or not this is option is required to be set
             const std::string & ht = "", ///< [in] The type to display in the help message
             const std::string & he = ""  ///< [in] The explanation to display in the help message
           );

   ///Parse the command line, updating the targets
   void parseCommandLine( int argc,  ///< [in] standard command line result specifying number of argumetns in argv
                          char ** argv, ///< [in] standard command line result containing the arguments.
                          const std::string & oneTarget = "" ///< [in] [optional] if not empty, then only this target is extracted by the parser.
                        );


   ///Read and parse a config/ini file, updating the targets
   /** \todo handle += here, by appending to the last value as if a vector.
     */
   int readConfig( const std::string & fname,     ///< [in] the config file name 
                   bool reportFileNotFound = true ///< [in] [optiona] control whether a file not found is reported.
                 );

   /// Check if a target has been set by the configuration
   /**
     * \returns true if the configuration set at least one value for this target
     * \returns false if no value was set.
     */
   bool isSet( const std::string & name,                               ///< [in] the target name 
               std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
             );
   
   /// Check if a target has been set by the configuration
   /**
     * \overload
     * 
     * \returns true if the configuration set at least one value for this target
     * \returns false if no value was set.
     */
   bool isSet(const std::string & name /**< [in] the target name */);

   ///Get the number of different values set for the specified config target
   /**
     * \returns the number of different values set for name.
     */
   int count( const std::string & name,                               ///< [in] the target name 
              std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
            );
   
   ///Get the number of different values set for the specified config target
   /**
     * \overload
     * 
     * \returns the number of different values set for name.
     */
   int count( const std::string & name /**< [in] the target name */);

   ///Get the command line verbosity count for this option.
   /** E.g., -v ==> 1, -vv ==> 2, -vvv ==> 3, etc.  Note that for this to work
     * properly, this must be of type mx::argType::True.
     * 
     * \returns the verbosity count.
     */
   int verbosity( const std::string & name,                               ///< [in] the target name 
                  std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
                );

   ///Get the command line verbosity count for this option.
   /** E.g., -v ==> 1, -vv ==> 2, -vvv ==> 3, etc.  Note that for this to work
     * properly, this must be of type mx::argType::True.
     * 
     * \overload
     * 
     * \returns the verbosity count.
     */
   int verbosity( const std::string & name /**< [in] the target name */);
   
   /// Get the i-th value of the target, converted to the specified type
   /** The supplied value is only altered if the config target was set, which preserves.
     * default values.
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename typeT>
   int get( typeT & v,                                              ///< [out] the variable to store the value in, unaltered if not set.
            const std::string & name,                               ///< [in] the config target name
            size_t i,                                               ///< [in] the number of config specification to get.
            std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
          );

   /// Get the i-th value of the target from the used set, converted to the specified type
   /** The supplied value is only altered if the config target was set, which preserves.
     * default values.
     * 
     * \overload
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename typeT>
   int get( typeT & v,                ///< [out] the variable to store the value in
            const std::string & name, ///< [in] the config target name
            size_t i                  ///< [in] the number of config specification to get.
          );
   
   /// Get the final value of the target, converted to the specified type
   /** The supplied value is only altered if the config target was set, which preserves.
     * default values.
     *
     * \overload
     *
     * \retval 0 on success.
     * \retval -1 on error.
     */
   template<typename typeT>
   int get( typeT & v,                                              ///< [out] the variable to store the value in
            const std::string & name,                               ///< [in] the config target name
            std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
          );
   
   /// Get the final value of the target from the used set, converted to the specified type
   /** The supplied value is only altered if the config target was set, which preserves.
     * default values.
     *
     * \overload
     *
     * \retval 0 on success.
     * \retval -1 on error.
     */
   template<typename typeT>
   int get( typeT & v, ///< [out] the variable to store the value in
            const std::string & name ///< [in] the config target name
          );

   /// Get the i-th value of the target, converted to the specified config target
   /** The vector is only populated if the config target was set.  If it is populated,
     * it is cleared first. Thus if a vector filled with default values is passed
     * in, it will only be overwritten if the user specified new values.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename typeT>
   int get( std::vector<typeT> & v,                                 ///< [out] the vector to populate
            const std::string & name,                               ///< [in] the config target name.
            size_t i,                                               ///< [in] the number of config specification to get.
            std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
          );
   
   /// Get the i-th value of the target, converted to the specified config target
   /** The vector is only populated if the config target was set.  If it is populated,
     * it is cleared first. Thus if a vector filled with default values is passed
     * in, it will only be overwritten if the user specified new values.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename typeT>
   int get( std::vector<typeT> & v, ///< [out] the vector to populate
            const std::string & name, ///< [in] the config target name.
            size_t i ///< [in] the number of config specification to get.
          );
   
   /// Get the i-th value of the target as a vector containing the specified type.
   /** The vector is only populated if the config target was set.  If it is populated,
     * it is cleared first. Thus if a vector filled with default values is passed
     * in, it will only be overwritten if the user specified new values.
     *
     * \retval 0 on success.
     * \retval -1 on error.
     */
   template<typename typeT>
   int get( std::vector<typeT> & v,                                 ///< [out] the vector to populate
            const std::string & name,                               ///< [in] the config target name.
            std::unordered_map<std::string, configTarget> & targets ///< [in] the map of config targets to use
          );

   /// Get the i-th value of the target in the used set, as a vector containing the specified type.
   /** The vector is only populated if the config target was set.  If it is populated,
     * it is cleared first. Thus if a vector filled with default values is passed
     * in, it will only be overwritten if the user specified new values.
     *
     * \retval 0 on success.
     * \retval -1 on error.
     */
   template<typename typeT>
   int get( std::vector<typeT> & v,   ///< [out] the vector to populate
            const std::string & name ///< [in] the config target name.
          );
   
   /// Access operator, configures a value by calling get.
   /**
     * \retval 0 on success
     * \retval -1 on error
     */ 
   template<typename typeT>
   int operator()( typeT & v,               ///< [out] the variable to populate (either scalar or vector), will be unaltered if not set.
                   const std::string & name ///< [in] the config target name.
                 );

   /// Configure a value from the unused map, using the iniFile key.
   /**
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename typeT>
   int configUnused( typeT & v,              ///< [out] the variable to populate (either scalar or vector), will be unaltered if not set.
                     const std::string & key ///< [in] the iniFile key for this target. 
                   );
   
   /// Configure a value from the unused map, using the section and keyword.
   /**
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename typeT>
   int configUnused( typeT & v,                   ///< [out] the variable to populate (either scalar or vector), will be unaltered if not set.
                     const std::string & section, ///< [in] the section name for this target
                     const std::string & keyword  ///< [in] the keyword for this target.
                   );
   
   /// Get the unique sections in the unused config targets.
   /**
     * \retval 0 on success
     * \retval -1 on error
     */
   int unusedSections( std::vector<std::string> & sections );
   
   /// Check if a target has been set in the unused configuration
   /**
     * \returns true if the unused configuration set at least one value for this target
     * \returns false if no value was set.
     */
   int isSetUnused( const std::string & name /**< [in] the target name */);
   
   
   /// Call an external logging function whenever a config value is accessed by get or operator().
   /** Only called if this is not a nullptr (the default), otherwise no logging or reporting is done.
     */
   void (*configLog)( const std::string & name,     ///< [in] The name of the config target.
                      const int & code,             ///< [in] The type code from mx::typeDescription
                      const std::string & valueStr, ///< [in] The value in its string form as found in the configuration 
                      const std::string & source    ///< [in] The source of the value, either default, command line, or a path.
                    ) {nullptr};
                  
};

inline
void appConfigurator::clear()
{
   m_targets.clear();
   clOnlyTargets.clear();
   nonOptions.clear();
   m_unusedConfigs.clear();
}

inline
void appConfigurator::add( const configTarget & tgt )
{

   //First check for duplicate name and command line only
   if(m_targets.count(tgt.name) > 0 && tgt.section == "" && tgt.keyword == "")
   {
      clOnlyTargets.push_back(tgt);
      clOnlyTargets.back().orderAdded = nAdded;
   }
   else
   {
      std::pair<targetIterator, bool> res = m_targets.insert({tgt.name, tgt});
      res.first->second.orderAdded = nAdded;
   }

   ++nAdded;

}

inline
void appConfigurator::add( const std::string &n,
                           const std::string &so,
                           const std::string &lo,
                           int clt,
                           const std::string & s,
                           const std::string & kw,
                           bool isReq,
                           const std::string & ht,
                           const std::string & he )
{
   add( configTarget(n,so,lo,clt, s, kw, isReq, ht, he) );
}

inline
void appConfigurator::parseCommandLine( int argc,
                                        char ** argv,
                                        const std::string & oneTarget
                                      )
{
   if(argc == 0) return;

   clOptions clOpts;

   //First load the options into the clOptions parser
   targetIterator it;
   for(it = m_targets.begin(); it != m_targets.end(); ++it)
   {
      if(it->second.shortOpt == "" && it->second.longOpt == "") continue; //No command line opts specified.

      clOpts.add(it->second.name,it->second.shortOpt.c_str(),it->second.longOpt.c_str(), it->second.clType);
   }

   //Then load the command-line-only options.
   clOnlyTargetIterator cloit;
   for(cloit = clOnlyTargets.begin(); cloit != clOnlyTargets.end(); ++cloit)
   {
      if(cloit->shortOpt == "" && cloit->longOpt == "") continue; //Nothing to add?

      clOpts.add(cloit->name,cloit->shortOpt.c_str(),cloit->longOpt.c_str(), cloit->clType);
   }



   //Now we parse
   clOpts.parse(argc, argv, &nonOptions);

   //If nothing more to do, get out
   if(clOpts.nOpts == 0)
   {
      return;
   }

   //And then load the results in the config target map.
   for(it = m_targets.begin(); it != m_targets.end(); ++it)
   {
      if(oneTarget != "" && it->second.name != oneTarget) continue;

      if(clOpts.optSet(it->second.name))
      {
         std::vector<std::string> args;

         clOpts.getAll(args, it->second.name);
         it->second.values.insert( it->second.values.end(), args.begin(), args.end());
         if(m_sources)
         {
            for(size_t n=0; n < args.size(); ++n) it->second.sources.push_back("command line");
         }
         
         it->second.verbosity = clOpts.count(it->second.name);
         it->second.set = true;
      }
   }

}

int appConfigurator::readConfig( const std::string & fname,
                                 bool reportFileNotFound
                               )
{
   //Handle empty string quietly
   if(fname == "") return 0;

   iniFile iF;

   ///\todo update error handling to include >0 (line numer of parse error) and -2 memory allocation error.
   int prv = iF.parse(fname); 
   
   if( prv == -1)
   {
      if(!reportFileNotFound) return -1;

      mxError("appConfigurator: ", MXE_FILENOTFOUND, "The file " + fname + " was not found");
      return -1;
   }
   
   if( prv == -2)
   {
      mxError("appConfigurator: ", MXE_ALLOCERR, "Memory allocation error in config file parser");
      return -1;
   }
   
   if( prv > 0)
   {
      mxError("appConfigurator: ", MXE_PARSEERR, "Parsing error in " + fname + " at line " + std::to_string(prv));
      return -1;
   }

   targetIterator it;

   for(it = m_targets.begin(); it != m_targets.end(); ++it)
   {
      if(iF.count(it->second.section, it->second.keyword) > 0)
      {
         it->second.values.push_back(iF(it->second.section, it->second.keyword));
         if(m_sources)
         {
            it->second.sources.push_back(fname);
         }
         it->second.set = true;
         
         iF.erase(it->second.section, it->second.keyword); //Erase it from the iniFile map.
      }
   }
   
   //here set aside non-deleted iF entries
   for( auto iFit = iF.names.begin(); iFit != iF.names.end(); ++iFit)
   {
      //Insert or update existing
      m_unusedConfigs[iFit->first].name = iFit->first;
      
      std::string sect, nam;
      iniFile::parseKey(sect, nam, iFit->first);
      m_unusedConfigs[iFit->first].section = sect;
      m_unusedConfigs[iFit->first].keyword = nam;
      
      m_unusedConfigs[iFit->first].values.push_back(iFit->second);
      
      if(m_sources) m_unusedConfigs[iFit->first].sources.push_back(fname);
      
      m_unusedConfigs[iFit->first].set = true;
   }
   
   return 0;
}

inline
bool appConfigurator::isSet( const std::string & name,
                             std::unordered_map<std::string, configTarget> & targets
                           )
{
   if(targets.count(name) == 0) return false;
   return targets[name].set;
}

inline
bool appConfigurator::isSet( const std::string & name )
{
   return isSet(name, m_targets);
}

inline
int appConfigurator::count( const std::string & name,
                            std::unordered_map<std::string, configTarget> & targets
                          )
{
   return targets[name].values.size();
}

inline
int appConfigurator::count( const std::string & name )
{
   return count(name, m_targets);
}

inline
int appConfigurator::verbosity( const std::string & name,
                                std::unordered_map<std::string, configTarget> & targets
                              )
{
   return targets[name].verbosity;
}

inline
int appConfigurator::verbosity( const std::string & name )
{
   return verbosity(name, m_targets);
}

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name,
                          size_t i,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   if(!isSet(name, targets)) 
   {
      if(configLog)
      {
         configLog(name, meta::typeDescription<typeT>::code(), ioutils::convertToString(v), "default"); 
      }
      
      return 0;
   }
   
   if( targets[name].values.size() <= i) return -1;

   v = ioutils::convertFromString<typeT>(targets[name].values[i]);

   //Log it here.
   if(configLog)
   {
      if(m_sources)
      {
         configLog(name, meta::typeDescription<typeT>::code(), targets[name].values[i], targets[name].sources[i]);
      }
      else
      {
         configLog(name, meta::typeDescription<typeT>::code(), targets[name].values[i], "");
      }
   }
   
   return 0;
}

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name,
                          size_t i
                        )
{
   return get(v, name, i, m_targets);
}

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   if(!isSet(name, targets)) 
   {
      if(configLog)
      {
         configLog(name, meta::typeDescription<typeT>::code(), ioutils::convertToString(v), "default");
      }
      
      return 0;
   }
   
   int i = targets[name].values.size() - 1;

   if(i < 0) return -1;

   return get(v, name, i, targets);
}

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name)
{
   return get(v, name, m_targets);
}

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          size_t i,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   if(!isSet(name, targets)) 
   {
      if(configLog)
      {
         configLog(name, meta::typeDescription<typeT>::code(), "[need a vector to string]", "default");
      }
      
      return 0;
   }
   
   if( targets[name].values.size() <= i) return -1;

   std::string s;
   
   
   s = ioutils::convertFromString<std::string>(targets[name].values[i]);
   
   //if( get<std::string>(s, name, i, targets) < 0) return -1;

   //Case that s was set to be empty.
   if(s.size() == 0)
   {
      if(configLog)
      {
         if(m_sources)
         {
            configLog(name, meta::typeDescription<typeT>::code(), "", targets[name].sources[i]);
         }
         else
         {
            configLog(name, meta::typeDescription<typeT>::code(), "", "");
         }
      }
      
      v.clear(); //We clear the vector passed as default
      return 0;
   }

   size_t st;
   size_t com;

   st = 0;
   
   while(::isspace(s[st]) && st < s.size()-1) ++st;
   
   com = s.find(',', st);

   v.clear();

   while(com != std::string::npos)
   {
      v.push_back( ioutils::convertFromString<typeT>(s.substr(st, com-st)) );
      st = com + 1;
      while(::isspace(s[st]) && st < s.size()-1) ++st;
      
      com = s.find(',', st);
   }
   v.push_back( ioutils::convertFromString<typeT>(s.substr(st, s.size()-st)));

   //Log it here.
   if(configLog)
   {
      if(m_sources)
      {
         configLog(name, meta::typeDescription<typeT>::code(), targets[name].values[i], targets[name].sources[i]);
      }
      else
      {
         configLog(name, meta::typeDescription<typeT>::code(), targets[name].values[i], "");
      }
   }
   
   return 0;
}

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          size_t i
                        )
{
   return get(v, name, i, m_targets);
}

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   if(!isSet(name, targets)) 
   {
      if(configLog)
      {
         configLog(name, meta::typeDescription<typeT>::code(), "[need a vector to string]", "default");
      }
      return 0;
   }
   int i = targets[name].values.size() - 1;

   if(i < 0) return -1;

   return get(v, name, i, targets);
}

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name
                        )
{
   return get(v, name, m_targets);
}

template<typename typeT>
int appConfigurator::operator()( typeT & v,
                                 const std::string & name
                               )
{
   return get(v, name, m_targets);
}

template<typename typeT>
int appConfigurator::configUnused( typeT & v,
                                   const std::string & key
                                 )
{
   return get(v, key, m_unusedConfigs);
}

template<typename typeT>
int appConfigurator::configUnused( typeT & v,
                                   const std::string & section,
                                   const std::string & keyword
                                 )
{
   return configUnused(v, iniFile::makeKey(section,keyword));
}

inline
int appConfigurator::unusedSections( std::vector<std::string> & sections )
{
   sections.clear();
   
   //Wind through all the targets
   for(auto it = m_unusedConfigs.begin(); it != m_unusedConfigs.end(); ++it)
   {
      std::string sect = it->second.section;
      
      bool add = true;

      //Check if this section is already in the vector -- is there a std::algorithms way to do this?
      for(size_t i=0;i < sections.size(); ++i)
      {
         if(sections[i] == sect)
         {
            add = false;
            break;
         }
      }
      
      //Add it if it wasn't found.
      if(add) sections.push_back(sect);
   }
   
   return 0;
}

inline
int appConfigurator::isSetUnused( const std::string & name )
{
   return isSet( name, m_unusedConfigs);
}

/// A simple config file writing function, useful for testing.
/** Write a config file to the path specified by `fname`.  Each of the parameters
  * is a vector, and each component is required of each vector.
  * 
  * Example:
  *\code
   writeConfigFile( "/tmp/test.conf", {"",     "",     "sect1", "sect1", "sect2", "sect2"},
                                      {"key0", "key1", "key2",  "key3",  "key4",  "key5"},
                                      {"val0", "val1", "val2",  "val3",  "val4",  "val5"} );
                                            
  *\endcode
  *results in the file `/tmp/test.conf' containing
  * \verbatim
  key0=val0
  key1=val1
  
  [sect1]
  key2=val2
  key3=val3
  
  [sect2]
  key4=val4
  key5=val5
  
  * \endverbatim
  * 
  * No error checking is done.
  */
inline
void writeConfigFile( const std::string & fname,                 ///< [in] the complete path for the output file.
                      const std::vector<std::string> & sections, ///< [in] sections, one per config value
                      const std::vector<std::string> & keywords, ///< [in] keywords, one per config value
                      const std::vector<std::string> & values    ///< [in] the values
                    )
{
   std::ofstream fout;
   
   fout.open(fname);
   
   std::string lastSection;
   for(size_t i=0; i< sections.size(); ++i)
   {
      std::string currSection = sections[i];
      
      if( currSection != lastSection && currSection != "")
      {
         fout << "\n[" << currSection << "]\n";
      }
         
      fout << keywords[i] << "=" << values[i] << "\n";
         
      lastSection = currSection;
   }
   
   fout.close();
   
   return;
}

} //namespace app
} //namespace mx

#endif // appConfigurator_hpp
