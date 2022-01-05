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

#ifndef app_appConfigurator_hpp
#define app_appConfigurator_hpp

#include <list>
#include <fstream>

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
  * The configuration files are TOML/ini-style, with sections.  That is
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
  * \todo should just swith to strict TOML
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
             bool isReq = false,          ///< [in] Whether or not this option is required to be set
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
    
   /// Get the number of unknown options found during config processing.
   int numUnknownOptions();
};

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name,
                          size_t i,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   targets[name].used = true; //This means this was checked.
   
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
//explicits
extern template
int appConfigurator::get<char>( char & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<signed char>( signed char & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<short>( short & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<int>( int & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long>( long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long long>( long long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<float>( float & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<double>( double & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long double>( long double & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);
#endif

extern template
int appConfigurator::get<bool>( bool & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<std::string>( std::string & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name,
                          size_t i
                        )
{
   return get(v, name, i, m_targets);
}

//explicits
extern template
int appConfigurator::get<char>( char & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<signed char>( signed char & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<short>( short & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<int>( int & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<long>( long & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<long long>( long long & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<float>( float & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<double>( double & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<long double>( long double & v, const std::string & name, size_t i);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name, size_t i);
#endif

extern template
int appConfigurator::get<bool>( bool & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<std::string>( std::string & v, const std::string & name, size_t i);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   targets[name].used = true; //This means this was checked.
   
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
//explicits:

extern template
int appConfigurator::get<char>( char & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<signed char>( signed char & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<short>( short & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<int>( int & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long>( long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long long>( long long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<float>( float & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<double>( double & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long double>( long double & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);
#endif

extern template
int appConfigurator::get<bool>( bool & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<std::string>( std::string & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name)
{
   return get(v, name, m_targets);
}
//explicits:

extern template
int appConfigurator::get<char>( char & v, const std::string & name);

extern template
int appConfigurator::get<char16_t>( char16_t & v, const std::string & name);

extern template
int appConfigurator::get<char32_t>( char32_t & v, const std::string & name);

extern template
int appConfigurator::get<wchar_t>( wchar_t & v, const std::string & name);

extern template
int appConfigurator::get<signed char>( signed char & v, const std::string & name);

extern template
int appConfigurator::get<short>( short & v, const std::string & name);

extern template
int appConfigurator::get<int>( int & v, const std::string & name);

extern template
int appConfigurator::get<long>( long & v, const std::string & name);

extern template
int appConfigurator::get<long long>( long long & v, const std::string & name);

extern template
int appConfigurator::get<unsigned char>( unsigned char & v, const std::string & name);

extern template
int appConfigurator::get<unsigned short>( unsigned short & v, const std::string & name);

extern template
int appConfigurator::get<unsigned int>( unsigned int & v, const std::string & name);

extern template
int appConfigurator::get<unsigned long>( unsigned long & v, const std::string & name);

extern template
int appConfigurator::get<unsigned long long>( unsigned long long & v, const std::string & name);

extern template
int appConfigurator::get<float>( float & v, const std::string & name);

extern template
int appConfigurator::get<double>( double & v, const std::string & name);

extern template
int appConfigurator::get<long double>( long double & v, const std::string & name);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( __float128 & v, const std::string & name);
#endif

extern template
int appConfigurator::get<bool>( bool & v, const std::string & name);

extern template
int appConfigurator::get<std::string>( std::string & v, const std::string & name);

//+++++++++++++++++++++++++++++++++++++
template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          size_t i,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   targets[name].used = true; //This means this was checked.
   
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
//explicits:

extern template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);
#endif

extern template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name, size_t i, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          size_t i
                        )
{
   return get(v, name, i, m_targets);
}

//explicits:

extern template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name, size_t i);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name, size_t i);
#endif

extern template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name, size_t i);

extern template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name, size_t i);

//+++++++++++++++++++++++++++++++++++++
template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          std::unordered_map<std::string, configTarget> & targets
                        )
{
   targets[name].used = true; //This means this was checked.
   
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
//explicits:

extern template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);
#endif

extern template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

extern template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name, std::unordered_map<std::string, configTarget> & targets);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name
                        )
{
   return get(v, name, m_targets);
}
//explicits:

extern template
int appConfigurator::get<char>( std::vector<char> & v, const std::string & name);

extern template
int appConfigurator::get<char16_t>( std::vector<char16_t> & v, const std::string & name);

extern template
int appConfigurator::get<char32_t>( std::vector<char32_t> & v, const std::string & name);

extern template
int appConfigurator::get<wchar_t>( std::vector<wchar_t> & v, const std::string & name);

extern template
int appConfigurator::get<signed char>( std::vector<signed char> & v, const std::string & name);

extern template
int appConfigurator::get<short>( std::vector<short> & v, const std::string & name);

extern template
int appConfigurator::get<int>( std::vector<int> & v, const std::string & name);

extern template
int appConfigurator::get<long>( std::vector<long> & v, const std::string & name);

extern template
int appConfigurator::get<long long>( std::vector<long long> & v, const std::string & name);

extern template
int appConfigurator::get<unsigned char>( std::vector<unsigned char> & v, const std::string & name);

extern template
int appConfigurator::get<unsigned short>( std::vector<unsigned short> & v, const std::string & name);

extern template
int appConfigurator::get<unsigned int>( std::vector<unsigned int> & v, const std::string & name);

extern template
int appConfigurator::get<unsigned long>( std::vector<unsigned long> & v, const std::string & name);

extern template
int appConfigurator::get<unsigned long long>( std::vector<unsigned long long> & v, const std::string & name);

extern template
int appConfigurator::get<float>( std::vector<float> & v, const std::string & name);

extern template
int appConfigurator::get<double>( std::vector<double> & v, const std::string & name);

extern template
int appConfigurator::get<long double>( std::vector<long double> & v, const std::string & name);

#ifdef HASQUAD
extern template
int appConfigurator::get<__float128>( std::vector<__float128> & v, const std::string & name);
#endif

extern template
int appConfigurator::get<bool>( std::vector<bool> & v, const std::string & name);

extern template
int appConfigurator::get<std::string>( std::vector<std::string> & v, const std::string & name);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::operator()( typeT & v,
                                 const std::string & name
                               )
{
   return get(v, name, m_targets);
}

//explicits:

extern template
int appConfigurator::operator()<char>( char & v, const std::string & name);

extern template
int appConfigurator::operator()<char16_t>( char16_t & v, const std::string & name);

extern template
int appConfigurator::operator()<char32_t>( char32_t & v, const std::string & name);

extern template
int appConfigurator::operator()<wchar_t>( wchar_t & v, const std::string & name);

extern template
int appConfigurator::operator()<signed char>( signed char & v, const std::string & name);

extern template
int appConfigurator::operator()<short>( short & v, const std::string & name);

extern template
int appConfigurator::operator()<int>( int & v, const std::string & name);

extern template
int appConfigurator::operator()<long>( long & v, const std::string & name);

extern template
int appConfigurator::operator()<long long>( long long & v, const std::string & name);

extern template
int appConfigurator::operator()<unsigned char>( unsigned char & v, const std::string & name);

extern template
int appConfigurator::operator()<unsigned short>( unsigned short & v, const std::string & name);

extern template
int appConfigurator::operator()<unsigned int>( unsigned int & v, const std::string & name);

extern template
int appConfigurator::operator()<unsigned long>( unsigned long & v, const std::string & name);

extern template
int appConfigurator::operator()<unsigned long long>( unsigned long long & v, const std::string & name);

extern template
int appConfigurator::operator()<float>( float & v, const std::string & name);

extern template
int appConfigurator::operator()<double>( double & v, const std::string & name);

extern template
int appConfigurator::operator()<long double>( long double & v, const std::string & name);

#ifdef HASQUAD
extern template
int appConfigurator::operator()<__float128>( __float128 & v, const std::string & name);
#endif

extern template
int appConfigurator::operator()<bool>( bool & v, const std::string & name);

extern template
int appConfigurator::operator()<std::string>( std::string & v, const std::string & name);

//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::configUnused( typeT & v,
                                   const std::string & key
                                 )
{
   return get(v, key, m_unusedConfigs);
}

extern template
int appConfigurator::configUnused<char>( char & v, const std::string & key);

extern template
int appConfigurator::configUnused<char16_t>( char16_t & v, const std::string & key);

extern template
int appConfigurator::configUnused<char32_t>( char32_t & v, const std::string & key);

extern template
int appConfigurator::configUnused<wchar_t>( wchar_t & v, const std::string & key);

extern template
int appConfigurator::configUnused<signed char>( signed char & v, const std::string & key);

extern template
int appConfigurator::configUnused<short>( short & v, const std::string & key);

extern template
int appConfigurator::configUnused<int>( int & v, const std::string & key);

extern template
int appConfigurator::configUnused<long>( long & v, const std::string & key);

extern template
int appConfigurator::configUnused<long long>( long long & v, const std::string & key);

extern template
int appConfigurator::configUnused<unsigned char>( unsigned char & v, const std::string & key);

extern template
int appConfigurator::configUnused<unsigned short>( unsigned short & v, const std::string & key);

extern template
int appConfigurator::configUnused<unsigned int>( unsigned int & v, const std::string & key);

extern template
int appConfigurator::configUnused<unsigned long>( unsigned long & v, const std::string & key);

extern template
int appConfigurator::configUnused<unsigned long long>( unsigned long long & v, const std::string & key);

extern template
int appConfigurator::configUnused<float>( float & v, const std::string & key);

extern template
int appConfigurator::configUnused<double>( double & v, const std::string & key);

extern template
int appConfigurator::configUnused<long double>( long double & v, const std::string & key);

#ifdef HASQUAD
extern template
int appConfigurator::configUnused<__float128>( __float128 & v, const std::string & key);
#endif

extern template
int appConfigurator::configUnused<bool>( bool & v, const std::string & key);

extern template
int appConfigurator::configUnused<std::string>( std::string & v, const std::string & key);
//+++++++++++++++++++++++++++++++++++++

template<typename typeT>
int appConfigurator::configUnused( typeT & v,
                                   const std::string & section,
                                   const std::string & keyword
                                 )
{
   return configUnused(v, iniFile::makeKey(section,keyword));
}

extern template
int appConfigurator::configUnused<char>( char & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<char16_t>( char16_t & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<char32_t>( char32_t & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<wchar_t>( wchar_t & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<signed char>( signed char & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<short>( short & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<int>( int & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<long>( long & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<long long>( long long & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<unsigned char>( unsigned char & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<unsigned short>( unsigned short & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<unsigned int>( unsigned int & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<unsigned long>( unsigned long & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<unsigned long long>( unsigned long long & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<float>( float & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<double>( double & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<long double>( long double & v, const std::string & section, const std::string & keyword);

#ifdef HASQUAD
extern template
int appConfigurator::configUnused<__float128>( __float128 & v, const std::string & section, const std::string & keyword);
#endif

extern template
int appConfigurator::configUnused<bool>( bool & v, const std::string & section, const std::string & keyword);

extern template
int appConfigurator::configUnused<std::string>( std::string & v, const std::string & section, const std::string & keyword);

//+++++++++++++++++++++++++++++++++++++

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
void writeConfigFile( const std::string & fname,                 ///< [in] the complete path for the output file.
                      const std::vector<std::string> & sections, ///< [in] sections, one per config value
                      const std::vector<std::string> & keywords, ///< [in] keywords, one per config value
                      const std::vector<std::string> & values    ///< [in] the values
                    );

} //namespace app
} //namespace mx

#endif // app_appConfigurator_hpp
