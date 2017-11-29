/** \file appConfigurator.hpp
 * \author Jared R. Males
 * \brief An application configuration manager
 * 
 * \ingroup mxApp_files
 *
 */

#ifndef __appConfigurator_hpp__
#define __appConfigurator_hpp__

#include "../stringUtils.hpp"
#include "../meta/trueFalseT.hpp"
#include "../meta/typeTraits.hpp"

#include "clOptions.hpp"
#include "iniFile.hpp"


#include <list>

namespace mx
{
   
/// A configuration target
/** Specifies the details of a configuration target, which is a value that can be set from the command line and/or a config file.
  * A target has a name used as a key for accessing it, and defines how it is set with short and long command line options and the section and
  * key for the config file.  Can also include help message details.
  * 
  * \ingroup mxApp
  */ 
struct configTarget
{
   std::string name; ///<The name of the target
   std::string shortOpt; ///< The command-line short option (e.g. "f" for -f)
   std::string longOpt; ///< The command-line long option (e.g. "file" for --file)
   int clType; ///< The command-line option type, argType::false, argType::true, argType::optional, argType::required
   std::string section; ///< The config file section name, can be empty ""
   std::string keyword; ///< The config file keyword, read in a "keyword=value" pair
   bool set; ///< true if the value has been set by the configuration, use to distinguish empty strings
   

   bool isRequired; ///< Whether or not this is option is required to be set.
   std::string helpType;         ///< The type to display in the help message.
   std::string helpExplanation;  ///< The explanation to display in the help message.
   
   
   std::vector<std::string> values; ///< holds the values in the order they are set by the configuration
   
   int orderAdded; ///< The order in which this was added.  Useful for displaying help messages.
   
   /// Default c'tor
   configTarget()
   {
   }

   /// Construct and set values
   configTarget( const std::string &n,  ///< [in] The name of the target
                 const std::string &so, ///< [in] The command-line short option (e.g. "f" for -f)
                 const std::string &lo, ///< [in] The command-line long option (e.g. "file" for --file)
                 int clt,               ///< [in] The command-line option type, argType::false, argType::true, argType::optional, argType::required
                 const std::string & s, ///< [in] The config file section name, can be empty ""
                 const std::string & kw, ///< [in] The config file keyword, read in a "keyword=value" pair
                 bool isReq = false, ///< [in] Whether or not this is option is required to be set
                 const std::string & ht = "",  ///< [in] The type to display in the help message
                 const std::string & he = "",  ///< [in] The explanation to display in the help message
                 int oa = 0 ///< [in] [optional] ///< the order in which this was added. 
               )
   {
      name = n;
      shortOpt = so;
      longOpt = lo;
      clType = clt;
      section = s;
      keyword = kw;
      
      isRequired = isReq;
      helpType = ht;
      helpExplanation = he;
      
      set = false;
   }
};

/// Class to manage a set of configurable values, and read their values from config/ini files and the command line.
/** \bug a config=value pair listed in a conf file twice seems to cause a failure, even if they are the same value.
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
   std::unordered_map<std::string, configTarget> targets;
   
   ///Targets which are only for the command line are stored separately in a list.
   std::list<configTarget> clOnlyTargets;
   
   
   std::vector<std::string> nonOptions;
      
   int nAdded;
   
   appConfigurator()
   {
      nAdded = 0;
   }
   
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
   void add( const std::string &n,  ///< [in] The name of the target
             const std::string &so, ///< [in] The command-line short option (e.g. "f" for -f)
             const std::string &lo, ///< [in] The command-line long option (e.g. "file" for --file)
             int clt,               ///< [in] The command-line option type, argType::false, argType::true, argType::optional, argType::required
             const std::string & s, ///< [in] The config file section name, can be empty ""
             const std::string & kw, ///< [in] The config file keyword, read in a "keyword=value" pair
             bool isReq = false, ///< [in] Whether or not this is option is required to be set
             const std::string & ht = "",  ///< [in] The type to display in the help message
             const std::string & he = ""  ///< [in] The explanation to display in the help message
           );
      
   ///Parse the command line, updating the targets
   void parseCommandLine( int argc,  ///< [in] standard command line result specifying number of argumetns in argv
                          char ** argv, ///< [in] standard command line result containing the arguments.
                          const std::string & oneTarget = "" ///< [in] [optional] if not empty, then only this target is extracted by the parser.
                        );
   

   ///Parse a config/ini file, updating the targets
   /** \todo handle += here, by appending to the last value as if a vector.
     */
   void readConfig(const std::string & fname /**< [in] the config file name */);
      
   /// Check if a target has been set by the configuration
   /**
     * \returns true if the configuration set at least one value for this target
     * \returns false if no value was set.
     */ 
   bool isSet(const std::string & name /**< [in] the target name */);
 
   ///Get the number of different values set for the specified config target
   /** 
     * \returns the number of different values set for name.
     */ 
   int count( const std::string & name /**< [in] the target name */);

   
   /// Get the i-th value of the target, converted to the specified type
   /** The supplied value is only altered if the config target was set, which preserves.  
     * default values.
     * 
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename typeT>
   int get( typeT & v, ///< [out] the variable to store the value in
            const std::string & name, ///< [in] the config target name
            int i ///< [in] the number of config specification to get.
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
   int get( typeT & v, ///< [out] the variable to store the value in
            const std::string & name ///< [in] the config target name
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
   int get( std::vector<typeT> & v,   ///< [out] the vector to populate
            const std::string & name ///< [in[ the config target name.
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
            int i ///< [in] the number of config specification to get.
          );
   
  
   
   template<typename typeT>
   int operator()( typeT & v,
                   const std::string & name
                 );

private:
   //Tag-dispatching: These overloads provide different behavior depending on whether typeT is a vector, effectively 
   // overloading get(name) based on return type.   
   template<typename typeT>
   typeT internalGet( const std::string & name, //[in] the target name
                      meta::trueFalseT<false> isNotVector //[in] the false type makes this the overload for non-std::vector typeT
                    );
   
   template<typename typeT>
   typeT internalGet( const std::string & name, //[in] the target name
                      meta::trueFalseT<true> isVector //[in] the true type makes this the overload for typeT == std::vector
                    );

public:   
   /// Get the final value of the target, returned as the specified type
   /** The default-constructed/initialized value is returned if \ref isSet() == false, so this should carefully be used
     * in conjunction with \ref isSet() to avoid overwriting default values.
     * 
     * This works for vectors as well (through tag dispatching).  In any case, you must explicitly specify typeT, as in 
     * \code
       float f = get<float>(name);
       std::vector<float> v = get<std::vector<float>>(name); 
       \endcode
     * 
     * 
     * \retval 0 on success.
     * \retval -1 on error.
     */
   template<typename typeT>
   typeT get( const std::string & name /**< [in] the config target name */ );
   
   
};

inline
void appConfigurator::clear()
{
   targets.clear();
   clOnlyTargets.clear();
}
   
inline
void appConfigurator::add( const configTarget & tgt )
{
   
   //First check for duplicate name and command line only
   if(targets.count(tgt.name) > 0 && tgt.section == "" && tgt.keyword == "")
   {
      clOnlyTargets.push_back(tgt);
      clOnlyTargets.back().orderAdded = nAdded;
   }
   else
   {
      std::pair<targetIterator, bool> res = targets.insert({tgt.name, tgt});
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
   for(it = targets.begin(); it != targets.end(); ++it)
   {
      if(it->second.shortOpt == "" && it->second.longOpt == "") continue; //No command line opts specified.
      
      if(oneTarget != "" && it->second.name != oneTarget) continue;
      
      clOpts.add(it->second.name,it->second.shortOpt.c_str(),it->second.longOpt.c_str(), it->second.clType);
   }
   
   //Then load the command-line-only options.
   clOnlyTargetIterator cloit;
   for(cloit = clOnlyTargets.begin(); cloit != clOnlyTargets.end(); ++cloit)
   {
      if(cloit->shortOpt == "" && cloit->longOpt == "") continue; //Nothing to add?
      
      if(oneTarget != "" && cloit->name != oneTarget) continue;
      
      clOpts.add(cloit->name,cloit->shortOpt.c_str(),cloit->longOpt.c_str(), cloit->clType);
   }

   //If we get here and nothing to do, get out
   if(clOpts.nOpts == 0) 
   {
      return;
   }

   //Now we parse
   clOpts.parse(argc, argv, &nonOptions);

   //And then load the results in the config target map.
   for(it = targets.begin(); it != targets.end(); ++it)
   {
      if(clOpts.optSet(it->second.name)) 
      {
         std::vector<std::string> args;

         clOpts.getAll(args, it->second.name);
         it->second.values.insert( it->second.values.end(), args.begin(), args.end());
         it->second.set = true;
      }
   }      
   
}   

void appConfigurator::readConfig(const std::string & fname)
{
   iniFile iF;
   
   iF.parse(fname);
   
   targetIterator it;
   
   for(it = targets.begin(); it != targets.end(); ++it)
   {
      if(iF.count(it->second.section, it->second.keyword) > 0) 
      {
         it->second.values.push_back(iF(it->second.section, it->second.keyword));
         it->second.set = true;
      }
   }
   
}
   
inline
bool appConfigurator::isSet(const std::string & name)
{
   if(targets.count(name) == 0) return false;
   return targets[name].set;
}
   
inline
int appConfigurator::count(const std::string & name)
{
   return targets[name].values.size();
}

template<typename typeT>
int appConfigurator::get( typeT & v,
                          const std::string & name, 
                          int i
                        )
{
   if(!isSet(name)) return 0;

   if( targets[name].values.size() <= i) return -1;
   
   v = convertFromString<typeT>(targets[name].values[i]);
   
   return 0;
}

template<typename typeT>
int appConfigurator::get( typeT & v, 
                          const std::string & name)
{
      
   if(!isSet(name)) return 0;
   
   int i = targets[name].values.size() - 1;
   
   if(i < 0) return -1;
   
   return get(v, name, i);
}

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name,
                          int i
                        )
{
   
   if(!isSet(name)) return 0;
   
   if( targets[name].values.size() <= i) return -1;
   
   std::string s;
   if( get<std::string>(s, name, i) < 0) return -1;
   
   //Case that s was set to be empty.   
   if(s.size() == 0)
   {
      v.clear();
      return 0;
   }
   
   int st;
   int com;

   st = 0;
   com = s.find(',', st);
   
   v.clear();
      
   while(com != std::string::npos)
   {
      v.push_back( convertFromString<typeT>(s.substr(st, com-st)) );
      st = com + 1;
      com = s.find(',', st);
   }
   v.push_back( convertFromString<typeT>(s.substr(st, s.size()-st)));
   
   return 0;
}

template<typename typeT>
int appConfigurator::get( std::vector<typeT> & v,
                          const std::string & name
                        )
{
   
   if(!isSet(name)) return 0;
   
   int i = targets[name].values.size() - 1;
   
   if(i < 0) return -1;
   
   return get(v, name, i);
}
   

template<typename typeT>
int appConfigurator::operator()( typeT & v,
                                 const std::string & name
                               )
{
   return get(v, name);
}    




template<typename typeT>
typeT appConfigurator::internalGet( const std::string & name,
                                    meta::trueFalseT<false> isNotVector
                                  )
{
   static_cast<void>(isNotVector);
   
   if(!isSet(name)) return typeT(); //this zero-initializes any fundamental types, calls default constructor for class-types.
   
   typeT v;
   get(v, name);
   
   return v;
}
   
template<typename typeT>
typeT appConfigurator::internalGet( const std::string & name,
                                    meta::trueFalseT<true> isVector
                                  )
{
   static_cast<void>(isVector);
   
   typeT v;
   
   get(v, name); //Will only be modified if isSet == true.
   
   return v;
}
   
template<typename typeT>
typeT appConfigurator::get( const std::string & name )
{
   return internalGet<typeT>(name, meta::trueFalseT< meta::is_std_vector<typeT>::value >());
}


} //namespace mx

#endif // __appConfigurator_hpp__
