/** \file clOptions.hpp
 * \author Jared R. Males
 * \brief A command line parser
 *
 * \ingroup mxApp_files
 */

#ifndef __clOptions_hpp__
#define __clOptions_hpp__


#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

#include "optionparser/optionparser.h"

namespace mx
{

   
/// Argument types
enum argType
{
   None,
   False,
   True,
   Optional,
   Required
};

static const char * falseStr = "false";
static const char * trueStr = "true";
static const char * blankStr = "";

static option::ArgStatus Arg_Required(const option::Option& option, bool msg)
{
   if (option.arg != 0)
   return option::ARG_OK;

   return option::ARG_ILLEGAL;
}
   

/// Command line options parser.
/** This is a wrapper for the <a href="https://sourceforge.net/projects/optionparser/">"The Lean Mean C++ Option Parser"</a> command line parser.  
  * 
  * \ingroup mxApp
  */ 
struct clOptions
{
   std::unordered_map<std::string, unsigned int> map;
   std::unordered_map<std::string, unsigned int> typeMap;
   
   typedef std::unordered_map<std::string, unsigned int>::iterator mapIterator;

   std::vector<option::Descriptor> descriptions;
   
   option::Option * options {nullptr};
   option::Option * buffer {nullptr};
   
   unsigned int nOpts {0};


   /// D'tor.  Deletes the options and buffer pointers.
   ~clOptions()
   {
      if(options) delete options;
      if(buffer) delete buffer;
   }

   /// Clear the memory held by this object.
   void clear()
   {
      if(options) delete options;
      options = 0;
      
      if(buffer) delete buffer;
      buffer = 0;
      
      map.clear();
      typeMap.clear();
      descriptions.clear();
   }
      
   /// Add a command line option target.
   void add( const std::string & optName, ///< [in] The name of the option, used as its key 
             const char * const shortOpt, ///< [in] The short option character 
             const char * const longOpt,  ///< [in] The long option keyword
             int argT                     ///< [in] The option type
           )
   {
      mapIterator it;
      it = map.find(optName);

      
      if(it == map.end())
      {
         if(argT == argType::Optional)
         {
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::Optional, ""});
         }
         else if(argT == argType::Required)
         {
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, Arg_Required, ""});
         }
         else
         {
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::None, ""});
         }
         
         map.insert({optName, nOpts});
         typeMap.insert({optName, argT});
         
         ++nOpts;
         return;
      }
      
//       if(argT == argType::Optional)
//             descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::Optional, ""});
//          else if(argT == argType::Required)
//             descriptions.push_back({nOpts, argT, shortOpt, longOpt, Arg_Required, ""});
//          else
//             descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::None, ""});
//      descriptions.push_back({it->second, argT, shortOpt, longOpt, option::Arg::Optional, ""});
   }

   ///Parse the command line
   void parse( int argc,   ///< [in] From main(argc, argv), the number of command line arguments
               char **argv, ///< in] From main(argc, argv), the command line arguments
               std::vector<std::string> * nonOptions = 0 ///< [out] [optional] the elements in argv which are not option or option-arguments.
             )
   {
      argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
      
      descriptions.push_back({0,0,0,0,0,0});
      
      option::Stats  stats(descriptions.data(), argc, argv);
      options = new option::Option[stats.options_max];
      buffer  = new option::Option[stats.buffer_max];
      
      option::Parser parse(false, descriptions.data(), argc, argv, options, buffer);
            
      if(nonOptions)
      {
         nonOptions->resize(parse.nonOptionsCount());
      
         for(int i=0;i<parse.nonOptionsCount(); ++i)
         {
            (*nonOptions)[i] = parse.nonOption(i);
         }
      }
   }
      
   ///Get the value of the option, if present.
   /**
     * For a true/false type, returns true or false as appropriate.
     * Returns an empty string if no argument for the key. 
     * Otherwise return the argument.
     * 
     * \returns 0 on error
     * \returns an empty string on no argument
     * \returns a string otherwise
     */
   const char * operator[](const std::string & key /**< [in] the key identifying the element */)
   {
      mapIterator it = map.find(key);

      if(it == map.end())
      {
         std::cerr << "oh no\n";
         return 0;
      }

      int typeDesc = typeMap[key]; 

      //If this option is not pressent, either return the opposite T/F condition or blank
      if(options[it->second].type() == argType::None) 
      {
         if(typeDesc == argType::False) return trueStr;
         if(typeDesc == argType::True) return falseStr;
         return blankStr;
      }
      
      if(typeDesc == argType::False || typeDesc == argType::True)
      {
         if(options[it->second].last()->type() == argType::False) return falseStr;
         else return trueStr;
      }

      if(options[it->second].arg == 0) return blankStr;
      
      return options[it->second].last()->arg;
   }
   
   ///Get the number of times the option was set on the command line.
   int count(const std::string & key /**< [in] the key identifying the element */)
   {
      mapIterator it = map.find(key);
      
      if(it == map.end()) return -1;
      
      return options[it->second].count();
   }
   
   ///Fill a vector of strings with the arguments supplied for key, last first. 
   void getAll( std::vector<std::string> & args, ///< [out] will be resized and populated with the arguments.  Will be empty if no arguments specified.
                const std::string & key ///< [in] the key identifying the element 
              )
   {
      mapIterator it = map.find(key);
      
      
      if(it == map.end())
      {
         std::cerr << "oh no\n";
         return;
      }

      int typeDesc = typeMap[key]; 

      //If this option is not present, either return the opposite T/F condition or blank
      if(options[it->second].type() == argType::None) 
      {
         
         if(typeDesc == argType::False) 
         {
            args.resize(1);
            args[0]= trueStr;
            return;
         }
         
         if(typeDesc == argType::True) 
         {
            args.resize(1);
            args[0]= falseStr;
            return;
         }
         
         args.clear();
         
         return;
      }
      
      if(typeDesc == argType::False || typeDesc == argType::True)
      {
         args.resize(1);
         
         if(options[it->second].last()->type() == argType::False) args[0] = falseStr;
         else args[0] = trueStr;         
         return;
      }
      
      if(options[it->second].arg == 0) 
      {
         args.clear();
         return;
      }
      
      int N = options[it->second].count();
      
      args.resize(N);
      
      int i=0;
      
      for (option::Option* opt = options[it->second]; opt != NULL && i < N; opt = opt->next())
      {
         args[i] = opt->arg;
         ++i;
      }
   }
      
   ///Test whether this option was set.
   bool optSet(const std::string & key /**< [in] the key identifying the element */)
   {
     mapIterator it = map.find(key);
     
     if(it == map.end()) return false; //Not found --> e.g. if neither command line short nor long option set.
     
     if( options[it->second].type() != argType::None) return true;
     return false;
   };
   
 
};


} //namespace mx


#endif //__clOptions_hpp__

