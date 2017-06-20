/** \file clOptions.hpp
 * \author Jared R. Males
 * \brief A command line parser
 *
 */

#ifndef __clOptions_hpp__
#define __clOptions_hpp__


#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

#include "optionparser.h"

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
struct clOptions
{
   std::unordered_map<std::string, unsigned int> map;
   std::unordered_map<std::string, unsigned int> typeMap;
   
   typedef std::unordered_map<std::string, unsigned int>::iterator mapIterator;

   std::vector<option::Descriptor> descriptions;
   
   option::Option *  options;
   option::Option* buffer;
   
   unsigned int nOpts;


   
   clOptions()
   {
      options = 0;
      buffer = 0;
      nOpts = 0;
   }

   ~clOptions()
   {
      if(options) delete options;
      if(buffer) delete buffer;
   }

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
      
   void add(const std::string & optName, const char * const shortOpt, const char * const longOpt, int argT)
   {
      mapIterator it;
      it = map.find(optName);

      
      if(it == map.end())
      {
         if(argT == argType::Optional)
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::Optional, ""});
         else if(argT == argType::Required)
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, Arg_Required, ""});
         else
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::None, ""});
         
         map.insert({optName, nOpts});
         typeMap.insert({optName, argT});
         
         ++nOpts;
         return;
      }
      
      if(argT == argType::Optional)
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::Optional, ""});
         else if(argT == argType::Required)
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, Arg_Required, ""});
         else
            descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::None, ""});
//      descriptions.push_back({it->second, argT, shortOpt, longOpt, option::Arg::Optional, ""});
   }

   void parse(int argc, char **argv, std::vector<std::string> * nonOptions = 0)
   {
      argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
      
      
      descriptions.push_back({0,0,0,0,0,0});
      
      option::Stats  stats(descriptions.data(), argc, argv);
      options = new option::Option[stats.options_max];
      buffer  = new option::Option[stats.buffer_max];
      
      option::Parser parse(descriptions.data(), argc, argv, options, buffer);
            
      if(nonOptions)
      {
         nonOptions->resize(parse.nonOptionsCount());
      
         for(int i=0;i<parse.nonOptionsCount(); ++i)
         {
            (*nonOptions)[i] = parse.nonOption(i);
         }
      }
   }
      
   const char * operator[](const std::string & key)
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
   
   int count(const std::string & key)
   {
      mapIterator it = map.find(key);
      return options[it->second].count();
   }
   
   ///Fill a vector of strings with the arguments supplied for key, last first. 
   void getAll(std::vector<std::string> & args, const std::string & key)
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
      
   bool optSet(const std::string & key)
   {
     mapIterator it = map.find(key);
     if(it == map.end()) return false; //Not found --> e.g. if neither command line short nor long option set.
     
     if( options[it->second].type() != argType::None) return true;
     return false;
   };
   
 
};


} //namespace mx


#endif //__clOptions_hpp__

