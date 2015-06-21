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
         descriptions.push_back({nOpts, argT, shortOpt, longOpt, option::Arg::Optional, ""});
         map.insert({optName, nOpts});
         typeMap.insert({optName, argT});
         
         ++nOpts;
         return;
      }
      
      descriptions.push_back({it->second, argT, shortOpt, longOpt, option::Arg::Optional, ""});
   }

   void parse(int argc, char **argv)
   {
      argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
      
      
      descriptions.push_back({0,0,0,0,0,0});
      
      option::Stats  stats(descriptions.data(), argc, argv);
      options = new option::Option[stats.options_max];
      buffer  = new option::Option[stats.buffer_max];
      
      option::Parser parse(descriptions.data(), argc, argv, options, buffer);
      
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
      
      return options[it->second].arg;
   }
   
   bool optSet(const std::string & key)
   {
     mapIterator it = map.find(key);
     
     if( options[it->second].type() != argType::None) return true;
     
     return false;
   };
   
 
};


} //namespace mx


#endif //__clOptions_hpp__

