/** \file clOptions.hpp
 * \author Jared R. Males
 * \brief A command line parser
 *
 * \ingroup mxApp_files
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

#ifndef clOptions_hpp__
#define clOptions_hpp__


#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

#include "../mxlib.hpp"

#include "optionparser/optionparser.h"

namespace mx
{
namespace app 
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

static option::ArgStatus Arg_Required(const option::Option& option, bool UNUSED(msg))
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
   
   unsigned int nOpts {0}; ///< The number of options added.

   /// D'tor.  Deletes the options and buffer pointers.
   ~clOptions();

   /// Clear the memory held by this object.
   void clear();
      
   /// Add a command line option target.
   void add( const std::string & optName, ///< [in] The name of the option, used as its key 
             const char * const shortOpt, ///< [in] The short option character 
             const char * const longOpt,  ///< [in] The long option keyword
             int argT                     ///< [in] The option type
           );

   ///Parse the command line
   void parse( int argc,                                 ///< [in] From main(argc, argv), the number of command line arguments
               char **argv,                              ///< in] From main(argc, argv), the command line arguments
               std::vector<std::string> * nonOptions = 0 ///< [out] [optional] the elements in argv which are not option or option-arguments.
             );
      
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
   const char * operator[](const std::string & key /**< [in] the key identifying the element */);
   
   ///Get the number of times the option was set on the command line.
   /**
     * \returns the number of times this option was found.
     */ 
   int count( const std::string & key /**< [in] the key identifying the element */ );
   
   ///Fill a vector of strings with the arguments supplied for key, last to first. 
   void getAll( std::vector<std::string> & args, ///< [out] will be resized and populated with the arguments.  Will be empty if no arguments specified.
                const std::string & key          ///< [in] the key identifying the element 
              );
   
   ///Test whether this option was set.
   /** 
     * \returns true if the options was set on the command line 
     * \returns false otherwise
     */ 
   bool optSet( const std::string & key /**< [in] the key identifying the element */ );
 
   /// Get the number of unknown options found by the parser
   /**
     * \returns the count of options found that do  no match an input option
     */ 
   int numUnknown();
   
   /// Get a vector of the unknown options found by the parser
   /**
     * \returns 0 on success (always)
     */ 
   int unknown( std::vector<std::string> & unk /**<[out] the vector to populate with the unknown options found*/);
};


inline
clOptions::~clOptions()
{
   if(options) delete[] options;
   if(buffer) delete[] buffer;
}

inline
void clOptions::clear()
{
   if(options) delete[] options;
   options = 0;
   
   if(buffer) delete[] buffer;
   buffer = 0;
   
   map.clear();
   typeMap.clear();
   descriptions.clear();
}
   
inline
void clOptions::add( const std::string & optName, 
                     const char * const shortOpt, 
                     const char * const longOpt,  
                     int argT                    
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
}

inline
void clOptions::parse( int argc,    
                       char **argv, 
                       std::vector<std::string> * nonOptions 
                     )
{
   argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
   
   //If not already done, we push the unknown catcher and the termination descriptor.
   if(descriptions.back().index != 0)
   {
      descriptions.push_back({nOpts, 0, "", "", option::Arg::None, ""}); //This is inserted to catch unknown options
      descriptions.push_back({0,0,0,0,0,0});
   }
   
   //Now allocate.
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
   
inline
const char * clOptions::operator[]( const std::string & key )
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

inline
int clOptions::count( const std::string & key )
{
   mapIterator it = map.find(key);
   
   if(it == map.end()) return -1;
   
   return options[it->second].count();
}
   
inline
void clOptions::getAll( std::vector<std::string> & args,
                        const std::string & key 
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
   
inline
bool clOptions::optSet( const std::string & key )
{
   mapIterator it = map.find(key);
  
   if(it == map.end()) return false; //Not found --> e.g. if neither command line short nor long option set.
  
   if( options[it->second].type() != argType::None) return true;
   return false;
};
   
inline   
int clOptions::numUnknown()
{
   //The dummy description to catch unknown options is inserted at position nOpts.
   return options[nOpts].count();
}
   
inline
int clOptions::unknown( std::vector<std::string> & unk )
{
   unk.clear();
   
   if(numUnknown() == 0) return 0;
   
   for (option::Option* opt = options[nOpts]; opt != NULL; opt = opt->next())
   {
      unk.push_back( opt->name );
   }
   
   return 0;
}
   
} //namespace app 
} //namespace mx


#endif //clOptions_hpp__

