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

#ifndef app_clOptions_hpp
#define app_clOptions_hpp


#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

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
   
} //namespace app 
} //namespace mx


#endif //app_clOptions_hpp

