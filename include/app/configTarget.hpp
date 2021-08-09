/** \file configTarget.hpp
 * \author Jared R. Males
 * \brief Targets for the configuration manager, and utiltities.
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

#ifndef configTarget_hpp
#define configTarget_hpp

#include <unordered_map>


#include "../ioutils/stringUtils.hpp"
#include "../meta/trueFalseT.hpp"
#include "../meta/typeTraits.hpp"
#include "../meta/typeDescription.hpp"

#include "../mxError.hpp"

namespace mx
{
namespace app 
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
   int clType {0}; ///< The command-line option type, argType::false, argType::true, argType::optional, argType::required
   std::string section; ///< The config file section name, can be empty ""
   std::string keyword; ///< The config file keyword, read in a "keyword=value" pair
   bool set {false}; ///< true if the value has been set by the configuration, use to distinguish empty strings


   bool isRequired {false}; ///< Whether or not this is option is required to be set.
   std::string helpType;         ///< The type to display in the help message.
   std::string helpExplanation;  ///< The explanation to display in the help message.


   std::vector<std::string> values; ///< holds the values in the order they are set by the configuration
   std::vector<std::string> sources; ///< holds the sources of the values (command line or config file path)
   
   int verbosity {0}; ///< Records the verbosity of command line options.  E.g. for -v:1, -vv:2, -vvv:3 etc. 
   int orderAdded {0}; ///< The order in which this was added.  Useful for displaying help messages.

   bool used {false};
   
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
               ) : name(n), shortOpt(so), longOpt(lo), clType (clt), section(s), keyword(kw), isRequired(isReq), helpType(ht), helpExplanation(he), orderAdded(oa)
   {
   }
};




} //namespace app 
} //namespace mx

#endif // appConfigurator_hpp
