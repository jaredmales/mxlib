/** \file application.hpp
  * \author Jared R. Males
  * \brief A class for managing application configuration and execution
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

#ifndef app_application_hpp
#define app_application_hpp


#include "appConfigurator.hpp"
#include "../ioutils/textTable.hpp"

namespace mx
{
namespace app 
{
   
   
/// A class for managing application configuration and execution
/** Derived classes should implement at a minimum
  *
  * \code
    virtual void setupConfig();
    virtual void loadConfig();
    virtual int execute();
   \endcode
  *
  * These are executed in the order shown by the call to \ref main(). 
  *
  * This class uses a cascaded configuration system.  The application configuration is built up from the following sources, in increasing order of precedence:
  * - A global configuration file
  * - A user configuration file (which may be specified relative to the users home directory)
  * - A local configuration file (in the pwd)
  * - A configuration file specified on the command line
  * - The command line
  *
  * At each step in the above order, values read from that configuration source override any previous settings. So the command line has the highest precedence.
  * 
  * The configuration is set up and accessed  using an object of type  \ref appConfigurator named `config`.
  * For specification of the configuration file syntax and command line arguments see \ref appConfigurator.
  *
  * The configuration should be set up in the \ref setupConfig function.  This is normally done using
  * the appConfigurator::add method, as in:
  * \code
    void derived_class::setupConfig()
    {
       config.add("name", "s", "long", argType::Required, "section", "keyword", false, "int", "help message");
    }
    \endcode
  * The configuration is then accessed using the `config` member's operator as in
  \code
    void derived_class::loadConfig()
    {
         int val = 0;
         config(val, "name");
    }
  \endcode
  * In these examples, a config target with name "name" is created.  If the user sets a value of 2
  * via the command line with `-s 2`, `--long=2` or in a configuration file with
  \verbatim
  [section]
  keyword=2
  \endverbatim
  then val will be set to 2.  Otherwise, `val` will remain 0.
  *
  * Standard target `-h`/`--help` with name "help" to trigger printing a help message,
  * and `-c`/`--config` with name "config" to pass a config file via the command line. This
  * behavior can be changed by overriding functions.
  *
  * \note After loadConfig() but before execute(), the containers in \ref config are de-allocated , so they can not be used inside execute.
  *
  * A standard help message is produced when requested by the `-h`/`--help` option.  This behavior can be changed
  * by overriding the \ref help method.
  *
  * \ingroup mxApp
  */
class application
{

protected:

   std::string invokedName; ///< The name used to invoke this application.

   std::string m_configPathGlobal;     ///< The path to the gobal configuration file.  Set in constructor to use.
   std::string m_configPathGlobal_env; ///< Environment variable to check for the global config path.  Set in constructor to use.
   std::string m_configPathUser;       ///< The path to the user's configuration file.  If the first character is not '/' or '~', this is added to HOME. Set in constructor to use.
   std::string m_configPathUser_env;   ///< Environment variable to check fo the user config path.  Set in constructor to use.
   std::string m_configPathLocal;      ///< The path to a local configuration file. Set in constructor to use.
   std::string m_configPathLocal_env;  ///< Environment variable to check for the local config path.  Set in constructor to use.
   
   bool m_requireConfigPathLocal {true}; ///< Flag controlling whether lack of a local  configuration file should be reported.
   
   std::string m_configPathCL;         ///< The path to a configuration file specified on the command line. 
   std::string m_configPathCLBase;     ///< A base path to add to the CL path.  Set in constructor to use.
   std::string m_configPathCLBase_env; ///< Environment variable to check for the CL base path.  Set in constructor to use.
   
   appConfigurator config; ///< The structure used for parsing and storing the configuration.

   bool m_preserveConfig {false}; ///< Flag controlling whether the configuration is cleared before execution.  Set in derived constructor.
   
   bool doHelp {false}; ///< Flag to control whether the help message is printed or not.

   int m_helpWidth {120} ; ///< The total text width available for the help message.
   int m_helpSOColWidth {2}; ///< The width of the short option (-o) column in the help message.
   int m_helpLOColWidth {25}; ///< The width of the long option (--option) column in the help message.
   int m_helpCFColWidth {25}; ///< The width of the config file option column in the help message.
   int m_helpTypeColWidth {15}; ///< The width of the argument type column in the help message.

   int m_argc; ///< Store argc for later use. E.g. in reReadConfig().
   char ** m_argv; ///< Store argv for later use. E.g. in reReadConfig().

public:
   //application();

   virtual ~application();

   ///The application main function.
   /** Call this from the true main function, passing the command line arguments to be processed.
     * This calls \ref setup(), then checks if the doHelp flag was set.  If so, it calls \ref help() and returns.
     * If doHelp is not set, it then clears the config structure, and then calls \ref execute().
     * 
     * The configuration is cleared before the call to execute, unless m_preserveConfig = true.
     * 
     * \returns 1 if help is executed.
     * \returns -1 on error.
     * \returns the value of \ref execute() otherwise.
     */
   int main( int argc, ///< [in] standard command line result specifying number of argumetns in argv
             char **argv ///< [in] standard command line result containing the arguments.
           );

   ///Set the global configuration path
   void setConfigPathGlobal(const std::string & s /**< [in] the new path */);

   ///Set the user configuration path
   /** If the provided path does not star with a '/' or '~', then it will appended to
     * the users HOME directory obtained from the environment.
     */ 
   void setConfigPathUser(const std::string & s /**< [in] the new path to the user config file*/);

   ///Set the local configuration path
   void setConfigPathLocal( const std::string & s /**< [in] the new path */);

protected:

   /** \name Required Overrides
     * These methods should be implemented in derived classes to use an mx::application in its default behavior.
     * @{
     */

   ///In derived classes this is where the config targets are added to \ref config.
   virtual void setupConfig();

   ///Override this function to extract the configured values and set up the application.
   virtual void loadConfig();

   ///This function is where the derived class should do its work.  
   /** The application will exit with the return value of execute.
     */
   virtual int execute();

   ///@}

   /** \name Optional Overrides
     * These methods do not need to be implemented in derived classes unless it is desired to change behavior.
     * @{
     */

   ///Sets up the application by executing the configuration steps
   /** This is the key method which defines the mx::application configuration system.
     * This will not normally need to be implemented by derived clasess --
     * only do so if you intend to change the configuration process!
     */
   virtual void setup( int argc, ///< [in] standard command line result specifying number of argumetns in argv
                       char ** argv ///< [in] standard command line result containing the arguments.
                     );


   ///Re-read the config stack.
   /** This would be used if some config targets can only be constructed after 
     * a first pass.  Note that all previously read values will be appended as if
     * entered twice, so you must be careful to only access new targets
     * after calling this.
     *
     * \returns 0 on success.
     * \returns -1 on error.
     */ 
   int reReadConfig();
   
   ///Set the default search paths for config files
   /** In general you should not need to redefine this.
     *
     * The comand line arguments are not used by the default version, but are parameters in
     * case derived classes need access when overriding.
     *
     */
   virtual void setDefaults( int argc, ///< [in] standard command line result specifying number of arguments in argv
                             char ** argv ///< [in] standard command line result containing the arguments.
                           );
   

   ///Set up the command-line config option in a standard way.
   /** This adds "-c --config" as command line options.
     * You can override this function to change this behavior, or simply clear config
     * at the beginning of \ref setupConfig().  If you do this, you should also override \ref loadStandardConfig().
     */
   virtual void setupStandardConfig();

   ///Set up the help an options in a standard way.
   /** This adds "-h and --help" as command line options.
     * You can override this function to change this behavior, or simply clear config
     * at the beginning of \ref setupConfig().  If you do this, you should also override \ref loadStandardHelp().
     */
   virtual void setupStandardHelp();

   ///Loads the values of "config".
   /** Override this function if you do not want to use this, or have different behavior.
     * See also \ref setupStandardConfig().
     */
   virtual void loadStandardConfig();

   ///Loads the value of "help" into doHelp.
   /** Override this function if you do not want to use this, or have different behavior.
     * See also \ref setupStandardConfig().
     */
   virtual void loadStandardHelp();

   /// Format a configTarget for the help message.
   virtual void optionHelp( configTarget & tgt,     ///< [in] The Target
                            ioutils::textTable & tt ///< [out] the textTable being populated
                          );

   /// Setup a basic configuration.  Can be used in an intermediate derived class to allow its children to use setupConfig.
   /** This is called just before setupConfig().
     */
   virtual void setupBasicConfig();

   /// Load a basic configuration.  Can be used in an intermediate derived class to allow its children to use loadConfig.
   /** This is called just before loadConfig().
     */
   virtual void loadBasicConfig();
   
   /// Check the config.  This is called at the end of setup, before the configuration is cleared.
   /** It is up to you to decide how to handle the outcome.  If a bad config results in printing help,
     * you can set the doHelp flag.
     */
   virtual void checkConfig();
   
   ///Print a formatted help message, based on the config target inputs.
   virtual void help();

   ///@}

};





} //namespace app
} //namespace mx

#endif // app_application_hpp

//Deprecations:
// Produce errors if older code is compiled which attempts to define-config the config:

#ifdef MX_APP_DEFAULT_configPathGlobal
   #error MX_APP_DEFAULT_configPathGlobal no longer works.  Set m_configPathGlobal in constructor.
#endif

#ifdef MX_APP_DEFAULT_configPathGlobal_env
   #error MX_APP_DEFAULT_configPathGlobal_env no longer works.  Set m_configPathGlobal_env in constructor.
#endif
      
#ifdef MX_APP_DEFAULT_configPathUser
   #error MX_APP_DEFAULT_configPathUser no longer works.  Set m_configPathUser in constructor.
#endif

#ifdef MX_APP_DEFAULT_configPathUser_env
   #error MX_APP_DEFAULT_configPathUser_env no longer works.  Set m_configPathUser_env in constructor.
#endif

#ifdef MX_APP_DEFAULT_configPathLocal
   #error MX_APP_DEFAULT_configPathLocal no longer works.  Set m_configPathLocal in constructor.
#endif

#ifdef MX_APP_DEFAULT_configPathLocal_env
   #error MX_APP_DEFAULT_configPathLocal_env no longer works.  Set m_configPathLocal_env in constructor.
#endif

 #ifdef MX_APP_DEFAULT_configPathCLBase_env 
   #error MX_APP_DEFAULT_configPathCLBase_env no longer works.  Set m_configPathCLBase_env in constructor.
#endif
