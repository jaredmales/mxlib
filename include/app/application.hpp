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

#ifndef application_hpp__
#define application_hpp__


#include "../mxlib.hpp"

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
  * - A user configuration file (specified relative to the users home directory)
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

   std::string configPathGlobal; ///< The path to the gobal configuration file.
   std::string configPathUser; ///< The path to the user's configuration file.
   std::string configPathLocal; ///< The path to a local configuration file.
   bool m_requireConfigPathLocal {true}; ///< Flag controlling whether lack of a configuration file should be reported.
   
   std::string m_configPathCL; ///< The path to a configuration file specified on the command line.
   std::string m_configPathCLBase; ///< A base path to add to the CL path.  Can be set by environment variable defined in MX_APP_DEFAULT_configPathCLBase_env.
   
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
   void setConfigPathUser(const std::string & s /**< [in] the new path */);

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
   virtual void setDefaults( int argc, ///< [in] standard command line result specifying number of argumetns in argv
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

   ///Print a formatted help message, based on the config target inputs.
   virtual void help();

   ///@}

};



inline
application::~application()
{
   return;
}

inline
int application::main( int argc,
                       char **argv
                     )
{
   m_argc = argc;
   m_argv = argv;
   
   setup(argc, argv);

   if(doHelp)
   {
      help();
      return 1;
   }

   if(!m_preserveConfig)
   {
      config.clear();
   }
   
   return execute();
}

inline
void application::setConfigPathGlobal(const std::string & s)
{
   configPathGlobal = s;
}

inline
void application::setConfigPathUser(const std::string & s)
{
   configPathUser = s;
}

inline
void application::setConfigPathLocal( const std::string & s )
{
   configPathLocal = s;
}

inline
void application::setupConfig() //virtual
{
   return;
}

inline
void application::loadConfig() //virtual
{
   return;
}

inline
int application::execute() //virtual
{
   return 0;
}



inline
void application::setup( int argc,
                         char ** argv
                       )
{
   invokedName = argv[0];

   setupStandardConfig();
   setupStandardHelp();

   setupBasicConfig();
   setupConfig();

   setDefaults(argc, argv);

   config.readConfig(configPathGlobal);
   config.readConfig(configPathUser);
   config.readConfig(configPathLocal, m_requireConfigPathLocal);

   //Parse CL just to get the CL config.
   config.parseCommandLine(argc, argv, "config");

   //And now get the value of it and parse it.
   loadStandardConfig();
   config.readConfig(m_configPathCL);

   //Now parse the command line for real.
   config.parseCommandLine(argc, argv);

   loadStandardHelp();

   loadBasicConfig();
   loadConfig();
}

inline
int application::reReadConfig()
{
   config.readConfig(configPathGlobal);
   config.readConfig(configPathUser);
   config.readConfig(configPathLocal);

   config.readConfig(m_configPathCL);

   //Now parse the command line for real.
   config.parseCommandLine(m_argc, m_argv);
   
   return 0;
}

inline
void application::setDefaults( int UNUSED(argc),
                               char ** UNUSED(argv)
                             ) //virtual
{
   std::string tmp;

   char * tmpstr;
   
   #ifdef MX_APP_DEFAULT_configPathGlobal
      configPathGlobal = MX_APP_DEFAULT_configPathGlobal;
   #endif
   #ifdef MX_APP_DEFAULT_configPathGlobal_env
      tmpstr = getenv(MX_APP_DEFAULT_configPathGlobal_env);
      if(tmpstr != 0) configPathGlobal = tmpstr;
   #endif


   #ifdef MX_APP_DEFAULT_configPathUser
      configPathUser = MX_APP_DEFAULT_configPathUser;
   #endif
   #ifdef MX_APP_DEFAULT_configPathUser_env
      tmpstr = getenv(MX_APP_DEFAULT_configPathUser_env);
      if(tmpstr != 0) configPathUser = tmpstr;
   #endif

   if(configPathUser != "")
   {
      tmp = getenv("HOME");
      tmp += "/" + configPathUser;
      configPathUser = tmp;
   }

   #ifdef MX_APP_DEFAULT_configPathLocal
      configPathLocal = MX_APP_DEFAULT_configPathLocal;
   #endif
   #ifdef MX_APP_DEFAULT_configPathLocal_env
      tmpstr = getenv(MX_APP_DEFAULT_configPathLocal_env);
      if(tmpstr != 0) configPathLocal = tmpstr;
   #endif

   #ifdef MX_APP_DEFAULT_configPathCLBase_env 
      tmpstr = getenv(MX_APP_DEFAULT_configPathCLBase_env);
      if(tmpstr != 0) m_configPathCLBase = tmpstr;
      if(m_configPathCLBase.size()>0)
         if(m_configPathCLBase[m_configPathCLBase.size()-1] != '/')
            m_configPathCLBase += '/';
   #endif
      
   return;

}

inline
void application::setupStandardConfig() //virtual
{
   config.add("config","c", "config",argType::Required, "", "config", false, "string", "A local config file");
}

inline
void application::setupStandardHelp() //virtual
{
   config.add("help", "h", "help", argType::True,  "", "", false, "none", "Print this message and exit");
}

inline
void application::loadStandardConfig() //virtual
{
   config(m_configPathCL, "config");
   m_configPathCL = m_configPathCLBase + m_configPathCL; 
}

inline
void application::loadStandardHelp() //virtual
{
   config(doHelp, "help");
}

inline
void application::setupBasicConfig() //virtual
{
   return;
}

inline
void application::loadBasicConfig() //virtual
{
   return;
}

inline
void application::optionHelp( configTarget & tgt,
                              ioutils::textTable & tt
                            ) //virtual
{
   std::string tmp;
   int row = tgt.orderAdded;

   if(tgt.shortOpt != "")
   {
      tmp = "-" + tgt.shortOpt;
      tt.addCell(row, 0, tmp);
   }

   if(tgt.longOpt != "")
   {
      tmp =  "--" + tgt.longOpt;
      tt.addCell(row, 1, tmp);
   }

   tmp = "";
   if(tgt.section != "")
   {
      tmp = tgt.section + ".";
   }

   if(tgt.keyword !="")
   {
      tmp += tgt.keyword;
      tt.addCell(row, 2, tmp);
   }

   if(tgt.helpType != "")
   {
      tmp = "<" + tgt.helpType + "> ";
      tt.addCell(row, 3, tmp);
   }

   tt.addCell(row, 4, tgt.helpExplanation);


}

inline
void application::help() //virtual
{
   appConfigurator::targetIterator targIt;
   appConfigurator::clOnlyTargetIterator clOnlyTargIt;

   ioutils::textTable tt;

   int otherColWidth = m_helpSOColWidth + m_helpLOColWidth + m_helpCFColWidth + m_helpTypeColWidth;

   tt.m_colWidths = {m_helpSOColWidth,m_helpLOColWidth,m_helpCFColWidth, m_helpTypeColWidth, m_helpWidth-otherColWidth-4-4};
   tt.m_lineStart = "    ";
   tt.m_colSep = " ";
   tt.m_rowSep = "";

   std::cerr << "usage: " << invokedName << " [options] \n";
   std::cerr << "\n";
   std::cerr << "  Required arguments:\n";

   for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt !=  config.clOnlyTargets.end(); ++clOnlyTargIt)
   {
      if( clOnlyTargIt->isRequired == true)
      {
         optionHelp(*clOnlyTargIt, tt);
      }
   }

   for( targIt = config.m_targets.begin(); targIt !=  config.m_targets.end(); ++targIt)
   {
      if( targIt->second.isRequired == true)
      {
         optionHelp(targIt->second, tt);
      }
   }

   tt.outPut(std::cerr);
   tt.m_rows.clear();

   //row = 0;
   std::cerr << "\n  Optional arguments:\n";

   for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt !=  config.clOnlyTargets.end(); ++clOnlyTargIt)
   {
      if( clOnlyTargIt->isRequired == false)
      {
         optionHelp(*clOnlyTargIt, tt);
      }
   }

   for( targIt = config.m_targets.begin(); targIt !=  config.m_targets.end(); ++targIt)
   {
      if( targIt->second.isRequired == false)
      {
         optionHelp(targIt->second, tt);
      }
   }

   tt.outPut(std::cerr);

}

} //namespace app
} //namespace mx

#endif // application_hpp__
