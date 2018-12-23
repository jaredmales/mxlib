/** \file application.hpp
  * \author Jared R. Males
  * \brief A class for managing application configuration and execution
  *
  */

#ifndef __application_hpp__
#define __application_hpp__


#include "../mxlib.hpp"

#include "appConfigurator.hpp"
#include "../ioutils/textTable.hpp"

namespace mx
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
  * - A user configuration file
  * - A local configuration file
  * - A configuration file specified on the command line
  * - The command line
  *
  * At each step in the above order, values read from that configuration source override any previous settings. So the command line has the highest precedence.
  *
  * The configuration files are ini-style, with sections.  That is
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
  *
  * The command line parser handles both short-opt ("-h -vArg -n Arg") and long-opt ("--help --value=Arg --number=Arg") styles.
  *
  * \note After loadConfig() but before execute(), the containers in \ref config are de-allocated , so they can not be used inside execute.
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
   std::string configPathCL; ///< The path to a configuration file specified on the command line.

   appConfigurator config; ///< The structure used for parsing and storing the configuration.

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

   ///Loads the values of "help" and "config".
   /** Override this function if you do not want to use these, or have different behavior.
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

   config.clear();

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
   config.readConfig(configPathLocal);

   //Parse CL just to get the CL config.
   config.parseCommandLine(argc, argv, "config");

   //And now get the value of it and parse it.
   loadStandardConfig();
   config.readConfig(configPathCL);

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

   config.readConfig(configPathCL);

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

   #ifdef MX_APP_DEFAULT_configPathGlobal
      configPathGlobal = MX_APP_DEFAULT_configPathGlobal;
   #endif
   #ifdef MX_APP_DEFAULT_configPathGlobal_env
      char * tmpstr;
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

   return;

}

inline
void application::setupStandardConfig() //virtual
{
   config.add("config","c", "config",mx::argType::Required, "", "config", false, "string", "A local config file");
}

inline
void application::setupStandardHelp() //virtual
{
   config.add("help", "h", "help", mx::argType::True,  "", "", false, "none", "Print this message and exit");
}

inline
void application::loadStandardConfig() //virtual
{
   config(configPathCL, "config");
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

  // int row = 0;
   for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt !=  config.clOnlyTargets.end(); ++clOnlyTargIt)
   {
      if( clOnlyTargIt->isRequired == true)
      {
         optionHelp(*clOnlyTargIt, tt);
         //++row;
      }
   }

   for( targIt = config.targets.begin(); targIt !=  config.targets.end(); ++targIt)
   {
      if( targIt->second.isRequired == true)
      {
         optionHelp(targIt->second, tt);
         //++row;
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
         //++row;
      }
   }

   for( targIt = config.targets.begin(); targIt !=  config.targets.end(); ++targIt)
   {
      if( targIt->second.isRequired == false)
      {
         optionHelp(targIt->second, tt);
         //++row;
      }
   }

   tt.outPut(std::cerr);

}

} //namespace mx

#endif // __application_hpp__
