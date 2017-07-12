/** \file application.hpp
  * \author Jared R. Males
  * \brief A class for managing application configuration and execution
  *
  */

#ifndef __application_hpp__
#define __application_hpp__


#include "appConfigurator.hpp"
#include "../textTable.hpp"

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
  * \note After loadConfig() but before execute(), the containers in \ref config are de-allocated , so they can not be used inside execute.
  * 
  * \ingroup mxApp
  */
class application
{
   
protected:
   
   std::string invokedName;
   
   std::string configPathGlobal;
   std::string configPathUser;
   std::string configPathLocal;
   std::string configPathCL;
   
   appConfigurator config;
   
   bool doHelp;

   int helpWidth;
   
public:
   application()
   {
      doHelp = false;
      helpWidth = 120;
      
      return;
   }
   
   ~application()
   {
      return;
   }
   
   ///The application main function.
   /** Call this from the true main function, passing the command line arguments to be processed.
     * This calls \ref setup(), then deletes the config structure, and then calls \ref execute(). 
     */
   int main( int argc, ///< [in] standard command line result specifying number of argumetns in argv 
             char **argv ///< [in] standard command line result containing the arguments.
           )
   {
      setup(argc, argv);
      
      if(doHelp)
      {
         help();
         return 0;
      }
      
      config.clear();
      
      return execute();
   }
   
   void setConfigPathGlobal(const std::string & s)
   {
      configPathGlobal = s;
   }
   
   void setConfigPathUser(const std::string & s)
   {
      configPathUser = s;
   }
   
   void setConfigPathLocal( const std::string & s )
   {
      configPathLocal = s;
   }
   
protected:   
   
   ///Sets up the application by executing the configuration steps
   /** This will not normally need to be implemented by derived clasess.
     * Only do so if you intende to change the configuraiton process. 
     */
   virtual void setup(int argc, char ** argv)
   {
      invokedName = argv[0];
      
      setDefaults();
      
      setupConfig();
      
      config.readConfig(configPathGlobal);
      config.readConfig(configPathUser);
      config.readConfig(configPathLocal);
      
      config.parseCommandLine(argc, argv);

      setConfigPathCL();
      
      config.readConfig(configPathCL);
      
      loadConfig();
   }


   ///Set the default search paths for config files
   /** In general you should not need to redefine this.
     */
   virtual void setDefaults()
   {
      char *tmpstr;
      std::string tmp;
      
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
         
      return;
      
   }
   
   virtual void setupConfig()
   {
      return;
   }
   
   virtual void setConfigPathCL()
   {
      configPathCL = "";
      return;
   }
   
   virtual void loadConfig()
   {
      return;
   }
   
   virtual void optionHelp( configTarget & tgt, textTable & tt, int row)
   {
      std::string tmp;
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
   
   ///Print a formatted help message, based on the config target inputs.
   virtual void help()
   {
      appConfigurator::targetIterator targIt;
      appConfigurator::clOnlyTargetIterator clOnlyTargIt;
      
      textTable tt;
      
      tt.colWidths = {2,15,25, 15, helpWidth-57-4-4};
      tt.lineStart = "    ";
      tt.colSep = " ";
      
      std::cerr << "usage: " << invokedName << " [options] \n";
      std::cerr << "\n";
      std::cerr << "  Required options:\n";
      
      int row = 0;
      for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt !=  config.clOnlyTargets.end(); ++clOnlyTargIt)
      {
         if( clOnlyTargIt->isRequired == true)
         {
            optionHelp(*clOnlyTargIt, tt, row);
            ++row;
         }
      }
            
      for( targIt = config.targets.begin(); targIt !=  config.targets.end(); ++targIt)
      {
         if( targIt->second.isRequired == true)
         {
            optionHelp(targIt->second, tt, row);
            ++row;
         }
      }

      tt.outPut(std::cerr);
      tt.rows.clear();
      
      row = 0;
      std::cerr << "\n  Optional options:\n";
      
      for( clOnlyTargIt = config.clOnlyTargets.begin(); clOnlyTargIt !=  config.clOnlyTargets.end(); ++clOnlyTargIt)
      {
         if( clOnlyTargIt->isRequired == false)
         {
            optionHelp(*clOnlyTargIt, tt, row);
            ++row;
         }
      }
            
      for( targIt = config.targets.begin(); targIt !=  config.targets.end(); ++targIt)
      {
         if( targIt->second.isRequired == false)
         {
            optionHelp(targIt->second, tt, row);
            ++row;
         }
      }
      
      tt.outPut(std::cerr);
      
   }
   
   ///This function is where the derived class should do its work.
   virtual int execute()
   {
      return 0;
   }
};
      
} //namespace mx

#endif // __application_hpp__

