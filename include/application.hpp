/** \file application.hpp
  * \author Jared R. Males
  * \brief A class for managing application configuration and execution
  *
  */

#ifndef __application_hpp__
#define __application_hpp__


#include "appConfigurator.hpp"

namespace mx
{

/// A class for managing application configuration and execution   
class application
{
   
protected:
   
   std::string invokedName;
   
   std::string configPathGlobal;
   std::string configPathUser;
   std::string configPathLocal;
   std::string configPathCL;
   
   appConfigurator config;

public:
   application()
   {
      return;
   }
   
   ~application()
   {
      return;
   }
   
   int main(int argc, char **argv)
   {
      setup(argc, argv);
      
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
   
   void setConfigPathLocal(const std::string & s)
   {
      configPathLocal = s;
   }
   
protected:   
   
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


   virtual void setDefaults()
   {
      std::string tmp;
      
      #ifdef MX_APP_DEFAULT_configPathGlobal
         configPathGlobal = MX_APP_DEFAULT_configPathGlobal;
      #endif
      #ifdef MX_APP_DEFAULT_configPathGlobal_env
         tmp = getenv(MX_APP_DEFAULT_configPathGlobal_env);
         if(tmp != "") configPathGlobal = tmp;
      #endif
         
         
      #ifdef MX_APP_DEFAULT_configPathUser
         configPathUser = MX_APP_DEFAULT_configPathUser;
      #endif
      #ifdef MX_APP_DEFAULT_configPathUser_env
         tmp = getenv(MX_APP_DEFAULT_configPathUser_env);
         if(tmp != "") configPathUser = tmp;
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
         tmp = getenv(MX_APP_DEFAULT_configPathLocal_env);
         if(tmp != "") configPathLocal = tmp;
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
   
   virtual int execute()
   {
      return 0;
   }
};
      
} //namespace mx

#endif // __application_hpp__

