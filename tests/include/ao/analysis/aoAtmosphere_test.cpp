/** \file fileAtmosphere_test.cpp
 */
#include "../../../catch2/catch.hpp"


#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/aoAtmosphere.hpp"

typedef double realT;

using namespace mx::app;
using namespace mx::AO::analysis;

/** Scenario: Loading aoAtmosphere config settings
  *   
  * Verify parsing of config
  *
  * \anchor tests_ao_analysis_aoAtmosphere_config
  */
SCENARIO( "Loading aoAtmosphere config settings", "[ao::analysis::aoAtmosphere]" ) 
{
   GIVEN("a valid config vile")
   {
      aoAtmosphere<realT> atm; //This will be cumulative

      WHEN("all normal settings")
      {
         appConfigurator config;
         

         writeConfigFile( "aoAtmosphere.conf",
                          {"atm",  "atm",    "atm",        "atm",            "atm",            "atm",   "atm", "atm",       "atm",            "atm"},
                          {"r_0",  "lam_0",  "L_0",        "l_0",            "layer_z",        "h_obs", "H",   "layer_Cn2", "layer_v_wind",   "layer_dir"},
                          {"0.25", "0.4e-9", "10,15,22.5", "0.1,0.01,0.001", "1000,2000,5000", "3002",  "2.3", "1,1,3.0",   "10.5,11.2,23.7", "0.1,0.6,1.6" }
                        );

         atm.setupConfig(config);
         config.readConfig("aoAtmosphere.conf");
         atm.loadConfig(config);

         REQUIRE(atm.r_0() == 0.25);
         REQUIRE(atm.lam_0() == 0.4e-9);

         REQUIRE(atm.L_0(0) == 10);
         REQUIRE(atm.L_0(1) == 15);
         REQUIRE(atm.L_0(2) == 22.5);

         REQUIRE(atm.l_0(0) == 0.1);
         REQUIRE(atm.l_0(1) == 0.01);
         REQUIRE(atm.l_0(2) == 0.001);

         REQUIRE(atm.layer_z(0) == 1000);
         REQUIRE(atm.layer_z(1) == 2000);
         REQUIRE(atm.layer_z(2) == 5000);

         REQUIRE(atm.h_obs() == 3002);

         REQUIRE(atm.H() == 2.3);

         REQUIRE(atm.layer_Cn2(0) == Approx(0.2));
         REQUIRE(atm.layer_Cn2(1) == Approx(0.2));
         REQUIRE(atm.layer_Cn2(2) == Approx(0.6));

         REQUIRE(atm.layer_v_wind(0) == Approx(10.5));
         REQUIRE(atm.layer_v_wind(1) == Approx(11.2));
         REQUIRE(atm.layer_v_wind(2) == Approx(23.7));

         REQUIRE(atm.layer_dir(0) == Approx(0.1));
         REQUIRE(atm.layer_dir(1) == Approx(0.6));
         REQUIRE(atm.layer_dir(2) == Approx(1.6));

         REQUIRE(atm.v_wind() == Approx(19.278644205193533));
         REQUIRE(atm.z_mean() == Approx(3886.4496373530733));

         REQUIRE(atm.nonKolmogorov() == false);
      }
      WHEN("setting v_wind and z_mean")
      {
         appConfigurator config;
         
         writeConfigFile( "aoAtmosphere.conf",
                          {"atm",            "atm",       "atm",            "atm",         "atm",    "atm" },
                          {"layer_z",        "layer_Cn2", "layer_v_wind",   "layer_dir",   "v_wind", "z_mean"  },
                          {"1000,2000,5000", "1,1,3.0",   "10.5,11.2,23.7", "0.1,0.6,1.6", "15.0",    "5001.0" }
                        );
         
         atm.setupConfig(config);
         config.readConfig("aoAtmosphere.conf");
         atm.loadConfig(config);

         REQUIRE(atm.v_wind() == Approx(15.0));
         REQUIRE(atm.z_mean() == Approx(5001.0));

         REQUIRE(atm.nonKolmogorov() == false);
      }
      WHEN("setting nonKolmogorov to false")
      {
         appConfigurator config;
    
         writeConfigFile( "aoAtmosphere.conf",
                          {"atm"},
                          {"nonKolmogorov"},
                          {"false"}
                         );   

         atm.setupConfig(config);
         config.readConfig("aoAtmosphere.conf");
         atm.loadConfig(config);

         REQUIRE(atm.nonKolmogorov() == false);
      }
      WHEN("setting nonKolmogorov to true")
      {
         appConfigurator config;
    
         writeConfigFile( "aoAtmosphere.conf",
                          {"atm"},
                          {"nonKolmogorov"},
                          {"true"}
                         );   

         atm.setupConfig(config);
         config.readConfig("aoAtmosphere.conf");
         atm.loadConfig(config);

         REQUIRE(atm.nonKolmogorov() == true);
      }

      WHEN("setting nonKolmogorov alpha")
      {
         atm.nonKolmogorov(false);
         REQUIRE(atm.nonKolmogorov() == false);

         appConfigurator config;
    
         writeConfigFile( "aoAtmosphere.conf",
                          {"atm"},
                          {"alpha"},
                          {"4.7"}
                         );   

         atm.setupConfig(config);
         config.readConfig("aoAtmosphere.conf");
         atm.loadConfig(config);

         REQUIRE(atm.alpha() == Approx(4.7));
         REQUIRE(atm.nonKolmogorov() == true);
      }

      WHEN("setting nonKolmogorov beta")
      {
         atm.nonKolmogorov(false);
         REQUIRE(atm.nonKolmogorov() == false);

         appConfigurator config;
    
         writeConfigFile( "aoAtmosphere.conf",
                          {"atm"},
                          {"beta"},
                          {"0.026"}
                         );   

         atm.setupConfig(config);
         config.readConfig("aoAtmosphere.conf");
         atm.loadConfig(config);

         REQUIRE(atm.beta() == Approx(0.026));
         REQUIRE(atm.nonKolmogorov() == true);
      }

   }
   GIVEN("command line options")
   {
      aoAtmosphere<realT> atm; //This will be cumulative

      WHEN("setting nonKolmogorov to true")
      {
         appConfigurator config;
         std::vector<std::string> argvs(2);

         argvs[0] = "test";
         argvs[1] = "--atm.nonKolmogorov=true";

         char *argv[2];
         argv[0] = (char *) argvs[0].data();
         argv[1] = (char *) argvs[1].data();

         REQUIRE(atm.nonKolmogorov() == false);
         atm.setupConfig(config);
         config.parseCommandLine(argvs.size(), argv);
         atm.loadConfig(config);
         REQUIRE(atm.nonKolmogorov() == true);
      }

      WHEN("setting nonKolmogorov to true")
      {
         appConfigurator config;
         std::vector<std::string> argvs(2);

         argvs[0] = "test";
         argvs[1] = "--atm.nonKolmogorov=false";

         char *argv[2];
         argv[0] = (char *) argvs[0].data();
         argv[1] = (char *) argvs[1].data();

         atm.nonKolmogorov(true);
         REQUIRE(atm.nonKolmogorov() == true);
         atm.setupConfig(config);
         config.parseCommandLine(argvs.size(), argv);
         atm.loadConfig(config);
         REQUIRE(atm.nonKolmogorov() == false);
      }
   }
}

         
