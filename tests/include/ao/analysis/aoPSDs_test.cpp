/** \file fileAtmosphere_test.cpp
 */
#include "../../../catch2/catch.hpp"


#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/aoPSDs.hpp"

typedef double realT;

using namespace mx::app;
using namespace mx::AO::analysis;

/** Scenario: Loading aoAtmosphere config settings
  *   
  * Verify parsing of config
  *
  * \anchor tests_ao_analysis_aoAtmosphere_config
  */
SCENARIO( "Loading aoPSD config settings", "[ao::analysis::aoPSD]" ) 
{
   GIVEN("a valid config file")
   {
      vonKarmanSpectrum<realT> psd; //This will be cumulative

      WHEN("all normal settings: phase")
      {
         appConfigurator config;

         writeConfigFile( "aoPSDs.conf",
                          {"psd",  "psd",       "psd",         "psd",            "psd"            },
                          {"D",    "subPiston", "subTipTilt",  "scintillation",  "component" },
                          {"5.5",  "true",      "true",        "true",           "phase" }
                        );

         psd.setupConfig(config);
         config.readConfig("aoPSDs.conf");
         psd.loadConfig(config);

         REQUIRE(psd.D() == 5.5);
         REQUIRE(psd.subPiston() == true);
         REQUIRE(psd.subTipTilt() == true);
         REQUIRE(psd.scintillation() == true);
         REQUIRE(psd.component() == PSDComponent::phase);

      }
      WHEN("all normal settings: amplitude")
      {
         appConfigurator config;

         writeConfigFile( "aoPSDs.conf",
                          {"psd",  "psd",       "psd",         "psd",            "psd"            },
                          {"D",    "subPiston", "subTipTilt",  "scintillation",  "component" },
                          {"3.5",  "false",      "false",        "false",           "amplitude" }
                        );

         psd.setupConfig(config);
         config.readConfig("aoPSDs.conf");
         psd.loadConfig(config);

         REQUIRE(psd.D() == 3.5);
         REQUIRE(psd.subPiston() == false);
         REQUIRE(psd.subTipTilt() == false);
         REQUIRE(psd.scintillation() == false);
         REQUIRE(psd.component() == PSDComponent::amplitude);

      }
      WHEN("all normal settings: dispPhase")
      {
         appConfigurator config;

         writeConfigFile( "aoPSDs.conf",
                          {"psd",  "psd",       "psd",         "psd",            "psd"            },
                          {"D",    "subPiston", "subTipTilt",  "scintillation",  "component" },
                          {"7.2",  "true",      "false",        "true",           "dispPhase" }
                        );

         psd.setupConfig(config);
         config.readConfig("aoPSDs.conf");
         psd.loadConfig(config);

         REQUIRE(psd.D() == 7.2);
         REQUIRE(psd.subPiston() == true);
         REQUIRE(psd.subTipTilt() == false);
         REQUIRE(psd.scintillation() == true);
         REQUIRE(psd.component() == PSDComponent::dispPhase);

      }
      WHEN("all normal settings: dispAmplitude")
      {
         appConfigurator config;

         writeConfigFile( "aoPSDs.conf",
                          {"psd",  "psd",       "psd",         "psd",            "psd"            },
                          {"D",    "subPiston", "subTipTilt",  "scintillation",  "component" },
                          {"25.4",  "false",      "true",        "false",           "dispAmplitude" }
                        );

         psd.setupConfig(config);
         config.readConfig("aoPSDs.conf");
         psd.loadConfig(config);

         REQUIRE(psd.D() == 25.4);
         REQUIRE(psd.subPiston() == false);
         REQUIRE(psd.subTipTilt() == true);
         REQUIRE(psd.scintillation() == false);
         REQUIRE(psd.component() == PSDComponent::dispAmplitude);

      }
   }
}

         
