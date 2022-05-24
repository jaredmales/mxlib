/** \file aoSystem_test.cpp
 */
#include "../../../catch2/catch.hpp"


#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/aoSystem.hpp"

typedef double realT;

using namespace mx::app;
using namespace mx::AO::analysis;

/** Scenario: Loading aoSystem config settings
  * 
  * Verify parsing of config
  * \anchor tests_ao_analysis_aoSystem_config
  */
SCENARIO( "Loading aoSystem config settings", "[ao::analysis::aoSystem]" ) 
{
   GIVEN("no config file")
   {
      aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aosys; //This will be cumulative
      WHEN("verifying the defaults")
      {
         
         REQUIRE(aosys.D() == 0.0);
         REQUIRE(aosys.d_min(0) == 0.0);
         REQUIRE(aosys.d_opt() == 1e-50);
         REQUIRE(aosys.optd() == false);
         REQUIRE(aosys.optd_delta() == 1.0);
         REQUIRE(aosys.wfsBeta()->_id == "Ideal WFS");
         REQUIRE(aosys.lam_wfs() == 0);
         REQUIRE(aosys.npix_wfs().size() == 1);
         REQUIRE(aosys.npix_wfs((size_t) 0) == 0);
         REQUIRE(aosys.ron_wfs().size() == 1);
         REQUIRE(aosys.ron_wfs((size_t) 0) == 0);
         REQUIRE(aosys.Fbg().size() == 1);
         REQUIRE(aosys.Fbg((size_t) 0) == 0);
         REQUIRE(aosys.minTauWFS().size() == 1);
         REQUIRE(aosys.minTauWFS((size_t) 0) == 0);
         REQUIRE(aosys.bin_npix() == true);
         REQUIRE(aosys.tauWFS() == 0);
         REQUIRE(aosys.deltaTau() == 0);
         REQUIRE(aosys.optTau() == true);
         REQUIRE(aosys.lam_sci() == 0);
         REQUIRE(aosys.zeta() == 0);
         REQUIRE(aosys.secZeta() == 1.0);
         REQUIRE(aosys.fit_mn_max() == 100);
         REQUIRE(aosys.circularLimit() == false);
         REQUIRE(aosys.spatialFilter_ku() == std::numeric_limits<realT>::max());
         REQUIRE(aosys.spatialFilter_kv() == std::numeric_limits<realT>::max());
         REQUIRE(aosys.ncp_wfe() == 0);
         REQUIRE(aosys.ncp_alpha() == 2.0);
         REQUIRE(aosys.F0() == 0.0);
         REQUIRE(aosys.starMag() == 0);
         REQUIRE(aosys.Fg() == 0);
         REQUIRE(aosys.Fg(2.5) == 0);
      }
   }

   GIVEN("a valid config file")
   {
      aoSystem<realT, mx::AO::analysis::vonKarmanSpectrum<realT>> aosys; //This will be cumulative

      WHEN("setting wfs to ideal")
      {
         appConfigurator config;
         

         writeConfigFile( "aoSystem.conf",
                          {"aosys"},
                          {"wfs"},
                          {"ideal"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.wfsBeta()->_id == "Ideal WFS");
      }
      WHEN("setting wfs to unmodPyWFS")
      {
         appConfigurator config;
         

         writeConfigFile( "aoSystem.conf",
                          {"aosys"},
                          {"wfs"},
                          {"unmodPyWFS"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");

         aosys.loadConfig(config);

         REQUIRE(aosys.wfsBeta()->_id == "Unmodulated Pyramid");
      }
      WHEN("setting wfs to asympModPyWFS")
      {
         appConfigurator config;
         

         writeConfigFile( "aoSystem.conf",
                          {"aosys"},
                          {"wfs"},
                          {"asympModPyWFS"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.wfsBeta()->_id == "Asymptotic Modulated Pyramid");
      }
      WHEN("setting wfs to SHWFS")
      {
         appConfigurator config;
         

         writeConfigFile( "aoSystem.conf",
                          {"aosys"},
                          {"wfs"},
                          {"SHWFS"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.wfsBeta()->_id == "Shack Hartmann");
      }
      WHEN("typical settings part 1") //These are broken into parts just to keep things <  1 line
      {
         appConfigurator config;
         

         writeConfigFile( "aoSystem.conf",
                          {"aosys","aosys","aosys","aosys","aosys"},
                          {"D","d_min","optd","optd_delta","lam_wfs"},
                          {"7.6","0.122","true","0.0037","800e-9"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.D() == Approx(7.6));
         REQUIRE(aosys.d_min(0) == Approx(0.122));
         REQUIRE(aosys.optd() == true);
         REQUIRE(aosys.optd_delta() == Approx(0.0037));
         REQUIRE(aosys.lam_wfs() == Approx(800e-9));
      }
      WHEN("typical settings part 2")
      {
         appConfigurator config;
         
         writeConfigFile( "aoSystem.conf",
                          {"aosys",    "aosys",   "aosys",      "aosys",         "aosys"},
                          {"npix_wfs", "ron_wfs",  "Fbg",       "minTauWFS",     "bin_npix"},
                          {"100,1001", "0.5,0.25", "0.01,0.03", "0.0001,0.0005", "false"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.npix_wfs(0) == 100);
         REQUIRE(aosys.npix_wfs(1) == 1001);
         REQUIRE(aosys.ron_wfs(0) == 0.5);
         REQUIRE(aosys.ron_wfs(1) == 0.25);
         REQUIRE(aosys.Fbg(0) == 0.01);
         REQUIRE(aosys.Fbg(1) == 0.03);
         REQUIRE(aosys.minTauWFS(0) == 0.0001);
         REQUIRE(aosys.minTauWFS(1) == 0.0005);
         REQUIRE(aosys.bin_npix() == false);
      }
      WHEN("typical settings part 3")
      {
         appConfigurator config;
         
         writeConfigFile( "aoSystem.conf",
                          {"aosys",  "aosys",    "aosys",  "aosys",   "aosys", "aosys",     "aosys",         "aosys",            "aosys"},
                          {"tauWFS", "deltaTau", "optTau", "lam_sci", "zeta", "fit_mn_max", "circularLimit", "spatialFilter_ku", "spatialFilter_kv"},
                          {"0.0007", "0.0001",   "false",  "802e-9",  "0.15", "12",         "true",          "50",               "53"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.tauWFS() == Approx(0.0007));
         REQUIRE(aosys.deltaTau() == Approx(0.0001));
         REQUIRE(aosys.optTau() == false);
         REQUIRE(aosys.lam_sci() == Approx(802e-9));
         REQUIRE(aosys.zeta() == Approx(0.15));
         REQUIRE(aosys.secZeta() == Approx(1/cos(0.15)));
         REQUIRE(aosys.fit_mn_max() == 12);
         REQUIRE(aosys.circularLimit() == true);
         REQUIRE(aosys.spatialFilter_ku() == 50);
         REQUIRE(aosys.spatialFilter_kv() == 53);
      }
      WHEN("typical settings part 4")
      {
         appConfigurator config;
         
         writeConfigFile( "aoSystem.conf",
                          {"aosys",   "aosys",     "aosys",   "aosys", "atm"}, //last one verifies that atm config gets processed
                          {"ncp_wfe", "ncp_alpha", "starMag", "F0",    "r_0"},
                          {"30e-9",   "2.34",      "8.2",     "5e9",   "0.25"}
                        );

         aosys.setupConfig(config);
         config.readConfig("aoSystem.conf");
         aosys.loadConfig(config);

         REQUIRE(aosys.ncp_wfe() == Approx(30e-9));
         REQUIRE(aosys.ncp_alpha() == Approx(2.34));
         REQUIRE(aosys.F0() == Approx(5e9));
         REQUIRE(aosys.starMag() == Approx(8.2));
         REQUIRE(aosys.Fg() == Approx(5e9*pow(10, -0.4*8.2)));
         REQUIRE(aosys.Fg(9.0) == Approx(5e9*pow(10, -0.4*9.0)));
         REQUIRE(aosys.atm.r_0() == Approx(0.25));
      }
   }
}

         
