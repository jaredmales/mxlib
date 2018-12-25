#define CATCH_CONFIG_MAIN
#include "../../catch2/catch.hpp"

#include <fstream>

#include "../../../include/app/appConfigurator.hpp"


void writeConfigFile( const std::string & fname,
                      const std::vector<std::string> & sections,
                      const std::vector<std::string> & names,
                      const std::vector<std::string> & values 
                    )
{
   std::ofstream fout;
   
   fout.open(fname);
   
   std::string lastSection;
   for(int i=0; i< sections.size(); ++i)
   {
      std::string currSection = sections[i];
      
      if( currSection != lastSection && currSection != "")
      {
         fout << "\n[" << currSection << "]\n";
      }
         
      fout << names[i] << "=" << values[i] << "\n";
         
      lastSection = currSection;
   }
   
   fout.close();
   
   return;
}

SCENARIO( "config file parsing", "[appConfigurator]" ) 
{
   GIVEN("a basic config file")
   {
      WHEN("no sections")
      {
         writeConfigFile( "/tmp/test.conf", {"",     "",     "",      "",      "",      ""},
                                            {"key0", "key1", "key2",  "key3",  "key4",  "key5"},
                                            {"val0", "val1", "val2",  "val3",  "val4",  "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "", "key1", false, "", "");
         config.add("key2", "", "", 0, "", "key2", false, "", "");
         config.add("key3", "", "", 0, "", "key3", false, "", "");
         config.add("key4", "", "", 0, "", "key4", false, "", "");
         config.add("key5", "", "", 0, "", "key5", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0");
         
         config(val, "key1");
         REQUIRE( val == "val1");
         
         config(val, "key2");
         REQUIRE( val == "val2");
         
         config(val, "key3");
         REQUIRE( val == "val3");
         
         config(val, "key4");
         REQUIRE( val == "val4");
         
         config(val, "key5");
         REQUIRE( val == "val5");
      }
      
      WHEN("sections, unique keys")
      {
         writeConfigFile( "/tmp/test.conf", {"",     "",     "sect1", "sect1", "sect2", "sect2"},
                                            {"key0", "key1", "key2",  "key3",  "key4",  "key5"},
                                            {"val0", "val1", "val2",  "val3",  "val4",  "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "", "key1", false, "", "");
         config.add("sect1.key2", "", "", 0, "sect1", "key2", false, "", "");
         config.add("sect1.key3", "", "", 0, "sect1", "key3", false, "", "");
         config.add("sect2.key4", "", "", 0, "sect2", "key4", false, "", "");
         config.add("sect2.key5", "", "", 0, "sect2", "key5", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0");
         
         config(val, "key1");
         REQUIRE( val == "val1");
         
         config(val, "sect1.key2");
         REQUIRE( val == "val2");
         
         config(val, "sect1.key3");
         REQUIRE( val == "val3");
         
         config(val, "sect2.key4");
         REQUIRE( val == "val4");
         
         config(val, "sect2.key5");
         REQUIRE( val == "val5");
      }
      
      WHEN("sections, repeated keys")
      {
         writeConfigFile( "/tmp/test.conf", {"",     "",     "sect1", "sect1", "sect2", "sect2"},
                                            {"key0", "key1", "key2",  "key3",  "key2",  "key3"},
                                            {"val0", "val1", "val2",  "val3",  "val4",  "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "", "key1", false, "", "");
         config.add("sect1.key2", "", "", 0, "sect1", "key2", false, "", "");
         config.add("sect1.key3", "", "", 0, "sect1", "key3", false, "", "");
         config.add("sect2.key2", "", "", 0, "sect2", "key2", false, "", "");
         config.add("sect2.key3", "", "", 0, "sect2", "key3", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0");
         
         config(val, "key1");
         REQUIRE( val == "val1");
         
         config(val, "sect1.key2");
         REQUIRE( val == "val2");
         
         config(val, "sect1.key3");
         REQUIRE( val == "val3");
         
         config(val, "sect2.key2");
         REQUIRE( val == "val4");
         
         config(val, "sect2.key3");
         REQUIRE( val == "val5");
      }
   }
   
   GIVEN("a config file with unused entries")
   {
      
      WHEN("no sections, unused entry in middle")
      {
         writeConfigFile( "/tmp/test.conf", {"",     "",     "",      "",      "",      ""},
                                            {"key0", "key1", "key2",  "key3",  "key4",  "key5"},
                                            {"val0", "val1", "val2",  "val3",  "val4",  "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "", "key1", false, "", "");
         config.add("key2", "", "", 0, "", "key2", false, "", "");
         config.add("key4", "", "", 0, "", "key4", false, "", "");
         config.add("key5", "", "", 0, "", "key5", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         //Check normal parsing
         config(val, "key0");
         REQUIRE( val == "val0");
         
         config(val, "key1");
         REQUIRE( val == "val1");
         
         config(val, "key2");
         REQUIRE( val == "val2");
         
         config(val, "key4");
         REQUIRE( val == "val4");
         
         config(val, "key5");
         REQUIRE( val == "val5");

         //Check that the unused one is available.
         config.configUnused(val, "", "key3");
         REQUIRE( val == "val3");         
      }
      
      WHEN("sections, repeated keys, unused sections")
      {
         writeConfigFile( "/tmp/test.conf", {"",     "",     "sect1", "sect1", "sect2", "sect2", "sect3"},
                                            {"key0", "key1", "key2",  "key3",  "key2",  "key3", "key4" },
                                            {"val0", "val1", "val2",  "val3",  "val4",  "val5", "val6"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "", "key1", false, "", "");
         config.add("sect2.key2", "", "", 0, "sect2", "key2", false, "", "");
         config.add("sect2.key3", "", "", 0, "sect2", "key3", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0");
         
         config(val, "key1");
         REQUIRE( val == "val1");
         
         config(val, "sect2.key2");
         REQUIRE( val == "val4");
         
         config(val, "sect2.key3");
         REQUIRE( val == "val5");
         
         config.configUnused(val, "sect1", "key2");
         REQUIRE( val == "val2");
         
         config.configUnused(val, "sect1", "key3");
         REQUIRE( val == "val3");
         
         config.configUnused(val, "sect3", "key4");
         REQUIRE( val == "val6");
         
         std::vector<std::string> sections;
         config.unusedSections(sections);
         REQUIRE( sections.size()==2);
         REQUIRE( (sections[0] == "sect1" || sections[1] == "sect1") );
         REQUIRE( (sections[0] == "sect3" || sections[1] == "sect3") );
         
      }
   }
   
#if 1
   GIVEN("a config file with repeated keys")
   {
      WHEN("no sections")
      {
         writeConfigFile( "/tmp/test.conf", {"",     "",     "",      "",      "",      ""},
                                            {"key0", "key0", "key2",  "key2",  "key4",  "key4"},
                                            {"val0", "val1", "val2",  "val3",  "val4",  "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key2", "", "", 0, "", "key2", false, "", "");
         config.add("key4", "", "", 0, "", "key4", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0val1");
         
         config(val, "key2");
         REQUIRE( val == "val2val3");
         
         config(val, "key4");
         REQUIRE( val == "val4val5");
         
      }
      
      WHEN("repeated sections and keys")
      {
         writeConfigFile( "/tmp/test7.conf", {"",     "sect1", "sect2", "sect1", "sect2", "sect2",  "sect3"},
                                            {"key0", "key1",  "key2",  "key1",  "key2",  "key2", "key3"},
                                            {"val0", "val1",  "val2",  "val3",  "val4",  "val4.1", "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "sect1", "key1", false, "", "");
         config.add("key2", "", "", 0, "sect2", "key2", false, "", "");
         config.add("key3", "", "", 0, "sect3", "key3", false, "", "");
         
         config.readConfig("/tmp/test7.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0");
         
         config(val, "key1");
         REQUIRE( val == "val1val3");
         
         config(val, "key2");
         REQUIRE( val == "val2val4val4.1");
         
         config(val, "key3");
         REQUIRE( val == "val5");
         
      }
      
      WHEN("multi-line keys")
      {
         writeConfigFile( "/tmp/test.conf", {"",             "sect1", "sect2", "sect1",         "sect2", "sect3"},
                                            {"key0",         "key1",  "key2",  "key1",          "key2",  "key3"},
                                            {"val0\n    val0.1", "val1",  "val2",  "val3\n val3.1",  "val4",  "val5"} );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "sect1", "key1", false, "", "");
         config.add("key2", "", "", 0, "sect2", "key2", false, "", "");
         config.add("key3", "", "", 0, "sect3", "key3", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::string val;
         
         config(val, "key0");
         REQUIRE( val == "val0val0.1");
         
         config(val, "key1");
         REQUIRE( val == "val1val3val3.1");
         
         config(val, "key2");
         REQUIRE( val == "val2val4");
         
         config(val, "key3");
         REQUIRE( val == "val5");
         
      }
   }
   
   GIVEN("a config file with vectors")
   {
      WHEN("no sections")
      {
         writeConfigFile( "/tmp/test.conf", {"",               "",                 ""  },
                                            {"key0",           "key1",             "key2" },
                                            {"val0,val1,val2", " val3, val4,    ", "val5" } );
      
         mx::app::appConfigurator config;
         config.add("key0", "", "", 0, "", "key0", false, "", "");
         config.add("key1", "", "", 0, "", "key1", false, "", "");
         config.add("key2", "", "", 0, "", "key2", false, "", "");
         
         config.readConfig("/tmp/test.conf");
         
         std::vector<std::string> vals;
         
         config(vals, "key0");
         REQUIRE( vals[0] == "val0");
         REQUIRE( vals[1] == "val1");
         REQUIRE( vals[2] == "val2");
         
         config(vals, "key1");
         REQUIRE( vals[0] == "val3");
         REQUIRE( vals[1] == "val4");
         REQUIRE( vals[2] == "");
         
         config(vals, "key2");
         REQUIRE( vals[0] == "val5");
      }
   }
#endif
}

