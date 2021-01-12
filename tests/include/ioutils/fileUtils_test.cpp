/** \file fileUtils_test.cpp
 */
#include "../../catch2/catch.hpp"


#define MX_NO_ERROR_REPORTS

#include "../../../include/ioutils/fileUtils.hpp"

/** Verify creation of sequential file names
  * 
  * \anchor tests_ioutils_fileUtils_getSequentialFilename
  */
SCENARIO( "creating sequential filenames", "[ioutils::fileUtils]" ) 
{
   GIVEN("a varying numbers of digits desired")
   {
      WHEN("default 4 digits, starting at 0")
      {
         std::string fname = mx::ioutils::getSequentialFilename("base", ".test");
         REQUIRE(fname == "base0000.test");
      }
      
      WHEN("default 4 digits, starting at 1")
      {
         std::string fname = mx::ioutils::getSequentialFilename("base", ".test", 1);
         REQUIRE(fname == "base0001.test");
      }
      
      WHEN("default 7 digits, starting at 0")
      {
         std::string fname = mx::ioutils::getSequentialFilename("base", ".test",0,7);
         REQUIRE(fname == "base0000000.test");
      }
      
      WHEN("default 7 digits, starting at 1")
      {
         std::string fname = mx::ioutils::getSequentialFilename("base", ".test", 1, 7);
         REQUIRE(fname == "base0000001.test");
      }
      
      WHEN("default 12 digits, starting at 0")
      {
         std::string fname = mx::ioutils::getSequentialFilename("base", ".test",0,12);
         REQUIRE(fname == "base000000000000.test");
      }
      
      WHEN("default 12 digits, starting at 1")
      {
         std::string fname = mx::ioutils::getSequentialFilename("base", ".test", 1, 12);
         REQUIRE(fname == "base000000000001.test");
      }
   }
}

         
