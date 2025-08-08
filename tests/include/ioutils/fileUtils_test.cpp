/** \file fileUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../include/mxException.hpp"
#include "../../../include/ioutils/fileUtils.hpp"

/** Verify creation of sequential file names
 *
 * \anchor tests_ioutils_fileUtils_getSequentialFilename
 */
SCENARIO( "creating sequential filenames", "[ioutils::fileUtils]" )
{
    GIVEN( "a varying numbers of digits desired" )
    {
        WHEN( "default 4 digits, starting at 0" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test" );
            REQUIRE( fname == "base0000.test" );
        }

        WHEN( "default 4 digits, starting at 1" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 1 );
            REQUIRE( fname == "base0001.test" );
        }

        WHEN( "default 7 digits, starting at 0" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 0, 7 );
            REQUIRE( fname == "base0000000.test" );
        }

        WHEN( "default 7 digits, starting at 1" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 1, 7 );
            REQUIRE( fname == "base0000001.test" );
        }

        WHEN( "default 12 digits, starting at 0" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 0, 12 );
            REQUIRE( fname == "base000000000000.test" );
        }

        WHEN( "default 12 digits, starting at 1" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 1, 12 );
            REQUIRE( fname == "base000000000001.test" );
        }
    }
}

#include <fstream>
#include <filesystem>

void createfiles( const std::string &basedir )
{
    std::filesystem::create_directories( basedir );

    std::ofstream fout;

    fout.open( basedir + "/file_txt_1.txt" );
    fout.close();
    fout.open( basedir + "/file_xtx_2.txt" );
    fout.close();
    fout.open( basedir + "/elif_txt_3.txt" );
    fout.close();
    fout.open( basedir + "/file_xtx_4.txt" );
    fout.close();
    fout.open( basedir + "/elif_txt_5.txt" );
    fout.close();
    fout.open( basedir + "/elif_xtx_1.xxx" );
    fout.close();
    fout.open( basedir + "/elif_txt_2.xxx" );
    fout.close();
    fout.open( basedir + "/file_xtx_3.xxx" );
    fout.close();
    fout.open( basedir + "/elif_txt_4.xxx" );
    fout.close();
    fout.open( basedir + "/file_xtx_5.xxx" );
    fout.close();
}

/** \test Scenario: Getting a list of files
 *
 * Tests that files are read according to the specification and sorted.  Also tests basic errors.
 *
 * \anchor tests_ioutils_fileUtils_getFileNames
 */
SCENARIO( "Getting a list of files", "[ioutils::fileUtils]" )
{
    std::string basedir = "/tmp/fileUtils_test/dir";
    createfiles( basedir );

    GIVEN( "a directory with files of various type and names" )
    {
        WHEN( "directory only" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "", "", "" );
            REQUIRE( fnames.size() == 10 );
        }

        WHEN( "directory only, overload 1" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "" );
            REQUIRE( fnames.size() == 10 );
        }

        WHEN( "directory only, overload 2" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir );
            REQUIRE( fnames.size() == 10 );
        }

        WHEN( "single extension with ." )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "", "", ".txt" );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_3.txt" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_5.txt" );
            REQUIRE( fnames[2] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_4.txt" );
        }

        WHEN( "single extension with ., overload" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, ".txt" );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_3.txt" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_5.txt" );
            REQUIRE( fnames[2] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_4.txt" );
        }

        WHEN( "single extension without ." )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "", "", "txt" );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_3.txt" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_5.txt" );
            REQUIRE( fnames[2] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_4.txt" );
        }

        WHEN( "different extension with ." )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "", "", ".xxx" );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_2.xxx" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_4.xxx" );
            REQUIRE( fnames[2] == basedir + "/elif_xtx_1.xxx" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_5.xxx" );
        }

        WHEN( "a prefix, no extension" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "file", "", "" );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_4.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_5.xxx" );
        }

        WHEN( "a prefix, with extension with ." )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "file", "", ".txt" );
            REQUIRE( fnames.size() == 3 );
            REQUIRE( fnames[0] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_4.txt" );
        }

        WHEN( "a substr alone" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "", "xtx", "" );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_xtx_1.xxx" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_4.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_5.xxx" );
        }

        WHEN( "a substr which is actually the prefix" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "", "file", "" );
            REQUIRE( fnames.size() == 0 );
        }

        WHEN( "a prefix and a substr which is actually the prefix" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "file", "file", "" );
            REQUIRE( fnames.size() == 0 );
        }

        WHEN( "a prefix, and a substr, no extension" )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "file", "xtx", "" );
            REQUIRE( fnames.size() == 4 );
            REQUIRE( fnames[0] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_4.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_5.xxx" );
        }

        WHEN( "a prefix, a substr, and an extension with ." )
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir, "file", "xtx", ".txt" );
            REQUIRE( fnames.size() == 2 );
            REQUIRE( fnames[0] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_4.txt" );
        }
    }

    GIVEN( "a directory which does not exist" )
    {
        bool caught = false;
        try
        {
            std::vector<std::string> fnames = mx::ioutils::getFileNames( basedir + "nf", "file", "xtx", ".txt" );
        }
        catch( const mx::err::notfound &e )
        {
            caught = true;
        }
        catch( ... )
        {
        }

        REQUIRE( caught == true );
    }

    GIVEN( "a directory which is a file" )
    {
        bool caught = false;
        try
        {
            std::vector<std::string> fnames =
                mx::ioutils::getFileNames( basedir + "/file_xtx_2.txt", "file", "xtx", ".txt" );
        }
        catch( const mx::err::invalidarg &e )
        {
            caught = true;
        }
        catch( ... )
        {
        }

        REQUIRE( caught == true );
    }
}
