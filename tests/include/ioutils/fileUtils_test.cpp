/** \file fileUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include <fstream>
#include "../../../include/ioutils/fileUtils.hpp"

#undef ioutils_fileUtils_hpp
#define MXLIBTEST_NAMESPACE MXLIBTEST_DIREXISTSIS_ISEXISTSERR_ns
#define MXLIBTEST_DIREXISTSIS_ISEXISTSERR
#include "../../../include/ioutils/fileUtils.hpp"
#undef MXLIBTEST_NAMESPACE
#undef MXLIBTEST_DIREXISTSIS_ISEXISTSERR

#undef ioutils_fileUtils_hpp
#define MXLIBTEST_NAMESPACE MXLIBTEST_DIREXISTSIS_ISDIRERR_ns
#define MXLIBTEST_DIREXISTSIS_ISDIRERR
#include "../../../include/ioutils/fileUtils.hpp"
#undef MXLIBTEST_NAMESPACE
#undef MXLIBTEST_DIREXISTSIS_ISDIRERR

namespace unitTest
{
namespace ioutilsTest
{
namespace fileUtilsTest
{
/// Checking existence and whether or not a path is a directory
/**
 * \ingroup fileUtils_unit_tests
 */
TEST_CASE( "Checking existence and whether or not a path is a directory", "[ioutils::fileUtils]" )
{
    SECTION( "directory exists and is a directory" )
    {
        std::filesystem::create_directories( "/tmp/fileUtils_test/dir_is_exists/" );

        mx::error_t errc;
        bool isdir = mx::ioutils::dir_exists_is( "/tmp/fileUtils_test/dir_is_exists/", errc );

        REQUIRE( errc == mx::error_t::noerror );
        REQUIRE( isdir == true );
    }

    SECTION( "directory does not exist" )
    {
        mx::error_t errc;
        bool isdir = mx::ioutils::dir_exists_is( "/tmp/fileUtils_test/dir_is_exists0/", errc );

        REQUIRE( errc == mx::error_t::noerror );
        REQUIRE( isdir == false );
    }

    SECTION( "directory is a file does not exist" )
    {
        std::filesystem::create_directories( "/tmp/fileUtils_test/dir_is_exists/" );
        std::ofstream fout;
        fout.open( "/tmp/fileUtils_test/dir_is_exists/test.txt" );
        fout << "test";
        fout.close();

        mx::error_t errc;
        bool isdir = mx::ioutils::dir_exists_is( "/tmp/fileUtils_test/dir_is_exists/test.txt", errc );

        REQUIRE( errc == mx::error_t::noerror );
        REQUIRE( isdir == false );
    }

    SECTION( "error from exists" )
    {
        mx::error_t errc;
        bool isdir =
            mx::ioutils::MXLIBTEST_DIREXISTSIS_ISEXISTSERR_ns::dir_exists_is( "/tmp/fileUtils_test/dir_is_exists/",
                                                                              errc );

        REQUIRE( errc == mx::error_t::eexist );
        REQUIRE( isdir == false );
    }

    SECTION( "error from is_directory" )
    {
        mx::error_t errc;
        bool isdir =
            mx::ioutils::MXLIBTEST_DIREXISTSIS_ISDIRERR_ns::dir_exists_is( "/tmp/fileUtils_test/dir_is_exists/", errc );

        REQUIRE( errc == mx::error_t::eacces );
        REQUIRE( isdir == false );
    }
}
/// creating sequential filenames
/**
 * \ingroup fileUtils_unit_tests
 */
TEST_CASE( "creating sequential filenames", "[ioutils::fileUtils]" )
{
    SECTION( "a varying numbers of digits desired" )
    {
        SECTION( "default 4 digits, starting at 0" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test" );
            REQUIRE( fname == "base0000.test" );
        }

        SECTION( "default 4 digits, starting at 1" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 1 );
            REQUIRE( fname == "base0001.test" );
        }

        SECTION( "default 7 digits, starting at 0" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 0, 7 );
            REQUIRE( fname == "base0000000.test" );
        }

        SECTION( "default 7 digits, starting at 1" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 1, 7 );
            REQUIRE( fname == "base0000001.test" );
        }

        SECTION( "default 12 digits, starting at 0" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 0, 12 );
            REQUIRE( fname == "base000000000000.test" );
        }

        SECTION( "default 12 digits, starting at 1" )
        {
            std::string fname = mx::ioutils::getSequentialFilename( "base", ".test", 1, 12 );
            REQUIRE( fname == "base000000000001.test" );
        }
    }
}

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

/// Getting a list of files
/**
 * Tests that files are read according to the specification and sorted.  Also tests basic errors.
 *
 *//// Getting a list of files

/// Getting a list of files
/**
 * \ingroup fileUtils_unit_tests
 */
TEST_CASE( "Getting a list of files", "[ioutils::fileUtils]" )
{
    std::string basedir = "/tmp/fileUtils_test/dir";
    createfiles( basedir );

    SECTION( "a directory with files of various type and names" )
    {
        SECTION( "directory only" )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "", "", "" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 10 );
        }

        SECTION( "single extension with ." )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "", "", ".txt" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_3.txt" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_5.txt" );
            REQUIRE( fnames[2] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_4.txt" );
        }

        SECTION( "single extension without ." )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "", "", "txt" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_3.txt" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_5.txt" );
            REQUIRE( fnames[2] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_4.txt" );
        }

        SECTION( "different extension with ." )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "", "", ".xxx" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_txt_2.xxx" );
            REQUIRE( fnames[1] == basedir + "/elif_txt_4.xxx" );
            REQUIRE( fnames[2] == basedir + "/elif_xtx_1.xxx" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_5.xxx" );
        }

        SECTION( "a prefix, no extension" )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "file", "", "" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_4.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_5.xxx" );
        }

        SECTION( "a prefix, with extension with ." )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "file", "", ".txt" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 3 );
            REQUIRE( fnames[0] == basedir + "/file_txt_1.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_4.txt" );
        }

        SECTION( "a substr alone" )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "", "xtx", "" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 5 );
            REQUIRE( fnames[0] == basedir + "/elif_xtx_1.xxx" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_4.txt" );
            REQUIRE( fnames[4] == basedir + "/file_xtx_5.xxx" );
        }

        SECTION( "a substr which is actually the prefix" )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "", "file", "" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 0 );
        }

        SECTION( "a prefix and a substr which is actually the prefix" )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "file", "file", "" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 0 );
        }

        SECTION( "a prefix, and a substr, no extension" )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "file", "xtx", "" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 4 );
            REQUIRE( fnames[0] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_3.xxx" );
            REQUIRE( fnames[2] == basedir + "/file_xtx_4.txt" );
            REQUIRE( fnames[3] == basedir + "/file_xtx_5.xxx" );
        }

        SECTION( "a prefix, a substr, and an extension with ." )
        {
            std::vector<std::string> fnames;
            mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir, "file", "xtx", ".txt" );

            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( fnames.size() == 2 );
            REQUIRE( fnames[0] == basedir + "/file_xtx_2.txt" );
            REQUIRE( fnames[1] == basedir + "/file_xtx_4.txt" );
        }
    }

    SECTION( "a directory which does not exist" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir + "nf", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::dirnotfound );
    }

    SECTION( "a directory which is a file" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc = mx::ioutils::getFileNames( fnames, basedir + "/file_xtx_2.txt", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::invalidarg );
    }

    SECTION( "a directory which does not exist, verbose = vv" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc = mx::ioutils::getFileNames<mx::verbose::vv>( fnames, basedir + "nf", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::dirnotfound );
    }

    SECTION( "a directory which is a file, verbose = vv" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc =
            mx::ioutils::getFileNames<mx::verbose::vv>( fnames, basedir + "/file_xtx_2.txt", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::invalidarg );
    }

    SECTION( "a directory which does not exist, verbose = v" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc = mx::ioutils::getFileNames<mx::verbose::v>( fnames, basedir + "nf", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::dirnotfound );
    }

    SECTION( "a directory which is a file, verbose = v" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc =
            mx::ioutils::getFileNames<mx::verbose::v>( fnames, basedir + "/file_xtx_2.txt", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::invalidarg );
    }

    SECTION( "a directory which does not exist, verbose = o" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc = mx::ioutils::getFileNames<mx::verbose::o>( fnames, basedir + "nf", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::dirnotfound );
    }

    SECTION( "a directory which is a file, verbose = o" )
    {
        std::vector<std::string> fnames;
        mx::error_t errc =
            mx::ioutils::getFileNames<mx::verbose::o>( fnames, basedir + "/file_xtx_2.txt", "file", "xtx", ".txt" );

        REQUIRE( errc == mx::error_t::invalidarg );
    }
}

} // namespace fileUtilsTest
} // namespace ioutilsTest
} // namespace unitTest
