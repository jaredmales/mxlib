/** \file readColumns_test.cpp
 */
#include "../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../include/ioutils/readColumns.hpp"

using namespace Catch::Matchers;

/** \test Scenario: Reading space delimited numeric data
 *
 * \anchor tests_ioutils_readColumns_space_numeric
 */
SCENARIO( "Reading space delimited numeric data", "[ioutils::readColumns]" )
{
    GIVEN( "a single column of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_single_col_float.dat";
        std::ofstream fout;
        fout.open(fname);
        fout << "#a commment\n";
        fout << "1.23\n";
        fout << "2.15\n";
        fout << "4.96 #with comment\n";
        fout << "5.23\n";
        fout << "1e-9"; //no newline on last line
        fout.close();

        WHEN( "reading as int" )
        {
            std::vector<int> data0;
            mx::error_t errc = mx::ioutils::readColumns(fname, data0);
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE(data0.size() == 5);
            REQUIRE(data0[0] == 1);
            REQUIRE(data0[1] == 2);
            REQUIRE(data0[2] == 4);
            REQUIRE(data0[3] == 5);
            REQUIRE(data0[4] == 1); //this is a consequence of stoi
        }

        WHEN( "reading as float" )
        {
            std::vector<float> data0;
            mx::error_t errc = mx::ioutils::readColumns(fname, data0);
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE(data0.size() == 5);
            REQUIRE_THAT(data0[0], WithinRel(1.23, 1e-5));
            REQUIRE_THAT(data0[1], WithinRel(2.15, 1e-5));
            REQUIRE_THAT(data0[2], WithinRel(4.96, 1e-5));
            REQUIRE_THAT(data0[3], WithinRel(5.23, 1e-5));
            REQUIRE_THAT(data0[4], WithinRel(1e-9, 1e-5));

        }
    }

    GIVEN( "a two columns of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_two_cols_float.dat";
        std::ofstream fout;
        fout.open(fname);
        fout << "#a commment\n";
        fout << "1.23 6.78\n";
        fout << "2.15 8.88\n";
        fout << "4.96 -2.33#with comment\n";
        fout << "\n"; //empty line
        fout << "     \n"; //blank line
        fout << "5.23 9.9e5\n";
        fout << "1e-9 -5.6e2\n";
        fout << "   #comment at end not on first char and no newline";
        fout.close();

        WHEN( "reading both as int" )
        {
            std::vector<int> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns(fname, data0, data1);
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE(data0.size() == 5);
            REQUIRE(data0[0] == 1);
            REQUIRE(data0[1] == 2);
            REQUIRE(data0[2] == 4);
            REQUIRE(data0[3] == 5);
            REQUIRE(data0[4] == 1); //this is a consequence of stoi

            REQUIRE(data1.size() == 5);
            REQUIRE(data1[0] == 6);
            REQUIRE(data1[1] == 8);
            REQUIRE(data1[2] == -2);
            REQUIRE(data1[3] == 9);//this is a consequence of stoi
            REQUIRE(data1[4] == -5); //this is a consequence of stoi
        }

        WHEN( "reading both as float" )
        {
            std::vector<float> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns(fname, data0, data1);
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE(data0.size() == 5);
            REQUIRE_THAT(data0[0], WithinRel(1.23, 1e-5));
            REQUIRE_THAT(data0[1], WithinRel(2.15, 1e-5));
            REQUIRE_THAT(data0[2], WithinRel(4.96, 1e-5));
            REQUIRE_THAT(data0[3], WithinRel(5.23, 1e-5));
            REQUIRE_THAT(data0[4], WithinRel(1e-9, 1e-5));

            REQUIRE(data1.size() == 5);
            REQUIRE_THAT(data1[0], WithinRel(6.78, 1e-5));
            REQUIRE_THAT(data1[1], WithinRel(8.88, 1e-5));
            REQUIRE_THAT(data1[2], WithinRel(-2.33, 1e-5));
            REQUIRE_THAT(data1[3], WithinRel(9.9e5, 1e-5));
            REQUIRE_THAT(data1[4], WithinRel(-5.6e2, 1e-5));
        }

    }
}

/** \test Scenario: Reading space delimited numeric data with errors
 *
 * \anchor tests_ioutils_readColumns_twocol_space_numeric_errors
 */
SCENARIO( "Reading space delimited numeric data with errors", "[ioutils::readColumns]" )
{
    GIVEN( "two columns of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_two_col_float_intoverflow.dat";
        std::ofstream fout;
        fout.open(fname);
        fout << "1.23 6.7\n";
        fout << "2.15 8.8\n";
        fout << "4.96 9.9\n";
        fout << "5.23 10.11\n";
        fout << "1e-9 55555555555555555555555555555555555\n";
        fout.close();

        WHEN( "reading as int" )
        {
            std::vector<int> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::v>(fname, data0, data1);
            REQUIRE( errc == mx::error_t::std_out_of_range );
        }

        /*WHEN( "reading as int64_t" )
        {
            std::vector<int64_t> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns(fname, data0, data1);
            REQUIRE( errc == mx::error_t::std_out_of_range );
        }*/
    }
}
