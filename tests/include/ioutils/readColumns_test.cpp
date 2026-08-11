/** \file readColumns_test.cpp
 * \brief Tests delimited-column input.
 */
#include "../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../include/ioutils/readColumns.hpp"

#undef __readColumns_hpp__
#define MXLIBTEST_NAMESPACE MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns
#define MXLIB_CATCH_ALL_EXCEPTIONS
#include "../../../include/ioutils/readColumns.hpp"
#undef MXLIBTEST_NAMESPACE
#undef MXLIB_CATCH_ALL_EXCEPTIONS

using namespace Catch::Matchers;

/** \cond */
namespace
{

enum class injectedReadColumnsFailure
{
    none,
    invalidArgument,
    outOfRange,
    badAlloc,
    standard,
    unknown,
    streamErrno,
    streamFallback
};

mx::ioutils::readColumnsDetail::operation failingOperation = mx::ioutils::readColumnsDetail::operation::afterReadLine;
injectedReadColumnsFailure injectedFailure = injectedReadColumnsFailure::none;

void injectReadColumnsFailure( mx::ioutils::readColumnsDetail::operation operation, std::ifstream &stream )
{
    if( operation != failingOperation )
    {
        return;
    }

    switch( injectedFailure )
    {
        case injectedReadColumnsFailure::badAlloc:
            throw std::bad_alloc();
        case injectedReadColumnsFailure::standard:
            throw std::runtime_error( "injected readColumns failure" );
        case injectedReadColumnsFailure::unknown:
            throw 1;
        case injectedReadColumnsFailure::streamErrno:
            errno = EIO;
            stream.setstate( operation == mx::ioutils::readColumnsDetail::operation::afterReadLine
                                 ? std::ios::badbit
                                 : std::ios::failbit );
            return;
        case injectedReadColumnsFailure::streamFallback:
            errno = 0;
            stream.setstate( operation == mx::ioutils::readColumnsDetail::operation::afterReadLine
                                 ? std::ios::badbit
                                 : std::ios::failbit );
            return;
        default:
            return;
    }
}

class readColumnsHookGuard
{
  public:
    readColumnsHookGuard() : m_saved( mx::ioutils::readColumnsDetail::operationHook() )
    {
        injectedFailure = injectedReadColumnsFailure::none;
        mx::ioutils::readColumnsDetail::operationHook() = injectReadColumnsFailure;
    }

    ~readColumnsHookGuard()
    {
        mx::ioutils::readColumnsDetail::operationHook() = m_saved;
    }

  private:
    mx::ioutils::readColumnsDetail::operationHookT m_saved;
};

struct throwingColumn
{
    using value_type = std::string;

    void push_back( const std::string & )
    {
        switch( injectedFailure )
        {
            case injectedReadColumnsFailure::invalidArgument:
                throw std::invalid_argument( "injected readcol failure" );
            case injectedReadColumnsFailure::outOfRange:
                throw std::out_of_range( "injected readcol failure" );
            case injectedReadColumnsFailure::badAlloc:
                throw std::bad_alloc();
            case injectedReadColumnsFailure::standard:
                throw std::runtime_error( "injected readcol failure" );
            case injectedReadColumnsFailure::unknown:
                throw 1;
            default:
                return;
        }
    }
};

} // namespace
/** \endcond */

namespace unitTest
{
namespace ioutilsTest
{
namespace readColumnsTest
{

/// Reading space delimited numeric data
/**
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Reading space delimited numeric data", "[ioutils::readColumns]" )
{
    SECTION( "a single column of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_single_col_float.dat";
        std::ofstream fout;
        fout.open( fname );
        fout << "#a commment\n";
        fout << "1.23\n";
        fout << "2.15\n";
        fout << "4.96 #with comment\n";
        fout << "5.23\n";
        fout << "1e-9"; // no newline on last line
        fout.close();

        SECTION( "reading as int" )
        {
            std::vector<int> data0;
            mx::error_t errc = mx::ioutils::readColumns( fname, data0 );
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( data0.size() == 5 );
            REQUIRE( data0[0] == 1 );
            REQUIRE( data0[1] == 2 );
            REQUIRE( data0[2] == 4 );
            REQUIRE( data0[3] == 5 );
            REQUIRE( data0[4] == 1 ); // this is a consequence of stoi
        }

        SECTION( "reading as float" )
        {
            std::vector<float> data0;
            mx::error_t errc = mx::ioutils::readColumns( fname, data0 );
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( data0.size() == 5 );
            REQUIRE_THAT( data0[0], WithinRel( 1.23, 1e-5 ) );
            REQUIRE_THAT( data0[1], WithinRel( 2.15, 1e-5 ) );
            REQUIRE_THAT( data0[2], WithinRel( 4.96, 1e-5 ) );
            REQUIRE_THAT( data0[3], WithinRel( 5.23, 1e-5 ) );
            REQUIRE_THAT( data0[4], WithinRel( 1e-9, 1e-5 ) );
        }
    }

    SECTION( "a two columns of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_two_cols_float.dat";
        std::ofstream fout;
        fout.open( fname );
        fout << "#a commment\n";
        fout << "1.23 6.78\n";
        fout << "2.15 8.88\n";
        fout << "4.96 -2.33#with comment\n";
        fout << "\n";      // empty line
        fout << "     \n"; // blank line
        fout << "5.23 9.9e5\n";
        fout << "1e-9 -5.6e2\n";
        fout << "   #comment at end not on first char and no newline";
        fout.close();

        SECTION( "reading both as int" )
        {
            std::vector<int> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns( fname, data0, data1 );
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( data0.size() == 5 );
            REQUIRE( data0[0] == 1 );
            REQUIRE( data0[1] == 2 );
            REQUIRE( data0[2] == 4 );
            REQUIRE( data0[3] == 5 );
            REQUIRE( data0[4] == 1 ); // this is a consequence of stoi

            REQUIRE( data1.size() == 5 );
            REQUIRE( data1[0] == 6 );
            REQUIRE( data1[1] == 8 );
            REQUIRE( data1[2] == -2 );
            REQUIRE( data1[3] == 9 );  // this is a consequence of stoi
            REQUIRE( data1[4] == -5 ); // this is a consequence of stoi
        }

        SECTION( "reading both as float" )
        {
            std::vector<float> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns( fname, data0, data1 );
            REQUIRE( errc == mx::error_t::noerror );
            REQUIRE( data0.size() == 5 );
            REQUIRE_THAT( data0[0], WithinRel( 1.23, 1e-5 ) );
            REQUIRE_THAT( data0[1], WithinRel( 2.15, 1e-5 ) );
            REQUIRE_THAT( data0[2], WithinRel( 4.96, 1e-5 ) );
            REQUIRE_THAT( data0[3], WithinRel( 5.23, 1e-5 ) );
            REQUIRE_THAT( data0[4], WithinRel( 1e-9, 1e-5 ) );

            REQUIRE( data1.size() == 5 );
            REQUIRE_THAT( data1[0], WithinRel( 6.78, 1e-5 ) );
            REQUIRE_THAT( data1[1], WithinRel( 8.88, 1e-5 ) );
            REQUIRE_THAT( data1[2], WithinRel( -2.33, 1e-5 ) );
            REQUIRE_THAT( data1[3], WithinRel( 9.9e5, 1e-5 ) );
            REQUIRE_THAT( data1[4], WithinRel( -5.6e2, 1e-5 ) );
        }
    }
}

/// Reading space delimited numeric data with errors
/**
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Reading space delimited numeric data with errors", "[ioutils::readColumns]" )
{
    SECTION( "a missing file" )
    {
        std::vector<std::string> filenames;
        std::vector<double> values;
        REQUIRE( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>(
                     "/tmp/readcol_test_missing_file.dat",
                     filenames,
                     values ) == mx::error_t::enoent );
        REQUIRE( filenames.empty() );
        REQUIRE( values.empty() );
    }

    SECTION( "two columns of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_two_col_float_intoverflow.dat";
        std::ofstream fout;
        fout.open( fname );
        fout << "1.23 6.7\n";
        fout << "2.15 8.8\n";
        fout << "4.96 9.9\n";
        fout << "5.23 10.11\n";
        fout << "1e-9 55555555555555555555555555555555555\n";
        fout.close();

        SECTION( "reading as int" )
        {
            std::vector<int> data0, data1;
            mx::error_t errc =
                mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::v>( fname, data0, data1 );
            REQUIRE( errc == mx::error_t::erange );
        }

        /*SECTION( "reading as int64_t" )
        {
            std::vector<int64_t> data0, data1;
            mx::error_t errc = mx::ioutils::readColumns(fname, data0, data1);
            REQUIRE( errc == mx::error_t::std_out_of_range );
        }*/
    }

#if 0 // for doxygen
    std::string fname;
    std::vector<int> data0, data1;
    mx::ioutils::readColumns( fname, data0, data1 );
#endif
}

/// Reading combined type data
/**
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Reading combined type data", "[ioutils::readColumns]" )
{
    SECTION( "a single column of floating point numbers" )
    {
        std::string fname = "/tmp/readcol_test_combos.dat";
        std::ofstream fout;
        fout.open( fname );
        fout << "test1 1.23 2.56\n";
        fout << "test2 2.15 8.93\n";
        fout << "test3 0 0\n";
        fout << "test4 6.7 1\n";
        fout << "test5 0 22.2\n";
        fout << "test6 1e-4 0\n";
        fout.close();

        std::vector<std::string> str0;
        std::vector<float> data1, data2;
        mx::error_t errc =
            mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( fname, str0, data1, data2 );

        REQUIRE( errc == mx::error_t::noerror );

        REQUIRE( str0[0] == "test1" );
        REQUIRE( str0[1] == "test2" );
        REQUIRE( str0[2] == "test3" );
        REQUIRE( str0[3] == "test4" );
        REQUIRE( str0[4] == "test5" );
        REQUIRE( str0[5] == "test6" );

        REQUIRE_THAT( data1[0], WithinRel( 1.23, 1e-7 ) );
        REQUIRE_THAT( data1[1], WithinRel( 2.15, 1e-7 ) );
        REQUIRE( data1[2] == 0 );
        REQUIRE_THAT( data1[3], WithinRel( 6.7, 1e-7 ) );
        REQUIRE( data1[4] == 0 );
        REQUIRE_THAT( data1[5], WithinRel( 1e-4, 1e-7 ) );

        REQUIRE_THAT( data2[0], WithinRel( 2.56, 1e-7 ) );
        REQUIRE_THAT( data2[1], WithinRel( 8.93, 1e-7 ) );
        REQUIRE( data2[2] == 0 );
        REQUIRE_THAT( data2[3], WithinRel( 1, 1e-7 ) );
        REQUIRE_THAT( data2[4], WithinRel( 22.2, 1e-7 ) );
        REQUIRE( data2[5] == 0 );
    }
}

/** \brief Verifies filename and double-valued metadata columns used by HCI quality and weight lists.
 *
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Reading filename and double metadata columns", "[ioutils::readColumns]" )
{
    const std::string filename{ "/tmp/readcol_test_hci_metadata.dat" };
    std::ofstream fout{ filename };
    fout << "# image quality or weight\n";
    fout << "target/frame001.fits 0.25\n";
    fout << "reference/frame002.fits 1.5 # retained inline comment\n";
    fout.close();

    std::vector<std::string> filenames{ "existing.fits" };
    std::vector<double> values{ -1.0 };
    REQUIRE( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename, filenames, values ) ==
             mx::error_t::noerror );

    REQUIRE( filenames ==
             std::vector<std::string>{ "existing.fits", "target/frame001.fits", "reference/frame002.fits" } );
    REQUIRE( values == std::vector<double>{ -1.0, 0.25, 1.5 } );
}

/** \brief Verifies direct column parsing handles empty, missing, terminated, and skipped fields.
 *
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Parsing column edge cases", "[ioutils::readColumns]" )
{
    SECTION( "empty input" )
    {
        std::vector<std::string> values;
        int column = 0;
        REQUIRE( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "", 0, column, values ) ==
                 mx::error_t::noerror );
        REQUIRE( values.empty() );
    }

    SECTION( "whitespace-only final field" )
    {
        std::vector<std::string> values;
        int column = 0;
        REQUIRE( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "   ", 3, column, values ) ==
                 mx::error_t::noerror );
        REQUIRE( values == std::vector<std::string>{ "" } );
    }

    SECTION( "invalid whitespace-only numeric final field" )
    {
        std::vector<double> values;
        int column = 0;
        REQUIRE( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "   ", 3, column, values ) ==
                 mx::error_t::invalidarg );
    }

    SECTION( "missing comma-delimited field" )
    {
        std::vector<int> values;
        int column = 0;
        REQUIRE( mx::ioutils::readcol<mx::ioutils::readColCommaDelim, mx::verbose::vv>( ",value", 6, column, values ) ==
                 mx::error_t::noerror );
        REQUIRE( values == std::vector<int>{ -99 } );
    }

    SECTION( "terminated direct field" )
    {
        std::vector<std::string> values;
        int column = 0;
        REQUIRE(
            mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "value\n", 6, column, values ) ==
            mx::error_t::noerror );
        REQUIRE( values == std::vector<std::string>{ "value" } );
    }

    SECTION( "skipped field" )
    {
        const std::string filename{ "/tmp/readcol_test_skip.dat" };
        std::ofstream fout{ filename };
        fout << "first ignored 2.5\n";
        fout.close();

        std::vector<std::string> names;
        mx::ioutils::skipCol skipped;
        std::vector<double> values;
        REQUIRE( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename,
                                                                                            names,
                                                                                            skipped,
                                                                                            values ) ==
                 mx::error_t::noerror );
        REQUIRE( names == std::vector<std::string>{ "first" } );
        REQUIRE( values == std::vector<double>{ 2.5 } );
    }
}

/** \brief Verifies per-column conversion exceptions follow the configured propagation policy.
 *
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Handling column conversion exceptions", "[ioutils::readColumns]" )
{
    throwingColumn values;
    int column = 3;

    SECTION( "always-translated conversion exceptions" )
    {
        injectedFailure = injectedReadColumnsFailure::invalidArgument;
        REQUIRE( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "value", 5, column, values ) ==
                 mx::error_t::std_invalid_argument );

        injectedFailure = injectedReadColumnsFailure::outOfRange;
        REQUIRE( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "value", 5, column, values ) ==
                 mx::error_t::std_out_of_range );
    }

    SECTION( "allocation failure" )
    {
        injectedFailure = injectedReadColumnsFailure::badAlloc;
        REQUIRE_THROWS_AS(
            ( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "value", 5, column, values ) ),
            mx::exception<mx::verbose::vv> );
        REQUIRE( mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readcol<
                     mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColSpaceDelim,
                     mx::verbose::vv>( "value", 5, column, values ) == mx::error_t::std_bad_alloc );
    }

    SECTION( "standard failure" )
    {
        injectedFailure = injectedReadColumnsFailure::standard;
        REQUIRE_THROWS_AS(
            ( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "value", 5, column, values ) ),
            mx::exception<mx::verbose::vv> );
        REQUIRE( mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readcol<
                     mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColSpaceDelim,
                     mx::verbose::vv>( "value", 5, column, values ) == mx::error_t::std_exception );
    }

    SECTION( "unknown failure" )
    {
        injectedFailure = injectedReadColumnsFailure::unknown;
        REQUIRE_THROWS_AS(
            ( mx::ioutils::readcol<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( "value", 5, column, values ) ),
            mx::exception<mx::verbose::vv> );
        REQUIRE( mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readcol<
                     mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColSpaceDelim,
                     mx::verbose::vv>( "value", 5, column, values ) == mx::error_t::exception );
    }
}

/** \brief Verifies line-read exceptions follow the configured propagation policy.
 *
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Handling line-read exceptions", "[ioutils::readColumns]" )
{
    const std::string filename{ "/tmp/readcol_test_line_exceptions.dat" };
    std::ofstream fout{ filename };
    fout << "value 1.5\n";
    fout.close();

    readColumnsHookGuard guard;
    failingOperation = mx::ioutils::readColumnsDetail::operation::afterReadLine;
    std::vector<std::string> names;
    std::vector<double> values;

    SECTION( "allocation failure" )
    {
        injectedFailure = injectedReadColumnsFailure::badAlloc;
        REQUIRE_THROWS_AS(
            ( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename, names, values ) ),
            mx::exception<mx::verbose::vv> );
        REQUIRE( mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColumns<
                     mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColSpaceDelim,
                     mx::verbose::vv>( filename, names, values ) == mx::error_t::std_bad_alloc );
    }

    SECTION( "standard failure" )
    {
        injectedFailure = injectedReadColumnsFailure::standard;
        REQUIRE_THROWS_AS(
            ( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename, names, values ) ),
            mx::exception<mx::verbose::vv> );
        REQUIRE( mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColumns<
                     mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColSpaceDelim,
                     mx::verbose::vv>( filename, names, values ) == mx::error_t::std_exception );
    }

    SECTION( "unknown failure" )
    {
        injectedFailure = injectedReadColumnsFailure::unknown;
        REQUIRE_THROWS_AS(
            ( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename, names, values ) ),
            mx::exception<mx::verbose::vv> );
        REQUIRE( mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColumns<
                     mx::ioutils::MXLIBTEST_CATCH_ALL_EXCEPTIONS_ns::readColSpaceDelim,
                     mx::verbose::vv>( filename, names, values ) == mx::error_t::exception );
    }
}

/** \brief Verifies open, read, and close stream errors retain errno or use the documented fallback code.
 *
 * \ingroup readColumns_unit_tests
 */
TEST_CASE( "Handling column file stream errors", "[ioutils::readColumns]" )
{
    const std::string filename{ "/tmp/readcol_test_stream_errors.dat" };
    std::ofstream fout{ filename };
    fout << "value 1.5\n";
    fout.close();

    readColumnsHookGuard guard;
    std::vector<std::string> names;
    std::vector<double> values;

    for( auto operation : { mx::ioutils::readColumnsDetail::operation::afterOpen,
                            mx::ioutils::readColumnsDetail::operation::afterReadLine,
                            mx::ioutils::readColumnsDetail::operation::afterClose } )
    {
        failingOperation = operation;

        injectedFailure = injectedReadColumnsFailure::streamErrno;
        REQUIRE( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename, names, values ) ==
                 mx::error_t::eio );

        injectedFailure = injectedReadColumnsFailure::streamFallback;
        const mx::error_t expected =
            operation == mx::ioutils::readColumnsDetail::operation::afterOpen       ? mx::error_t::fileoerr
            : operation == mx::ioutils::readColumnsDetail::operation::afterReadLine ? mx::error_t::filererr
                                                                                    : mx::error_t::filecerr;
        REQUIRE( mx::ioutils::readColumns<mx::ioutils::readColSpaceDelim, mx::verbose::vv>( filename, names, values ) ==
                 expected );
    }
}

} // namespace readColumnsTest
} // namespace ioutilsTest
} // namespace unitTest
