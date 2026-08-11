/** \file fitsHeader_test.cpp
 * \brief Tests FITS-header value extraction.
 */
#include "../../../catch2/catch.hpp"

using namespace Catch::Matchers;

#include "../../../../include/ioutils/fits/fitsHeader.hpp"
using namespace mx::fits;

#include "../../../../include/improc/eigenImage.hpp"

namespace mx
{
namespace unitTest
{
namespace fitsTest
{
namespace fitsHeaderTest
{

/// Test functionality of headersToValues
/**
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Test functionality of headersToValues", "[ioutils::fits::fitsHeader]" )
{
// to force doxygen linking
#ifdef MXLIB_DOXYGEN_PROTECTED_REF
#endif

    SECTION( "a vector of floating point values, no errors" )
    {
        std::vector<fitsHeader<mx::verbose::d>> fhs( 5 );

        fhs[0].append<float>( "KEYTEST", 1, "test comment" );
        fhs[1].append<float>( "KEYTEST", 2, "test comment" );
        fhs[2].append<float>( "KEYTEST", 3, "test comment" );
        fhs[3].append<float>( "KEYTEST", 4, "test comment" );
        fhs[4].append<float>( "KEYTEST", 5, "test comment" );

        std::vector<float> vals;
        std::vector<size_t> bad;
        mx::error_t errc = headersToValues( vals, bad, fhs, "KEYTEST" );

        REQUIRE( !errc );

        REQUIRE( vals[0] == 1 );
        REQUIRE( vals[1] == 2 );
        REQUIRE( vals[2] == 3 );
        REQUIRE( vals[3] == 4 );
        REQUIRE( vals[4] == 5 );
    }

    SECTION( "a vector of floating point values, no errors" )
    {
        std::vector<fitsHeader<mx::verbose::d>> fhs( 5 );

        fhs[0].append<float>( "KEYTEST", 1, "test comment" );
        fhs[1].append<float>( "KEYTEST", 2, "test comment" );
        // fhs[2] is not set
        fhs[3].append<float>( "KEYTEST", 4, "test comment" );
        fhs[4].append<float>( "KEYTEST", 5, "test comment" );

        std::vector<float> vals;
        std::vector<size_t> bad;
        mx::error_t errc = headersToValues( vals, bad, fhs, "KEYTEST" );

        REQUIRE( !!errc );

        REQUIRE( vals[0] == 1 );
        REQUIRE( vals[1] == 2 );
        REQUIRE( vals[2] == std::numeric_limits<float>::max() );
        REQUIRE( vals[3] == 4 );
        REQUIRE( vals[4] == 5 );
    }
}

/** \brief Verifies metadata, repeated provenance cards, and copied headers used by HCI output writers.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Appending HCI-style FITS provenance", "[ioutils::fits::fitsHeader][hciReduce]" )
{
    fitsHeader<mx::verbose::vv> header;
    REQUIRE( header.append( "MJD", 60000.25, "coadd midpoint" ) == mx::error_t::noerror );
    REQUIRE( header.append( "QFILE", "quality.dat", "quality file" ) == mx::error_t::noerror );
    REQUIRE( header.append( "", fitsCommentType(), "HCI parameters" ) == mx::error_t::noerror );
    REQUIRE( header.append( "", fitsCommentType(), "HCI parameters continued" ) == mx::error_t::noerror );
    REQUIRE( header.append( "HISTORY", fitsHistoryType(), "coadded frame001.fits" ) == mx::error_t::noerror );
    REQUIRE( header.append( "HISTORY", fitsHistoryType(), "coadded frame002.fits" ) == mx::error_t::noerror );

    fitsHeader<mx::verbose::vv> source;
    REQUIRE( source.append( "ORIGIN", "hciReduce", "producer" ) == mx::error_t::noerror );
    REQUIRE( header.append( source ) == mx::error_t::noerror );

    REQUIRE( header.size() == 7 );
    REQUIRE( header["MJD"].value<double>() == 60000.25 );
    REQUIRE( header["QFILE"].String() == "quality.dat" );
    REQUIRE( header["ORIGIN"].String() == "hciReduce" );
    REQUIRE( header.append( "MJD", 60001.0, "duplicate" ) == mx::error_t::invalidarg );
    REQUIRE( header.size() == 7 );
}

/** \brief Verifies HCI coadd updates of midpoint, start, end, and delta metadata cards.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Updating HCI coadd time and keyword cards", "[ioutils::fits::fitsHeader][hciReduce]" )
{
    fitsHeader<mx::verbose::vv> header;
    REQUIRE( header.append( "MJD-OBS", "2000-01-01T00:00:00.000000000Z", "coadd midpoint" ) ==
             mx::error_t::noerror );
    REQUIRE( header.append( "ANGLE", 10.0, "derotation angle" ) == mx::error_t::noerror );

    REQUIRE( header["MJD-OBS"].value( "2000-01-01T12:00:00.000000000Z" ) == mx::error_t::noerror );
    REQUIRE( header["START MJD-OBS"].value( "2000-01-01T00:00:00.000000000Z" ) == mx::error_t::noerror );
    REQUIRE( header["END MJD-OBS"].value( "2000-01-01T12:00:00.000000000Z" ) == mx::error_t::noerror );
    REQUIRE( header["DELTA MJD-OBS"].value( 0.5 ) == mx::error_t::noerror );

    REQUIRE( header["ANGLE"].value( 15.0 ) == mx::error_t::noerror );
    REQUIRE( header["START ANGLE"].value( 10.0 ) == mx::error_t::noerror );
    REQUIRE( header["END ANGLE"].value( 20.0 ) == mx::error_t::noerror );
    REQUIRE( header["DELTA ANGLE"].value( 10.0 ) == mx::error_t::noerror );

    REQUIRE( header["MJD-OBS"].String() == "2000-01-01T12:00:00.000000000Z" );
    REQUIRE( header["START MJD-OBS"].String() == "2000-01-01T00:00:00.000000000Z" );
    REQUIRE( header["END MJD-OBS"].String() == "2000-01-01T12:00:00.000000000Z" );
    REQUIRE( header["DELTA MJD-OBS"].value<double>() == 0.5 );
    REQUIRE( header["ANGLE"].value<double>() == 15.0 );
    REQUIRE( header["START ANGLE"].value<double>() == 10.0 );
    REQUIRE( header["END ANGLE"].value<double>() == 20.0 );
    REQUIRE( header["DELTA ANGLE"].value<double>() == 10.0 );
}

/** \brief Verifies FITS HISTORY provenance for HCI final-image outputs.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Appending HCI Git status provenance", "[ioutils::fits::fitsHeader][hciReduce]" )
{
    fitsHeader<mx::verbose::vv> cleanHeader;
    fitsHeaderGitStatus( cleanHeader, "mxlib_uncomp", "abc123", 0 );

    std::vector<std::string> cleanHistory;
    for( auto &card : cleanHeader )
    {
        if( card.type() == fitsType<fitsHistoryType>() )
        {
            cleanHistory.push_back( card.comment() );
        }
    }

    REQUIRE( cleanHistory == std::vector<std::string>{ "Git status for mxlib_uncomp:", "   sha1=abc123" } );

    fitsHeader<mx::verbose::vv> modifiedHeader;
    fitsHeaderGitStatus( modifiedHeader, "mxlib_uncomp", "abc123", 1 );

    std::vector<std::string> modifiedHistory;
    for( auto &card : modifiedHeader )
    {
        if( card.type() == fitsType<fitsHistoryType>() )
        {
            modifiedHistory.push_back( card.comment() );
        }
    }

    REQUIRE( modifiedHistory ==
             std::vector<std::string>{ "Git status for mxlib_uncomp:", "   sha1=abc123, modified" } );
}

} // namespace fitsHeaderTest
} // namespace fitsTest
} // namespace unitTest
} // namespace mx
