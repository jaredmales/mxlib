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

/// \cond
template <class headerT>
std::vector<std::string> keywords( headerT &header )
{
    std::vector<std::string> result;
    for( auto &card : header )
    {
        result.push_back( card.keyword() );
    }
    return result;
}

enum class mutationFailureKind
{
    returnedError,
    badAllocation,
    standardException
};

fitsHeaderDetail::mutation failingMutation = fitsHeaderDetail::mutation::appendList;
mutationFailureKind mutationFailure = mutationFailureKind::returnedError;

mx::error_t failMutation( fitsHeaderDetail::mutation operation )
{
    if( operation != failingMutation )
    {
        return mx::error_t::noerror;
    }

    if( mutationFailure == mutationFailureKind::badAllocation )
    {
        throw std::bad_alloc();
    }
    if( mutationFailure == mutationFailureKind::standardException )
    {
        throw std::runtime_error( "fitsHeader mutation failure" );
    }
    return mx::error_t::error;
}

class mutationHookGuard
{
  public:
    mutationHookGuard() : m_saved( fitsHeaderDetail::mutationHook() )
    {
        fitsHeaderDetail::mutationHook() = failMutation;
    }

    ~mutationHookGuard()
    {
        fitsHeaderDetail::mutationHook() = m_saved;
    }

  private:
    fitsHeaderDetail::mutationHookT m_saved;
};
/// \endcond

/** \brief Verifies empty-state, append, lookup, copy, assignment, and clearing behavior.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Managing the fitsHeader lifecycle", "[ioutils::fits::fitsHeader]" )
{
    fitsHeader<> header;
    REQUIRE( header.empty() );
    REQUIRE( header.size() == 0 );
    REQUIRE( header.begin() == header.end() );
    REQUIRE( header.iterator( "MISSING" ) == header.end() );
    REQUIRE( header.count( "MISSING" ) == 0 );

    REQUIRE( header.append( "UNKNOWN" ) == mx::error_t::noerror );
    REQUIRE( header.append( "INTEGER", 42 ) == mx::error_t::noerror );
    REQUIRE( header.append( "STRING", "value", "a comment" ) == mx::error_t::noerror );
    REQUIRE( header.append( fitsHeaderCard<>( "FLOAT", 2.5F, "float comment" ) ) == mx::error_t::noerror );
    REQUIRE_FALSE( header.empty() );
    REQUIRE( header.size() == 4 );
    REQUIRE( header.iterator( "INTEGER" ) != header.end() );
    REQUIRE( header.iterator( "INTEGER" )->keyword() == "INTEGER" );
    REQUIRE( header.iterator( "INTEGER" )->type() == fitsType<int>() );
    REQUIRE( header.iterator( "INTEGER" )->valueGood() );
    REQUIRE( header["INTEGER"].value<int>() == 42 );
    REQUIRE( header["STRING"].String() == "value" );
    REQUIRE( header["STRING"].comment() == "a comment" );
    REQUIRE( header["FLOAT"].value<float>() == 2.5F );
    REQUIRE( header.append( "INTEGER", 43 ) == mx::error_t::invalidarg );
    REQUIRE( header.size() == 4 );

    fitsHeader<> copied( header );
    REQUIRE( copied.size() == 4 );
    REQUIRE( copied["INTEGER"].value<int>() == 42 );
    REQUIRE( copied["INTEGER"].value( 7 ) == mx::error_t::noerror );
    REQUIRE( header["INTEGER"].value<int>() == 42 );

    fitsHeader<> assigned;
    REQUIRE( assigned.append( "OLD", 1 ) == mx::error_t::noerror );
    assigned = header;
    REQUIRE( assigned.size() == 4 );
    REQUIRE( assigned.count( "OLD" ) == 0 );
    REQUIRE( assigned["STRING"].String() == "value" );
    assigned = assigned;
    REQUIRE( assigned.size() == 4 );
    REQUIRE( assigned["INTEGER"].value<int>() == 42 );

    REQUIRE( assigned.clear() == mx::error_t::noerror );
    REQUIRE( assigned.clear() == mx::error_t::noerror );
    REQUIRE( assigned.empty() );
    REQUIRE( assigned.iterator( "INTEGER" ) == assigned.end() );
}

/** \brief Verifies mutable and const fitsHeader keyword access, including missing-card sentinels.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Accessing fitsHeader cards by keyword", "[ioutils::fits::fitsHeader]" )
{
    fitsHeader<> header;
    REQUIRE( header.append( "EXISTING", 3, "existing card" ) == mx::error_t::noerror );
    REQUIRE( header["EXISTING"].value<int>() == 3 );

    fitsHeaderCard<> &created = header["CREATED"];
    REQUIRE( created.keyword() == "CREATED" );
    REQUIRE( created.type() == fitsType<fitsUnknownType>() );
    REQUIRE( header.size() == 2 );
    REQUIRE( header.count( "CREATED" ) == 1 );

    const fitsHeader<> &constHeader = header;
    REQUIRE( constHeader["EXISTING"].keyword() == "EXISTING" );
    const fitsHeaderCard<> &missing = constHeader["MISSING"];
    REQUIRE( missing.keyword().empty() );
    const_cast<fitsHeaderCard<> &>( missing ).keyword( "changed sentinel" );
    REQUIRE( constHeader["ALSO MISSING"].keyword().empty() );
    REQUIRE( header.size() == 2 );
}

/** \brief Verifies fitsHeader erasure by keyword and exact list iterator.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Erasing fitsHeader cards", "[ioutils::fits::fitsHeader]" )
{
    fitsHeader<> header;
    REQUIRE( header.append( "FIRST", 1 ) == mx::error_t::noerror );
    REQUIRE( header.append( "COMMENT", fitsCommentType(), "same" ) == mx::error_t::noerror );
    REQUIRE( header.append( "COMMENT", fitsCommentType(), "same" ) == mx::error_t::noerror );
    REQUIRE( header.append( "HISTORY", fitsHistoryType(), "event" ) == mx::error_t::noerror );
    REQUIRE( header.append( "LAST", 2 ) == mx::error_t::noerror );

    REQUIRE( header.count( "COMMENT" ) == 2 );
    REQUIRE( header.erase( "FIRST" ) == mx::error_t::noerror );
    REQUIRE( header.count( "FIRST" ) == 0 );
    REQUIRE( header.erase( "FIRST" ) == mx::error_t::notfound );
    REQUIRE( header.erase( "COMMENT" ) == mx::error_t::invalidarg );
    REQUIRE( header.erase( "HISTORY" ) == mx::error_t::invalidarg );
    REQUIRE( header.erase( header.end() ) == mx::error_t::invalidarg );

    auto comment = header.begin();
    REQUIRE( comment->keyword() == "COMMENT" );
    REQUIRE( header.erase( comment ) == mx::error_t::noerror );
    REQUIRE( header.count( "COMMENT" ) == 1 );
    REQUIRE( header.size() == 3 );

    auto history = header.iterator( "HISTORY" );
    REQUIRE( history != header.end() );
    REQUIRE( header.erase( history ) == mx::error_t::noerror );
    REQUIRE( header.count( "HISTORY" ) == 0 );
    REQUIRE( keywords( header ) == std::vector<std::string>{ "COMMENT", "LAST" } );
}

/** \brief Verifies ordered fitsHeader insertion and duplicate-key rejection.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Inserting ordered fitsHeader cards", "[ioutils::fits::fitsHeader]" )
{
    fitsHeader<> header;
    REQUIRE( header.append( "SECOND", 2 ) == mx::error_t::noerror );
    REQUIRE( header.append( "FOURTH", 4 ) == mx::error_t::noerror );

    REQUIRE( header.insert_before( header.begin(), "FIRST", 1, "first" ) == mx::error_t::noerror );
    REQUIRE( header.insert_after( header.iterator( "SECOND" ), "THIRD", 3, "third" ) == mx::error_t::noerror );
    REQUIRE( header.insert_before( header.end(), fitsHeaderCard<>( "FIFTH", 5, "" ) ) == mx::error_t::noerror );
    REQUIRE( keywords( header ) == std::vector<std::string>{ "FIRST", "SECOND", "THIRD", "FOURTH", "FIFTH" } );
    REQUIRE( header["FIRST"].comment() == "first" );
    REQUIRE( header["THIRD"].comment() == "third" );

    REQUIRE( header.insert_before( header.begin(), "FIRST", 10 ) == mx::error_t::invalidarg );
    REQUIRE( header.insert_after( header.iterator( "SECOND" ), "THIRD", 30 ) == mx::error_t::invalidarg );
    REQUIRE( header.insert_after( header.end(), "INVALID", 0 ) == mx::error_t::invalidarg );
    REQUIRE( header.count( "INVALID" ) == 0 );

    REQUIRE( header.insert_before( header.begin(), fitsHeaderCard<>( "COMMENT", fitsCommentType(), "before" ) ) ==
             mx::error_t::noerror );
    REQUIRE( header.insert_after( header.begin(), fitsHeaderCard<>( "COMMENT", fitsCommentType(), "after" ) ) ==
             mx::error_t::noerror );
    REQUIRE( header.count( "COMMENT" ) == 2 );
}

/** \brief Verifies removal of standard top-of-header FITS cards and boilerplate comments.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Removing standard FITS header cards", "[ioutils::fits::fitsHeader]" )
{
    fitsHeader<> header;
    for( const std::string &keyword :
         { "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND", "BZERO", "BSCALE", "LONGSTRN" } )
    {
        REQUIRE( header.append( keyword, 1 ) == mx::error_t::noerror );
    }
    REQUIRE( header.append( "COMMENT", fitsCommentType(), "FITS (Flexible Image Transport System)" ) ==
             mx::error_t::noerror );
    REQUIRE( header.append( "COMMENT", fitsCommentType(), "Astronomy and Astrophysics' convention" ) ==
             mx::error_t::noerror );
    REQUIRE( header.append( "COMMENT", fitsCommentType(), "preserve this comment" ) == mx::error_t::noerror );
    REQUIRE( header.append( "HISTORY", fitsHistoryType(), "preserve this history" ) == mx::error_t::noerror );
    REQUIRE( header.append( "OBJECT", "target" ) == mx::error_t::noerror );

    REQUIRE( header.eraseStandardTop() == mx::error_t::noerror );
    REQUIRE( keywords( header ) == std::vector<std::string>{ "COMMENT", "HISTORY", "OBJECT" } );
    REQUIRE( header.begin()->comment() == "preserve this comment" );
}

/** \brief Verifies fitsHeader concatenation, duplicate propagation, self-append safety, and CONTINUE validation.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Appending complete FITS headers", "[ioutils::fits::fitsHeader]" )
{
    fitsHeader<> header;
    fitsHeader<> source;
    REQUIRE( header.append( "FIRST", 1 ) == mx::error_t::noerror );
    REQUIRE( source.append( "SECOND", 2 ) == mx::error_t::noerror );
    REQUIRE( source.append( "THIRD", 3 ) == mx::error_t::noerror );
    REQUIRE( header.append( source ) == mx::error_t::noerror );
    REQUIRE( keywords( header ) == std::vector<std::string>{ "FIRST", "SECOND", "THIRD" } );

    fitsHeader<> duplicate;
    REQUIRE( duplicate.append( "THIRD", 30 ) == mx::error_t::noerror );
    REQUIRE( duplicate.append( "FOURTH", 4 ) == mx::error_t::noerror );
    REQUIRE( header.append( duplicate ) == mx::error_t::invalidarg );
    REQUIRE( header.count( "FOURTH" ) == 0 );

    REQUIRE( header.append( header ) == mx::error_t::invalidarg );
    REQUIRE( header.size() == 3 );

    fitsHeader<> empty;
    fitsHeaderCard<> continuation( "CONTINUE", "'continued text'", fitsType<fitsContinueType>() );
    REQUIRE( empty.append( continuation ) == mx::error_t::invalidarg );
    REQUIRE( empty.empty() );

    fitsHeader<> continued;
    REQUIRE( continued.append( "LONG", "first&", "initial comment" ) == mx::error_t::noerror );
    fitsHeaderCard<> finalContinuation( "CONTINUE", "", fitsType<fitsContinueType>(), "'second' / final comment" );
    REQUIRE( continued.append( finalContinuation ) == mx::error_t::noerror );
    REQUIRE( continued.size() == 1 );
    REQUIRE( continued["LONG"].String() == "firstsecond" );
    REQUIRE( continued["LONG"].comment() == "final comment" );
}

/** \brief Verifies fitsHeader propagates allocation failures and rolls back partially completed mutations.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Propagating fitsHeader mutation failures", "[ioutils::fits::fitsHeader]" )
{
    mutationHookGuard guard;

    failingMutation = fitsHeaderDetail::mutation::appendList;
    mutationFailure = mutationFailureKind::returnedError;
    fitsHeader<> lookupHeader;
    fitsHeaderCard<> &missing = lookupHeader["FAILED"];
    REQUIRE( missing.keyword().empty() );
    REQUIRE( lookupHeader.empty() );

    for( mutationFailureKind failure : { mutationFailureKind::returnedError,
                                         mutationFailureKind::badAllocation,
                                         mutationFailureKind::standardException } )
    {
        for( fitsHeaderDetail::mutation operation :
             { fitsHeaderDetail::mutation::appendList, fitsHeaderDetail::mutation::appendMap } )
        {
            fitsHeader<> header;
            failingMutation = operation;
            mutationFailure = failure;
            if( failure == mutationFailureKind::returnedError )
            {
                REQUIRE( header.append( "CARD", 1 ) == mx::error_t::error );
            }
            else
            {
                REQUIRE_THROWS( header.append( "CARD", 1 ) );
            }
            REQUIRE( header.empty() );
            REQUIRE( header.count( "CARD" ) == 0 );
        }
    }

    for( mutationFailureKind failure : { mutationFailureKind::returnedError,
                                         mutationFailureKind::badAllocation,
                                         mutationFailureKind::standardException } )
    {
        for( fitsHeaderDetail::mutation operation :
             { fitsHeaderDetail::mutation::insertBeforeList, fitsHeaderDetail::mutation::insertBeforeMap } )
        {
            fitsHeaderDetail::resetMutationHook();
            fitsHeader<> header;
            REQUIRE( header.append( "ANCHOR", 1 ) == mx::error_t::noerror );
            fitsHeaderDetail::mutationHook() = failMutation;
            failingMutation = operation;
            mutationFailure = failure;
            if( failure == mutationFailureKind::returnedError )
            {
                REQUIRE( header.insert_before( header.begin(), "CARD", 2 ) == mx::error_t::error );
            }
            else
            {
                REQUIRE_THROWS( header.insert_before( header.begin(), "CARD", 2 ) );
            }
            REQUIRE( keywords( header ) == std::vector<std::string>{ "ANCHOR" } );
            REQUIRE( header.count( "CARD" ) == 0 );
        }
    }

    for( mutationFailureKind failure : { mutationFailureKind::returnedError,
                                         mutationFailureKind::badAllocation,
                                         mutationFailureKind::standardException } )
    {
        for( fitsHeaderDetail::mutation operation :
             { fitsHeaderDetail::mutation::insertAfterList, fitsHeaderDetail::mutation::insertAfterMap } )
        {
            fitsHeaderDetail::resetMutationHook();
            fitsHeader<> header;
            REQUIRE( header.append( "ANCHOR", 1 ) == mx::error_t::noerror );
            fitsHeaderDetail::mutationHook() = failMutation;
            failingMutation = operation;
            mutationFailure = failure;
            if( failure == mutationFailureKind::returnedError )
            {
                REQUIRE( header.insert_after( header.begin(), "CARD", 2 ) == mx::error_t::error );
            }
            else
            {
                REQUIRE_THROWS( header.insert_after( header.begin(), "CARD", 2 ) );
            }
            REQUIRE( keywords( header ) == std::vector<std::string>{ "ANCHOR" } );
            REQUIRE( header.count( "CARD" ) == 0 );
        }
    }
}

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
    REQUIRE( header.append( "MJD-OBS", "2000-01-01T00:00:00.000000000Z", "coadd midpoint" ) == mx::error_t::noerror );
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

/** \brief Verifies the scalar and string FITS metadata types emitted by HCI stdFitsHeader.
 *
 * \ingroup fitsHeader_unit_tests
 */
TEST_CASE( "Appending HCI reduction parameter metadata", "[ioutils::fits::fitsHeader][hciReduce]" )
{
    fitsHeader<mx::verbose::vv> header;
    char modes[] = "1,5,10";

    REQUIRE( header.append<int>( "NUMIMS", 5, "images processed" ) == mx::error_t::noerror );
    REQUIRE( header.append<float>( "QTHRESH", 0.75F, "quality threshold" ) == mx::error_t::noerror );
    REQUIRE( header.append<bool>( "PREPROC MASK", true, "apply preprocessing mask" ) == mx::error_t::noerror );
    REQUIRE( header.append<std::string>( "COADMTHD", "mean", "coadd method" ) == mx::error_t::noerror );
    REQUIRE( header.append<char *>( "NMODES", modes, "number of modes" ) == mx::error_t::noerror );

    REQUIRE( header["NUMIMS"].value<int>() == 5 );
    REQUIRE_THAT( header["QTHRESH"].value<float>(), WithinAbs( 0.75F, 1e-6F ) );
    REQUIRE( header["PREPROC MASK"].type() == fitsType<bool>() );
    REQUIRE( header["PREPROC MASK"].valueGood() );
    REQUIRE( header["COADMTHD"].String() == "mean" );
    REQUIRE( header["NMODES"].String() == "1,5,10" );
}

} // namespace fitsHeaderTest
} // namespace fitsTest
} // namespace unitTest
} // namespace mx
