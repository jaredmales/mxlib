/** \file fitsHeaderCard_test.cpp
 * \brief Tests FITS header-card parsing and formatting.
 */
#include "../../../catch2/catch.hpp"

#include <filesystem>
#include <unistd.h>

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ioutils/fits/fitsHeaderCard.hpp"
using namespace mx::fits;

namespace mx
{
namespace unitTest
{
namespace fitsTest
{
namespace fitsHeaderCardTest
{

/** \cond */
template <class verboseT = mx::verbose::d>
class testHeaderCard : public fitsHeaderCard<verboseT>
{
  public:
    using fitsHeaderCard<verboseT>::fitsHeaderCard;

    void forceTypeWithValue( int type )
    {
        this->value( 1 );
        this->m_type = type;
        this->m_valueGood = true;
        this->m_valueStrGood = false;
    }
};

class temporaryFitsFile
{
  public:
    temporaryFitsFile()
    {
        m_path = std::filesystem::temp_directory_path() /
                 ( "mxlib-fitsHeaderCard-" + std::to_string( static_cast<long long>( getpid() ) ) + ".fits" );
        std::filesystem::remove( m_path );

        int status = 0;
        std::string createPath = "!" + m_path.string();
        fits_create_file( &m_file, createPath.c_str(), &status );
        if( status != 0 )
        {
            throw std::runtime_error( "could not create temporary FITS file" );
        }

        long axes[1] = { 1 };
        fits_create_img( m_file, BYTE_IMG, 1, axes, &status );
        if( status != 0 )
        {
            throw std::runtime_error( "could not create temporary FITS image" );
        }
    }

    ~temporaryFitsFile()
    {
        if( m_file != nullptr )
        {
            int status = 0;
            fits_close_file( m_file, &status );
        }
        std::error_code error;
        std::filesystem::remove( m_path, error );
    }

    fitsfile *get()
    {
        return m_file;
    }

  private:
    std::filesystem::path m_path;
    fitsfile *m_file{ nullptr };
};
/** \endcond */

/** \brief Verifies construction, copying, assignment, and scalar access for every supported numeric card type.
 *
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Round-tripping native FITS header card values", "[ioutils::fits::fitsHeaderCard]" )
{
    fitsHeaderCard<> empty;
    REQUIRE( empty.keyword().empty() );
    REQUIRE( empty.type() == fitsType<fitsUnknownType>() );
    REQUIRE_FALSE( empty.valueGood() );
    REQUIRE_FALSE( empty.valueStrGood() );

    fitsHeaderCard<> character( "CHAR", static_cast<char>( 12 ), "char" );
    fitsHeaderCard<> unsignedCharacter( "UCHAR", static_cast<unsigned char>( 13 ), "uchar" );
    fitsHeaderCard<> shortInteger( "SHORT", static_cast<short>( 14 ), "short" );
    fitsHeaderCard<> unsignedShort( "USHORT", static_cast<unsigned short>( 15 ), "ushort" );
    fitsHeaderCard<> integer( "INT", 16, "int" );
    fitsHeaderCard<> unsignedInteger( "UINT", 17U, "uint" );
    fitsHeaderCard<> longInteger( "LONG", 18L, "long" );
    fitsHeaderCard<> unsignedLong( "ULONG", 19UL, "ulong" );
    fitsHeaderCard<> longLongInteger( "LONGLONG", 20LL, "long long" );
    fitsHeaderCard<> unsignedLongLong( "ULONGLONG", 21ULL, "unsigned long long" );
    fitsHeaderCard<> single( "FLOAT", 22.5F, "float" );
    fitsHeaderCard<> real( "DOUBLE", 23.5, "double" );
    fitsHeaderCard<> complexSingle( "COMPLEXF", std::complex<float>( 1.0F, -2.0F ), "complex float" );
    fitsHeaderCard<> complexReal( "COMPLEXD", std::complex<double>( 3.0, -4.0 ), "complex double" );
    fitsHeaderCard<> extended( "LONGDOUBLE", static_cast<long double>( 5.0 ), "long double" );
    fitsHeaderCard<> complexExtended(
        "COMPLEXLD",
        std::complex<long double>( static_cast<long double>( 6.0 ), static_cast<long double>( -7.0 ) ),
        "complex long double" );
#ifdef HAS_QUAD
    fitsHeaderCard<> quad( "QUAD", static_cast<__float128>( 8.0 ), "quad" );
    fitsHeaderCard<> complexQuad(
        "COMPLEXQ",
        std::complex<__float128>( static_cast<__float128>( 9.0 ), static_cast<__float128>( -10.0 ) ),
        "complex quad" );
#endif

    mx::error_t error = mx::error_t::error;
    REQUIRE( character.Char( &error ) == 12 );
    REQUIRE( error == mx::error_t::noerror );
    REQUIRE( unsignedCharacter.UChar() == 13 );
    REQUIRE( shortInteger.Short() == 14 );
    REQUIRE( unsignedShort.UShort() == 15 );
    REQUIRE( integer.Int() == 16 );
    REQUIRE( unsignedInteger.UInt() == 17U );
    REQUIRE( longInteger.Long() == 18L );
    REQUIRE( unsignedLong.ULong() == 19UL );
    REQUIRE( longLongInteger.LongLong() == 20LL );
    REQUIRE( unsignedLongLong.ULongLong() == 21ULL );
    REQUIRE( single.Float() == 22.5F );
    REQUIRE( real.Double() == 23.5 );

    REQUIRE( character.valueStr() == "12" );
    REQUIRE( unsignedCharacter.valueStr() == "13" );
    REQUIRE( shortInteger.valueStr() == "14" );
    REQUIRE( unsignedShort.valueStr() == "15" );
    REQUIRE( integer.valueStr() == "16" );
    REQUIRE( unsignedInteger.valueStr() == "17" );
    REQUIRE( longInteger.valueStr() == "18" );
    REQUIRE( unsignedLong.valueStr() == "19" );
    REQUIRE( longLongInteger.valueStr() == "20" );
    REQUIRE( unsignedLongLong.valueStr() == "21" );
    REQUIRE( single.valueStr() == "22.5" );
    REQUIRE( real.valueStr() == "23.5" );
    REQUIRE( complexSingle.valueStr() == "(1,-2)" );
    REQUIRE( complexReal.valueStr() == "(3,-4)" );
    REQUIRE( extended.valueGood() );
    REQUIRE( complexExtended.valueGood() );
#ifdef HAS_QUAD
    REQUIRE( quad.valueGood() );
    REQUIRE( complexQuad.valueGood() );
#endif

    fitsHeaderCard<> copied( real );
    REQUIRE( copied.keyword() == "DOUBLE" );
    REQUIRE( copied.Double() == 23.5 );
    REQUIRE( copied.comment() == "double" );

    fitsHeaderCard<> assigned( "OLD", 1, "old" );
    assigned = unsignedLongLong;
    REQUIRE( assigned.keyword() == "ULONGLONG" );
    REQUIRE( assigned.ULongLong() == 21ULL );
    REQUIRE( assigned.comment() == "unsigned long long" );
    assigned = assigned;
    REQUIRE( assigned.ULongLong() == 21ULL );

    REQUIRE( assigned.keyword( "RENAMED" ) == mx::error_t::noerror );
    REQUIRE( assigned.comment( "updated" ) == mx::error_t::noerror );
    REQUIRE( assigned.keyword() == "RENAMED" );
    REQUIRE( assigned.comment() == "updated" );
}

/** \brief Verifies parsing a stored FITS value string into every supported numeric scalar type.
 *
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Parsing numeric FITS header card strings", "[ioutils::fits::fitsHeaderCard]" )
{
    fitsHeaderCard<> character( "CHAR", "12", fitsType<char>() );
    fitsHeaderCard<> unsignedCharacter( "UCHAR", "13", fitsType<unsigned char>() );
    fitsHeaderCard<> shortInteger( "SHORT", "14", fitsType<short>() );
    fitsHeaderCard<> unsignedShort( "USHORT", "15", fitsType<unsigned short>() );
    fitsHeaderCard<> integer( "INT", "16", fitsType<int>() );
    fitsHeaderCard<> unsignedInteger( "UINT", "17", fitsType<unsigned int>() );
    fitsHeaderCard<> longInteger( "LONG", "18", fitsType<long>() );
    fitsHeaderCard<> unsignedLong( "ULONG", "19", fitsType<unsigned long>() );
    fitsHeaderCard<> longLongInteger( "LONGLONG", "20", fitsType<long long>() );
    fitsHeaderCard<> unsignedLongLong( "ULONGLONG", "21", fitsType<unsigned long long>() );
    fitsHeaderCard<> single( "FLOAT", "22.5", fitsType<float>() );
    fitsHeaderCard<> real( "DOUBLE", "23.5", fitsType<double>() );

    REQUIRE( character.Char() == 12 );
    REQUIRE( unsignedCharacter.UChar() == 13 );
    REQUIRE( shortInteger.Short() == 14 );
    REQUIRE( unsignedShort.UShort() == 15 );
    REQUIRE( integer.Int() == 16 );
    REQUIRE( unsignedInteger.UInt() == 17U );
    REQUIRE( longInteger.Long() == 18L );
    REQUIRE( unsignedLong.ULong() == 19UL );
    REQUIRE( longLongInteger.LongLong() == 20LL );
    REQUIRE( unsignedLongLong.ULongLong() == 21ULL );
    REQUIRE( single.Float() == 22.5F );
    REQUIRE( real.Double() == 23.5 );

    fitsHeaderCard<> invalid( "BAD", "not-a-number", fitsType<int>() );
    mx::error_t error = mx::error_t::noerror;
    REQUIRE( invalid.Int( &error ) == 0 );
    REQUIRE( error != mx::error_t::noerror );
    REQUIRE_FALSE( invalid.valueGood() );
}

/** \brief Verifies numeric fitsHeaderCard conversions between every supported scalar source and target type.
 *
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Converting numeric FITS header card types", "[ioutils::fits::fitsHeaderCard]" )
{
    fitsHeaderCard<> character( "CHAR", static_cast<char>( 7 ), "" );
    fitsHeaderCard<> unsignedCharacter( "UCHAR", static_cast<unsigned char>( 7 ), "" );
    fitsHeaderCard<> shortInteger( "SHORT", static_cast<short>( 7 ), "" );
    fitsHeaderCard<> unsignedShort( "USHORT", static_cast<unsigned short>( 7 ), "" );
    fitsHeaderCard<> integerSource( "INT", 7, "" );
    fitsHeaderCard<> unsignedInteger( "UINT", 7U, "" );
    fitsHeaderCard<> longInteger( "LONG", 7L, "" );
    fitsHeaderCard<> unsignedLong( "ULONG", 7UL, "" );
    fitsHeaderCard<> longLongInteger( "LONGLONG", 7LL, "" );
    fitsHeaderCard<> unsignedLongLong( "ULONGLONG", 7ULL, "" );
    fitsHeaderCard<> single( "FLOAT", 7.0F, "" );
    fitsHeaderCard<> real( "DOUBLE", 7.0, "" );

    REQUIRE( character.Double() == 7.0 );
    REQUIRE( unsignedCharacter.Double() == 7.0 );
    REQUIRE( shortInteger.Double() == 7.0 );
    REQUIRE( unsignedShort.Double() == 7.0 );
    REQUIRE( integerSource.Double() == 7.0 );
    REQUIRE( unsignedInteger.Double() == 7.0 );
    REQUIRE( longInteger.Double() == 7.0 );
    REQUIRE( unsignedLong.Double() == 7.0 );
    REQUIRE( longLongInteger.Double() == 7.0 );
    REQUIRE( unsignedLongLong.Double() == 7.0 );
    REQUIRE( single.Double() == 7.0 );
    REQUIRE( real.Int() == 7 );

    fitsHeaderCard<> target( "TARGET", 9, "" );
    REQUIRE( target.type( fitsType<unsigned char>() ) == mx::error_t::noerror );
    REQUIRE( target.UChar() == 9 );
    REQUIRE( target.type( fitsType<char>() ) == mx::error_t::noerror );
    REQUIRE( target.Char() == 9 );
    REQUIRE( target.type( fitsType<short>() ) == mx::error_t::noerror );
    REQUIRE( target.Short() == 9 );
    REQUIRE( target.type( fitsType<unsigned short>() ) == mx::error_t::noerror );
    REQUIRE( target.UShort() == 9 );
    REQUIRE( target.type( fitsType<int>() ) == mx::error_t::noerror );
    REQUIRE( target.Int() == 9 );
    REQUIRE( target.type( fitsType<unsigned int>() ) == mx::error_t::noerror );
    REQUIRE( target.UInt() == 9U );
    REQUIRE( target.type( fitsType<long>() ) == mx::error_t::noerror );
    REQUIRE( target.Long() == 9L );
    REQUIRE( target.type( fitsType<unsigned long>() ) == mx::error_t::noerror );
    REQUIRE( target.ULong() == 9UL );
    REQUIRE( target.type( fitsType<long long>() ) == mx::error_t::noerror );
    REQUIRE( target.LongLong() == 9LL );
    REQUIRE( target.type( fitsType<unsigned long long>() ) == mx::error_t::noerror );
    REQUIRE( target.ULongLong() == 9ULL );
    REQUIRE( target.type( fitsType<float>() ) == mx::error_t::noerror );
    REQUIRE( target.Float() == 9.0F );
    REQUIRE( target.type( fitsType<double>() ) == mx::error_t::noerror );
    REQUIRE( target.Double() == 9.0 );
    REQUIRE( target.type( fitsType<double>() ) == mx::error_t::noerror );

    fitsHeaderCard<> textTarget( "TEXT", 11, "" );
    REQUIRE( textTarget.type( fitsType<std::string>() ) == mx::error_t::noerror );
    REQUIRE( textTarget.String() == "11" );

    fitsHeaderCard<> unparsed( "UNPARSED", "12", fitsType<int>() );
    REQUIRE_FALSE( unparsed.valueGood() );
    REQUIRE( unparsed.type( fitsType<double>() ) == mx::error_t::noerror );
    REQUIRE( unparsed.type() == fitsType<double>() );
    REQUIRE_FALSE( unparsed.valueGood() );

    fitsHeaderCard<> unsetType( "UNSETTYPE" );
    REQUIRE( unsetType.type( fitsType<float>() ) == mx::error_t::noerror );
    REQUIRE( unsetType.type() == fitsType<float>() );
}

/** \brief Verifies fitsHeaderCard reports unsupported source and target conversions without undefined values.
 *
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Rejecting unsupported FITS header card conversions", "[ioutils::fits::fitsHeaderCard]" )
{
    mx::error_t error = mx::error_t::noerror;
    testHeaderCard<> source;

    source.forceTypeWithValue( fitsType<std::complex<float>>() );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::notimpl );
    source.forceTypeWithValue( fitsType<std::complex<double>>() );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::notimpl );
    source.forceTypeWithValue( fitsType<fitsCommentType>() );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::invalidarg );
    source.forceTypeWithValue( fitsType<fitsHistoryType>() );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::invalidarg );
    source.forceTypeWithValue( fitsType<fitsContinueType>() );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::invalidarg );
    source.forceTypeWithValue( fitsType<std::string>() );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::invalidarg );
    source.forceTypeWithValue( -12345 );
    REQUIRE( source.Int( &error ) == 0 );
    REQUIRE( error == mx::error_t::invalidarg );

    source.forceTypeWithValue( fitsType<std::complex<float>>() );
    REQUIRE( source.type( fitsType<int>() ) == mx::error_t::notimpl );

    fitsHeaderCard<> target( "TARGET", 1, "" );
    REQUIRE( target.type( fitsType<std::complex<float>>() ) == mx::error_t::notimpl );
    REQUIRE( target.type( fitsType<std::complex<double>>() ) == mx::error_t::notimpl );
    REQUIRE( target.type( fitsType<fitsCommentType>() ) == mx::error_t::invalidarg );
    REQUIRE( target.type( fitsType<fitsHistoryType>() ) == mx::error_t::invalidarg );
    REQUIRE( target.type( fitsType<fitsContinueType>() ) == mx::error_t::invalidarg );
    REQUIRE( target.type( -12345 ) == mx::error_t::invalidarg );

    testHeaderCard<> special;
    for( int type : { fitsType<fitsCommentType>(), fitsType<fitsHistoryType>(), fitsType<fitsContinueType>() } )
    {
        special.forceTypeWithValue( type );
        REQUIRE( special.String( &error ).empty() );
        REQUIRE( error == mx::error_t::noerror );
    }

    special.forceTypeWithValue( fitsType<std::string>() );
    REQUIRE( special.String( &error ).empty() );
    REQUIRE( error == mx::error_t::noerror );
    REQUIRE( special.valueStrGood() );

    special.forceTypeWithValue( -12345 );
    REQUIRE( special.String( &error ).empty() );
    REQUIRE( error == mx::error_t::invalidarg );

    fitsHeaderCard<> unset( "UNSET" );
    REQUIRE( unset.String( &error ).empty() );
    REQUIRE( error == mx::error_t::paramnotset );

    fitsHeaderCard<> continuation( "CONTINUE", "", "'continued'" );
    fitsHeaderCard<> numeric( "NUMBER", 1, "" );
    REQUIRE( numeric.appendContinue( continuation ) == mx::error_t::invalidarg );
    testHeaderCard<> missingString;
    missingString.forceTypeWithValue( fitsType<std::string>() );
    REQUIRE( missingString.appendContinue( continuation ) == mx::error_t::invalidconfig );
}

/** \brief Verifies writing every supported fitsHeaderCard representation through CFITSIO.
 *
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Writing FITS header card types", "[ioutils::fits::fitsHeaderCard]" )
{
    temporaryFitsFile file;

    fitsHeaderCard<> boolean( "BOOLEAN", true, "bool" );
    fitsHeaderCard<> character( "CHAR", static_cast<char>( 1 ), "char" );
    fitsHeaderCard<> unsignedCharacter( "UCHAR", static_cast<unsigned char>( 2 ), "uchar" );
    fitsHeaderCard<> shortInteger( "SHORT", static_cast<short>( 3 ), "short" );
    fitsHeaderCard<> unsignedShort( "USHORT", static_cast<unsigned short>( 4 ), "ushort" );
    fitsHeaderCard<> integer( "INT", 5, "int" );
    fitsHeaderCard<> unsignedInteger( "UINT", 6U, "uint" );
    fitsHeaderCard<> longInteger( "LONG", 7L, "long" );
    fitsHeaderCard<> unsignedLong( "ULONG", 8UL, "ulong" );
    fitsHeaderCard<> longLongInteger( "LONGLONG", 9LL, "long long" );
    fitsHeaderCard<> unsignedLongLong( "ULONGLONG", 10ULL, "unsigned long long" );
    fitsHeaderCard<> single( "FLOAT", 11.5F, "float" );
    fitsHeaderCard<> real( "DOUBLE", 12.5, "double" );
    fitsHeaderCard<> text( "STRING", "text", "string" );
    fitsHeaderCard<> serialized( "SERIAL", "13", fitsType<int>(), "serialized" );
    fitsHeaderCard<> comment( "COMMENT", fitsCommentType(), "comment text" );
    fitsHeaderCard<> history( "HISTORY", fitsHistoryType(), "history text" );

    REQUIRE( boolean.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( character.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( unsignedCharacter.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( shortInteger.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( unsignedShort.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( integer.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( unsignedInteger.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( longInteger.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( unsignedLong.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( longLongInteger.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( unsignedLongLong.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( single.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( real.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( text.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( serialized.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( comment.write( file.get() ) == mx::error_t::noerror );
    REQUIRE( history.write( file.get() ) == mx::error_t::noerror );

    fitsHeaderCard<> invalid( "INVALID", -12345 );
    REQUIRE( invalid.write( file.get() ) == mx::error_t::invalidarg );
    fitsHeaderCard<> complexSingle( "COMPLEX", std::complex<float>( 1.0F, 2.0F ), "complex" );
    REQUIRE( complexSingle.write( file.get() ) == mx::error_t::invalidarg );
}

/// Verify conversion of types
/**
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "fitsHeaderCard setting types", "[ioutils::fits::fitsHeaderCard]" )
{
    SECTION( "a fitsHeaderCard constructed with char type" )
    {
        SECTION( "setting a char from a string" )
        {
            fitsHeaderCard fhc( "KEYWORD", "39", fitsType<char>(), "this comment" );

            REQUIRE( fhc.keyword() == "KEYWORD" );
            std::string s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.type() == fitsType<char>() );
            REQUIRE( fhc.valueGood() == false );
            REQUIRE( fhc.valueStrGood() == true );
            REQUIRE( fhc.comment() == "this comment" );

            char c = fhc.value<char>();
            REQUIRE( c == 39 );
            REQUIRE( fhc.valueGood() == true );

            int i = fhc.value<int>();
            REQUIRE( i == 39 );
            REQUIRE( fhc.type() == fitsType<char>() );

            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == true );

            fhc.type( fitsType<int>() );
            REQUIRE( fhc.type() == fitsType<int>() );
            REQUIRE( fhc.Int() == 39 );
            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == false );

            s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.valueStrGood() == true );
        }
        SECTION( "setting a char from a char" )
        {
            fitsHeaderCard fhc( "KEYWORD", static_cast<char>( 39 ), "this comment" );

            REQUIRE( fhc.keyword() == "KEYWORD" );
            REQUIRE( fhc.type() == fitsType<char>() );
            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == false );
            REQUIRE( fhc.comment() == "this comment" );

            char c = fhc.value<char>();
            REQUIRE( c == 39 );
            REQUIRE( fhc.valueGood() == true );

            // Test that valueStr stayed false and then read it
            REQUIRE( fhc.valueStrGood() == false );
            std::string s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.valueStrGood() == true );

            int i = fhc.value<int>();
            REQUIRE( i == 39 );
            REQUIRE( fhc.type() == fitsType<char>() );

            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == true );

            mx::error_t errc = fhc.type( fitsType<int>() );
            REQUIRE( !errc );

            REQUIRE( fhc.type() == fitsType<int>() );
            errc = mx::error_t::error;
            REQUIRE( fhc.Int( &errc ) == 39 );
            REQUIRE( errc == mx::error_t::noerror );

            REQUIRE( fhc.valueGood() == true );
            REQUIRE( fhc.valueStrGood() == false );

            s = fhc.valueStr();
            REQUIRE( s == "39" );
            REQUIRE( fhc.valueStrGood() == true );
        }
    }
}

/// Removing white space around string values
/**
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Removing white space around string values", "[ioutils::fits::fitsHeaderCard]" )
{
    // clang-format off
    #ifdef MXLIBTEST_DOXYGEN_REF
        fitsHeaderCard fhc;
        mx::error_t errc;
        fhc.value( mx::meta::tagT<std::string>(), errc );
    #endif
    // clang-format on

    SECTION( "typical case" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'simple         '", "comment" );
        REQUIRE( fhc.String() == "simple         " );
    }

    SECTION( "no ' and no space" )
    {
        fitsHeaderCard fhc( "KEYTEST", "simple", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "no spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'simple'", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "no '" )
    {
        fitsHeaderCard fhc( "KEYTEST", "simple         ", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "space at beginning" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'   simple         '", "comment" );
        REQUIRE( fhc.String() == "   simple         " );
    }

    SECTION( "spaces at beginning, no spaces at end" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'     simple'", "comment" );
        REQUIRE( fhc.String() == "     simple" );
    }

    SECTION( "spaces at beginning, no '" )
    {
        fitsHeaderCard fhc( "KEYTEST", "   simple         ", "comment" );
        REQUIRE( fhc.String() == "simple" );
    }

    SECTION( "empty" )
    {
        fitsHeaderCard fhc( "KEYTEST", "", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "one '" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "two ''" )
    {
        fitsHeaderCard fhc( "KEYTEST", "''", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "two '', spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'  '", "comment" );
        REQUIRE( fhc.String() == "  " );
    }

    SECTION( "spaces only" )
    {
        fitsHeaderCard fhc( "KEYTEST", "    ", "comment" );
        REQUIRE( fhc.String() == "" );
    }

    SECTION( "' part of value" )
    {
        fitsHeaderCard fhc( "KEYTEST", "z'", "comment" );
        REQUIRE( fhc.String() == "z'" );
    }

    SECTION( "' part of value at end with spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'z'   '", "comment" );
        REQUIRE( fhc.String() == "z'   " );
    }

    SECTION( "' part of value at beginning with spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "''z   '", "comment" );
        REQUIRE( fhc.String() == "'z   " );
    }

    SECTION( "spaces before and after '' with ' in value surrounded by spaces" )
    {
        fitsHeaderCard fhc( "KEYTEST", "   '  'z   '  ", "comment" );
        REQUIRE( fhc.String() == "  'z   " );
    }
}

/// CONTINUE-ing a card
/**
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "CONTINUE-ing a card", "[ioutils::fits::fitsHeaderCard]" )
{
    SECTION( "normal CONTINUE, trailing space" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'one,two,three,&'", "" );
        fitsHeaderCard fhcc( "CONTINUE", "", "'four,five,six' / the comment" );

        REQUIRE( fhc.appendContinue( fhcc ) == error_t::noerror );

        REQUIRE( fhc.keyword() == "KEYTEST" );
        REQUIRE( fhc.String() == "one,two,three,four,five,six" );
        REQUIRE( fhc.comment() == "the comment" );
    }

    SECTION( "normal CONTINUE, trailing space" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'one two three &'", "" );
        fitsHeaderCard fhcc( "CONTINUE", "", "'four five six' / the comment" );

        REQUIRE( fhc.appendContinue( fhcc ) == error_t::noerror );

        REQUIRE( fhc.keyword() == "KEYTEST" );
        REQUIRE( fhc.String() == "one two three four five six" );
        REQUIRE( fhc.comment() == "the comment" );
    }

    SECTION( "normal CONTINUE, leading space" )
    {
        fitsHeaderCard fhc( "KEYTEST", "'one two three&'", "" );
        fitsHeaderCard fhcc( "CONTINUE", "", "' four five six' / the comment" );

        REQUIRE( fhc.appendContinue( fhcc ) == error_t::noerror );

        REQUIRE( fhc.keyword() == "KEYTEST" );
        REQUIRE( fhc.String() == "one two three four five six" );
        REQUIRE( fhc.comment() == "the comment" );
    }
}

/** \brief Verifies the typed card updates used for HCI coadd provenance.
 *
 * \ingroup fitsHeaderCard_unit_tests
 */
TEST_CASE( "Updating HCI coadd header cards", "[ioutils::fits::fitsHeaderCard][hciReduce]" )
{
    fitsHeaderCard date( "MJD-OBS", "2000-01-01T00:00:00.000000000Z", "mean observation time" );

    REQUIRE( date.value( "2000-01-01T12:00:00.000000000Z" ) == error_t::noerror );
    REQUIRE( date.String() == "2000-01-01T12:00:00.000000000Z" );
    REQUIRE( date.comment() == "mean observation time" );

    fitsHeaderCard angle( "ANGLE", 10.0, "derotation angle" );
    REQUIRE( angle.value( 15.5 ) == error_t::noerror );
    REQUIRE( angle.value<double>() == Approx( 15.5 ) );
    REQUIRE( angle.comment( "mean derotation angle" ) == error_t::noerror );
    REQUIRE( angle.comment() == "mean derotation angle" );
}

} // namespace fitsHeaderCardTest
} // namespace fitsTest
} // namespace unitTest
} // namespace mx
