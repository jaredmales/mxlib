/** \file fitsFile_test.cpp
 * \brief Tests FITS file I/O.
 */
#include "../../../catch2/catch.hpp"

#include <filesystem>

using namespace Catch::Matchers;

#include "../../../../include/ioutils/fits/fitsFile.hpp"
using namespace mx::fits;

#include "../../../../include/improc/eigenImage.hpp"
#include "../../../../include/improc/eigenCube.hpp"

namespace mx
{
namespace unitTest
{
namespace fitsTest
{
namespace fitsFileTest
{

// Test harnass for fitsFile
template <typename dataT>
class fitsFile_test : public fitsFile<dataT>
{
  public:
    typedef fitsFile<dataT>::pixarrT pixarrTT;

    // Interface to protected calcPixarrs
    void calcPixarrs( pixarrTT &pixarrs )
    {
        fitsFile<dataT>::calcPixarrs( pixarrs );
    }

    // Interface to set protected axis dimensions
    void setnax( int naxis, long naxes0, long naxes1, long naxes2 )
    {
        fitsFile<dataT>::m_naxis = naxis;

        if( naxis > 0 )
        {
            fitsFile<dataT>::m_naxes[0] = naxes0;
        }

        if( naxis > 1 )
        {
            fitsFile<dataT>::m_naxes[1] = naxes1;
        }

        if( naxis > 2 )
        {
            fitsFile<dataT>::m_naxes[2] = naxes2;
        }
    }

    mx::error_t writeTestFrame( const std::string &fname, int n )
    {
        mx::improc::eigenImage<float> im( 1024, 1024 );

        for( int cc = 0; cc < im.cols(); ++cc )
        {
            for( int rr = 0; rr < im.rows(); ++rr )
            {
                im( rr, cc ) = cc * im.rows() + rr;
            }
        }

        fitsHeader fh;
        fh.append( std::format( "FLOAT{}", n ), static_cast<float>( 1.1 ), "testing 1.1" );
        fh.append( std::format( "INT{}", n ), static_cast<int>( 2 ), "testing 1.1" );

        return this->write( std::format( "{}{}.fits", fname, n ), im, fh );
    }
};

/** \cond */
class temporaryDirectory
{
  public:
    temporaryDirectory( const std::string &name ) : m_path( std::filesystem::temp_directory_path() / name )
    {
        std::filesystem::remove_all( m_path );
        std::filesystem::create_directories( m_path );
    }

    ~temporaryDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all( m_path, error );
    }

    const std::filesystem::path &path() const
    {
        return m_path;
    }

  private:
    std::filesystem::path m_path;
};
/** \endcond */

/// Calculating subimage sizes
/** Verify calculation of subimage sizes
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Calculating subimage sizes", "[ioutils::fits::fitsFile]" )
{
// document what protected members are tested here
#ifdef MXLIB_DOXYGEN_PROTECTED_REF
    fitsFile<float> ff;
    fitsFile<float>::pixarrTT pixarrs;
    ff.calcPixarrs( pixarrs );
    ff.naxis();
    ff.naxes( 0 );
    ff.setReadSize( 10, 9, 7, 10 ) ff.setCubeReadSize( 5, 3 );
#endif

    SECTION( "a 1D image to read" )
    {
        SECTION( "reading the whole image with default setup" )
        {
            fitsFile_test<float> fft;

            fft.setnax( 1, 64, 0, 0 );

            REQUIRE( fft.naxis() == 1 );
            REQUIRE( fft.naxes( 0 ) == 64 );
            REQUIRE( fft.naxes( 1 ) == -1 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );

            REQUIRE( pixarrs.fpix[0] == 1 );
            REQUIRE( pixarrs.lpix[0] == 64 );
            REQUIRE( pixarrs.inc[0] == 1 );
        }

        SECTION( "reading a subimage" )
        {
            fitsFile_test<float> fft;

            fft.setnax( 1, 64, 0, 0 );

            REQUIRE( fft.naxis() == 1 );
            REQUIRE( fft.naxes( 0 ) == 64 );
            REQUIRE( fft.naxes( 1 ) == -1 );

            fft.setReadSize( 10, 0, 7, 0 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );

            REQUIRE( pixarrs.fpix[0] == 11 );
            REQUIRE( pixarrs.lpix[0] == 17 );
            REQUIRE( pixarrs.inc[0] == 1 );
        }
    }

    SECTION( "a 2D image to read" )
    {
        SECTION( "reading the whole image with default setup" )
        {
            fitsFile_test<float> fft;

            fft.setnax( 2, 64, 64, 0 );

            REQUIRE( fft.naxis() == 2 );
            REQUIRE( fft.naxes( 0 ) == 64 );
            REQUIRE( fft.naxes( 1 ) == 64 );
            REQUIRE( fft.naxes( 2 ) == -1 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );

            REQUIRE( pixarrs.fpix[0] == 1 );
            REQUIRE( pixarrs.fpix[1] == 1 );
            REQUIRE( pixarrs.lpix[0] == 64 );
            REQUIRE( pixarrs.lpix[1] == 64 );
            REQUIRE( pixarrs.inc[0] == 1 );
            REQUIRE( pixarrs.inc[1] == 1 );
        }

        SECTION( "reading a subimage" )
        {
            fitsFile_test<float> fft;

            fft.setnax( 2, 64, 64, 0 );

            REQUIRE( fft.naxis() == 2 );
            REQUIRE( fft.naxes( 0 ) == 64 );
            REQUIRE( fft.naxes( 1 ) == 64 );
            REQUIRE( fft.naxes( 2 ) == -1 );

            fft.setReadSize( 10, 9, 7, 10 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );

            REQUIRE( pixarrs.fpix[0] == 11 );
            REQUIRE( pixarrs.fpix[1] == 10 );
            REQUIRE( pixarrs.lpix[0] == 17 );
            REQUIRE( pixarrs.lpix[1] == 19 );
            REQUIRE( pixarrs.inc[0] == 1 );
            REQUIRE( pixarrs.inc[1] == 1 );
        }
    }

    SECTION( "a 3D image to read" )
    {
        WHEN( "reading the whole image with default setup" )
        {
            fitsFile_test<float> fft;

            fft.setnax( 3, 64, 64, 64 );

            REQUIRE( fft.naxis() == 3 );
            REQUIRE( fft.naxes( 0 ) == 64 );
            REQUIRE( fft.naxes( 1 ) == 64 );
            REQUIRE( fft.naxes( 2 ) == 64 );
            REQUIRE( fft.naxes( 3 ) == -1 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );

            REQUIRE( pixarrs.fpix[0] == 1 );
            REQUIRE( pixarrs.fpix[1] == 1 );
            REQUIRE( pixarrs.fpix[2] == 1 );

            REQUIRE( pixarrs.lpix[0] == 64 );
            REQUIRE( pixarrs.lpix[1] == 64 );
            REQUIRE( pixarrs.lpix[2] == 64 );
            REQUIRE( pixarrs.inc[0] == 1 );
            REQUIRE( pixarrs.inc[1] == 1 );
            REQUIRE( pixarrs.inc[2] == 1 );
        }

        SECTION( "reading a subimage" )
        {
            fitsFile_test<float> fft;

            fft.setnax( 3, 64, 64, 64 );

            REQUIRE( fft.naxis() == 3 );
            REQUIRE( fft.naxes( 0 ) == 64 );
            REQUIRE( fft.naxes( 1 ) == 64 );
            REQUIRE( fft.naxes( 2 ) == 64 );
            REQUIRE( fft.naxes( 3 ) == -1 );

            fft.setReadSize( 10, 9, 7, 10 );
            fft.setCubeReadSize( 5, 3 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );

            REQUIRE( pixarrs.fpix[0] == 11 );
            REQUIRE( pixarrs.fpix[1] == 10 );
            REQUIRE( pixarrs.fpix[2] == 6 );
            REQUIRE( pixarrs.lpix[0] == 17 );
            REQUIRE( pixarrs.lpix[1] == 19 );
            REQUIRE( pixarrs.lpix[2] == 8 );
            REQUIRE( pixarrs.inc[0] == 1 );
            REQUIRE( pixarrs.inc[1] == 1 );
            REQUIRE( pixarrs.inc[2] == 1 );
        }
    }
}

/// Basic writing and reading
/**
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Basic writing and reading", "[ioutils::fits::fitsFile]" )
{
    // clang-format off
    // document what fitsFile_test is doing
    #ifdef MXLIB_DOXYGEN_PROTECTED_REF
        fitsFile<float> ff;
        fitsHeader fh;
        fh.append( "dummy", 2, "dummy" );
        mx::improc::eigenImage<float> im;
        ff.write( "dummy", im, fh );
    #endif
    // clang-format on

    SECTION( "write a frame and read it back in" )
    {
        fitsFile_test<float> fft;

        REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 1 ) == mx::error_t::noerror );

        fitsFile<float> ff;
        fitsHeader fh;
        mx::improc::eigenImage<float> im;

        REQUIRE( ff.read( im, fh, "/tmp/fitsFile_test1.fits" ) == mx::error_t::noerror );

        REQUIRE_THAT( fh["FLOAT1"].Float(), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( fh["INT1"].Int() == 2 );
        REQUIRE( im.rows() == 1024 );
        REQUIRE( im.cols() == 1024 );
    }
}

/// Cube writing and reading
/**
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Cube writing and reading", "[ioutils::fits::fitsFile]" )
{
    // clang-format off
    // document what fitsFile_test is doing
    #ifdef MXLIB_DOXYGEN_PROTECTED_REF
        fitsFile<float> ff;
        fitsHeader fh;
        fh.append( "dummy", 2, "dummy" );
        mx::improc::eigenImage<float> im;
        ff.write( "dummy", im, fh );
    #endif
    // clang-format on

    SECTION( "write 5 frames and read them back in as a cube" )
    {
        fitsFile_test<float> fft;

        REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 1 ) == mx::error_t::noerror );
        REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 2 ) == mx::error_t::noerror );
        REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 3 ) == mx::error_t::noerror );
        REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 4 ) == mx::error_t::noerror );
        REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 5 ) == mx::error_t::noerror );

        std::vector<std::string> flist( { "/tmp/fitsFile_test1.fits",
                                          "/tmp/fitsFile_test2.fits",
                                          "/tmp/fitsFile_test3.fits",
                                          "/tmp/fitsFile_test4.fits",
                                          "/tmp/fitsFile_test5.fits" } );

        fitsFile<float> ff;
        mx::improc::eigenCube<float> imc;

        std::vector<mx::fits::fitsHeader<mx::verbose::d>> headers( flist.size() );
        for( size_t index = 0; index < headers.size(); ++index )
        {
            headers[index].append( std::format( "FLOAT{}", index + 1 ) );
        }

        REQUIRE( ff.read( imc, headers, flist ) == mx::error_t::noerror );

        REQUIRE( imc.rows() == 1024 );
        REQUIRE( imc.cols() == 1024 );
        REQUIRE( imc.planes() == 5 );
        for( size_t index = 0; index < headers.size(); ++index )
        {
            REQUIRE_THAT( headers[index][std::format( "FLOAT{}", index + 1 )].Float(), WithinAbs( 1.1, 1e-5 ) );
        }

        std::vector<float> rawImages( 1024 * 1024 * flist.size() );
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> rawHeaders( flist.size() );
        for( size_t index = 0; index < rawHeaders.size(); ++index )
        {
            rawHeaders[index].append( std::format( "INT{}", index + 1 ) );
        }

        REQUIRE( ff.read( rawImages.data(), rawHeaders, flist ) == mx::error_t::noerror );
        REQUIRE( rawImages.front() == 0.0F );
        REQUIRE( rawImages[1024 * 1024 - 1] == 1048575.0F );
        REQUIRE( rawImages[1024 * 1024] == 0.0F );
        for( size_t index = 0; index < rawHeaders.size(); ++index )
        {
            REQUIRE( rawHeaders[index][std::format( "INT{}", index + 1 )].Int() == 2 );
        }
    }
}

/// Validate HCI bulk image/header read arguments and later-file failures.
/**
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Validating HCI batch image and header reads", "[ioutils::fits::fitsFile][hciReduce]" )
{
    fitsFile_test<float> fixture;
    REQUIRE( fixture.writeTestFrame( "/tmp/fitsFile-batch-validation", 1 ) == mx::error_t::noerror );

    const std::vector<std::string> oneFile{ "/tmp/fitsFile-batch-validation1.fits" };
    fitsFile<float> fitsFile;

    SECTION( "empty file lists are rejected" )
    {
        mx::improc::eigenCube<float> cube;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> headers;
        const std::vector<std::string> noFiles;

        REQUIRE( fitsFile.read( cube, headers, noFiles ) == mx::error_t::invalidarg );
    }

    SECTION( "header vectors must contain one header per input image" )
    {
        std::vector<float> data( 1024 * 1024 );
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> headers;
        REQUIRE( fitsFile.read( data.data(), headers, oneFile ) == mx::error_t::invalidarg );

        mx::improc::eigenCube<float> cube;
        REQUIRE( fitsFile.read( cube, headers, oneFile ) == mx::error_t::invalidarg );
    }

    SECTION( "a missing later input file is returned as an error" )
    {
        mx::improc::eigenCube<float> cube;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> headers( 2 );
        const std::vector<std::string> files{ oneFile.front(), "/tmp/fitsFile-batch-validation-missing.fits" };

        REQUIRE( fitsFile.read( cube, headers, files ) != mx::error_t::noerror );
    }
}

/// Verify hciReduce's float Eigen image and cube writes round-trip with optional FITS headers.
/** Exercises mx::fits::fitsFile::write Eigen-array overloads with mx::verbose::vv. */
/**
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Verbose float Eigen FITS writes round-trip", "[ioutils::fits::fitsFile][hciReduce]" )
{
    temporaryDirectory testDirectory( "mxlib-fitsFile-write-test" );
    fitsFile<float, verbose::vv> fitsFile;

    mx::improc::eigenImage<float> image( 2, 3 );
    image << 1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F;

    const std::filesystem::path imagePath = testDirectory.path() / "image.fits";
    REQUIRE( fitsFile.write( imagePath.string(), image ) == mx::error_t::noerror );

    mx::improc::eigenImage<float> readImage;
    REQUIRE( fitsFile.read( readImage, imagePath.string() ) == mx::error_t::noerror );
    REQUIRE( readImage.isApprox( image ) );

    REQUIRE( fitsFile.open( imagePath.string() ) == mx::error_t::noerror );
    REQUIRE( fitsFile.write( imagePath.string(), image ) == mx::error_t::noerror );

    fitsHeader<verbose::vv> header;
    header.append( "TESTKEY", 17, "test header value" );

    const std::filesystem::path headerImagePath = testDirectory.path() / "image-with-header.fits";
    REQUIRE( fitsFile.write( headerImagePath.string(), image, header ) == mx::error_t::noerror );
    const std::filesystem::path unavailableImagePath = testDirectory.path() / "missing" / "image.fits";
    REQUIRE( fitsFile.write( unavailableImagePath.string(), image ) != mx::error_t::noerror );
    REQUIRE( fitsFile.write( unavailableImagePath.string(), image, header ) != mx::error_t::noerror );

    fitsHeader<verbose::vv> readHeader;
    REQUIRE( fitsFile.read( readImage, readHeader, headerImagePath.string() ) == mx::error_t::noerror );
    REQUIRE( readImage.isApprox( image ) );
    REQUIRE( readHeader["TESTKEY"].value<int>() == 17 );

    mx::improc::eigenCube<float> cube( 2, 3, 2 );
    cube.image( 0 ) = image;
    cube.image( 1 ) = image + 10.0F;

    const std::filesystem::path cubePath = testDirectory.path() / "cube.fits";
    REQUIRE( fitsFile.write( cubePath.string(), cube ) == mx::error_t::noerror );

    mx::improc::eigenCube<float> readCube;
    REQUIRE( fitsFile.read( readCube, cubePath.string() ) == mx::error_t::noerror );
    REQUIRE( readCube.rows() == cube.rows() );
    REQUIRE( readCube.cols() == cube.cols() );
    REQUIRE( readCube.planes() == cube.planes() );
    REQUIRE( readCube.image( 0 ).isApprox( cube.image( 0 ) ) );
    REQUIRE( readCube.image( 1 ).isApprox( cube.image( 1 ) ) );

    const std::filesystem::path headerCubePath = testDirectory.path() / "cube-with-header.fits";
    REQUIRE( fitsFile.write( headerCubePath.string(), cube, header ) == mx::error_t::noerror );
    const std::filesystem::path unavailableCubePath = testDirectory.path() / "missing" / "cube.fits";
    REQUIRE( fitsFile.write( unavailableCubePath.string(), cube, header ) != mx::error_t::noerror );

    REQUIRE( fitsFile.read( readCube, readHeader, headerCubePath.string() ) == mx::error_t::noerror );
    REQUIRE( readCube.image( 0 ).isApprox( cube.image( 0 ) ) );
    REQUIRE( readCube.image( 1 ).isApprox( cube.image( 1 ) ) );
    REQUIRE( readHeader["TESTKEY"].value<int>() == 17 );
}

/// Reading headers
/**
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Reading headers", "[ioutils::fits::fitsFile]" )
{
// document what fitsFile_test is doing
#ifdef MXLIB_DOXYGEN_PROTECTED_REF
    fitsFile<float> ff;
    fitsHeader fh;
    fh.append( "dummy", 2, "dummy" );
    mx::improc::eigenImage<float> im;
    ff.write( "dummy", im, fh );
#endif

    fitsFile_test<float> fft;

    std::vector<std::string> fnames(
        { "/tmp/fitsFile_test1.fits", "/tmp/fitsFile_test2.fits", "/tmp/fitsFile_test3.fits" } );
    REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 1 ) == mx::error_t::noerror );
    REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 2 ) == mx::error_t::noerror );
    REQUIRE( fft.writeTestFrame( "/tmp/fitsFile_test", 3 ) == mx::error_t::noerror );

    SECTION( "write three frames and read their headers as a vector, not allocated" )
    {

        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        REQUIRE( ff.readHeader( fh, fnames ) == mx::error_t::noerror );
        REQUIRE( fh.size() == fnames.size() );

        REQUIRE_THAT( fh[0]["FLOAT1"].Float(), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( fh[0]["INT1"].Int() == 2 );
        REQUIRE_THAT( fh[1]["FLOAT2"].Float(), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( fh[1]["INT2"].Int() == 2 );
        REQUIRE_THAT( fh[2]["FLOAT3"].Float(), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( fh[2]["INT3"].Int() == 2 );
    }

    SECTION( "write three frames and read their headers as a vector, allocated" )
    {

        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        fh.resize( fnames.size() );
        REQUIRE( ff.readHeader( fh, fnames ) == mx::error_t::noerror );

        mx::error_t errc;
        errc = mx::error_t::error;
        REQUIRE_THAT( fh[0]["FLOAT1"].Float( &errc ), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( errc == mx::error_t::noerror );
        errc = mx::error_t::error;
        REQUIRE( fh[0]["INT1"].Int( &errc ) == 2 );
        REQUIRE( errc == mx::error_t::noerror );
        errc = mx::error_t::error;
        REQUIRE_THAT( fh[1]["FLOAT2"].Float( &errc ), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( errc == mx::error_t::noerror );
        errc = mx::error_t::error;
        REQUIRE( fh[1]["INT2"].Int( &errc ) == 2 );
        REQUIRE( errc == mx::error_t::noerror );
        errc = mx::error_t::error;
        REQUIRE_THAT( fh[2]["FLOAT3"].Float( &errc ), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( errc == mx::error_t::noerror );
        errc = mx::error_t::error;
        REQUIRE( fh[2]["INT3"].Int( &errc ) == 2 );
        REQUIRE( errc == mx::error_t::noerror );
    }

    SECTION( "write three frames and read their headers as a vector, allocated to wrong size" )
    {
        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        fh.resize( fnames.size() + 1 );
        REQUIRE( ff.readHeader( fh, fnames ) == mx::error_t::invalidarg );
    }

    SECTION( "write three frames and read their headers as a vector, extracting only 1 keyword" )
    {

        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        fh.resize( fnames.size() );
        fh[0].append( "FLOAT1" );
        fh[1].append( "FLOAT2" );
        fh[2].append( "FLOAT3" );

        REQUIRE( ff.readHeader( fh, fnames ) == mx::error_t::noerror );

        mx::error_t errc;

        REQUIRE_THAT( fh[0]["FLOAT1"].Float( &errc ), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( errc == mx::error_t::noerror );
        REQUIRE( fh[0]["FLOAT1"].valueGood() == true );

        REQUIRE( fh[0]["INT1"].Int( &errc ) == 0 );
        REQUIRE( errc == mx::error_t::invalidarg );
        REQUIRE( fh[0]["INT1"].valueGood() == false );

        REQUIRE_THAT( fh[1]["FLOAT2"].Float( &errc ), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( errc == mx::error_t::noerror );
        REQUIRE( fh[1]["FLOAT2"].valueGood() == true );

        REQUIRE( fh[1]["INT2"].Int( &errc ) == 0 );
        REQUIRE( errc == mx::error_t::invalidarg );
        REQUIRE( fh[1]["INT2"].valueGood() == false );

        REQUIRE_THAT( fh[2]["FLOAT3"].Float( &errc ), WithinAbs( 1.1, 1e-5 ) );
        REQUIRE( errc == mx::error_t::noerror );
        REQUIRE( fh[2]["FLOAT3"].valueGood() == true );

        REQUIRE( fh[2]["INT3"].Int( &errc ) == 0 );
        REQUIRE( errc == mx::error_t::invalidarg );
        REQUIRE( fh[2]["INT3"].valueGood() == false );
    }
}

} // namespace fitsFileTest
} // namespace fitsTest
} // namespace unitTest
} // namespace mx
