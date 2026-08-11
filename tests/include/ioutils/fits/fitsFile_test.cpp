/** \file fitsFile_test.cpp
 * \brief Tests FITS file I/O.
 */
#include "../../../catch2/catch.hpp"

#include <algorithm>
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
    mx::error_t calcPixarrs( pixarrTT &pixarrs )
    {
        return fitsFile<dataT>::calcPixarrs( pixarrs );
    }

    mx::error_t allocatePixarrs( pixarrTT &pixarrs, int naxis )
    {
        return pixarrs.allocate( naxis );
    }

    void noComment( bool noComment )
    {
        fitsFile<dataT>::m_noComment = noComment;
    }

    bool isOpen() const
    {
        return fitsFile<dataT>::m_isOpen;
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

template <int exceptionKind>
void throwResizeException()
{
    if constexpr( exceptionKind == 0 )
    {
        throw std::bad_alloc();
    }
    else if constexpr( exceptionKind == 1 )
    {
        throw std::runtime_error( "resize failure" );
    }
    else
    {
        throw 1;
    }
}

template <int exceptionKind>
struct throwingCube
{
    using Scalar = float;

    void resize( int, int, int )
    {
        throwResizeException<exceptionKind>();
    }
};

template <int exceptionKind>
struct throwingImage
{
    using Scalar = float;

    void resize( int, int )
    {
        throwResizeException<exceptionKind>();
    }
};

class fitsFileOpsGuard
{
  public:
    fitsFileOpsGuard() : m_saved( fitsFileDetail::fitsFileCfitsioOpsInstance() )
    {
    }

    ~fitsFileOpsGuard()
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance() = m_saved;
    }

  private:
    fitsFileDetail::fitsFileCfitsioOps m_saved;
};

int failOpen( fitsfile **, const char *, int, int *status )
{
    *status = FILE_NOT_OPENED;
    return *status;
}

int failGetImageDim( fitsfile *, int *, int *status )
{
    *status = BAD_NAXIS;
    return *status;
}

int failGetImageSize( fitsfile *, int, long *, int *status )
{
    *status = BAD_DIMEN;
    return *status;
}

int failClose( fitsfile *, int *status )
{
    *status = FILE_NOT_CLOSED;
    return *status;
}

int failGetHeaderSpace( fitsfile *, int *, int *, int *status )
{
    *status = READ_ERROR;
    return *status;
}

int getHeaderSpaceCall = 0;
int failGetHeaderSpaceOnCall = 1;

int controlledGetHeaderSpace( fitsfile *file, int *keysExist, int *moreKeys, int *status )
{
    ++getHeaderSpaceCall;
    if( getHeaderSpaceCall == failGetHeaderSpaceOnCall )
    {
        *status = READ_ERROR;
        return *status;
    }

    return ffghsp( file, keysExist, moreKeys, status );
}

int failReadKey( fitsfile *, int, char *, char *, char *, int *status )
{
    *status = READ_ERROR;
    return *status;
}

int failCreateFile( fitsfile **, const char *, int *status )
{
    *status = FILE_NOT_CREATED;
    return *status;
}

int failCreateImage( fitsfile *file, int bitpix, int naxis, long *naxes, int *status )
{
    ffcrim( file, bitpix, naxis, naxes, status );
    if( *status != 0 )
    {
        return *status;
    }
    *status = BAD_NAXIS;
    return *status;
}

int failWritePixels( fitsfile *, int, long *, LONGLONG, void *, int *status )
{
    *status = WRITE_ERROR;
    return *status;
}

int readSubsetCall = 0;
int failReadSubsetOnCall = 1;
int readSubsetStatus = READ_ERROR;

int controlledReadSubset( fitsfile *file,
                          int dataType,
                          long *firstPixel,
                          long *lastPixel,
                          long *increment,
                          void *nullValue,
                          void *data,
                          int *anyNull,
                          int *status )
{
    ++readSubsetCall;
    if( readSubsetCall == failReadSubsetOnCall )
    {
        *status = readSubsetStatus;
        return *status;
    }

    return ffgsv( file, dataType, firstPixel, lastPixel, increment, nullValue, data, anyNull, status );
}

enum class allocationFailure
{
    none,
    badAllocation,
    standardException,
    unknownException,
    nullPointer
};

allocationFailure allocationFailureKind = allocationFailure::none;
int allocationCall = 0;
int failAllocationOnCall = 1;

long *controlledAllocateLongs( size_t count )
{
    ++allocationCall;
    if( allocationCall == failAllocationOnCall )
    {
        if( allocationFailureKind == allocationFailure::badAllocation )
        {
            throw std::bad_alloc();
        }
        if( allocationFailureKind == allocationFailure::standardException )
        {
            throw std::runtime_error( "coordinate allocation failure" );
        }
        if( allocationFailureKind == allocationFailure::unknownException )
        {
            throw 1;
        }
        if( allocationFailureKind == allocationFailure::nullPointer )
        {
            return nullptr;
        }
    }

    return new long[count];
}
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

        SECTION( "reading complete images from selected cube planes" )
        {
            fitsFile_test<float> fft;
            fft.setnax( 3, 64, 32, 8 );
            fft.setCubeReadSize( 2, 3 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );
            REQUIRE( pixarrs.fpix[0] == 1 );
            REQUIRE( pixarrs.lpix[0] == 64 );
            REQUIRE( pixarrs.fpix[1] == 1 );
            REQUIRE( pixarrs.lpix[1] == 32 );
            REQUIRE( pixarrs.fpix[2] == 3 );
            REQUIRE( pixarrs.lpix[2] == 5 );
        }

        SECTION( "reading subimages from every cube plane" )
        {
            fitsFile_test<float> fft;
            fft.setnax( 3, 64, 32, 8 );
            fft.setReadSize( 10, 5, 7, 4 );

            fitsFile_test<float>::pixarrTT pixarrs;
            fft.calcPixarrs( pixarrs );
            REQUIRE( pixarrs.fpix[0] == 11 );
            REQUIRE( pixarrs.lpix[0] == 17 );
            REQUIRE( pixarrs.fpix[1] == 6 );
            REQUIRE( pixarrs.lpix[1] == 9 );
            REQUIRE( pixarrs.fpix[2] == 1 );
            REQUIRE( pixarrs.lpix[2] == 8 );
        }
    }
}

/** \brief Verifies fitsFile construction, state transitions, dimensions, and read-window resets.
 *
 * Exercises mx::fits::fitsFile constructors, fileName, open, close, getDimensions, getSize, setReadSize,
 * setCubeReadSize, and repeated protected pixel-array allocation through the test fixture.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Managing FITS file state and dimensions", "[ioutils::fits::fitsFile]" )
{
    temporaryDirectory testDirectory( "mxlib-fitsFile-state-test" );
    const std::filesystem::path imagePath = testDirectory.path() / "image.fits";

    mx::improc::eigenImage<float> image( 2, 3 );
    image << 1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F;
    fitsFile<float> writer;
    REQUIRE( writer.write( imagePath.string(), image ) == mx::error_t::noerror );

    mx::error_t error = mx::error_t::error;
    fitsFile<float> defaultWithError( error );
    REQUIRE( error == mx::error_t::noerror );
    REQUIRE( defaultWithError.fileName().empty() );
    REQUIRE( defaultWithError.open() == mx::error_t::invalidconfig );
    REQUIRE( defaultWithError.close() == mx::error_t::noerror );
    REQUIRE( defaultWithError.getDimensions( error ) == -1 );
    REQUIRE( error == mx::error_t::invalidconfig );
    REQUIRE( defaultWithError.getSize( error ) == -1 );
    REQUIRE( error == mx::error_t::invalidconfig );
    REQUIRE( defaultWithError.getSize( 0, error ) == -1 );
    REQUIRE( error == mx::error_t::invalidconfig );

    fitsFile<float> named( imagePath.string(), error );
    REQUIRE( error == mx::error_t::noerror );
    REQUIRE( named.fileName() == imagePath.string() );

    fitsFile<float> namedClosed( imagePath.string(), false );
    REQUIRE( namedClosed.fileName() == imagePath.string() );
    REQUIRE( namedClosed.open() == mx::error_t::noerror );
    REQUIRE( namedClosed.open() == mx::error_t::noerror );
    REQUIRE( namedClosed.getDimensions( error ) == 2 );
    REQUIRE( namedClosed.getSize( error ) == 6 );
    REQUIRE( namedClosed.getSize( 0, error ) == 2 );
    REQUIRE( namedClosed.getSize( 1, error ) == 3 );
    REQUIRE( namedClosed.getSize( 2, error ) == -1 );
    REQUIRE( error == mx::error_t::invalidarg );

    namedClosed.setReadSize( 0, 0, 1, 2 );
    REQUIRE( namedClosed.getSize( error ) == 2 );
    REQUIRE( namedClosed.getSize( 0, error ) == 1 );
    REQUIRE( namedClosed.getSize( 1, error ) == 2 );
    namedClosed.setReadSize();
    REQUIRE( namedClosed.getSize( error ) == 6 );
    namedClosed.setCubeReadSize( 0, 1 );
    namedClosed.setCubeReadSize();

    REQUIRE( namedClosed.fileName( imagePath.string(), false ) == mx::error_t::noerror );
    REQUIRE( namedClosed.close() == mx::error_t::noerror );

    fitsFile<float> opened( imagePath.string(), true, error );
    REQUIRE( error == mx::error_t::noerror );
    REQUIRE( opened.getDimensions( error ) == 2 );

    fitsFile<float> missing( ( testDirectory.path() / "missing.fits" ).string(), true );
    REQUIRE( missing.naxis() == 0 );

    fitsFile_test<float> fixture;
    fitsFile_test<float>::pixarrTT pixelArrays;
    REQUIRE( fixture.allocatePixarrs( pixelArrays, 0 ) == mx::error_t::paramnotset );
    REQUIRE( fixture.allocatePixarrs( pixelArrays, 2 ) == mx::error_t::noerror );
    REQUIRE( fixture.allocatePixarrs( pixelArrays, 3 ) == mx::error_t::noerror );
}

/** \brief Verifies Eigen image and cube resize adapters preserve each exception category.
 *
 * Exercises mx::fits::eigenArrResize for bad allocation, standard, and non-standard exceptions from both supported
 * resize signatures.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Propagating FITS destination resize exceptions", "[ioutils::fits::fitsFile]" )
{
    throwingCube<0> badAllocationCube;
    throwingCube<1> standardExceptionCube;
    throwingCube<2> unknownExceptionCube;
    eigenArrResize<throwingCube<0>, verbose::vv, true> badAllocationCubeResize;
    eigenArrResize<throwingCube<1>, verbose::vv, true> standardExceptionCubeResize;
    eigenArrResize<throwingCube<2>, verbose::vv, true> unknownExceptionCubeResize;

    REQUIRE_THROWS_AS( badAllocationCubeResize.resize( badAllocationCube, 1, 1, 1 ), std::bad_alloc );
    REQUIRE_THROWS_AS( standardExceptionCubeResize.resize( standardExceptionCube, 1, 1, 1 ), std::runtime_error );
    REQUIRE_THROWS_AS( unknownExceptionCubeResize.resize( unknownExceptionCube, 1, 1, 1 ), int );

    throwingImage<0> badAllocationImage;
    throwingImage<1> standardExceptionImage;
    throwingImage<2> unknownExceptionImage;
    eigenArrResize<throwingImage<0>, verbose::vv> badAllocationImageResize;
    eigenArrResize<throwingImage<1>, verbose::vv> standardExceptionImageResize;
    eigenArrResize<throwingImage<2>, verbose::vv> unknownExceptionImageResize;

    REQUIRE_THROWS_AS( badAllocationImageResize.resize( badAllocationImage, 1, 1, 1 ), std::bad_alloc );
    REQUIRE_THROWS_AS( standardExceptionImageResize.resize( standardExceptionImage, 1, 1, 1 ), std::runtime_error );
    REQUIRE_THROWS_AS( unknownExceptionImageResize.resize( unknownExceptionImage, 1, 1, 1 ), int );
}

/** \brief Verifies fitsFile reports CFITSIO open, dimension, size, and close failures and remains recoverable.
 *
 * Exercises mx::fits::fitsFile::open and mx::fits::fitsFile::close with deterministic CFITSIO failures.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Propagating FITS file lifecycle failures", "[ioutils::fits::fitsFile]" )
{
    temporaryDirectory testDirectory( "mxlib-fitsFile-lifecycle-failure-test" );
    const std::filesystem::path imagePath = testDirectory.path() / "image.fits";
    mx::improc::eigenImage<float> image( 2, 2 );
    image.setOnes();
    fitsFile<float> writer;
    REQUIRE( writer.write( imagePath.string(), image ) == mx::error_t::noerror );

    SECTION( "open failure" )
    {
        fitsFileOpsGuard operationsGuard;
        fitsFileDetail::fitsFileCfitsioOpsInstance().openFile = failOpen;
        fitsFile_test<float> reader;
        REQUIRE( reader.open( imagePath.string() ) == mx::error_t::fits_file_not_opened );
        REQUIRE_FALSE( reader.isOpen() );
    }

    SECTION( "image-axis count failure" )
    {
        fitsFileOpsGuard operationsGuard;
        fitsFileDetail::fitsFileCfitsioOpsInstance().getImageDim = failGetImageDim;
        fitsFile_test<float> reader;
        REQUIRE( reader.open( imagePath.string() ) == mx::error_t::fits_bad_naxis );
        REQUIRE_FALSE( reader.isOpen() );
    }

    SECTION( "image-axis size failure" )
    {
        fitsFileOpsGuard operationsGuard;
        fitsFileDetail::fitsFileCfitsioOpsInstance().getImageSize = failGetImageSize;
        fitsFile_test<float> reader;
        REQUIRE( reader.open( imagePath.string() ) == mx::error_t::fits_bad_dimen );
        REQUIRE_FALSE( reader.isOpen() );
    }

    SECTION( "close failure" )
    {
        fitsFile_test<float> reader;
        REQUIRE( reader.open( imagePath.string() ) == mx::error_t::noerror );
        {
            fitsFileOpsGuard operationsGuard;
            fitsFileDetail::fitsFileCfitsioOpsInstance().closeFile = failClose;
            REQUIRE( reader.close() == mx::error_t::fits_file_not_closed );
            REQUIRE( reader.isOpen() );
        }
        REQUIRE( reader.close() == mx::error_t::noerror );
        REQUIRE_FALSE( reader.isOpen() );
    }

    fitsFileDetail::resetFitsFileCfitsioOps();
}

/** \brief Verifies fitsFile preserves allocation failures and rejects null pixel-coordinate arrays.
 *
 * Exercises the protected pixel-coordinate allocation used by mx::fits::fitsFile::read.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Propagating FITS pixel-coordinate allocation failures", "[ioutils::fits::fitsFile]" )
{
    fitsFileOpsGuard operationsGuard;
    fitsFileDetail::fitsFileCfitsioOpsInstance().allocateLongs = controlledAllocateLongs;
    allocationCall = 0;
    failAllocationOnCall = 1;

    fitsFile_test<float> fixture;
    fitsFile_test<float>::pixarrTT pixelArrays;

    SECTION( "bad allocation" )
    {
        allocationFailureKind = allocationFailure::badAllocation;
        REQUIRE_THROWS_AS( fixture.allocatePixarrs( pixelArrays, 2 ), std::bad_alloc );
    }

    SECTION( "standard exception" )
    {
        allocationFailureKind = allocationFailure::standardException;
        REQUIRE_THROWS_AS( fixture.allocatePixarrs( pixelArrays, 2 ), std::runtime_error );
    }

    SECTION( "unknown exception" )
    {
        allocationFailureKind = allocationFailure::unknownException;
        REQUIRE_THROWS_AS( fixture.allocatePixarrs( pixelArrays, 2 ), int );
    }

    SECTION( "null first-pixel array" )
    {
        allocationFailureKind = allocationFailure::nullPointer;
        failAllocationOnCall = 1;
        REQUIRE( fixture.allocatePixarrs( pixelArrays, 2 ) == mx::error_t::allocerr );
    }

    SECTION( "null last-pixel array" )
    {
        allocationFailureKind = allocationFailure::nullPointer;
        failAllocationOnCall = 2;
        REQUIRE( fixture.allocatePixarrs( pixelArrays, 2 ) == mx::error_t::allocerr );
    }

    SECTION( "null increment array" )
    {
        allocationFailureKind = allocationFailure::nullPointer;
        failAllocationOnCall = 3;
        REQUIRE( fixture.allocatePixarrs( pixelArrays, 2 ) == mx::error_t::allocerr );
    }

    SECTION( "allocation error propagates through read-coordinate calculation" )
    {
        allocationFailureKind = allocationFailure::nullPointer;
        failAllocationOnCall = 1;
        fixture.setnax( 2, 2, 2, 0 );
        REQUIRE( fixture.calcPixarrs( pixelArrays ) == mx::error_t::allocerr );
    }
}

/** \brief Verifies raw, Eigen, cube, and header reads propagate deterministic CFITSIO failures.
 *
 * Exercises mx::fits::fitsFile::read and mx::fits::fitsFile::readHeader error paths.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Propagating FITS read failures", "[ioutils::fits::fitsFile]" )
{
    temporaryDirectory testDirectory( "mxlib-fitsFile-read-failure-test" );
    const std::filesystem::path imagePath = testDirectory.path() / "image.fits";
    mx::improc::eigenImage<float> image( 2, 2 );
    image << 1.0F, 2.0F, 3.0F, 4.0F;
    fitsHeader<> writeHeader;
    REQUIRE( writeHeader.append( "TESTKEY", 42, "test value" ) == mx::error_t::noerror );
    fitsFile<float> writer;
    REQUIRE( writer.write( imagePath.string(), image, writeHeader ) == mx::error_t::noerror );

    fitsFileOpsGuard operationsGuard;

    SECTION( "raw and Eigen subset failures" )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset = controlledReadSubset;
        readSubsetCall = 0;
        failReadSubsetOnCall = 1;
        readSubsetStatus = READ_ERROR;

        std::vector<float> raw( 4 );
        fitsFile<float> rawReader( imagePath.string(), false );
        REQUIRE( rawReader.read( raw.data() ) == mx::error_t::fits_read_error );

        readSubsetCall = 0;
        mx::improc::eigenImage<float> readImage;
        fitsFile<float> eigenReader;
        REQUIRE( eigenReader.read( readImage, imagePath.string() ) == mx::error_t::fits_read_error );

        readSubsetCall = 0;
        fitsHeader<> readHeader;
        REQUIRE( eigenReader.read( readImage, readHeader, imagePath.string() ) == mx::error_t::fits_read_error );
    }

    SECTION( "end-of-file status is accepted" )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset = controlledReadSubset;
        readSubsetCall = 0;
        failReadSubsetOnCall = 1;
        readSubsetStatus = END_OF_FILE;
        std::vector<float> raw( 4 );
        fitsFile<float> reader;
        REQUIRE( reader.read( raw.data(), imagePath.string() ) == mx::error_t::noerror );
    }

    SECTION( "first and later cube plane failures" )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset = controlledReadSubset;
        const std::vector<std::string> files{ imagePath.string(), imagePath.string() };
        mx::improc::eigenCube<float> cube;
        fitsFile<float> reader;

        readSubsetCall = 0;
        failReadSubsetOnCall = 1;
        readSubsetStatus = READ_ERROR;
        REQUIRE( reader.read( cube, files ) == mx::error_t::fits_read_error );

        readSubsetCall = 0;
        failReadSubsetOnCall = 2;
        REQUIRE( reader.read( cube, files ) == mx::error_t::fits_read_error );
    }

    SECTION( "cube coordinate allocation failure" )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().allocateLongs = controlledAllocateLongs;
        allocationCall = 0;
        failAllocationOnCall = 1;
        allocationFailureKind = allocationFailure::nullPointer;
        mx::improc::eigenCube<float> cube;
        fitsFile<float> reader;
        REQUIRE( reader.read( cube, std::vector<std::string>{ imagePath.string() } ) == mx::error_t::allocerr );
    }

    SECTION( "first and later cube header failures" )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().getHeaderSpace = controlledGetHeaderSpace;
        const std::vector<std::string> files{ imagePath.string(), imagePath.string() };
        std::vector<fitsHeader<>> headers( files.size() );
        mx::improc::eigenCube<float> cube;
        fitsFile<float> reader;

        getHeaderSpaceCall = 0;
        failGetHeaderSpaceOnCall = 1;
        REQUIRE( reader.read( cube, headers, files ) == mx::error_t::fits_read_error );

        getHeaderSpaceCall = 0;
        failGetHeaderSpaceOnCall = 2;
        REQUIRE( reader.read( cube, headers, files ) == mx::error_t::fits_read_error );
    }

    SECTION( "header-space and header-card failures" )
    {
        fitsHeader<> readHeader;
        fitsFile<float> reader;
        fitsFileDetail::fitsFileCfitsioOpsInstance().getHeaderSpace = failGetHeaderSpace;
        REQUIRE( reader.readHeader( readHeader, imagePath.string() ) == mx::error_t::fits_read_error );

        fitsFileDetail::fitsFileCfitsioOpsInstance().getHeaderSpace = ffghsp;
        fitsFileDetail::fitsFileCfitsioOpsInstance().readKey = failReadKey;
        REQUIRE( reader.readHeader( readHeader, imagePath.string() ) == mx::error_t::fits_read_error );
    }

    SECTION( "header failure after a successful Eigen read" )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().getHeaderSpace = failGetHeaderSpace;
        mx::improc::eigenImage<float> readImage;
        fitsHeader<> readHeader;
        fitsFile<float> reader;
        REQUIRE( reader.read( readImage, readHeader, imagePath.string() ) == mx::error_t::fits_read_error );
    }

    SECTION( "header open failure" )
    {
        fitsHeader<> readHeader;
        fitsFile<float> reader;
        REQUIRE( reader.fileName( imagePath.string(), false ) == mx::error_t::noerror );
        fitsFileDetail::fitsFileCfitsioOpsInstance().openFile = failOpen;
        REQUIRE( reader.readHeader( readHeader ) == mx::error_t::fits_file_not_opened );
    }
}

/** \brief Verifies FITS writes propagate creation, image, pixel, header, and close failures.
 *
 * Exercises mx::fits::fitsFile::write error paths.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Propagating FITS write failures", "[ioutils::fits::fitsFile]" )
{
    temporaryDirectory testDirectory( "mxlib-fitsFile-write-failure-test" );
    mx::improc::eigenImage<float> image( 2, 2 );
    image.setOnes();

    SECTION( "file creation failure" )
    {
        fitsFileOpsGuard operationsGuard;
        fitsFileDetail::fitsFileCfitsioOpsInstance().createFile = failCreateFile;
        fitsFile_test<float> writer;
        REQUIRE( writer.write( ( testDirectory.path() / "create.fits" ).string(), image ) ==
                 mx::error_t::fits_file_not_created );
        REQUIRE_FALSE( writer.isOpen() );
    }

    SECTION( "image creation failure" )
    {
        fitsFile_test<float> writer;
        {
            fitsFileOpsGuard operationsGuard;
            fitsFileDetail::fitsFileCfitsioOpsInstance().createImage = failCreateImage;
            REQUIRE( writer.write( ( testDirectory.path() / "image.fits" ).string(), image ) ==
                     mx::error_t::fits_bad_naxis );
            REQUIRE( writer.isOpen() );
        }
        REQUIRE( writer.close() == mx::error_t::noerror );
    }

    SECTION( "pixel write failure" )
    {
        fitsFile_test<float> writer;
        {
            fitsFileOpsGuard operationsGuard;
            fitsFileDetail::fitsFileCfitsioOpsInstance().writePixels = failWritePixels;
            REQUIRE( writer.write( ( testDirectory.path() / "pixels.fits" ).string(), image ) ==
                     mx::error_t::fits_write_error );
            REQUIRE( writer.isOpen() );
        }
        REQUIRE( writer.close() == mx::error_t::noerror );
    }

    SECTION( "close failure" )
    {
        fitsFile_test<float> writer;
        {
            fitsFileOpsGuard operationsGuard;
            fitsFileDetail::fitsFileCfitsioOpsInstance().closeFile = failClose;
            REQUIRE( writer.write( ( testDirectory.path() / "close.fits" ).string(), image ) ==
                     mx::error_t::fits_file_not_closed );
            REQUIRE( writer.isOpen() );
        }
        REQUIRE( writer.close() == mx::error_t::noerror );
    }

    SECTION( "header-card failure" )
    {
        fitsHeader<> header;
        REQUIRE( header.append( "UNTYPED" ) == mx::error_t::noerror );
        fitsFile_test<float> writer;
        REQUIRE( writer.write( ( testDirectory.path() / "header.fits" ).string(), image, header ) !=
                 mx::error_t::noerror );
        if( writer.isOpen() )
        {
            REQUIRE( writer.close() == mx::error_t::noerror );
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

/** \brief Verifies raw-array overloads, one-dimensional images, subimages, and headerless batch reads.
 *
 * Exercises mx::fits::fitsFile::read through its raw pointer overloads, Eigen one-dimensional resize path,
 * configured two-dimensional subset, and raw multi-file overload without headers.
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Reading raw FITS arrays and configured subsets", "[ioutils::fits::fitsFile]" )
{
    temporaryDirectory testDirectory( "mxlib-fitsFile-raw-read-test" );
    fitsFile<float> fitsFile;

    const std::filesystem::path vectorPath = testDirectory.path() / "vector.fits";
    std::vector<float> vectorData{ 1.0F, 2.0F, 3.0F, 4.0F };
    REQUIRE( fitsFile.write( vectorPath.string(), vectorData.data(), 4, 0, 0 ) == mx::error_t::noerror );

    mx::error_t vectorSizeError = mx::error_t::error;
    REQUIRE( fitsFile.open( vectorPath.string() ) == mx::error_t::noerror );
    REQUIRE( fitsFile.getSize( vectorSizeError ) == 4 );
    REQUIRE( vectorSizeError == mx::error_t::noerror );
    fitsFile.setReadSize( 1, 0, 2, 0 );
    REQUIRE( fitsFile.getSize( vectorSizeError ) == 2 );
    REQUIRE( fitsFile.getSize( 0, vectorSizeError ) == 2 );
    fitsFile.setReadSize();
    REQUIRE( fitsFile.close() == mx::error_t::noerror );

    mx::improc::eigenImage<float> vectorImage;
    REQUIRE( fitsFile.read( vectorImage, vectorPath.string() ) == mx::error_t::noerror );
    REQUIRE( vectorImage.rows() == 4 );
    REQUIRE( vectorImage.cols() == 1 );
    REQUIRE( vectorImage( 3, 0 ) == 4.0F );

    mx::improc::eigenImage<float> image( 2, 3 );
    image << 1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F;
    fitsHeader<> header;
    REQUIRE( header.append( "TESTKEY", 23, "raw read header" ) == mx::error_t::noerror );

    const std::filesystem::path firstPath = testDirectory.path() / "first.fits";
    const std::filesystem::path secondPath = testDirectory.path() / "second.fits";
    REQUIRE( fitsFile.write( firstPath.string(), image.data(), image.rows(), image.cols(), 1, header ) ==
             mx::error_t::noerror );
    REQUIRE( fitsFile.write( secondPath.string(), image.data(), image.rows(), image.cols(), 1 ) ==
             mx::error_t::noerror );

    const std::filesystem::path directPath = testDirectory.path() / "direct.fits";
    REQUIRE( fitsFile.fileName( directPath.string(), false ) == mx::error_t::noerror );
    REQUIRE( fitsFile.write( image.data(), image.rows(), image.cols(), 1 ) == mx::error_t::noerror );
    REQUIRE( fitsFile.open() == mx::error_t::noerror );
    REQUIRE( fitsFile.write( image.data(), image.rows(), image.cols(), 1, header ) == mx::error_t::noerror );

    std::vector<float> raw( 6 );
    REQUIRE( fitsFile.fileName( firstPath.string(), false ) == mx::error_t::noerror );
    REQUIRE( fitsFile.read( raw.data() ) == mx::error_t::noerror );
    REQUIRE( std::equal( image.data(), image.data() + image.size(), raw.data() ) );

    fitsHeader<> readHeader;
    readHeader.append( "TESTKEY" );
    REQUIRE( fitsFile.read( raw.data(), readHeader ) == mx::error_t::noerror );
    REQUIRE( readHeader["TESTKEY"].value<int>() == 23 );
    REQUIRE( fitsFile.read( raw.data(), secondPath.string() ) == mx::error_t::noerror );
    REQUIRE( fitsFile.read( raw.data(), readHeader, firstPath.string() ) == mx::error_t::noerror );
    REQUIRE( readHeader["TESTKEY"].value<int>() == 23 );

    fitsFile.setReadSize( 1, 0, 1, 2 );
    mx::improc::eigenImage<float> subset;
    REQUIRE( fitsFile.read( subset, firstPath.string() ) == mx::error_t::noerror );
    REQUIRE( subset.rows() == 1 );
    REQUIRE( subset.cols() == 2 );
    REQUIRE( subset( 0, 0 ) == image( 1, 0 ) );
    REQUIRE( subset( 0, 1 ) == image( 1, 1 ) );
    fitsFile.setReadSize();

    std::vector<float> batch( 12 );
    const std::vector<std::string> files{ firstPath.string(), secondPath.string() };
    REQUIRE( fitsFile.read( batch.data(), files ) == mx::error_t::noerror );
    REQUIRE( std::equal( image.data(), image.data() + image.size(), batch.data() ) );
    REQUIRE( std::equal( image.data(), image.data() + image.size(), batch.data() + image.size() ) );

    const std::vector<std::string> noFiles;
    REQUIRE( fitsFile.read( batch.data(), noFiles ) == mx::error_t::invalidarg );
    std::vector<fitsHeader<>> noHeaders;
    REQUIRE( fitsFile.read( batch.data(), noHeaders, noFiles ) == mx::error_t::invalidarg );

    mx::fits::fitsFile<float> lazyReader( firstPath.string(), false );
    mx::improc::eigenImage<float> lazyImage;
    REQUIRE( lazyReader.read( lazyImage ) == mx::error_t::noerror );
    REQUIRE( lazyImage.isApprox( image ) );

    const std::string missingPath = ( testDirectory.path() / "missing.fits" ).string();
    REQUIRE( fitsFile.read( lazyImage, missingPath ) != mx::error_t::noerror );
    REQUIRE( fitsFile.read( lazyImage, readHeader, missingPath ) != mx::error_t::noerror );

    mx::improc::eigenCube<float> imageCube;
    REQUIRE( fitsFile.read( imageCube, std::vector<std::string>{ missingPath } ) != mx::error_t::noerror );

    fitsFile.setReadSize( 100, 100, 1, 1 );
    REQUIRE( fitsFile.read( imageCube, std::vector<std::string>{ firstPath.string() } ) != mx::error_t::noerror );
    fitsFile.setReadSize();

    mx::improc::eigenImage<float> smallImage( 1, 1 );
    smallImage.setConstant( 9.0F );
    const std::filesystem::path smallPath = testDirectory.path() / "small.fits";
    REQUIRE( fitsFile.write( smallPath.string(), smallImage ) == mx::error_t::noerror );
    REQUIRE( fitsFile.read( imageCube, std::vector<std::string>{ firstPath.string(), smallPath.string() } ) !=
             mx::error_t::noerror );

    fitsFile_test<float> noCommentReader;
    noCommentReader.noComment( true );
    fitsHeader<> headerWithoutComments;
    headerWithoutComments.append( "TESTKEY" );
    REQUIRE( noCommentReader.readHeader( headerWithoutComments, firstPath.string() ) == mx::error_t::noerror );
    REQUIRE( headerWithoutComments["TESTKEY"].value<int>() == 23 );
    REQUIRE( headerWithoutComments["TESTKEY"].comment().empty() );
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
    header.append( "PREPROC MASK", true, "apply preprocessing mask" );
    header.append( "QTHRESH", 0.75F, "quality threshold" );
    header.append( "", fitsCommentType(), "HCI preprocessing parameters" );
    header.append( "HISTORY", fitsHistoryType(), "coadded target000001.fits" );

    const std::filesystem::path headerImagePath = testDirectory.path() / "image-with-header.fits";
    REQUIRE( fitsFile.write( headerImagePath.string(), image, header ) == mx::error_t::noerror );
    const std::filesystem::path unavailableImagePath = testDirectory.path() / "missing" / "image.fits";
    REQUIRE( fitsFile.write( unavailableImagePath.string(), image ) != mx::error_t::noerror );
    REQUIRE( fitsFile.write( unavailableImagePath.string(), image, header ) != mx::error_t::noerror );

    fitsHeader<verbose::vv> readHeader;
    REQUIRE( fitsFile.read( readImage, readHeader, headerImagePath.string() ) == mx::error_t::noerror );
    REQUIRE( readImage.isApprox( image ) );
    REQUIRE( readHeader["TESTKEY"].value<int>() == 17 );
    REQUIRE( readHeader["PREPROC MASK"].String() == "1" );
    REQUIRE_THAT( readHeader["QTHRESH"].value<float>(), WithinAbs( 0.75F, 1e-6F ) );

    bool foundHistory = false;
    for( auto &card : readHeader )
    {
        if( card.type() == fitsType<fitsHistoryType>() && card.comment() == "coadded target000001.fits" )
        {
            foundHistory = true;
        }
    }
    REQUIRE( foundHistory );

    const std::filesystem::path rawImagePath = testDirectory.path() / "raw-image-with-header.fits";
    REQUIRE( fitsFile.write( rawImagePath.string(), image.data(), image.rows(), image.cols(), 1, header ) ==
             mx::error_t::noerror );
    REQUIRE( fitsFile.read( readImage, readHeader, rawImagePath.string() ) == mx::error_t::noerror );
    REQUIRE( readImage.isApprox( image ) );
    REQUIRE( readHeader["TESTKEY"].value<int>() == 17 );
    REQUIRE( readHeader["PREPROC MASK"].String() == "1" );

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
