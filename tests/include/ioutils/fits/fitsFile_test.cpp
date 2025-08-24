/** \file fitsFile_test.cpp
 */
#include "../../../catch2/catch.hpp"

using namespace Catch::Matchers;

#include "../../../../include/ioutils/fits/fitsFile.hpp"
using namespace mx::fits;

#include "../../../../include/improc/eigenImage.hpp"

namespace unitTest
{
namespace ioutilsTest
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

    mx::error_t writeTestFrame(const std::string & fname, int n)
    {
        mx::improc::eigenImage<float> im(1024,1024);

        for(int cc = 0; cc < im.cols(); ++cc)
        {
            for(int rr = 0; rr < im.rows(); ++rr)
            {
                im(rr,cc) = cc*im.rows() + rr;
            }
        }

        fitsHeader fh;
        fh.append(std::format("FLOAT{}",n), static_cast<float>(1.1), "testing 1.1");
        fh.append(std::format("INT{}",n), static_cast<int>(2), "testing 1.1");

        return this->write(std::format("{}{}.fits",fname,n), im, fh);
    }
};


/// Calculating subimage sizes
/** Verify calculation of subimage sizes
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Calculating subimage sizes", "[ioutils::fits::fitsFile]" )
{
    //document what protected members are tested here
    #ifdef MXLIB_DOXYGEN_PROTECTED_REF
    fitsFile<float> ff;
    fitsFile<float>::pixarrTT pixarrs;
    ff.calcPixarrs( pixarrs );
    ff.naxis();
    ff.naxes(0);
    ff.setReadSize( 10, 9, 7, 10 )
    ff.setCubeReadSize( 5, 3 );
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
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Basic writing and reading", "[ioutils::fits::fitsFile]" )
{
    //document what fitsFile_test is doing
    #ifdef MXLIB_DOXYGEN_PROTECTED_REF
    fitsFile<float> ff;
    fitsHeader fh;
    fh.append("dummy", 2, "dummy");
    mx::improc::eigenImage<float> im;
    ff.write("dummy", im, fh);
    #endif

    SECTION("write a frame and read it back in")
    {
        fitsFile_test<float> fft;

        REQUIRE(fft.writeTestFrame("/tmp/fitsFile_test", 1) == mx::error_t::noerror);

        fitsFile<float> ff;
        fitsHeader fh;
        mx::improc::eigenImage<float> im;

        REQUIRE(ff.read(im, fh, "/tmp/fitsFile_test1.fits") == mx::error_t::noerror);

        REQUIRE_THAT(fh["FLOAT1"].Float(), WithinAbs(1.1, 1e-5));
        REQUIRE(fh["INT1"].Int() == 2);
        REQUIRE(im.rows() == 1024);
        REQUIRE(im.cols() == 1024);


    }

}

/// Reading headers
/**
 *
 * \ingroup fitsFile_unit_tests
 */
TEST_CASE( "Reading headers", "[ioutils::fits::fitsFile]" )
{
    //document what fitsFile_test is doing
    #ifdef MXLIB_DOXYGEN_PROTECTED_REF
    fitsFile<float> ff;
    fitsHeader fh;
    fh.append("dummy", 2, "dummy");
    mx::improc::eigenImage<float> im;
    ff.write("dummy", im, fh);
    #endif

    fitsFile_test<float> fft;

    std::vector<std::string> fnames({"/tmp/fitsFile_test1.fits","/tmp/fitsFile_test2.fits","/tmp/fitsFile_test3.fits"});
    REQUIRE(fft.writeTestFrame("/tmp/fitsFile_test",1) == mx::error_t::noerror);
    REQUIRE(fft.writeTestFrame("/tmp/fitsFile_test",2) == mx::error_t::noerror);
    REQUIRE(fft.writeTestFrame("/tmp/fitsFile_test",3) == mx::error_t::noerror);

    SECTION("write three frames and read their headers as a vector, not allocated")
    {

        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        REQUIRE(ff.readHeader(fh, fnames) == mx::error_t::noerror);
        REQUIRE(fh.size() == fnames.size());

        REQUIRE_THAT(fh[0]["FLOAT1"].Float(), WithinAbs(1.1, 1e-5));
        REQUIRE(fh[0]["INT1"].Int() == 2);
        REQUIRE_THAT(fh[1]["FLOAT2"].Float(), WithinAbs(1.1, 1e-5));
        REQUIRE(fh[1]["INT2"].Int() == 2);
        REQUIRE_THAT(fh[2]["FLOAT3"].Float(), WithinAbs(1.1, 1e-5));
        REQUIRE(fh[2]["INT3"].Int() == 2);
    }

    SECTION("write three frames and read their headers as a vector, allocated")
    {

        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        fh.resize(fnames.size());
        REQUIRE(ff.readHeader(fh, fnames) == mx::error_t::noerror);

        mx::error_t errc;
        errc = mx::error_t::error;
        REQUIRE_THAT(fh[0]["FLOAT1"].Float(&errc), WithinAbs(1.1, 1e-5));
        REQUIRE(errc == mx::error_t::noerror);
        errc = mx::error_t::error;
        REQUIRE(fh[0]["INT1"].Int(&errc) == 2);
        REQUIRE(errc == mx::error_t::noerror);
        errc = mx::error_t::error;
        REQUIRE_THAT(fh[1]["FLOAT2"].Float(&errc), WithinAbs(1.1, 1e-5));
        REQUIRE(errc == mx::error_t::noerror);
        errc = mx::error_t::error;
        REQUIRE(fh[1]["INT2"].Int(&errc) == 2);
        REQUIRE(errc == mx::error_t::noerror);
        errc = mx::error_t::error;
        REQUIRE_THAT(fh[2]["FLOAT3"].Float(&errc), WithinAbs(1.1, 1e-5));
        REQUIRE(errc == mx::error_t::noerror);
        errc = mx::error_t::error;
        REQUIRE(fh[2]["INT3"].Int(&errc) == 2);
        REQUIRE(errc == mx::error_t::noerror);
    }

    SECTION("write three frames and read their headers as a vector, allocated to wrong size")
    {
        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        fh.resize(fnames.size()+1);
        REQUIRE(ff.readHeader(fh, fnames) == mx::error_t::invalidarg);
    }

    SECTION("write three frames and read their headers as a vector, extracting only 1 keyword")
    {

        fitsFile<float> ff;
        std::vector<mx::fits::fitsHeader<mx::verbose::d>> fh;
        mx::improc::eigenImage<float> im;

        fh.resize(fnames.size());
        fh[0].append("FLOAT1");
        fh[1].append("FLOAT2");
        fh[2].append("FLOAT3");


        REQUIRE(ff.readHeader(fh, fnames) == mx::error_t::noerror);

        mx::error_t errc;

        REQUIRE_THAT(fh[0]["FLOAT1"].Float(&errc), WithinAbs(1.1, 1e-5));
        REQUIRE(errc == mx::error_t::noerror);
        REQUIRE(fh[0]["FLOAT1"].valueGood() == true);

        REQUIRE(fh[0]["INT1"].Int(&errc) == 0);
        REQUIRE(errc == mx::error_t::invalidarg);
        REQUIRE(fh[0]["INT1"].valueGood() == false);

        REQUIRE_THAT(fh[1]["FLOAT2"].Float(&errc), WithinAbs(1.1, 1e-5));
        REQUIRE(errc == mx::error_t::noerror);
        REQUIRE(fh[1]["FLOAT2"].valueGood() == true);

        REQUIRE(fh[1]["INT2"].Int(&errc) == 0);
        REQUIRE(errc == mx::error_t::invalidarg);
        REQUIRE(fh[1]["INT2"].valueGood() == false);

        REQUIRE_THAT(fh[2]["FLOAT3"].Float(&errc), WithinAbs(1.1, 1e-5));
        REQUIRE(errc == mx::error_t::noerror);
        REQUIRE(fh[2]["FLOAT3"].valueGood() == true);

        REQUIRE(fh[2]["INT3"].Int(&errc) == 0);
        REQUIRE(errc == mx::error_t::invalidarg);
        REQUIRE(fh[2]["INT3"].valueGood() == false);
    }

}

} // namespace fitsFile
} // namespace fits
} // namespace ioutils
} // namespace test
