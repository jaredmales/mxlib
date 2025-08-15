/** \file fitsFile_test.cpp
 */
#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ioutils/fits/fitsFile.hpp"
using namespace mx::fits;

template <typename dataT>
class fitsFile_test : public fitsFile<dataT>
{
  public:
    typedef fitsFile<dataT>::pixarrT pixarrTT;

    void calcPixarrs( pixarrTT & pixarrs )
    {
        fitsFile<dataT>::calcPixarrs( pixarrs );
    }

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
};

/** Verify calculation of subimage sizes
 *
 * \anchor tests_ioutils_fits_fitsFile_subimage_sizes
 */
SCENARIO( "fitsFile calculating subimage sizes", "[ioutils::fits::fitsFile]" )
{
    GIVEN( "a 1D image to read" )
    {
        WHEN( "reading the whole image with default setup" )
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

        WHEN( "reading a subimage" )
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

    GIVEN( "a 2D image to read" )
    {
        WHEN( "reading the whole image with default setup" )
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

        WHEN( "reading a subimage" )
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

    GIVEN( "a 3D image to read" )
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

        WHEN( "reading a subimage" )
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
