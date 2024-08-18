/** \file psdFilter_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/signalWindows.hpp"
#include "../../../include/improc/eigenImage.hpp"
#include "../../../include/improc/milkImage.hpp"

/** Scenario: creating 2D Rectangular Tukey Windows
 *
 * Verify creation of a 2D rectangular Tukey window
 *
 * \anchor tests_sigproc_signalWindows_Tukey2DSquare
 */
SCENARIO( "creating 2D Rectangular Tukey Windows", "[sigproc::signalWindows::tukey2dSquare]" )
{
    GIVEN( "a centered square array" )
    {
        WHEN( "256x256, alpha=0" )
        {

            mx::improc::eigenImage<float> win;
            win.resize( 256, 256 );

            mx::sigproc::window::tukey2dSquare<float>( win.data(),
                                                       win.rows(),
                                                       win.cols(),
                                                       256,
                                                       256,
                                                       0.0,
                                                       0.5 * ( win.rows() - 1.0 ),
                                                       0.5 * ( win.cols() - 1.0 ) );

            REQUIRE( win.sum() == 256 * 256 );
        }

        WHEN( "256x256, alpha=1" )
        {

            mx::improc::eigenImage<float> win;
            win.resize( 256, 256 );

            mx::sigproc::window::tukey2dSquare<float>( win.data(),
                                                       win.rows(),
                                                       win.cols(),
                                                       256,
                                                       256,
                                                       1.0,
                                                       0.5 * ( win.rows() - 1.0 ),
                                                       0.5 * ( win.cols() - 1.0 ) );

            std::vector<float> win1( 256 );
            mx::sigproc::window::tukey<float>( win1, 1.0 );

            REQUIRE( win( 0, 0 ) == win1[0] * win1[0] );
            REQUIRE( win( 10, 15 ) == win1[10] * win1[15] );
        }
        WHEN( "256x256, alpha=0.5" )
        {

            mx::improc::eigenImage<float> win;
            win.resize( 256, 256 );

            mx::sigproc::window::tukey2dSquare<float>( win.data(),
                                                       win.rows(),
                                                       win.cols(),
                                                       256,
                                                       256,
                                                       0.5,
                                                       0.5 * ( win.rows() - 1.0 ),
                                                       0.5 * ( win.cols() - 1.0 ) );

            std::vector<float> win1( 256 );
            mx::sigproc::window::tukey<float>( win1, 0.5 );

            mx::improc::milkImage<float> mwin;
            mwin.create( "win", win );

            REQUIRE( win( 0, 0 ) == win1[0] * win1[0] );
            REQUIRE( win( 10, 15 ) == win1[10] * win1[15] );
        }
    }
}
