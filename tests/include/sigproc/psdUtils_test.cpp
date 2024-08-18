/** \file psdUtils_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/psdFilter.hpp"
#include "../../../include/sigproc/psdUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"
#include "../../../include/math/randomT.hpp"
#include "../../../include/math/vectorUtils.hpp"

/** Scenario: calculating variance from a 1D PSD
 *
 * Verify calculations of psdVar1sided, psdVar2sided, and psdVar.
 *
 * \anchor tests_sigproc_psdUtils_psdVar_1D
 */
SCENARIO( "calculating variance from a 1D PSD", "[sigproc::psdUtils]" )
{
    GIVEN( "a 1 sided PSD," )
    {
        WHEN( "a flat PSD, midpoint rule" )
        {
            std::vector<double> f( 5 ), psd( 5 );
            for( size_t n = 0; n < f.size(); ++n )
            {
                f[n] = n;
                psd[n] = 1;
            }

            REQUIRE( mx::sigproc::psdVar( f, psd, 1.0 ) == 5 );
        }
        WHEN( "a flat PSD, trapezoid rule by default" )
        {
            std::vector<double> f( 5 ), psd( 5 );
            for( size_t n = 0; n < f.size(); ++n )
            {
                f[n] = n;
                psd[n] = 1;
            }

            REQUIRE( mx::sigproc::psdVar( f, psd ) == 4 );
        }
    }
    GIVEN( "a 2 sided PSD," )
    {
        WHEN( "a flat PSD, midpoint rule" )
        {
            std::vector<double> f( 6 ), psd( 6 );
            f[0] = 0;
            f[1] = 1;
            f[2] = 2;
            f[3] = 3;
            f[4] = -2;
            f[5] = -1;

            for( size_t n = 0; n < f.size(); ++n )
            {
                psd[n] = 1;
            }
            psd[3] = 2; // This one gets twice the power

            REQUIRE( mx::sigproc::psdVar( f, psd, 1.0 ) == 7 );
        }
        WHEN( "a flat PSD, trapezoid rule by default" )
        {
            std::vector<double> f( 6 ), psd( 6 );
            f[0] = 0;
            f[1] = 1;
            f[2] = 2;
            f[3] = 3;
            f[4] = -2;
            f[5] = -1;

            for( size_t n = 0; n < f.size(); ++n )
            {
                psd[n] = 1;
            }
            psd[3] = 2; // This one gets twice the power

            REQUIRE( mx::sigproc::psdVar( f, psd ) == 6 );
        }
    }
}

/** Verify scaling and normalization of augment1SidedPSD
 *
 * \anchor tests_sigproc_psdUtils_augment1SidedPSD
 */
SCENARIO( "augmenting a 1 sided PSD", "[sigproc::psdUtils]" )
{
    GIVEN( "a 1 sided PSD, with a 0 freq value" )
    {
        WHEN( "1/f^2" )
        {
            std::vector<double> f( 5 ), psd( 5 );

            for( size_t n = 0; n < psd.size(); ++n )
            {
                f[n] = n;
                psd[n] = pow( f[n], -2. );
            }
            psd[0] = psd[1];

            mx::sigproc::normPSD( psd, f, 1.0 );
            // Now have a 1/f^2 PSD with total 1-sided variance of 1.0.
            REQUIRE( fabs( mx::sigproc::psdVar( f, psd, 1.0 ) - 1.0 ) < 1e-10 ); // handles epsilon

            // proceed to augment:
            std::vector<double> f2s, psd2s;
            mx::sigproc::augment1SidedPSDFreq( f2s, f );
            mx::sigproc::augment1SidedPSD( psd2s, psd );

            REQUIRE( f2s.size() == 8 );
            REQUIRE( psd2s.size() == 8 );

            REQUIRE( f2s[0] == 0 );
            REQUIRE( f2s[1] == 1 );
            REQUIRE( f2s[2] == 2 );
            REQUIRE( f2s[3] == 3 );
            REQUIRE( f2s[4] == 4 );
            REQUIRE( f2s[5] == -3 );
            REQUIRE( f2s[6] == -2 );
            REQUIRE( f2s[7] == -1 );

            // Should now have 1.0 in bin 0, 0.5 in all other bins.
            REQUIRE( psd2s[0] == psd[0] );
            REQUIRE( psd2s[1] == 0.5 * psd[1] );
            REQUIRE( psd2s[2] == 0.5 * psd[2] );
            REQUIRE( psd2s[3] == 0.5 * psd[3] );
            REQUIRE( psd2s[4] == psd[4] );
            REQUIRE( psd2s[5] == 0.5 * psd[3] );
            REQUIRE( psd2s[6] == 0.5 * psd[2] );
            REQUIRE( psd2s[7] == 0.5 * psd[1] );

            // handle machine precision
            REQUIRE( fabs( mx::sigproc::psdVar( f2s, psd2s, 1.0 ) - 1.0 ) < 1e-10 );
        }
    }
}

/** Verify creation of a 1D frequency grid
 *
 * \anchor tests_sigproc_psdUtils_frequencyGrid_1D
 */
SCENARIO( "creating a 1D frequency grid", "[sigproc::psdUtils]" )
{
    GIVEN( "2 sided FFT-order frequency grid" )
    {
        WHEN( "dt = 1" )
        {
            std::vector<double> f( 10 );

            REQUIRE( mx::sigproc::frequencyGrid( f, 1.0 ) == 0 );

            REQUIRE( fabs( f[0] - 0 ) < 1e-10 );
            REQUIRE( fabs( f[1] - 0.1 ) < 1e-10 );
            REQUIRE( fabs( f[2] - 0.2 ) < 1e-10 );
            REQUIRE( fabs( f[3] - 0.3 ) < 1e-10 );
            REQUIRE( fabs( f[4] - 0.4 ) < 1e-10 );
            REQUIRE( fabs( f[5] - 0.5 ) < 1e-10 );
            REQUIRE( fabs( f[6] - -0.4 ) < 1e-10 );
            REQUIRE( fabs( f[7] - -0.3 ) < 1e-10 );
            REQUIRE( fabs( f[8] - -0.2 ) < 1e-10 );
            REQUIRE( fabs( f[9] - -0.1 ) < 1e-10 );
        }

        WHEN( "dt = 2" )
        {
            std::vector<double> f( 10 );

            REQUIRE( mx::sigproc::frequencyGrid( f, 2.5 ) == 0 );

            REQUIRE( fabs( f[0] - 0 ) < 1e-10 );
            REQUIRE( fabs( f[1] - 0.04 ) < 1e-10 );
            REQUIRE( fabs( f[2] - 0.08 ) < 1e-10 );
            REQUIRE( fabs( f[3] - 0.12 ) < 1e-10 );
            REQUIRE( fabs( f[4] - 0.16 ) < 1e-10 );
            REQUIRE( fabs( f[5] - 0.2 ) < 1e-10 );
            REQUIRE( fabs( f[6] - -0.16 ) < 1e-10 );
            REQUIRE( fabs( f[7] - -0.12 ) < 1e-10 );
            REQUIRE( fabs( f[8] - -0.08 ) < 1e-10 );
            REQUIRE( fabs( f[9] - -0.04 ) < 1e-10 );
        }
    }

    GIVEN( "1 sided frequency grid" )
    {
        WHEN( "dt = 1, odd size" )
        {
            std::vector<double> f( 5 );

            REQUIRE( mx::sigproc::frequencyGrid( f, 1.0, false ) == 0 );

            REQUIRE( fabs( f[0] - 0.1 ) < 1e-10 );
            REQUIRE( fabs( f[1] - 0.2 ) < 1e-10 );
            REQUIRE( fabs( f[2] - 0.3 ) < 1e-10 );
            REQUIRE( fabs( f[3] - 0.4 ) < 1e-10 );
            REQUIRE( fabs( f[4] - 0.5 ) < 1e-10 );
        }

        WHEN( "dt = 1, even size" )
        {
            std::vector<double> f( 6 );

            REQUIRE( mx::sigproc::frequencyGrid( f, 1.0, false ) == 0 );

            REQUIRE( fabs( f[0] - 0.0 ) < 1e-10 );
            REQUIRE( fabs( f[1] - 0.1 ) < 1e-10 );
            REQUIRE( fabs( f[2] - 0.2 ) < 1e-10 );
            REQUIRE( fabs( f[3] - 0.3 ) < 1e-10 );
            REQUIRE( fabs( f[4] - 0.4 ) < 1e-10 );
            REQUIRE( fabs( f[5] - 0.5 ) < 1e-10 );
        }
    }
}
