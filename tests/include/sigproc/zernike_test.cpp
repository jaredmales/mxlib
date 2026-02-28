/** \file zernike_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/zernike.hpp"/// testing noll_nm


/// testing noll_nm
/**
 * \ingroup zernike_unit_tests
 */
TEST_CASE( "testing noll_nm", "[sigproc::zernike]" )
{

    SECTION( "a j value" )
    {
        SECTION( "j==0" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 0 );
            REQUIRE( rv == -1 );
        }

        SECTION( "j==1" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 1 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 0 );
            REQUIRE( m == 0 );
        }

        SECTION( "j==2" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 2 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 1 );
            REQUIRE( m == 1 );
        }

        SECTION( "j==3" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 3 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 1 );
            REQUIRE( m == -1 );
        }

        SECTION( "j==4" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 4 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 2 );
            REQUIRE( m == 0 );
        }

        SECTION( "j==5" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 5 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 2 );
            REQUIRE( m == -2 );
        }

        SECTION( "j==6" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 6 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 2 );
            REQUIRE( m == 2 );
        }

        SECTION( "j==7" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 7 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 3 );
            REQUIRE( m == -1 );
        }

        SECTION( "j==8" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 8 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 3 );
            REQUIRE( m == 1 );
        }

        SECTION( "j==9" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 9 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 3 );
            REQUIRE( m == -3 );
        }

        SECTION( "j==10" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 10 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 3 );
            REQUIRE( m == 3 );
        }

        SECTION( "j==11" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 11 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 4 );
            REQUIRE( m == 0 );
        }

        SECTION( "j==12" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 12 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 4 );
            REQUIRE( m == 2 );
        }

        SECTION( "j==13" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 13 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 4 );
            REQUIRE( m == -2 );
        }

        SECTION( "j==14" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 14 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 4 );
            REQUIRE( m == 4 );
        }

        SECTION( "j==15" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 15 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 4 );
            REQUIRE( m == -4 );
        }

        SECTION( "j==16" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 16 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 5 );
            REQUIRE( m == 1 );
        }

        SECTION( "j==17" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 17 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 5 );
            REQUIRE( m == -1 );
        }

        SECTION( "j==18" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 18 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 5 );
            REQUIRE( m == 3 );
        }

        SECTION( "j==19" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 19 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 5 );
            REQUIRE( m == -3 );
        }

        SECTION( "j==20" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 20 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 5 );
            REQUIRE( m == 5 );
        }

        SECTION( "j==21" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 21 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 5 );
            REQUIRE( m == -5 );
        }

        SECTION( "j==22" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 22 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 6 );
            REQUIRE( m == 0 );
        }

        SECTION( "j==23" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 23 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 6 );
            REQUIRE( m == -2 );
        }

        SECTION( "j==24" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 24 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 6 );
            REQUIRE( m == 2 );
        }

        SECTION( "j==25" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 25 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 6 );
            REQUIRE( m == -4 );
        }

        SECTION( "j==26" )
        {
            int m, n;
            int rv = mx::sigproc::noll_nm( n, m, 26 );
            REQUIRE( rv == 0 );
            REQUIRE( n == 6 );
            REQUIRE( m == 4 );
        }
    }
}/// testing zernikeQNorm


/// testing zernikeQNorm
/**
 * \ingroup zernike_unit_tests
 */
TEST_CASE( "testing zernikeQNorm", "[sigproc::zernike]" )
{
    SECTION( "an array" )
    {
        SECTION( "j==1" )
        {
            Eigen::Array<double, -1, -1> arr, k, phi;
            arr.resize( 32, 32 );
            k.resize( 32, 32 );
            phi.resize( 32, 32 );

            for( int i = 0; i < 32; ++i )
            {
                for( int j = 0; j < 32; ++j )
                {
                    double kx = i - 5;
                    double ky = j - 15;
                    k( i, j ) = sqrt( kx * kx + ky * ky );
                    phi( i, j ) = atan( ky / kx );
                }
            }
            int rv = mx::sigproc::zernikeQNorm( arr, k, phi, 1 );
            REQUIRE( rv == 0 );
        }
    }
}
