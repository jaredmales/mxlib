/** \file randomT_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/math/randomT.hpp"

/** Verify compilation and basic operation of randomT with std::distributions.
 * Basic tests include verification that seeding works.  Note that there is the very slight
 * possibility that normal operation could return two consecutive variates with the same value
 * which would cause some of the below tests to fail.  Will probably never happen . . .
 *
 * \anchor tests_math_randomT_basic
 */
SCENARIO( "Verify compilation and basic operation of randomT with std::distributions", "[math::randomT]" )
{
    GIVEN( "a uniform distribution is desired" )
    {
        WHEN( "two double random numbers from same generator requested" )
        {
            mx::math::randomT<double, std::mt19937_64, std::uniform_real_distribution<double>> uniDist;

            double r1 = uniDist;
            double r2 = uniDist;

            REQUIRE( r1 != r2 );
        }

        WHEN( "two double random numbers requested from different generators with same seed" )
        {
            mx::math::randomT<double, std::mt19937_64, std::uniform_real_distribution<double>> uniDist1( false );
            mx::math::randomT<double, std::mt19937_64, std::uniform_real_distribution<double>> uniDist2( false );

            uniDist1.seed( 10 );
            uniDist2.seed( 10 );
            double r1 = uniDist1;
            double r2 = uniDist2;

            REQUIRE( r1 == r2 );
        }

        WHEN( "two double random numbers requested from different generators with different seed" )
        {
            mx::math::randomT<double, std::mt19937_64, std::uniform_real_distribution<double>> uniDist1( false );
            mx::math::randomT<double, std::mt19937_64, std::uniform_real_distribution<double>> uniDist2( false );

            uniDist1.seed( 10 );
            uniDist2.seed( 11 );
            double r1 = uniDist1;
            double r2 = uniDist2;

            REQUIRE( r1 != r2 );
        }
    }

    GIVEN( "a normal distribution is desired" )
    {
        WHEN( "two double random numbers from same generator requested" )
        {
            mx::math::randomT<double, std::mt19937_64, std::normal_distribution<double>> normDist;

            double r1 = normDist;
            double r2 = normDist;

            REQUIRE( r1 != r2 );
        }

        WHEN( "two double random numbers requested from different generators with same seed" )
        {
            mx::math::randomT<double, std::mt19937_64, std::normal_distribution<double>> normDist1( false );
            mx::math::randomT<double, std::mt19937_64, std::normal_distribution<double>> normDist2( false );

            normDist1.seed( 10 );
            normDist2.seed( 10 );
            double r1 = normDist1;
            double r2 = normDist2;

            REQUIRE( r1 == r2 );
        }

        WHEN( "two double random numbers requested from different generators with different seed" )
        {
            mx::math::randomT<double, std::mt19937_64, std::normal_distribution<double>> normDist1( false );
            mx::math::randomT<double, std::mt19937_64, std::normal_distribution<double>> normDist2( false );

            normDist1.seed( 10 );
            normDist2.seed( 11 );
            double r1 = normDist1;
            double r2 = normDist2;

            REQUIRE( r1 != r2 );
        }
    }
}

/** Verify compilation and basic operation of randomT with the Lapace distribution
 * Basic tests include verification that seeding works.  Note that there is the very slight
 * possibility that normal operation could return two consecutive variates with the same value
 * which would cause some of the below tests to fail.  Will probably never happen . . .
 *
 * \anchor tests_math_randomT_basic_laplace
 */
SCENARIO( "Verify compilation and basic operation of randomT with the Lapace distribution",
          "[math::laplace_distribution]" )
{
    GIVEN( "a laplace distribution is desired" )
    {
        WHEN( "two double random numbers from same generator requested" )
        {
            mx::math::randomT<double, std::mt19937_64, mx::math::laplace_distribution<double>> lapDist;

            double r1 = lapDist;
            double r2 = lapDist;

            REQUIRE( r1 != r2 );
        }

        WHEN( "two double random numbers requested from different generators with same seed" )
        {
            mx::math::randomT<double, std::mt19937_64, mx::math::laplace_distribution<double>> lapDist1( false );
            mx::math::randomT<double, std::mt19937_64, mx::math::laplace_distribution<double>> lapDist2( false );

            lapDist1.seed( 10 );
            lapDist2.seed( 10 );
            double r1 = lapDist1;
            double r2 = lapDist2;

            REQUIRE( r1 == r2 );
        }

        WHEN( "two double random numbers requested from different generators with different seed" )
        {
            mx::math::randomT<double, std::mt19937_64, mx::math::laplace_distribution<double>> lapDist1( false );
            mx::math::randomT<double, std::mt19937_64, mx::math::laplace_distribution<double>> lapDist2( false );

            lapDist1.seed( 10 );
            lapDist2.seed( 11 );
            double r1 = lapDist1;
            double r2 = lapDist2;

            REQUIRE( r1 != r2 );
        }
    }
}
