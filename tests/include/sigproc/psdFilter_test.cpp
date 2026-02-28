/** \file psdFilter_test.cpp
 */
#include "../../catch2/catch.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#define MX_NO_ERROR_REPORTS

#include "../../../include/sigproc/psdFilter.hpp"
#include "../../../include/sigproc/psdUtils.hpp"
#include "../../../include/improc/eigenCube.hpp"
#include "../../../include/math/randomT.hpp"
#include "../../../include/math/vectorUtils.hpp"

#ifdef MXLIB_BUILD_COVERAGE
constexpr int psdFilterTrials = 100;
constexpr double psdFilterTol = 0.09;
#else
constexpr int psdFilterTrials = 10000;
constexpr double psdFilterTol = 0.02;
#endif/// compiling psdFilter

namespace
{
class TrialProgressLogger
{
  public:
    TrialProgressLogger( const char *label, int totalTrials )
        : m_label( label ),
          m_totalTrials( totalTrials ),
          m_reportInterval( std::max( 1, totalTrials / 4 ) ),
          m_start( std::chrono::steady_clock::now() )
    {
        std::cerr << "[psdFilter_test] " << m_label << ": starting " << m_totalTrials << " trials" << std::endl;
    }

    void update( int trialIndex )
    {
        const int trialCount = trialIndex + 1;
        if( trialCount % m_reportInterval != 0 && trialCount != m_totalTrials )
            return;

        std::cerr << "[psdFilter_test] " << m_label << ": " << trialCount << "/" << m_totalTrials
                  << " trials complete, elapsed " << elapsedSeconds() << " s" << std::endl;
    }

  private:
    double elapsedSeconds() const
    {
        return std::chrono::duration_cast<std::chrono::duration<double>>( std::chrono::steady_clock::now() - m_start )
            .count();
    }

    const char *m_label{ nullptr };
    int m_totalTrials{ 0 };
    int m_reportInterval{ 1 };
    std::chrono::steady_clock::time_point m_start;
};
} // namespace


/// compiling psdFilter
/**
 * \ingroup psdFilter_unit_tests
 */
TEST_CASE( "compiling psdFilter", "[sigproc::psdFilter]" )
{
    SECTION( "a psdFilter, sqrt pointer" )
    {
        SECTION( "rank==1" )
        {
            mx::sigproc::psdFilter<double, 1> psdF;

            std::vector<double> psdSqrt( 1024, 1 );

            int rv = psdF.psdSqrt( &psdSqrt, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 1024 );
            REQUIRE( psdF.cols() == 1 );
            REQUIRE( psdF.planes() == 1 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
        SECTION( "rank==2" )
        {
            mx::sigproc::psdFilter<double, 2> psdF;

            Eigen::Array<double, -1, -1> psdSqrt;
            psdSqrt.resize( 256, 256 );
            psdSqrt.setConstant( 1 );

            int rv = psdF.psdSqrt( &psdSqrt, 1, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 256 );
            REQUIRE( psdF.cols() == 256 );
            REQUIRE( psdF.planes() == 1 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
        SECTION( "rank==3" )
        {
            mx::sigproc::psdFilter<double, 3> psdF;

            mx::improc::eigenCube<double> psdSqrt;
            psdSqrt.resize( 128, 128, 256 );
            // psdSqrt.setConstant(1);

            int rv = psdF.psdSqrt( &psdSqrt, 1, 1, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 128 );
            REQUIRE( psdF.cols() == 128 );
            REQUIRE( psdF.planes() == 256 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
    }

    SECTION( "a psdFilter, sqrt reference" )
    {
        SECTION( "rank==1" )
        {
            mx::sigproc::psdFilter<double, 1> psdF;

            std::vector<double> psdSqrt( 1024, 1 );

            int rv = psdF.psdSqrt( psdSqrt, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 1024 );
            REQUIRE( psdF.cols() == 1 );
            REQUIRE( psdF.planes() == 1 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
        SECTION( "rank==2" )
        {
            mx::sigproc::psdFilter<double, 2> psdF;

            Eigen::Array<double, -1, -1> psdSqrt;
            psdSqrt.resize( 256, 256 );
            psdSqrt.setConstant( 1 );

            int rv = psdF.psdSqrt( psdSqrt, 1, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 256 );
            REQUIRE( psdF.cols() == 256 );
            REQUIRE( psdF.planes() == 1 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
        SECTION( "rank==3" )
        {
            mx::sigproc::psdFilter<double, 3> psdF;

            mx::improc::eigenCube<double> psdSqrt;
            psdSqrt.resize( 128, 128, 256 );
            // psdSqrt.setConstant(1);

            int rv = psdF.psdSqrt( psdSqrt, 1, 1, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 128 );
            REQUIRE( psdF.cols() == 128 );
            REQUIRE( psdF.planes() == 256 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
    }

    SECTION( "a psdFilter, psd reference" )
    {
        SECTION( "rank==1" )
        {
            mx::sigproc::psdFilter<double, 1> psdF;

            std::vector<double> psd( 1024, 1 );

            int rv = psdF.psd( psd, 1.0 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 1024 );
            REQUIRE( psdF.cols() == 1 );
            REQUIRE( psdF.planes() == 1 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
        SECTION( "rank==2" )
        {
            mx::sigproc::psdFilter<double, 2> psdF;

            Eigen::Array<double, -1, -1> psd;
            psd.resize( 256, 256 );
            psd.setConstant( 1 );

            int rv = psdF.psd( psd, 1.0, 1.0 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 256 );
            REQUIRE( psdF.cols() == 256 );
            REQUIRE( psdF.planes() == 1 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
        SECTION( "rank==3" )
        {
            mx::sigproc::psdFilter<double, 3> psdF;

            mx::improc::eigenCube<double> psd;
            psd.resize( 128, 128, 256 );
            // psdSqrt.setConstant(1);

            int rv = psdF.psd( psd, 1, 1, 1 );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 128 );
            REQUIRE( psdF.cols() == 128 );
            REQUIRE( psdF.planes() == 256 );

            psdF.clear();
            REQUIRE( psdF.rows() == 0 );
            REQUIRE( psdF.cols() == 0 );
            REQUIRE( psdF.planes() == 0 );
        }
    }
}

/** Verify filtering and noise normalization
 * Conducts random noise tests, verifying that the resultant rms is within 2% of expected value on average over many
 * trials. Results are usually better than 1%, but 2% makes sure we don't get false failures.
 *
 *//// filtering with psdFilter

/// filtering with psdFilter
/**
 * \ingroup psdFilter_unit_tests
 */
TEST_CASE( "filtering with psdFilter", "[sigproc::psdFilter]" )
{
    SECTION( "a rank 1 psd" )
    {
        SECTION( "alpha=-2.5, df Nyquist matched to array size, var=1" )
        {
            mx::sigproc::psdFilter<double, 1> psdF;

            std::vector<double> f( 2049 ), psd( 2049 );

            for( size_t n = 0; n < psd.size(); ++n )
                f[n] = n * 1.0 / 4096.;

            for( size_t n = 1; n < psd.size(); ++n )
                psd[n] = pow( f[n], -2.5 );
            psd[0] = psd[1];

            mx::sigproc::normPSD( psd, f, 1.0, -1e5, 1e5 );

            std::vector<double> f2s, psd2s;
            mx::sigproc::augment1SidedPSDFreq( f2s, f );
            mx::sigproc::augment1SidedPSD( psd2s, psd );

            double df = f2s[1] - f2s[0];
            // mx::sigproc::normPSD(psd2s, f2s, 1.0, -1e5, 1e5);
            int rv = psdF.psd( psd2s, df );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 2. * psd.size() - 2 );
            REQUIRE( psdF.cols() == 1 );
            REQUIRE( psdF.planes() == 1 );

            std::vector<double> noise( psdF.rows() );

            mx::math::normDistT<double> normVar;

            double avgRms = 0;
            TrialProgressLogger trialLogger( "rank1 alpha=-2.5 var=1", psdFilterTrials );

            for( int k = 0; k < psdFilterTrials; ++k )
            {
                for( size_t n = 0; n < noise.size(); ++n )
                    noise[n] = normVar;
                psdF( noise );
                avgRms += ( mx::math::vectorVariance( noise, 0.0 ) );
                trialLogger.update( k );
            }

            avgRms = sqrt( avgRms / psdFilterTrials );

            REQUIRE_THAT( avgRms, Catch::Matchers::WithinAbs( 1.0, psdFilterTol ) );
        }
        SECTION( "alpha=-1.5, df arbitrary, var = 2.2" )
        {
            mx::sigproc::psdFilter<double, 1> psdF;

            std::vector<double> f( 1025 ), psd( 1025 );

            for( size_t n = 0; n < psd.size(); ++n )
                f[n] = n * 1.0 / 7000.;

            for( size_t n = 1; n < psd.size(); ++n )
                psd[n] = pow( f[n], -1.5 );
            psd[0] = psd[1];

            std::vector<double> f2s, psd2s;
            mx::sigproc::augment1SidedPSDFreq( f2s, f );
            mx::sigproc::augment1SidedPSD( psd2s, psd );

            double df = f2s[1] - f2s[0];
            mx::sigproc::normPSD( psd2s, f2s, 2.2, -1e5, 1e5 );

            int rv = psdF.psd( psd2s, df );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == 2. * psd.size() - 2 );
            REQUIRE( psdF.cols() == 1 );
            REQUIRE( psdF.planes() == 1 );

            std::vector<double> noise( psdF.rows() );

            mx::math::normDistT<double> normVar;

            double avgRms = 0;
            TrialProgressLogger trialLogger( "rank1 alpha=-1.5 var=2.2", psdFilterTrials );

            for( int k = 0; k < psdFilterTrials; ++k )
            {
                for( size_t n = 0; n < noise.size(); ++n )
                    noise[n] = normVar;
                psdF( noise );
                avgRms += ( mx::math::vectorVariance( noise, 0.0 ) );
                trialLogger.update( k );
            }

            avgRms = sqrt( avgRms / psdFilterTrials );

            REQUIRE_THAT( avgRms, Catch::Matchers::WithinAbs( sqrt( 2.2 ), psdFilterTol * sqrt( 2.2 ) ) );
        }
    }
    SECTION( "a rank 2 psd" )
    {
        SECTION( "alpha=-2.5, dk Nyquist matched to array size, var=1" )
        {
            mx::sigproc::psdFilter<double, 2> psdF;

            Eigen::Array<double, -1, -1> k, psd;

            k.resize( 64, 64 );
            psd.resize( 64, 64 );

            mx::sigproc::frequencyGrid( k, 1. / 128. );
            for( int cc = 0; cc < psd.cols(); ++cc )
            {
                for( int rr = 0; rr < psd.rows(); ++rr )
                {
                    if( k( rr, cc ) == 0 )
                        psd( rr, cc ) = 0;
                    else
                        psd( rr, cc ) = pow( k( rr, cc ), -2.5 );
                }
            }

            double dk = k( 0, 1 ) - k( 0, 0 );

            mx::sigproc::normPSD( psd, k, 1.0 );

            int rv = psdF.psd( psd, dk, dk );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == psd.rows() );
            REQUIRE( psdF.cols() == psd.cols() );
            REQUIRE( psdF.planes() == 1 );

            Eigen::Array<double, -1, -1> noise( psdF.rows(), psdF.cols() );

            mx::math::normDistT<double> normVar;

            double avgRms = 0;
            TrialProgressLogger trialLogger( "rank2 alpha=-2.5 var=1", psdFilterTrials );

            for( int k = 0; k < psdFilterTrials; ++k )
            {
                for( int cc = 0; cc < psd.cols(); ++cc )
                {
                    for( int rr = 0; rr < psd.rows(); ++rr )
                    {
                        noise( rr, cc ) = normVar;
                    }
                }

                psdF( noise );
                avgRms += noise.square().sum(); //(mx::math::vectorVariance(noise,0.0));
                trialLogger.update( k );
            }

            avgRms = sqrt( avgRms / ( psd.rows() * psd.cols() ) / psdFilterTrials );

            REQUIRE_THAT( avgRms, Catch::Matchers::WithinAbs( 1.0, psdFilterTol ) );
        }
        SECTION( "alpha=-1.5, dk arb, var=2.2" )
        {
            mx::sigproc::psdFilter<double, 2> psdF;

            Eigen::Array<double, -1, -1> k, psd;

            k.resize( 64, 64 );
            psd.resize( 64, 64 );

            mx::sigproc::frequencyGrid( k, 1. / 302. );
            for( int cc = 0; cc < psd.cols(); ++cc )
            {
                for( int rr = 0; rr < psd.rows(); ++rr )
                {
                    if( k( rr, cc ) == 0 )
                        psd( rr, cc ) = 0;
                    else
                        psd( rr, cc ) = pow( k( rr, cc ), -1.5 );
                }
            }

            double dk = k( 0, 1 ) - k( 0, 0 );

            mx::sigproc::normPSD( psd, k, 2.2 );

            int rv = psdF.psd( psd, dk, dk );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == psd.rows() );
            REQUIRE( psdF.cols() == psd.cols() );
            REQUIRE( psdF.planes() == 1 );

            Eigen::Array<double, -1, -1> noise( psdF.rows(), psdF.cols() );

            mx::math::normDistT<double> normVar;

            double avgRms = 0;
            TrialProgressLogger trialLogger( "rank2 alpha=-1.5 var=2.2", psdFilterTrials );

            for( int k = 0; k < psdFilterTrials; ++k )
            {
                for( int cc = 0; cc < psd.cols(); ++cc )
                {
                    for( int rr = 0; rr < psd.rows(); ++rr )
                    {
                        noise( rr, cc ) = normVar;
                    }
                }

                psdF( noise );
                avgRms += noise.square().sum(); //(mx::math::vectorVariance(noise,0.0));
                trialLogger.update( k );
            }

            avgRms = sqrt( avgRms / ( psd.rows() * psd.cols() ) / psdFilterTrials );

            REQUIRE_THAT( avgRms, Catch::Matchers::WithinAbs( sqrt( 2.2 ), psdFilterTol * sqrt( 2.2 ) ) );
        }
    }
    SECTION( "a rank 3 psd" )
    {
        SECTION( "k-alpha=-2.5, f-alph=-2.5, dk Nyquist matched to array size, df Nyquist matched to array size, var=1" )
        {
            mx::sigproc::psdFilter<double, 3> psdF;

            Eigen::Array<double, -1, -1> k, psdk;
            std::vector<double> f, f2s, psd2s;

            mx::improc::eigenCube<double> psd;

            k.resize( 32, 32 );
            f.resize( 33 );

            mx::sigproc::frequencyGrid( k, 1. / 64. );
            psdk.resize( k.rows(), k.cols() );
            for( int cc = 0; cc < psdk.cols(); ++cc )
            {
                for( int rr = 0; rr < psdk.rows(); ++rr )
                {
                    if( k( rr, cc ) == 0 )
                        psdk( rr, cc ) = 0;
                    else
                        psdk( rr, cc ) = pow( k( rr, cc ), -2.5 );
                }
            }
            mx::sigproc::normPSD( psdk, k, 1.0 );

            for( size_t n = 0; n < f.size(); ++n )
                f[n] = n * 1.0 / 64.;
            mx::sigproc::augment1SidedPSDFreq( f2s, f );
            psd2s.resize( f2s.size() );
            for( size_t n = 0; n < psd2s.size(); ++n )
                psd2s[n] = pow( fabs( f2s[n] ), -2.5 );
            psd2s[0] = psd2s[1];

            psd.resize( k.rows(), k.cols(), f2s.size() );

            for( int cc = 0; cc < psd.cols(); ++cc )
            {
                for( int rr = 0; rr < psd.rows(); ++rr )
                {
                    if( k( rr, cc ) == 0 )
                        psd.pixel( rr, cc ).setZero();
                    else
                    {
                        double p = psdk( rr, cc );
                        mx::sigproc::normPSD( psd2s, f2s, p, -1e5, 1e5 );

                        for( int pp = 0; pp < psd.planes(); ++pp )
                            psd.image( pp )( rr, cc ) = psd2s[pp];
                    }
                }
            }

            double dk = k( 0, 1 ) - k( 0, 0 );
            double df = f[1] - f[0];

            int rv = psdF.psd( psd, dk, dk, df );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == psd.rows() );
            REQUIRE( psdF.cols() == psd.cols() );
            REQUIRE( psdF.planes() == psd.planes() );

            mx::improc::eigenCube<double> noise( psdF.rows(), psdF.cols(), psdF.planes() );

            mx::math::normDistT<double> normVar;

            double avgRms = 0;
            TrialProgressLogger trialLogger( "rank3 kalpha=-2.5 falpha=-2.5 var=1", psdFilterTrials );

            for( int k = 0; k < psdFilterTrials; ++k )
            {
                for( int pp = 0; pp < psd.planes(); ++pp )
                {
                    for( int cc = 0; cc < psd.cols(); ++cc )
                    {
                        for( int rr = 0; rr < psd.rows(); ++rr )
                        {
                            noise.image( pp )( rr, cc ) = normVar;
                        }
                    }
                }

                psdF( noise );
                for( int pp = 0; pp < noise.planes(); ++pp )
                    avgRms += noise.image( pp ).square().sum();
                trialLogger.update( k );
            }

            avgRms = sqrt( avgRms / ( psd.rows() * psd.cols() * psd.planes() ) / psdFilterTrials );

            REQUIRE_THAT( avgRms, Catch::Matchers::WithinAbs( 1.0, psdFilterTol ) );
        }
        SECTION( "k-alpha=-3.5, f-alph=-1.5, dk arb, df arb, var=2" )
        {
            mx::sigproc::psdFilter<double, 3> psdF;

            Eigen::Array<double, -1, -1> k, psdk;
            std::vector<double> f, f2s, psd2s;

            mx::improc::eigenCube<double> psd;

            k.resize( 32, 32 );
            f.resize( 33 );

            mx::sigproc::frequencyGrid( k, 1. / 640. );
            psdk.resize( k.rows(), k.cols() );
            for( int cc = 0; cc < psdk.cols(); ++cc )
            {
                for( int rr = 0; rr < psdk.rows(); ++rr )
                {
                    if( k( rr, cc ) == 0 )
                        psdk( rr, cc ) = 0;
                    else
                        psdk( rr, cc ) = pow( k( rr, cc ), -3.5 );
                }
            }
            mx::sigproc::normPSD( psdk, k, 2.0 );

            for( size_t n = 0; n < f.size(); ++n )
                f[n] = n * 1.0 / 78.;
            mx::sigproc::augment1SidedPSDFreq( f2s, f );
            psd2s.resize( f2s.size() );
            for( size_t n = 0; n < psd2s.size(); ++n )
                psd2s[n] = pow( fabs( f2s[n] ), -1.5 );
            psd2s[0] = psd2s[1];

            psd.resize( k.rows(), k.cols(), f2s.size() );

            for( int cc = 0; cc < psd.cols(); ++cc )
            {
                for( int rr = 0; rr < psd.rows(); ++rr )
                {
                    if( k( rr, cc ) == 0 )
                        psd.pixel( rr, cc ).setZero();
                    else
                    {
                        double p = psdk( rr, cc );
                        mx::sigproc::normPSD( psd2s, f2s, p, -1e5, 1e5 );

                        for( int pp = 0; pp < psd.planes(); ++pp )
                            psd.image( pp )( rr, cc ) = psd2s[pp];
                    }
                }
            }

            double dk = k( 0, 1 ) - k( 0, 0 );
            double df = f[1] - f[0];

            int rv = psdF.psd( psd, dk, dk, df );

            REQUIRE( rv == 0 );
            REQUIRE( psdF.rows() == psd.rows() );
            REQUIRE( psdF.cols() == psd.cols() );
            REQUIRE( psdF.planes() == psd.planes() );

            mx::improc::eigenCube<double> noise( psdF.rows(), psdF.cols(), psdF.planes() );

            mx::math::normDistT<double> normVar;

            double avgRms = 0;
            TrialProgressLogger trialLogger( "rank3 kalpha=-3.5 falpha=-1.5 var=2", psdFilterTrials );

            for( int k = 0; k < psdFilterTrials; ++k )
            {
                for( int pp = 0; pp < psd.planes(); ++pp )
                {
                    for( int cc = 0; cc < psd.cols(); ++cc )
                    {
                        for( int rr = 0; rr < psd.rows(); ++rr )
                        {
                            noise.image( pp )( rr, cc ) = normVar;
                        }
                    }
                }

                psdF( noise );
                for( int pp = 0; pp < noise.planes(); ++pp )
                    avgRms += noise.image( pp ).square().sum();
                trialLogger.update( k );
            }

            avgRms = sqrt( avgRms / ( psd.rows() * psd.cols() * psd.planes() ) / psdFilterTrials );

            REQUIRE_THAT( avgRms, Catch::Matchers::WithinAbs( sqrt( 2.0 ), psdFilterTol * sqrt( 2 ) ) );
        }
    }
}
