/** \file templateCublas_test.cpp
 */
#include "../../../catch2/catch.hpp"

#include "../../../../include/math/cuda/cudaPtr.hpp"
#include "../../../../include/math/cuda/templateCublas.hpp"

/** Scenario: scaling a vector with cublas
 * Tests cublasTscal, as well as basic cudaPtr operations.
 *
 * \anchor test_math_templateCublas_scal
 */
SCENARIO( "scaling a vector with cublas", "[math::cuda::templateCublas]" )
{
    GIVEN( "a vector" )
    {
        WHEN( "type is single precision real" )
        {
            std::vector<float> hx;
            mx::cuda::cudaPtr<float> dx;
            float alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = n;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTscal( handle, dx.size(), &alpha, dx(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dx.download( hx.data() );

            REQUIRE( hx[0] == 0 );
            REQUIRE( hx[1] == 2 );
            REQUIRE( hx[2] == 4 );
            REQUIRE( hx[3] == 6 );
            REQUIRE( hx[4] == 8 );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "type is double precision real" )
        {
            std::vector<double> hx;
            mx::cuda::cudaPtr<double> dx;
            double alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = n;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTscal( handle, dx.size(), &alpha, dx(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dx.download( hx.data() );

            REQUIRE( hx[0] == 0 );
            REQUIRE( hx[1] == 2 );
            REQUIRE( hx[2] == 4 );
            REQUIRE( hx[3] == 6 );
            REQUIRE( hx[4] == 8 );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "type is single precision complex" )
        {
            std::vector<std::complex<float>> hx;
            mx::cuda::cudaPtr<std::complex<float>> dx;
            std::complex<float> alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<float>( n, n );

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTscal( handle, dx.size(), (cuComplex *)&alpha, dx(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dx.download( hx.data() );

            REQUIRE( hx[0] == std::complex<float>( 0, 0 ) );
            REQUIRE( hx[1] == std::complex<float>( 2, 2 ) );
            REQUIRE( hx[2] == std::complex<float>( 4, 4 ) );
            REQUIRE( hx[3] == std::complex<float>( 6, 6 ) );
            REQUIRE( hx[4] == std::complex<float>( 8, 8 ) );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "type is double precision complex" )
        {
            std::vector<std::complex<double>> hx;
            mx::cuda::cudaPtr<std::complex<double>> dx;
            std::complex<double> alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<double>( n, n );

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTscal( handle, dx.size(), (cuDoubleComplex *)&alpha, dx(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dx.download( hx.data() );

            REQUIRE( hx[0] == std::complex<double>( 0, 0 ) );
            REQUIRE( hx[1] == std::complex<double>( 2, 2 ) );
            REQUIRE( hx[2] == std::complex<double>( 4, 4 ) );
            REQUIRE( hx[3] == std::complex<double>( 6, 6 ) );
            REQUIRE( hx[4] == std::complex<double>( 8, 8 ) );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }
    }
}

/** Scenario: scaling and accumulating a vector with cublas
 * Tests cublasTaxpy, as well as basic cudaPtr operations.
 *
 * \anchor test_math_templateCublas_axpy
 */
SCENARIO( "scaling and accumulating a vector with cublas", "[math::cuda::templateCublas]" )
{
    GIVEN( "a vector" )
    {
        WHEN( "type is single precision real" )
        {
            std::vector<float> hx, hy;
            mx::cuda::cudaPtr<float> dx, dy;
            float alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = n;

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = 1;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.upload( hy.data(), hy.size() );
            REQUIRE( dy.size() == hy.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTaxpy( handle, dx.size(), &alpha, dx(), 1, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 1 );
            REQUIRE( hy[1] == 3 );
            REQUIRE( hy[2] == 5 );
            REQUIRE( hy[3] == 7 );
            REQUIRE( hy[4] == 9 );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "type is double precision real" )
        {
            std::vector<double> hx, hy;
            mx::cuda::cudaPtr<double> dx, dy;
            double alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = n;

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = 1;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.upload( hy.data(), hy.size() );
            REQUIRE( dy.size() == hy.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTaxpy( handle, dx.size(), &alpha, dx(), 1, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 1 );
            REQUIRE( hy[1] == 3 );
            REQUIRE( hy[2] == 5 );
            REQUIRE( hy[3] == 7 );
            REQUIRE( hy[4] == 9 );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "type is single precision complex" )
        {
            std::vector<std::complex<float>> hx, hy;
            mx::cuda::cudaPtr<std::complex<float>> dx, dy;
            std::complex<float> alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<float>( n, n );

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = std::complex<float>( 1, 1 );

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.upload( hy.data(), hy.size() );
            REQUIRE( dy.size() == hy.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTaxpy( handle, dx.size(), (cuComplex *)&alpha, dx(), 1, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == std::complex<float>( 1, 1 ) );
            REQUIRE( hy[1] == std::complex<float>( 3, 3 ) );
            REQUIRE( hy[2] == std::complex<float>( 5, 5 ) );
            REQUIRE( hy[3] == std::complex<float>( 7, 7 ) );
            REQUIRE( hy[4] == std::complex<float>( 9, 9 ) );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "type is double precision complex" )
        {
            std::vector<std::complex<double>> hx, hy;
            mx::cuda::cudaPtr<std::complex<double>> dx, dy;
            std::complex<double> alpha = 2;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<double>( n, n );

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = std::complex<double>( 1, 1 );

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.upload( hy.data(), hy.size() );
            REQUIRE( dy.size() == hy.size() );

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTaxpy( handle, dx.size(), (cuDoubleComplex *)&alpha, dx(), 1, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == std::complex<double>( 1, 1 ) );
            REQUIRE( hy[1] == std::complex<double>( 3, 3 ) );
            REQUIRE( hy[2] == std::complex<double>( 5, 5 ) );
            REQUIRE( hy[3] == std::complex<double>( 7, 7 ) );
            REQUIRE( hy[4] == std::complex<double>( 9, 9 ) );

            stat = cublasDestroy( handle );

            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }
    }
}

/** Scenario: multiplying two vectors element by element
 * Tests mx::cuda::elementwiseXxY, as well as basic cudaPtr operations.
 *
 * \anchor test_math_templateCublas_elementwiseXxY
 */
SCENARIO( "multiplying two vectors element by element", "[math::cuda::templateCublas]" )
{
    GIVEN( "a vector" )
    {
        WHEN( "both types are single precision real" )
        {
            std::vector<float> hx, hy;
            mx::cuda::cudaPtr<float> dx, dy;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = n;

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = 2 * n;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.resize( hy.size() );
            dy.upload( hy.data() );
            REQUIRE( dy.size() == hy.size() );

            cudaError_t rv = mx::cuda::elementwiseXxY( dx(), dy(), dx.size() );
            REQUIRE( rv == cudaSuccess );

            dx.download( hx.data() );

            REQUIRE( hx[0] == 0 );
            REQUIRE( hx[1] == 2 );
            REQUIRE( hx[2] == 8 );
            REQUIRE( hx[3] == 18 );
            REQUIRE( hx[4] == 32 );
        }

        WHEN( "type1 is complex-float, and type2 is float" )
        {
            std::vector<std::complex<float>> hx;
            mx::cuda::cudaPtr<std::complex<float>> dx;

            std::vector<float> hy;
            mx::cuda::cudaPtr<float> dy;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<float>( n, n );

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = 2 * n;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.resize( hy.size() );
            dy.upload( hy.data() );
            REQUIRE( dy.size() == hy.size() );

            cudaError_t rv = mx::cuda::elementwiseXxY( dx(), dy(), dx.size() );
            REQUIRE( rv == cudaSuccess );

            dx.download( hx.data() );

            REQUIRE( hx[0] == std::complex<float>( 0, 0 ) );
            REQUIRE( hx[1] == std::complex<float>( 2, 2 ) );
            REQUIRE( hx[2] == std::complex<float>( 8, 8 ) );
            REQUIRE( hx[3] == std::complex<float>( 18, 18 ) );
            REQUIRE( hx[4] == std::complex<float>( 32, 32 ) );
        }

        WHEN( "type1 is complex-float, and type2 is complex-float" )
        {
            std::vector<std::complex<float>> hx;
            mx::cuda::cudaPtr<std::complex<float>> dx;

            std::vector<std::complex<float>> hy;
            mx::cuda::cudaPtr<std::complex<float>> dy;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<float>( n, n );

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = std::complex<float>( 0, 2 * n );

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.resize( hy.size() );
            dy.upload( hy.data() );
            REQUIRE( dy.size() == hy.size() );

            cudaError_t rv = mx::cuda::elementwiseXxY( dx(), dy(), dx.size() );
            REQUIRE( rv == cudaSuccess );

            dx.download( hx.data() );

            REQUIRE( hx[0] == std::complex<float>( 0, 0 ) );
            REQUIRE( hx[1] == std::complex<float>( -2, 2 ) ); //(1,1) * (0,2) = (0 + 2i + 0i -2) = (-2,2)
            REQUIRE( hx[2] == std::complex<float>( -8, 8 ) );
            REQUIRE( hx[3] == std::complex<float>( -18, 18 ) );
            REQUIRE( hx[4] == std::complex<float>( -32, 32 ) );
        }

        WHEN( "both types are double precision real" )
        {
            std::vector<double> hx, hy;
            mx::cuda::cudaPtr<double> dx, dy;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = n;

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = 2 * n;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.resize( hy.size() );
            dy.upload( hy.data() );
            REQUIRE( dy.size() == hy.size() );

            cudaError_t rv = mx::cuda::elementwiseXxY( dx(), dy(), dx.size() );
            REQUIRE( rv == cudaSuccess );

            dx.download( hx.data() );

            REQUIRE( hx[0] == 0 );
            REQUIRE( hx[1] == 2 );
            REQUIRE( hx[2] == 8 );
            REQUIRE( hx[3] == 18 );
            REQUIRE( hx[4] == 32 );
        }

        WHEN( "type1 is complex-double, and type2 is double" )
        {
            std::vector<std::complex<double>> hx;
            mx::cuda::cudaPtr<std::complex<double>> dx;

            std::vector<double> hy;
            mx::cuda::cudaPtr<double> dy;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<double>( n, n );

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = 2 * n;

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.resize( hy.size() );
            dy.upload( hy.data() );
            REQUIRE( dy.size() == hy.size() );

            cudaError_t rv = mx::cuda::elementwiseXxY( dx(), dy(), dx.size() );
            REQUIRE( rv == cudaSuccess );

            dx.download( hx.data() );

            REQUIRE( hx[0] == std::complex<double>( 0, 0 ) );
            REQUIRE( hx[1] == std::complex<double>( 2, 2 ) );
            REQUIRE( hx[2] == std::complex<double>( 8, 8 ) );
            REQUIRE( hx[3] == std::complex<double>( 18, 18 ) );
            REQUIRE( hx[4] == std::complex<double>( 32, 32 ) );
        }

        WHEN( "type1 is complex-double, and type2 is complex-double" )
        {
            std::vector<std::complex<double>> hx;
            mx::cuda::cudaPtr<std::complex<double>> dx;

            std::vector<std::complex<double>> hy;
            mx::cuda::cudaPtr<std::complex<double>> dy;

            hx.resize( 5 );
            for( size_t n = 0; n < hx.size(); ++n )
                hx[n] = std::complex<double>( n, n );

            hy.resize( 5 );
            for( size_t n = 0; n < hy.size(); ++n )
                hy[n] = std::complex<double>( 1, 2 * n );

            dx.upload( hx.data(), hx.size() );
            REQUIRE( dx.size() == hx.size() );

            dy.resize( hy.size() );
            dy.upload( hy.data() );
            REQUIRE( dy.size() == hy.size() );

            cudaError_t rv = mx::cuda::elementwiseXxY( dx(), dy(), dx.size() );
            REQUIRE( rv == cudaSuccess );

            dx.download( hx.data() );

            REQUIRE( hx[0] == std::complex<double>( 0, 0 ) );  //(0,0) * (0,0) = (0,0)
            REQUIRE( hx[1] == std::complex<double>( -1, 3 ) ); //(1,1) * (1,2) = (1 + 2i + i -2) = (-1,3)
            REQUIRE( hx[2] == std::complex<double>( -6, 10 ) );
            REQUIRE( hx[3] == std::complex<double>( -15, 21 ) );
            REQUIRE( hx[4] == std::complex<double>( -28, 36 ) );
        }
    }
}

/** Scenario: multiplying a vector by a matrix
 * Tests mx::cuda::cublasTgemv, as well as basic cudaPtr operations.
 *
 * \anchor test_math_templateCublas_cublasTgemv_inc
 */
SCENARIO( "multiplying a vector by a matrix giving increments", "[math::cuda::templateCublas]" )
{
    GIVEN( "a 2x2 matrix, float" )
    {
        WHEN( "float precision, beta is 0" )
        {
            std::vector<float> hA; // This will actually be a vector
            mx::cuda::cudaPtr<float> dA;

            std::vector<float> hx; // This will actually be a vector
            mx::cuda::cudaPtr<float> dx;

            std::vector<float> hy;
            mx::cuda::cudaPtr<float> dy;

            /* Column major order:
               1 3
               2 4
            */
            hA.resize( 4 );
            hA[0] = 1;
            hA[1] = 2;
            hA[2] = 3;
            hA[3] = 4;

            dA.resize( 2, 2 );
            dA.upload( hA.data() );

            hx.resize( 2 );
            hx[0] = 1;
            hx[1] = 2;

            dx.upload( hx.data(), hx.size() );

            hy.resize( 2 );

            dy.resize( 2 );
            dy.initialize();

            float alpha = 1;
            float beta = 0;

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTgemv( handle, CUBLAS_OP_N, 2, 2, &alpha, dA(), 2, dx(), 1, &beta, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 7 );
            REQUIRE( hy[1] == 10 );

            stat = cublasDestroy( handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "float precision, beta is 1, but y is all 0" )
        {
            std::vector<float> hA; // This will actually be a vector
            mx::cuda::cudaPtr<float> dA;

            std::vector<float> hx; // This will actually be a vector
            mx::cuda::cudaPtr<float> dx;

            std::vector<float> hy;
            mx::cuda::cudaPtr<float> dy;

            /* Column major order:
               1 3
               2 4
            */
            hA.resize( 4 );
            hA[0] = 1;
            hA[1] = 2;
            hA[2] = 3;
            hA[3] = 4;

            dA.resize( 2, 2 );
            dA.upload( hA.data() );

            hx.resize( 2 );
            hx[0] = 1;
            hx[1] = 2;

            dx.upload( hx.data(), hx.size() );

            hy.resize( 2 );

            dy.resize( 2 );
            dy.initialize();

            float alpha = 1;
            float beta = 1;

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTgemv( handle, CUBLAS_OP_N, 2, 2, &alpha, dA(), 2, dx(), 1, &beta, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 7 );
            REQUIRE( hy[1] == 10 );

            stat = cublasDestroy( handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }
        WHEN( "float precision, beta is 1, y is [1,2]" )
        {
            std::vector<float> hA; // This will actually be a vector
            mx::cuda::cudaPtr<float> dA;

            std::vector<float> hx; // This will actually be a vector
            mx::cuda::cudaPtr<float> dx;

            std::vector<float> hy;
            mx::cuda::cudaPtr<float> dy;

            /* Column major order:
               1 3
               2 4
            */
            hA.resize( 4 );
            hA[0] = 1;
            hA[1] = 2;
            hA[2] = 3;
            hA[3] = 4;

            dA.resize( 2, 2 );
            dA.upload( hA.data() );

            hx.resize( 2 );
            hx[0] = 1;
            hx[1] = 2;

            dx.upload( hx.data(), hx.size() );

            hy.resize( 2 );
            hy[0] = 1;
            hy[1] = 2;

            dy.resize( 2 );
            dy.upload( hx.data() );

            float alpha = 1;
            float beta = 1;

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTgemv( handle, CUBLAS_OP_N, 2, 2, &alpha, dA(), 2, dx(), 1, &beta, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 8 );
            REQUIRE( hy[1] == 12 );

            stat = cublasDestroy( handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }
    }

    GIVEN( "a 2x2 matrix, double" )
    {
        WHEN( "double precision, beta is 0" )
        {
            std::vector<double> hA; // This will actually be a vector
            mx::cuda::cudaPtr<double> dA;

            std::vector<double> hx; // This will actually be a vector
            mx::cuda::cudaPtr<double> dx;

            std::vector<double> hy;
            mx::cuda::cudaPtr<double> dy;

            /* Column major order:
               1 3
               2 4
            */
            hA.resize( 4 );
            hA[0] = 1;
            hA[1] = 2;
            hA[2] = 3;
            hA[3] = 4;

            dA.resize( 2, 2 );
            dA.upload( hA.data() );

            hx.resize( 2 );
            hx[0] = 1;
            hx[1] = 2;

            dx.upload( hx.data(), hx.size() );

            hy.resize( 2 );

            dy.resize( 2 );
            dy.initialize();

            double alpha = 1;
            double beta = 0;

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTgemv( handle, CUBLAS_OP_N, 2, 2, &alpha, dA(), 2, dx(), 1, &beta, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 7 );
            REQUIRE( hy[1] == 10 );

            stat = cublasDestroy( handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }

        WHEN( "double precision, beta is 1, but y is all 0" )
        {
            std::vector<double> hA; // This will actually be a vector
            mx::cuda::cudaPtr<double> dA;

            std::vector<double> hx; // This will actually be a vector
            mx::cuda::cudaPtr<double> dx;

            std::vector<double> hy;
            mx::cuda::cudaPtr<double> dy;

            /* Column major order:
               1 3
               2 4
            */
            hA.resize( 4 );
            hA[0] = 1;
            hA[1] = 2;
            hA[2] = 3;
            hA[3] = 4;

            dA.resize( 2, 2 );
            dA.upload( hA.data() );

            hx.resize( 2 );
            hx[0] = 1;
            hx[1] = 2;

            dx.upload( hx.data(), hx.size() );

            hy.resize( 2 );

            dy.resize( 2 );
            dy.initialize();

            double alpha = 1;
            double beta = 1;

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTgemv( handle, CUBLAS_OP_N, 2, 2, &alpha, dA(), 2, dx(), 1, &beta, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 7 );
            REQUIRE( hy[1] == 10 );

            stat = cublasDestroy( handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }
        WHEN( "double precision, beta is 1, y is [1,2]" )
        {
            std::vector<double> hA; // This will actually be a vector
            mx::cuda::cudaPtr<double> dA;

            std::vector<double> hx; // This will actually be a vector
            mx::cuda::cudaPtr<double> dx;

            std::vector<double> hy;
            mx::cuda::cudaPtr<double> dy;

            /* Column major order:
               1 3
               2 4
            */
            hA.resize( 4 );
            hA[0] = 1;
            hA[1] = 2;
            hA[2] = 3;
            hA[3] = 4;

            dA.resize( 2, 2 );
            dA.upload( hA.data() );

            hx.resize( 2 );
            hx[0] = 1;
            hx[1] = 2;

            dx.upload( hx.data(), hx.size() );

            hy.resize( 2 );
            hy[0] = 1;
            hy[1] = 2;

            dy.resize( 2 );
            dy.upload( hx.data() );

            double alpha = 1;
            double beta = 1;

            cublasHandle_t handle;
            cublasStatus_t stat;
            stat = cublasCreate( &handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            stat = mx::cuda::cublasTgemv( handle, CUBLAS_OP_N, 2, 2, &alpha, dA(), 2, dx(), 1, &beta, dy(), 1 );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );

            dy.download( hy.data() );

            REQUIRE( hy[0] == 8 );
            REQUIRE( hy[1] == 12 );

            stat = cublasDestroy( handle );
            REQUIRE( stat == CUBLAS_STATUS_SUCCESS );
        }
    }
}
