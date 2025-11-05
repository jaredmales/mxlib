/** \file syevdxT.hpp
 * \author Jared R. Males
 * \brief An interface to cuSOLVER Xsyevdx
 * \ingroup cuda_files
 *
 */

//***********************************************************************//
// Copyright 2025 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef math_syevdxT_hpp
#define math_syevdxT_hpp

#ifdef MXLIB_CUDA

#include <Eigen/Dense>

#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <cuda_runtime.h>
#include <cusolverDn.h>

// #include <cuComplex.h>
// #include <cuda_runtime_api.h>
// #include <cublas_api.h>
// #include <library_types.h>

// #include "/home/jrmales/Source/CUDALibrarySamples/cuSOLVER/utils/cusolver_utils.h"


namespace mx
{
namespace cuda
{

#define MXCUDA_EXCEPTION( ec, explan )                                                                                 \
    std::string msg = "Cuda Error: (";                                                                                 \
    msg += cudaGetErrorName( ec );                                                                                     \
    msg += ") ";                                                                                                       \
    msg += cudaGetErrorString( ec );                                                                                   \
    msg += "\nContext: " explan;                                                                                       \
    msg += ".\nAt line " + std::to_string( __LINE__ );                                                                 \
    msg += " in " __FILE__;                                                                                            \
    throw std::runtime_error( msg );

#define MXCUSOLVER_EXCEPTION( ec, explan )                                                                             \
    std::string msg = "cusolver Error: (";                                                                             \
    msg += std::to_string( static_cast<int>( ec ) );                                                                   \
    msg += ") ";                                                                                                       \
    msg += "\nContext: " explan;                                                                                       \
    msg += ".\nAt line " + std::to_string( __LINE__ );                                                                 \
    msg += " in " __FILE__;                                                                                            \
    throw std::runtime_error( msg );

template <typename floatT>
struct cusolver_traits;

template <>
struct cusolver_traits<float>
{
    static constexpr cudaDataType cuda_data_type = CUDA_R_32F;
};

template <>
struct cusolver_traits<double>
{
    static constexpr cudaDataType cuda_data_type = CUDA_R_64F;
};

template <typename floatT>
struct syevdxT
{

    cusolverDnHandle_t m_cusolverH{ nullptr };
    cusolverDnParams_t m_params{ nullptr };
    cudaStream_t m_stream{ nullptr };

    int m_allocations{ 0 };

    floatT *m_dev_A{ nullptr };
    floatT *m_dev_W{ nullptr };
    int *m_dev_info{ nullptr };

    cusolverEigRange_t m_range;
    cublasFillMode_t m_uplo;
    int64_t m_n{ 0 };
    int64_t m_lda{ 0 };

    floatT m_vu{ 0 };
    floatT m_vl{ 0 };
    int64_t m_il;
    int64_t m_iu;

    void *m_dev_work{ nullptr };  /* device workspace */
    size_t m_dev_wsBytes{ 0 };
    void *m_host_work{ nullptr }; /* device workspace */
    size_t m_host_wsBytes{ 0 };

    syevdxT();

    syevdxT( cusolverDnHandle_t cusolverH, cusolverDnParams_t params, cudaStream_t stream = nullptr );

    ~syevdxT();

    void free( bool null = true );

    void setup( cusolverDnHandle_t cusolverH, cusolverDnParams_t params, cudaStream_t stream = nullptr );

    void allocate( int64_t n, int64_t nVecs, cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER );

    int execute( int64_t &numEig, floatT *A );

    template <typename eigenT>
    int calcEigenVecs(eigenT &evecs,          /**< [out] on exit contains the eigen vectors*/
                      eigenT &evals,          /**< [out] on exit contains the eigen vectors*/
                      eigenT &cv,             /**< [in] a lower-triangle (in the Lapack sense) square
                                                             covariance matrix.*/
                      bool normalize = false, /**< [in] [opt] flag specifying whether or not to
                                                                        normalize the eigenvectors.*/
                      bool check = false      /**< [in] [opt] flag specifying whether or not to
                                                                         check the eigenvalues/vectors for
                                                                         validity. Requires normalize=true.*/ );
};

template <typename floatT>
syevdxT<floatT>::syevdxT()
{
}

template <typename floatT>
syevdxT<floatT>::syevdxT( cusolverDnHandle_t cusolverH, cusolverDnParams_t params, cudaStream_t stream )
    : m_cusolverH{ cusolverH }, m_params{ params }, m_stream{ stream }
{
}

template <typename floatT>
syevdxT<floatT>::~syevdxT()
{
    free( false );
}

template <typename floatT>
void syevdxT<floatT>::free( bool null )
{
    if( m_dev_A != nullptr )
    {
        cudaError_t ec = cudaFree( m_dev_A );
        if( ec != cudaSuccess )
        {
            MXCUDA_EXCEPTION( ec, "syevdxT::free: free-ing device memory for matrix" );
        }
    }

    if( m_dev_W != nullptr )
    {
        cudaError_t ec = cudaFree( m_dev_W );
        if( ec != cudaSuccess )
        {
            MXCUDA_EXCEPTION( ec, "syevdxT::free: free-ing device memory for eigen values" );
        }
    }

    if( m_dev_info != nullptr )
    {
        cudaError_t ec = cudaFree( m_dev_info );
        if( ec != cudaSuccess )
        {
            MXCUDA_EXCEPTION( ec, "syevdxT::free: free-ing device memory for info" );
        }
    }

    if( m_dev_work != nullptr )
    {
        cudaError_t ec = cudaFree( m_dev_work );
        if( ec != cudaSuccess )
        {
            MXCUDA_EXCEPTION( ec, "syevdxT::free: free-ing device memory for work" );
        }
    }

    if( m_host_work != nullptr )
    {
        ::free( m_host_work );
    }

    if( null )
    {
        m_dev_A = nullptr;
        m_dev_W = nullptr;
        m_dev_info = nullptr;
        m_dev_work = nullptr;
        m_host_work = nullptr;
    }
}

template <typename floatT>
void syevdxT<floatT>::setup( cusolverDnHandle_t cusolverH, cusolverDnParams_t params, cudaStream_t stream )
{
    m_cusolverH = cusolverH;
    m_params = params;
    m_stream = stream;
}

template <typename floatT>
void syevdxT<floatT>::allocate( int64_t n, int64_t nVecs, cublasFillMode_t uplo )

{
    if( m_cusolverH == nullptr )
    {
        throw std::runtime_error( "syevdx::allocate: m_cusolverH is null" );
    }

    if( m_params == nullptr )
    {
        throw std::runtime_error( "syevdx::allocate: m_params is null" );
    }

    cusolverEigRange_t range;
    int64_t il, iu;
    if( nVecs == 0 || nVecs == n )
    {
        range = CUSOLVER_EIG_RANGE_ALL;
        il = 0;
        iu = 0;
    }
    else
    {
        range = CUSOLVER_EIG_RANGE_I;
        if( nVecs > n )
        {
            nVecs = n;
        }
        il = n - nVecs + 1;
        iu = n;
    }

    if( n == m_n && range == m_range && il == m_il && iu == m_iu && uplo == m_uplo && m_dev_A && m_dev_W &&
        m_dev_info && m_dev_work )
    {
        return;
    }

    m_n = n;
    m_lda = n;

    m_range = range;
    m_il = il;
    m_iu = iu;

    m_uplo = uplo;

    free( false );

    cudaError_t ec = cudaMalloc( reinterpret_cast<void **>( &m_dev_A ), sizeof( floatT ) * ( n * n ) );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::allocate: allocating device memory for matrix" );
    }

    ec = cudaMalloc( reinterpret_cast<void **>( &m_dev_W ), sizeof( floatT ) * n );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::allocate: allocating device memory for eigen values" );
    }

    ec = cudaMalloc( reinterpret_cast<void **>( &m_dev_info ), sizeof( int ) );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::allocate: allocating device memory for info" );
    }

    int64_t h_meig;

    cusolverStatus_t csec = cusolverDnXsyevdx_bufferSize( m_cusolverH,
                                                          m_params,
                                                          CUSOLVER_EIG_MODE_VECTOR,
                                                          m_range,
                                                          m_uplo,
                                                          m_n,
                                                          cusolver_traits<floatT>::cuda_data_type,
                                                          m_dev_A,
                                                          m_lda,
                                                          &m_vl,
                                                          &m_vu,
                                                          m_il,
                                                          m_iu,
                                                          &h_meig,
                                                          cusolver_traits<floatT>::cuda_data_type,
                                                          m_dev_W,
                                                          cusolver_traits<floatT>::cuda_data_type,
                                                          &m_dev_wsBytes,
                                                          &m_host_wsBytes );
    if( csec != CUSOLVER_STATUS_SUCCESS )
    {
        MXCUSOLVER_EXCEPTION( csec, "syevdxT::allocate: call to cusolverDnXsyevdx_bufferSize" );
    }

    ec = cudaMalloc( reinterpret_cast<void **>( &m_dev_work ), m_dev_wsBytes );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::allocate: allocating device memory for work" );
    }

    if( m_host_wsBytes > 0 )
    {
        m_host_work = reinterpret_cast<void *>( malloc( m_host_wsBytes ) );
        if( m_host_work == nullptr )
        {
            mxThrowException( err::allocerr, "syevdxT::allocate:", "allocating host memory for work" );
        }
    }
    else
    {
        m_host_work = nullptr; // not set by free
    }

    ++m_allocations;
}

template <typename floatT>
int syevdxT<floatT>::execute( int64_t &numEig, floatT *A )
{
    if( m_cusolverH == nullptr )
    {
        throw std::runtime_error( "syevdx::execute: m_cusolverH is null" );
    }

    if( m_params == nullptr )
    {
        throw std::runtime_error( "syevdx::execute: m_params is null" );
    }

    cudaError_t ec = cudaMemcpyAsync( m_dev_A, A, sizeof( floatT ) * m_n * m_n, cudaMemcpyHostToDevice, m_stream );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::execute: copying matrix to device" );
    }

    cusolverStatus_t csec = cusolverDnXsyevdx( m_cusolverH,
                                               m_params,
                                               CUSOLVER_EIG_MODE_VECTOR,
                                               m_range,
                                               m_uplo,
                                               m_n,
                                               cusolver_traits<floatT>::cuda_data_type,
                                               m_dev_A,
                                               m_lda,
                                               &m_vl,
                                               &m_vu,
                                               m_il,
                                               m_iu,
                                               &numEig,
                                               cusolver_traits<floatT>::cuda_data_type,
                                               m_dev_W,
                                               cusolver_traits<floatT>::cuda_data_type,
                                               m_dev_work,
                                               m_dev_wsBytes,
                                               m_host_work,
                                               m_host_wsBytes,
                                               m_dev_info );

    if( csec != CUSOLVER_STATUS_SUCCESS )
    {
        MXCUSOLVER_EXCEPTION( csec, "syevdxT::execute: call to cusolverDnXsyevdx" );
    }

    int info = 0;
    ec = cudaMemcpyAsync( &info, m_dev_info, sizeof( int ), cudaMemcpyDeviceToHost, m_stream );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::execute: copying info from device" );
    }

    ec = cudaStreamSynchronize( m_stream );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::execute: synchronizing" );
    }

    return info;
}

template <typename floatT>
template <typename eigenT>
int syevdxT<floatT>::calcEigenVecs( eigenT &evecs, eigenT &evals, eigenT &cv, bool normalize, bool check )
{
    if( m_cusolverH == nullptr )
    {
        throw std::runtime_error( "syevdx::calcEigenVecs: m_cusolverH is null" );
    }

    if( m_params == nullptr )
    {
        throw std::runtime_error( "syevdx::calcEigenVecs: m_params is null" );
    }

    int64_t nVecs;

    int info = execute( nVecs, cv.data() );

    evecs.resize( m_n, nVecs );
    cudaError_t ec =
        cudaMemcpyAsync( evecs.data(), m_dev_A, sizeof( floatT ) * m_n * nVecs, cudaMemcpyDeviceToHost, m_stream );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::calcEigenVecs: copying eigenvectors from device" );
    }
    evals.resize( 1, nVecs );
    ec = cudaMemcpyAsync( evals.data(), m_dev_W, sizeof( floatT ) * nVecs, cudaMemcpyDeviceToHost, m_stream );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::calcEigenVecs: copying eigenvalues from device" );
    }

    ec = cudaStreamSynchronize( m_stream );
    if( ec != cudaSuccess )
    {
        MXCUDA_EXCEPTION( ec, "syevdxT::calcEigenVecs: synchronizing" );
    }

    if( normalize )
    {
        // Normalize the eigenvectors
        if( !check )
        {
            for( int i = 0; i < nVecs; ++i )
            {
                evecs.col( i ) = evecs.col( i ) / sqrt( evals( i ) );
            }
        }
        else // here we check for invalid results and 0 things out
        {
            for( int i = 0; i < nVecs; ++i )
            {
                if( evals( i ) == 0 )
                {
                    std::cerr << "got 0 eigenvalue (# " << i << ")\n";
                    evecs.col( i ) *= 0;
                }
                else if( evals( i ) < 0 )
                {
                    std::cerr << "got < 0 eigenvalue (# " << i << ")\n";
                    evecs.col( i ) *= 0;
                }
                else if( !std::isfinite( evals( i ) ) )
                {
                    std::cerr << "got not-normal eigenvalue (# " << i << ")\n";
                    evecs.col( i ) *= 0;
                }
                else
                {
                    evecs.col( i ) = evecs.col( i ) / sqrt( evals( i ) );
                }

                for( int r = 0; r < evecs.rows(); ++r )
                {
                    if( !std::isfinite( evecs.col( i )( r ) ) )
                    {
                        std::cerr << "got not-normal eigenvector entry (# " << i << "," << r
                                  << ") = " << evecs.col( i )( r ) << "\n";
                        evecs.col( i ) *= 0;
                        continue;
                    }
                }
            }
        }
    }

    return info;
}

} // namespace cuda
} // namespace mx

#endif //#MXLIB_CUDA

#endif // math_syevdxT_hpp
