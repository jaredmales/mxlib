/** \file imagingUtils.cu
 * \author Jared R. Males
 * \brief Implementation of GPU versions of imaging utilities
 * \ingroup wfp_files
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

#include <cuda_runtime.h>
#include <complex>
#include "math/cuda/templateCuda.hpp"
namespace mx
{
namespace wfp
{

template <typename realT, typename complexT>
__global__ void extractIntensityImageAccum_kern(
    realT *im, int im_cols, int imX0, int imXsz, int imY0, int imYsz, complexT *wf, int wf_cols, int wfX0, int wfY0 )
{
    // clang-format off
    #ifdef __CUDACC__ // clang-format on


    realT *im_data;
    complexT *wf_data;

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for( int j = threadID; j < imXsz; j += numThreads )
    {
        im_data = &im[imX0 + ( imY0 + j ) * im_cols];
        wf_data = &wf[wfX0 + ( wfY0 + j ) * wf_cols];
        for( int i = 0; i < imYsz; ++i )
        {
            realT a = reinterpret_cast<realT *>(&wf_data[i])[0];
            realT b = reinterpret_cast<realT *>(&wf_data[i])[1];

            im_data[i] += a*a + b*b;
        }
    }

    // clang-format off
    #endif //__CUDACC__
    // clang-format on
}

template <typename realT, typename complexT>
cudaError_t extractIntensityImageAccum_gpu_impl(
    realT * im, int im_cols, int imX0, int imXsz, int imY0, int imYsz, complexT * wf, int wf_cols, int wfX0, int wfY0 )
{

    cudaError_t rv = cudaSuccess;

    // clang-format off
    #ifdef __CUDACC__ // clang-format on

    rv = cudaGetLastError();

    extractIntensityImageAccum_kern<realT, complexT>
        <<<( imXsz + 255 ) / 256, 256>>>( im, im_cols, imX0, imXsz, imY0, imYsz, wf, wf_cols, wfX0, wfY0 );

    rv = cudaGetLastError();

    // clang-format off
    #endif // clang-format on


    return rv;
}

template <typename realT, typename complexT>
cudaError_t extractIntensityImageAccum_gpu(realT * im, int im_cols, int imX0, int imXsz, int imY0, int imYsz, complexT* wf, int wf_cols, int wfX0, int wfY0 );

template <>
cudaError_t extractIntensityImageAccum_gpu<float, cuComplex>(float * im, int im_cols, int imX0, int imXsz, int imY0, int imYsz, cuComplex* wf, int wf_cols, int wfX0, int wfY0 )
{
    return extractIntensityImageAccum_gpu_impl(im, im_cols, imX0, imXsz, imY0, imYsz, wf, wf_cols, wfX0, wfY0 );
}

template <>
cudaError_t extractIntensityImageAccum_gpu<double, cuDoubleComplex>(double * im, int im_cols, int imX0, int imXsz, int imY0, int imYsz, cuDoubleComplex* wf, int wf_cols, int wfX0, int wfY0 )
{
    return extractIntensityImageAccum_gpu_impl(im, im_cols, imX0, imXsz, imY0, imYsz, wf, wf_cols, wfX0, wfY0 );
}

} // namespace wfp
} // namespace mx
