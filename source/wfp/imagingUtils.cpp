/** \file imagingUtils.cpp
 * \author Jared R. Males
 * \brief Implementation of imaging utilities
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

#include "wfp/imagingUtils.hpp"

#include "improc/eigenImage.hpp"
using namespace mx::improc;

#include "math/cuda/cudaPtr.hpp"
using namespace mx::cuda;

//#include <iostream>
#include <Eigen/Dense>


namespace mx
{
namespace wfp
{


// CPU version
template <typename realImageT, typename complexImageT>
void extractIntensityImageAccum_cpu(
    realImageT &im, int imX0, int imXsz, int imY0, int imYsz, complexImageT &wf, int wfX0, int wfY0 )
{
    int im_rows = im.cols();

    int wf_rows = wf.cols();

    typename realImageT::Scalar *im_data;
    typename complexImageT::Scalar *wf_data;

    for( int j = 0; j < imXsz; ++j )
    {
        im_data = &im.data()[imX0 + ( imY0 + j ) * im_rows];
        wf_data = &wf.data()[wfX0 + ( wfY0 + j ) * wf_rows];
        for( int i = 0; i < imYsz; ++i )
        {
            im_data[i] += norm( wf_data[i] );
        }
    }
}

template <>
void extractIntensityImageAccum<eigenImage<float>, eigenImage<std::complex<float>>, 0>(
    eigenImage<float> &im, int imX0, int imXsz, int imY0, int imYsz, eigenImage<std::complex<float>> &wf, int wfX0, int wfY0 )
{
    return extractIntensityImageAccum_cpu(im, imX0, imXsz, imY0, imYsz, wf, wfX0, wfY0);
}

template <>
void extractIntensityImageAccum<eigenImage<double>, eigenImage<std::complex<double>>, 0>(
    eigenImage<double> &im, int imX0, int imXsz, int imY0, int imYsz, eigenImage<std::complex<double>> &wf, int wfX0, int wfY0 )
{
    return extractIntensityImageAccum_cpu(im, imX0, imXsz, imY0, imYsz, wf, wfX0, wfY0);
}

template <typename realT, typename complexT>
cudaError_t extractIntensityImageAccum_gpu(
    realT *im, int im_cols, int imX0, int imXsz, int imY0, int imYsz, complexT *wf, int wf_cols, int wfX0, int wfY0 );


template <>
void extractIntensityImageAccum<cudaPtr<float>, cudaPtr<std::complex<float>>, 1>(
    cudaPtr<float> &im,
    int imX0,
    int imXsz,
    int imY0,
    int imYsz,
    cudaPtr<std::complex<float>> & wf,
    int wfX0,
    int wfY0 )
{
    extractIntensityImageAccum_gpu( im.data(),im.cols(), imX0, imXsz, imY0, imYsz, wf.data(), wf.cols(), wfX0, wfY0 );
}


template <>
void extractIntensityImageAccum<cudaPtr<double>, cudaPtr<std::complex<double>>, 1>(
    cudaPtr<double> &im,
    int imX0,
    int imXsz,
    int imY0,
    int imYsz,
    cudaPtr<std::complex<double>> &wf,
    int wfX0,
    int wfY0 )
{
    extractIntensityImageAccum_gpu( im.data(), im.cols(), imX0, imXsz, imY0, imYsz, wf.data(), wf.cols(), wfX0, wfY0 );
}


} // namespace wfp
} // namespace mx
