/** \file fftT.hpp
 * \brief The Fast Fourier Transform interface
 * \ingroup ft_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2015-2025 Jared R. Males (jaredmales@gmail.com)
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

#ifndef fftT_hpp
#define fftT_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include "../cuda/templateCufft.hpp"

#include "fftT.hpp"

#include "ftTypes.hpp"

namespace mx
{
namespace math
{
namespace ft
{

/// Fast Fourier Transforms using the \ref cuda "CUDA library"
/** The \ref cufft_templates "CUFFT Templates" type resolution system is used to allow the compiler
 * to access the right plan and types for the transforms based on inputT and outputT.
 *
 *
 * \tparam _inputT is the input type of the transform, can be either real or complex
 * \tparam _outputT is the output type of the transform, can be either real or complex
 * \tparam _rank is the rank of the transform. Limited to 2 for now.
 *
 * \ingroup fft
 */
template <typename _inputT, typename _outputT, size_t _rank>
class fftT<_inputT, _outputT, _rank, 1>
{

  public:
    typedef _inputT inputT;
    typedef _outputT outputT;

    static const size_t rank = _rank;

    typedef typename fftwTypeSpec<inputT, outputT>::realT realT;

    typedef typename fftwTypeSpec<inputT, outputT>::complexT complexT;

    typedef Eigen::Array<inputT, -1, -1> eigenArrayInT;

    typedef Eigen::Array<outputT, -1, -1> eigenArrayOutT;

    typedef typename cufftHandle planT;

  protected:
    dir m_direction{ dir::forward };       ///< Direction of this FFT, either dir::forward (default) or dir::backward
    int m_cufftDirection{ CUFFT_FORWARD }; ///< Direction to pass to CUFFT routine.  Kept synchronized with m_direction.

    int m_szX{ 0 };                        ///< Size of the x dimension
    int m_szY{ 0 };                        ///< Size of the y dimension
    int m_szZ{ 0 };                        ///< size of the z dimension

    planT m_plan{ nullptr };               ///< The cufft handle object.

  public:
    /// Default c'tor
    fftT();

    /// Constructor for rank 1 FFT.
    template <int crank = _rank>
    fftT( int nx,                  ///< [in] the desired size of the FFT
          dir ndir = dir::forward, ///< [in] [optional] direction of this FFT, either dir::forward (default) or
                                   ///< dir::backward
          bool inPlace = false,    /**< [in] [optional] whether or not this is an in-place transform.
                                                        Default is false, out-of-place.*/
          typename std::enable_if<crank == 1>::type * = 0 );

    /// Constructor for rank 2 FFT.
    template <int crank = _rank>
    fftT( int nx,                  ///< [in] the desired x size of the FFT
          int ny,                  ///< [in] the desired y size of the FFT
          dir ndir = dir::forward, /**< [in] [optional] direction of this FFT, either dir::forward
                                                        (default) or dir::backward */
          bool inPlace = false,    /**< [in] [optional] whether or not this is an in-place transform.
                                                        Default is false, out-of-place. */
          typename std::enable_if<crank == 2>::type * = 0 );

    /// Constructor for rank 3 FFT.
    template <int crank = _rank>
    fftT( int nx,                  ///< [in] the desired x size of the FFT
          int ny,                  ///< [in] the desired y size of the FFT
          int nz,                  ///< [in] the desired z size of the FFT
          dir ndir = dir::forward, /**< [in] [optional] direction of this FFT, either dir::forward (default)
                                                        or dir::backward*/
          bool inPlace = false,    /**< [in] [optional] whether or not this is an in-place
                                                        transform.  Default is false, out-of-place. */
          typename std::enable_if<crank == 3>::type * = 0 );

    /// Destructor
    ~fftT();

    /// Destroy (de-allocate) the plan
    void destroyPlan();

    /// Get the direction of this FFT
    /** The direction is either dir::forward or dir::backward.
     *
     * \returns the current value of m_direction.
     */
    ft::dir direction();

    /// Planning routine for rank 1 transforms.
    template <int crank = _rank>
    void plan( int nx,                      ///< [in] the desired size of the FFT
               ft::dir ndir = dir::forward, /**< [in] [optional] direction of this FFT, either dir::forward (default)
                                                              or dir::backward */
               bool inPlace = false,        /**< [in] [optional] whether or not this is an in-place transform.
                                                                 Default is false, out-of-place. */
               typename std::enable_if<crank == 1>::type * = 0 );

    /// Planning routine for rank 2 transforms.
    template <int crank = _rank>
    void plan( int nx,                      ///< [in] the desired x size of the FFT
               int ny,                      ///< [in] the desired y size of the FFT
               ft::dir ndir = dir::forward, /**< [in] [optional] direction of this FFT, either dir::forward (default)
                                                              or dir::backward */
               bool inPlace = false,        /**< [in] [optional] whether or not this is an in-place transform.
                                                                 Default is false, out-of-place. */
               typename std::enable_if<crank == 2>::type * = 0 );

    /// Planning routine for rank 3 transforms.
    template <int crank = _rank>
    void plan( int nx,                      ///< [in] the desired x size of the FFT
               int ny,                      ///< [in] the desired y size of the FFT
               int nz,                      ///< [in] the desired z size of the FFT
               ft::dir ndir = dir::forward, /**< [in] [optional] direction of this FFT, either dir::forward
                                                                (default) or dir::backward */
               bool inPlace = false,        /**< [in] [optional] whether or not this is an in-place transform.
                                                                 Default is false, out-of-place. */
               typename std::enable_if<crank == 3>::type * = 0 );

    /// Conduct the FFT
    void operator()( outputT *out, ///< [out] the output of the FFT, must be pre-allocated
                     inputT *in    ///< [in] the input to the FFT
    ) const;

    /// Conduct the MFT
    void operator()( eigenArrayOutT &out, /**< [out] the output of the DFT */
                     eigenArrayInT &in /**< [in] the input to the DFT */ ) const;
};

template <typename inputT, typename outputT, size_t rank>
fftT<inputT, outputT, rank, 0>::fftT()
{
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
fftT<inputT, outputT, rank, 0>::fftT( int nx, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 1>::type * )
{
    static_assert( crank == 2, "only rank 2 is currently supported for cuda fftT" );
    // m_direction = ndir;

    // plan( nx, ndir, inPlace );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
fftT<inputT, outputT, rank, 0>::fftT(
    int nx, int ny, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 2>::type * )
{
    plan( nx, ny, ndir, inPlace );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
fftT<inputT, outputT, rank, 0>::fftT(
    int nx, int ny, int nz, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 3>::type * )
{
    static_assert( crank == 2, "only rank 2 is currently supported for cuda fftT" );

    // m_direction = ndir;

    // plan( nx, ny, nz, ndir, inPlace );
}

template <typename inputT, typename outputT, size_t rank>
fftT<inputT, outputT, rank, 0>::~fftT()
{
    destroyPlan();
}

template <typename inputT, typename outputT, size_t rank>
void fftT<inputT, outputT, rank, 0>::destroyPlan()
{
    if( m_plan )
        fftw_destroy_plan<realT>( m_plan );

    m_plan = 0;

    m_szX = 0;
    m_szY = 0;
}

template <typename inputT, typename outputT, size_t rank>
ft::dir fftT<inputT, outputT, rank, 0>::direction()
{
    return m_direction;
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
void fftT<inputT, outputT, rank, 0>::plan( int nx,
                                           ft::dir ndir,
                                           bool inPlace,
                                           typename std::enable_if<crank == 1>::type * )
{
    static_assert( crank == 2, "only rank 2 is currently supported for cuda fftT" );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
void fftT<inputT, outputT, rank, 0>::plan(
    int nx, int ny, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 2>::type * )
{
    m_direction = ndir;

    if( m_direction == dir::backward )
    {
        m_cufftDirection = CUFFT_INVERSE;
    }
    else
    {
        m_cufftDirection = CUFFT_FORWARD;
    }

    mx::cuda::cufftPlan2d<inputT, outputT>( m_plan, nx, ny );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
void fftT<inputT, outputT, rank, 0>::plan(
    int nx, int ny, int nz, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 3>::type * )
{
    static_assert( crank == 2, "only rank 2 is currently supported for cuda fftT" );
}

template <typename inputT, typename outputT, size_t rank>
void fftT<inputT, outputT, rank, 0>::operator()( outputT *out, inputT *in ) const
{
    mx::cuda::cufftExec<inputT, outputT>( m_plan, in, out, m_cufftDirection );
}

template <typename inputT, typename outputT, size_t rank>
void fftT<inputT, outputT, rank, 0>::operator()( eigenArrayOutT &out, eigenArrayInT &in ) const
{
    mx::cuda::cufftExec<inputT, outputT>( m_plan, in.data(), out.data(), m_cufftDirection );
}

} // namespace ft
} // namespace math
} // namespace mx

#endif // fftT_hpp
