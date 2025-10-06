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

#include <fftw3.h>

#include "fftwTemplates.hpp"
#include "../../meta/trueFalseT.hpp"

#include "ftTypes.hpp"

namespace mx
{
namespace math
{
namespace ft
{

template <typename _inputT, typename _outputT, size_t _rank, int _cudaGPU = 0>
class fftT;

template <int rank>
std::vector<int> fftwDimVec( int szX, int szY, int szZ );

/// Fast Fourier Transforms using the \ref fftw "FFTW library"
/** The \ref fftw_templates "FFTW Templates" type resolution system is used to allow the compiler
 * to access the right plan and types for the transforms based on inputT and outputT.
 *
 * Planning is done with "measure" to achieve optimum performance at the cost of extra time during
 * the planning step.  This can be mitigated using "wisdom", managed by \ref fftwEnvironment.
 *
 * Calls to the FFTW plan functions are protected by `\#pragma omp critical` directives
 * unless MX_FFTW_NOOMP is defined prior to the first inclusion of this file.  This means that
 * you do not need to use the `fftw_make_planner_thread_safe` facility in FFTW.
 *
 * \note I recommend not using fftw_make_planner_thread_safe or specifying number of threads in
 *       \ref fftwEnvironment.
 *
 * \todo add execute interface with fftw like signature
 * \todo add plan interface where user passes in pointers (to avoid allocations)
 *
 * \tparam _inputT is the input type of the transform, can be either real or complex
 * \tparam _outputT is the output type of the transform, can be either real or complex
 * \tparam _rank is the rank of the transform. Limited to 1, 2, or 3.
 *
 * \ingroup fft
 */
template <typename _inputT, typename _outputT, size_t _rank>
class fftT<_inputT, _outputT, _rank, 0>
{

  public:
    typedef _inputT inputT;
    typedef _outputT outputT;

    static const size_t rank = _rank;

    typedef typename fftwTypeSpec<inputT, outputT>::realT realT;

    typedef typename fftwTypeSpec<inputT, outputT>::complexT complexT;

    typedef Eigen::Array<inputT, -1, -1> eigenArrayInT;

    typedef Eigen::Array<outputT, -1, -1> eigenArrayOutT;

    typedef typename fftwPlanSpec<realT>::planT planT;

  protected:
    dir m_direction{ dir::forward }; ///< Direction of this FFT, either dir::forward (default) or dir::backward

    int m_szX{ 0 };                  ///< Size of the x dimension
    int m_szY{ 0 };                  ///< Size of the y dimension
    int m_szZ{ 0 };                  ///< size of the z dimension

    planT m_plan{ nullptr };         ///< The FFTW plan object.  This is a pointer, allocated by FFTW library calls.

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

    /// Call the FFTW planning routine for an out-of-place transform.
    void doPlan( const meta::trueFalseT<false> &inPlace );

    /// Call the FFTW planning routine for an in-place transform.
    void doPlan( const meta::trueFalseT<true> &inPlace );

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
    m_direction = ndir;

    plan( nx, ndir, inPlace );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
fftT<inputT, outputT, rank, 0>::fftT(
    int nx, int ny, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 2>::type * )
{
    m_direction = ndir;

    plan( nx, ny, ndir, inPlace );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
fftT<inputT, outputT, rank, 0>::fftT(
    int nx, int ny, int nz, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 3>::type * )
{
    m_direction = ndir;

    plan( nx, ny, nz, ndir, inPlace );
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
void fftT<inputT, outputT, rank, 0>::doPlan( const meta::trueFalseT<false> &inPlace )
{
    (void)inPlace;

    inputT *forplan1;
    outputT *forplan2;

    int sz;

    if( rank == 1 )
        sz = m_szX;
    if( rank == 2 )
        sz = m_szX * m_szY;
    if( rank == 3 )
        sz = m_szX * m_szY * m_szZ;

    forplan1 = fftw_malloc<inputT>( sz );
    forplan2 = fftw_malloc<outputT>( sz );

    int pdir = FFTW_FORWARD;
    if( m_direction == dir::backward )
    {
        pdir = FFTW_BACKWARD;
    }

#ifndef MX_FFTW_NOOMP
    #pragma omp critical
    { // scope for pragma
#endif
        m_plan = fftw_plan_dft<inputT, outputT>( fftwDimVec<rank>( m_szX, m_szY, m_szZ ),
                                                 forplan1,
                                                 forplan2,
                                                 pdir,
                                                 FFTW_MEASURE );
#ifndef MX_FFTW_NOOMP
    }
#endif

    fftw_free<inputT>( forplan1 );
    fftw_free<outputT>( forplan2 );
}

template <typename inputT, typename outputT, size_t rank>
void fftT<inputT, outputT, rank, 0>::doPlan( const meta::trueFalseT<true> &inPlace )
{
    (void)inPlace;

    complexT *forplan;

    int sz;

    if( rank == 1 )
        sz = m_szX;
    if( rank == 2 )
        sz = m_szX * m_szY;
    if( rank == 3 )
        sz = m_szX * m_szY * m_szZ;

    forplan = fftw_malloc<complexT>( sz );

    int pdir = FFTW_FORWARD;
    if( m_direction == dir::backward )
    {
        pdir = FFTW_BACKWARD;
    }

#ifndef MX_FFTW_NOOMP
    #pragma omp critical
    { // scope for pragma
#endif
        m_plan = fftw_plan_dft<inputT, outputT>( fftwDimVec<rank>( m_szX, m_szY, m_szZ ),
                                                 reinterpret_cast<inputT *>( forplan ),
                                                 reinterpret_cast<outputT *>( forplan ),
                                                 pdir,
                                                 FFTW_MEASURE );
#ifndef MX_FFTW_NOOMP
    }
#endif

    fftw_free<inputT>( reinterpret_cast<inputT *>( forplan ) );
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
void fftT<inputT, outputT, rank, 0>::plan( int nx,
                                           ft::dir ndir,
                                           bool inPlace,
                                           typename std::enable_if<crank == 1>::type * )
{
    if( m_szX == nx && m_direction == ndir && m_plan )
    {
        return;
    }

    destroyPlan();

    m_direction = ndir;

    m_szX = nx;
    m_szY = 0;
    m_szZ = 0;

    if( inPlace == false )
    {
        doPlan( meta::trueFalseT<false>() );
    }
    else
    {
        doPlan( meta::trueFalseT<true>() );
    }
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
void fftT<inputT, outputT, rank, 0>::plan(
    int nx, int ny, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 2>::type * )
{
    if( m_szX == nx && m_szY == ny && m_direction == ndir && m_plan )
    {
        return;
    }

    destroyPlan();

    m_direction = ndir;

    m_szX = nx;
    m_szY = ny;
    m_szZ = 0;

    if( inPlace == false )
    {
        doPlan( meta::trueFalseT<false>() );
    }
    else
    {
        doPlan( meta::trueFalseT<true>() );
    }
}

template <typename inputT, typename outputT, size_t rank>
template <int crank>
void fftT<inputT, outputT, rank, 0>::plan(
    int nx, int ny, int nz, ft::dir ndir, bool inPlace, typename std::enable_if<crank == 3>::type * )
{
    if( m_szX == nx && m_szY == ny && m_szZ == nz && m_direction == ndir && m_plan )
    {
        return;
    }

    destroyPlan();

    m_direction = ndir;

    m_szX = nx;
    m_szY = ny;
    m_szZ = nz;

    if( inPlace == false )
    {
        doPlan( meta::trueFalseT<false>() );
    }
    else
    {
        doPlan( meta::trueFalseT<true>() );
    }
}

template <typename inputT, typename outputT, size_t rank>
void fftT<inputT, outputT, rank, 0>::operator()( outputT *out, inputT *in ) const
{
    fftw_execute_dft<inputT, outputT>( m_plan, in, out );
}

template <typename inputT, typename outputT, size_t rank>
void fftT<inputT, outputT, rank, 0>::operator()( eigenArrayOutT &out, eigenArrayInT &in ) const
{
    fftw_execute_dft<inputT, outputT>( m_plan, in.data(), out.data() );
}

template <>
std::vector<int> fftwDimVec<1>( int szX, int szY, int szZ );

template <>
std::vector<int> fftwDimVec<2>( int szX, int szY, int szZ );

template <>
std::vector<int> fftwDimVec<3>( int szX, int szY, int szZ );

#ifdef MXLIB_CUDA

#include "fftTcuda.hpp"

#endif

} // namespace ft
} // namespace math
} // namespace mx

#endif // fftT_hpp
