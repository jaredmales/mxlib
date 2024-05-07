/** \file fft.hpp
  * \brief The fast Fourier transform interface
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef fft_hpp
#define fft_hpp


#include <fftw3.h>

#include <complex>

#include "fftwTemplates.hpp"
#include "../../meta/trueFalseT.hpp"

namespace mx
{
namespace math
{
namespace fft 
{
   
#define MXFFT_FORWARD  (FFTW_FORWARD)
#define MXFFT_BACKWARD (FFTW_BACKWARD)



template<typename _inputT, typename _outputT, size_t _rank, int _cudaGPU=0> 
class fftT;

template<int rank>
std::vector<int> fftwDimVec( int szX, int szY, int szZ);

/// Fast Fourier Transforms using the FFTW library interface
/** The fftwTemplates type resolution system is used to allow the compiler
  * to access the right plan and types for the transforms based on inputT and outputT.
  * 
  * Calls the FFTW plan functions are protected by '\#pragma omp critical' directives
  * unless MX_FFTW_NOOMP is define prior to the first inclusion of this file.
  * 
  * \todo add execute interface with fftw like signature
  * \todo add plan interface where user passes in pointers (to avoid allocations)
  * 
  * \tparam inputT is the input type of the transform, can be either real or complex
  * \tparam outputT is the output type of the transform, can be either real or complex
  * 
  * \ingroup fft
  */ 
template<typename _inputT, typename _outputT, size_t _rank>
class fftT< _inputT, _outputT, _rank, 0>
{  
   typedef _inputT inputT;
   typedef _outputT outputT;
   
   typedef typename fftwTypeSpec<inputT, outputT>::complexT complexT;
   
   static const size_t rank = _rank;
   
   typedef typename fftwTypeSpec<inputT, outputT>::realT realT;
   
   typedef typename fftwPlanSpec<realT>::planT planT;
   
   
protected:
   int m_dir {MXFFT_FORWARD}; ///< Direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
   
   int m_szX {0}; ///< Size of the x dimension
   int m_szY {0}; ///< Size of the y dimension
   int m_szZ {0}; ///< size of the z dimension
   
   planT m_plan {nullptr}; ///< The FFTW plan object.  This is a pointer, allocated by FFTW library calls.
   
public:
   
   /// Default c'tor
   fftT();

   ///Constructor for rank 1 FFT.
   template<int crank = _rank>
   fftT( int nx,                   ///< [in] the desired size of the FFT
         int ndir = MXFFT_FORWARD, ///< [in] [optional] direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
         bool inPlace = false,     ///< [in] [optional] whether or not this is an in-place transform.  Default is false, out-of-place.
         typename std::enable_if<crank==1>::type* = 0 
       );

   ///Constructor for rank 2 FFT.
   template<int crank = _rank>
   fftT( int nx,                   ///< [in] the desired x size of the FFT
         int ny,                   ///< [in] the desired y size of the FFT
         int ndir = MXFFT_FORWARD, ///< [in] [optional] direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
         bool inPlace = false,     ///< [in] [optional] whether or not this is an in-place transform.  Default is false, out-of-place.
         typename std::enable_if<crank==2>::type* = 0 );
   
   ///Constructor for rank 3 FFT.
   template<int crank = _rank>
   fftT( int nx,                   ///< [in] the desired x size of the FFT
         int ny,                   ///< [in] the desired y size of the FFT
         int nz,                   ///< [in] the desired z size of the FFT
         int ndir = MXFFT_FORWARD, ///< [in] [optional] direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
         bool inPlace = false,     ///< [in] [optional] whether or not this is an in-place transform.  Default is false, out-of-place.
         typename std::enable_if<crank==3>::type* = 0 );
   
   ///Destructor
   ~fftT();
   
   ///Destroy (de-allocate) the plan
   void destroyPlan();

   ///Get the direction of this FFT 
   /** The direction is either MXFFT_FORWARD or MXFFT_BACKWARD.
     * 
     * \returns the current value of m_dir.
     */ 
   int dir();
       
   /// Call the FFTW planning routine for an out-of-place transform.
   void doPlan(const meta::trueFalseT<false> & inPlace);
   
   /// Call the FFTW planning routine for an in-place transform.
   void doPlan(const meta::trueFalseT<true> & inPlace);
   
   /// Planning routine for rank 1 transforms.
   template<int crank = _rank>
   void plan( int nx,                 ///< [in] the desired size of the FFT
              int ndir=MXFFT_FORWARD, ///< [in] [optional] direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
              bool inPlace=false,     ///< [in] [optional] whether or not this is an in-place transform.  Default is false, out-of-place.
              typename std::enable_if<crank==1>::type* = 0 
            );
   
   /// Planning routine for rank 2 transforms.
   template<int crank = _rank>
   void plan( int nx,                 ///< [in] the desired x size of the FFT
              int ny,                 ///< [in] the desired y size of the FFT
              int ndir=MXFFT_FORWARD, ///< [in] [optional] direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
              bool inPlace=false,     ///< [in] [optional] whether or not this is an in-place transform.  Default is false, out-of-place.
              typename std::enable_if<crank==2>::type* = 0 
            );

   /// Planning routine for rank 3 transforms.
   template<int crank = _rank>
   void plan( int nx,                 ///< [in] the desired x size of the FFT
              int ny,                 ///< [in] the desired y size of the FFT
              int nz,                 ///< [in] the desired z size of the FFT
              int ndir=MXFFT_FORWARD, ///< [in] [optional] direction of this FFT, either MXFFT_FORWARD (default) or MXFFT_BACKWARD
              bool inPlace=false,     ///< [in] [optional] whether or not this is an in-place transform.  Default is false, out-of-place.
              typename std::enable_if<crank==3>::type* = 0 
            );
   
   /// Conduct the FFT
   void operator()( outputT *out, ///< [out] the output of the FFT, must be pre-allocated
                    inputT * in   ///< [in] the input to the FFT
                  ) const;
   
};

template<typename inputT, typename outputT, size_t rank>
fftT<inputT,outputT,rank,0>::fftT()
{
}

template<typename inputT, typename outputT, size_t rank>
template<int crank>
fftT<inputT,outputT,rank,0>::fftT( int nx, 
                                   int ndir,
                                   bool inPlace,
                                   typename std::enable_if<crank==1>::type* 
                                 )
{
   m_dir = ndir;
   
   plan(nx, ndir, inPlace);
}
   
template<typename inputT, typename outputT, size_t rank>
template<int crank>
fftT<inputT,outputT,rank,0>::fftT( int nx, 
                                   int ny, 
                                   int ndir,
                                   bool inPlace,
                                   typename std::enable_if<crank==2>::type* 
                                 )
{
   m_dir = ndir;
   
   plan(nx, ny, ndir, inPlace);
}

template<typename inputT, typename outputT, size_t rank>
template<int crank>
fftT<inputT,outputT,rank,0>::fftT( int nx, 
                                   int ny,
                                   int nz,
                                   int ndir,
                                   bool inPlace,
                                   typename std::enable_if<crank==3>::type* 
                                 )
{
   m_dir = ndir;
   
   plan(nx, ny, nz, ndir, inPlace);
}

template<typename inputT, typename outputT, size_t rank>
fftT<inputT,outputT,rank,0>::~fftT()
{
   destroyPlan();
}

template<typename inputT, typename outputT, size_t rank>
void fftT<inputT,outputT,rank,0>::destroyPlan()
{
   if(m_plan) fftw_destroy_plan<realT>(m_plan);
   
   m_plan = 0;
   
   m_szX = 0;
   m_szY = 0;

}

template<typename inputT, typename outputT, size_t rank>
int fftT<inputT,outputT,rank,0>::dir()
{
   return m_dir;
}
 
template<typename inputT, typename outputT, size_t rank> 
void fftT<inputT,outputT,rank,0>::doPlan(const meta::trueFalseT<false> & inPlace)
{
   (void) inPlace;
   
   inputT * forplan1;
   outputT * forplan2;
   
   int sz;
   
   if(rank == 1) sz = m_szX;
   if(rank == 2) sz = m_szX*m_szY;
   if(rank == 3) sz = m_szX*m_szY*m_szZ;
   
   forplan1 = fftw_malloc<inputT>(sz);
   forplan2 = fftw_malloc<outputT>(sz);
   
   int pdir = FFTW_FORWARD;
   if(m_dir == MXFFT_BACKWARD) pdir = FFTW_BACKWARD;
   
   #ifndef MX_FFTW_NOOMP
   #pragma omp critical
   #endif
   {//scope for pragma
      m_plan = fftw_plan_dft<inputT, outputT>( fftwDimVec<rank>(m_szX, m_szY, m_szZ), forplan1, forplan2,  pdir, FFTW_MEASURE);
   }

   fftw_free<inputT>(forplan1);
   fftw_free<outputT>(forplan2);
   
}

template<typename inputT, typename outputT, size_t rank>
void fftT<inputT,outputT,rank,0>::doPlan(const meta::trueFalseT<true> & inPlace)
{
   (void) inPlace;
   
   complexT * forplan;

   int sz;
   
   if(rank == 1) sz = m_szX;
   if(rank == 2) sz = m_szX*m_szY;
   if(rank == 3) sz = m_szX*m_szY*m_szZ;
   
   forplan = fftw_malloc<complexT>(sz);
   
   int pdir = FFTW_FORWARD;
   if(m_dir == MXFFT_BACKWARD) pdir = FFTW_BACKWARD;

   #ifndef MX_FFTW_NOOMP
   #pragma omp critical
   #endif
   {//scope for pragma
      m_plan = fftw_plan_dft<inputT, outputT>( fftwDimVec<rank>(m_szX, m_szY, m_szZ),  reinterpret_cast<inputT*>(forplan), reinterpret_cast<outputT*>(forplan),  pdir, FFTW_MEASURE);
   }

   fftw_free<inputT>(reinterpret_cast<inputT*>(forplan));
}

template<typename inputT, typename outputT, size_t rank>
template<int crank>
void fftT<inputT,outputT,rank,0>::plan( int nx, 
                                        int ndir, 
                                        bool inPlace,
                                        typename std::enable_if<crank==1>::type* 
                                      )
{
   if(m_szX == nx && m_dir == ndir && m_plan)
   {
      return;
   }
   
   destroyPlan();
   
   m_dir = ndir;
   
   m_szX = nx;
   m_szY = 0;
   m_szZ = 0;
   
   if(inPlace == false)
   {
      doPlan(meta::trueFalseT<false>());
   }
   else
   {
      doPlan(meta::trueFalseT<true>());
   }
}


template<typename inputT, typename outputT, size_t rank>
template<int crank>
void fftT<inputT,outputT,rank,0>::plan( int nx, 
                                        int ny, 
                                        int ndir, 
                                        bool inPlace,
                                        typename std::enable_if<crank==2>::type* 
                                      )
{
   if(m_szX == nx && m_szY == ny  && m_dir == ndir && m_plan)
   {
      return;
   }
   
   destroyPlan();
   
   m_dir = ndir;
   
   m_szX = nx;
   m_szY = ny;
   m_szZ = 0;
   
   if(inPlace == false)
   {
      doPlan(meta::trueFalseT<false>());
   }
   else
   {
      doPlan(meta::trueFalseT<true>());
   }
}

template<typename inputT, typename outputT, size_t rank>
template<int crank>
void fftT<inputT,outputT,rank,0>::plan( int nx, 
                                        int ny,
                                        int nz,
                                        int ndir, 
                                        bool inPlace,
                                        typename std::enable_if<crank==3>::type* 
                                      )
{
   if(m_szX == nx && m_szY == ny && m_szZ == nz && m_dir == ndir && m_plan)
   {
      return;
   }
   
   destroyPlan();
   
   m_dir = ndir;
   
   m_szX = nx;
   m_szY = ny;
   m_szZ = nz;
   
   if(inPlace == false)
   {
      doPlan(meta::trueFalseT<false>());
   }
   else
   {
      doPlan(meta::trueFalseT<true>());
   }
}

template<typename inputT, typename outputT, size_t rank>
void fftT<inputT,outputT,rank,0>::operator()( outputT *out, 
                                              inputT * in
                                            ) const
{
   fftw_execute_dft<inputT,outputT>(m_plan, in, out);
} 
 
 
 
 
 
template<>
std::vector<int> fftwDimVec<1>( int szX, 
                                int szY, 
                                int szZ
                              );

template<>
std::vector<int> fftwDimVec<2>( int szX, 
                                int szY, 
                                int szZ
                              );

template<>
std::vector<int> fftwDimVec<3>( int szX, 
                                int szY, 
                                int szZ
                              );

   
   
   
#ifdef MX_CUDA

template<>
class fft<std::complex<float>, std::complex<float>, 2, 1>
{  
protected:
   int m_szX;
   int m_szY;
    
   cufftHandle m_plan;
   
public:
   
   fft()
   {
      m_szX = 0;
      m_szY = 0;
   }
   
   fft(int nx, int ny)
   {
      m_szX = 0;
      m_szY = 0;
      
      plan(nx, ny);
   }
   
   void plan(int nx, int ny)
   {
      if(m_szX == nx && m_szY == ny)
      {
         return;
      }
      
      m_szX = nx;
      m_szY = ny;
         
      cufftPlan2d(&m_plan, m_szX, m_szY, CUFFT_C2C);
   }
   
   void fwdFFT(std::complex<float> * in, std::complex<float> * out)
   {
      cufftExecC2C(m_plan, (cufftComplex *) in, (cufftComplex *) out, CUFFT_FORWARD);
   }
};

#endif

}//namespace ffit 
}//namespace math
}//namespace mx

#endif // fft_hpp

