
#ifndef __fft_hpp__
#define __fft_hpp__


#include <fftw3.h>

#include <complex>


#include "fftwTemplates.hpp"
#include "../meta/trueFalseT.hpp"

namespace mx
{

#define MXFFT_FORWARD  (FFTW_FORWARD)
#define MXFFT_BACKWARD (FFTW_BACKWARD)



template<typename _inputT, typename _outputT, size_t _dim, int _cudaGPU=0> 
class fftT;

template<int dim>
std::vector<int> fftwDimVec( int szX, int szY, int szZ);

/// Fast Fourier Transforms using the FFTW library
/**
  * \todo document fftT 
  * \todo add execute interface with fftw like signature
  * \todo add plan interface where user passes in pointers (to avoid allocations)
  * 
  */ 
template<typename _inputT, typename _outputT, size_t _dim>
class fftT< _inputT, _outputT, _dim, 0>
{  
   typedef _inputT inputT;
   typedef _outputT outputT;
   
   typedef typename fftwTypeSpec<inputT, outputT>::complexT complexT;
   
   static const size_t dim = _dim;
   
   typedef typename fftwTypeSpec<inputT, outputT>::realT realT;
   
   typedef typename fftwPlanSpec<realT>::planT planT;
   
   
protected:
   int _dir;
   
   int _szX;
   int _szY;
   int _szZ;
   
   planT _plan;
   
public:
   
   
   fftT()
   {
      _szX = 0;
      _szY = 0;
      _szZ = 0;
      
      _plan = 0;
      
      _dir = MXFFT_FORWARD;
   }

   ///Constructor for rank 1 FFT.
   template<int cdim = _dim>
   fftT( int nx, 
         int ndir = MXFFT_FORWARD,
         bool inPlace = false,
         typename std::enable_if<cdim==1>::type* = 0 )
   {
      _plan = 0;
      
      _szX = 0;
      _szY = 0;
      _szZ = 0;
      
      _dir = ndir;
      
      plan(nx, 0, ndir, inPlace);
   }

   ///Constructor for rank 2 FFT.
   template<int cdim = _dim>
   fftT( int nx, 
         int ny, 
         int ndir = MXFFT_FORWARD,
         bool inPlace = false,
         typename std::enable_if<cdim==2>::type* = 0 )
   {
      _plan = 0;
      
      _szX = 0;
      _szY = 0;
      _szZ = 0;
      
      _dir = ndir;
      
      plan(nx, ny, ndir, inPlace);
   }
   
   ~fftT()
   {
      destroy_plan();
   }
   
   void destroy_plan()
   {
      if(_plan) fftw_destroy_plan<realT>(_plan);
      
      _plan = 0;
      
      _szX = 0;
      _szY = 0;

   }
   
   int dir()
   {
      return _dir;
   }
         
   void doPlan(const meta::trueFalseT<false> & inPlace)
   {
      inputT * forplan1;
      outputT * forplan2;
      
      int sz;
      
      if(_dim == 1) sz = _szX;
      if(_dim == 2) sz = _szX*_szY;
      
      forplan1 = fftw_malloc<inputT>(sz);
      forplan2 = fftw_malloc<outputT>(sz);
      
      int pdir = FFTW_FORWARD;
      if(_dir == MXFFT_BACKWARD) pdir = FFTW_BACKWARD;
      
#pragma omp critical
      _plan = fftw_plan_dft<inputT, outputT>( fftwDimVec<dim>(_szX, _szY, _szZ), forplan1, forplan2,  pdir, FFTW_MEASURE);
      
      fftw_free<inputT>(forplan1);
      fftw_free<outputT>(forplan2);
      
   }
   
   void doPlan(const meta::trueFalseT<true> & inPlace)
   {
      complexT * forplan;

      int sz;
      
      if(_dim == 1) sz = _szX;
      if(_dim == 2) sz = _szX*_szY;
      
      forplan = fftw_malloc<complexT>(sz);
      
      int pdir = FFTW_FORWARD;
      if(_dir == MXFFT_BACKWARD) pdir = FFTW_BACKWARD;

#pragma omp critical      
      _plan = fftw_plan_dft<inputT, outputT>( fftwDimVec<dim>(_szX, _szY, _szZ),  reinterpret_cast<inputT*>(forplan), reinterpret_cast<outputT*>(forplan),  pdir, FFTW_MEASURE);
            
      fftw_free<inputT>(forplan);
   }
   
   void plan(int nx, int ny, int ndir=MXFFT_FORWARD, bool inPlace=false)
   {
      if(_szX == nx && _szY == ny  && _dir == ndir && _plan)
      {
         return;
      }
      
      destroy_plan();
      
      _dir = ndir;
      
      _szX = nx;
      _szY = ny;
      
      if(inPlace == false)
      {
         doPlan(meta::trueFalseT<false>());
      }
      else
      {
         doPlan(meta::trueFalseT<true>());
      }
   }

   void operator()( outputT *out, inputT * in)
   {
      fftw_execute_dft<inputT,outputT>(_plan, in, out);
   }
   
};

template<>
std::vector<int> fftwDimVec<1>( int szX, int szY, int szZ)
{
   std::vector<int> v({szX});
   return v;
}

template<>
std::vector<int> fftwDimVec<2>( int szX, int szY, int szZ)
{
   std::vector<int> v({szX, szY});
   return v;
}

template<>
std::vector<int> fftwDimVec<3>( int szX, int szY, int szZ)
{
   std::vector<int> v({szX, szY, szZ});
   return v;
}




   
   
   
#ifdef MX_CUDA

template<>
class fft<std::complex<float>, std::complex<float>, 2, 1>
{  
protected:
   int _szX;
   int _szY;
    
   cufftHandle _plan;
   
public:
   
   fft()
   {
      _szX = 0;
      _szY = 0;
   }
   
   fft(int nx, int ny)
   {
      _szX = 0;
      _szY = 0;
      
      plan(nx, ny);
   }
   
   void plan(int nx, int ny)
   {
      if(_szX == nx && _szY == ny)
      {
         return;
      }
      
      _szX = nx;
      _szY = ny;
         
      cufftPlan2d(&_plan, _szX, _szY, CUFFT_C2C);
   }
   
   void fwdFFT(std::complex<float> * in, std::complex<float> * out)
   {
      cufftExecC2C(_plan, (cufftComplex *) in, (cufftComplex *) out, CUFFT_FORWARD);
   }
};

#endif

}//namespace mx

#endif // __fft_hpp__

