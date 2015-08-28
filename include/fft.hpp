
#ifndef __fft_hpp__
#define __fft_hpp__


#include <fftw3.h>



namespace mx
{

#define MXFFT_FORWARD  (FFTW_FORWARD)
#define MXFFT_BACKWARD (FFTW_BACKWARD)
   
template<typename inT, typename outT, size_t dim, int cudaGPU=0> 
class fftT;

typedef fftT<std::complex<float>, std::complex<float>, 2, 0> fft_cf_cf_2d;
typedef fftT<std::complex<double>, std::complex<double>, 2, 0> fft_cd_cd_2d;

//2D complex-float,  non-cudaGPU
template<>
class fftT<std::complex<float>, std::complex<float>, 2, 0>
{  
protected:
   int _dir;
   
   int _szX;
   int _szY;
    
   fftwf_plan _plan;
   
public:
   
   typedef std::complex<float> inputT;
   typedef std::complex<float> outputT;
   
   fftT()
   {
      _szX = 0;
      _szY = 0;
      _plan = 0;
      
      _dir = MXFFT_FORWARD;
   }
   
   fftT(int ndir)
   {
      _plan = 0;

      _szX = 0;
      _szY = 0;      
      _dir = ndir;
      
   }
   
   fftT(int nx, int ny, int ndir = MXFFT_FORWARD)
   {
      _plan = 0;

      _szX = 0;
      _szY = 0;      
      _dir = ndir;
      
      plan(nx, ny, _dir);
   }
   
   ~fftT()
   {
      destroy_plan();
   }
   
   void destroy_plan()
   {
      if(_plan) fftwf_destroy_plan(_plan);
      
      _plan = 0;
      
      _szX = 0;
      _szY = 0;
   }
      
   void plan(int nx, int ny, int ndir)
   {
      if(_szX == nx && _szY == ny && _dir == ndir && _plan)
      {
         return;
      }   
    
      destroy_plan();
      
      _dir = ndir;
      
      plan(nx, ny);
   }
   
   void plan(int nx, int ny)
   {
      if(_szX == nx && _szY == ny && _plan)
      {
         return;
      }

      destroy_plan();
            
      _szX = nx;
      _szY = ny;
      
      inputT * forplan1;
      outputT * forplan2;

      forplan1 = (inputT *) fftw_malloc(sizeof(inputT)*_szX*_szY);
      forplan2 = (outputT *) fftw_malloc(sizeof(outputT)*_szX*_szY);
      
      int pdir = FFTW_FORWARD;
      if(_dir == MXFFT_BACKWARD) pdir = FFTW_BACKWARD;
      
      _plan = fftwf_plan_dft_2d(_szX, _szY, reinterpret_cast<fftwf_complex*>(forplan1), reinterpret_cast<fftwf_complex*>(forplan2),  pdir, FFTW_MEASURE);
      
      fftw_free(forplan1);
      fftw_free(forplan2);
   }
   
   void fft(inputT * in, outputT * out)
   {
      fftwf_execute_dft(_plan, reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
   }
   
      
};

//2D complex-double, non-cudaGPU
template<>
class fftT<std::complex<double>, std::complex<double>, 2, 0>
{  
protected:
   int _dir;
   
   int _szX;
   int _szY;
    
   fftw_plan _plan;
   
public:
   
   typedef std::complex<double> inputT;
   typedef std::complex<double> outputT;
   
   fftT()
   {
      _szX = 0;
      _szY = 0;
      _plan = 0;
      
      _dir = MXFFT_FORWARD;
   }

   fftT(int ndir)
   {
      _plan = 0;

      _szX = 0;
      _szY = 0;      
      _dir = ndir;
      
   }
   
   fftT(int nx, int ny, int ndir = MXFFT_FORWARD)
   {
      _plan = 0;
      
      _szX = 0;
      _szY = 0;
      _dir = ndir;
      
      plan(nx, ny, ndir);
   }
   
   ~fftT()
   {
      destroy_plan();
   }
   
   void destroy_plan()
   {
      if(_plan) fftw_destroy_plan(_plan);
      
      _plan = 0;
      
      _szX = 0;
      _szY = 0;

   }
   
   int dir()
   {
      return _dir;
   }
   
   void plan(int nx, int ny, int ndir)
   {
      if(_szX == nx && _szY == ny && _dir == ndir && _plan)
      {
         return;
      }
      
      destroy_plan();
      
      _dir = ndir;
      
      plan(nx, ny);
   }

      
   void plan(int nx, int ny)
   {
      if(_szX == nx && _szY == ny && _plan)
      {
         return;
      }
      
      destroy_plan();
      
      _szX = nx;
      _szY = ny;
      

      inputT * forplan1;
      outputT * forplan2;

      forplan1 = (inputT *) fftw_malloc(sizeof(inputT)*_szX*_szY);
      forplan2 = (outputT *) fftw_malloc(sizeof(outputT)*_szX*_szY);
      
      int pdir = FFTW_FORWARD;
      if(_dir == MXFFT_BACKWARD) pdir = FFTW_BACKWARD;

      _plan = fftw_plan_dft_2d(_szX, _szY, reinterpret_cast<fftw_complex*>(forplan1), reinterpret_cast<fftw_complex*>(forplan2),  pdir, FFTW_MEASURE);
      
      fftw_free(forplan1);
      fftw_free(forplan2);
   }
   
   void fft(inputT * in, outputT * out)
   {
      fftw_execute_dft(_plan, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
   }
   
};

   
   
   
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

