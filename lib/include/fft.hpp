
#ifndef __fft_hpp__
#define __fft_hpp__


#include <fftw3.h>

namespace mx
{

template<class inT, class outT>
fftw_plan fftw_plan_simple(size_t N, inT * in, outT* out, int dir)
{
   std::cerr << "Must specialize the 1D fft plan\n";
}

template<class inT, class outT>
fftw_plan fftw_plan_simple(size_t N1, size_t N2, inT * in, outT* out, int dir)
{
   fftw_plan p;
   std::cerr << "Must specialize the 2D fft plan\n";
   return p;
}

template<>
fftw_plan fftw_plan_simple<double, std::complex<double> >( size_t N, 
                                                           double * in, 
                                                           std::complex<double> * out, 
                                                           int dir )
{
   fftw_plan p = fftw_plan_dft_r2c_1d(N, in, reinterpret_cast<fftw_complex*>(out),  FFTW_ESTIMATE);
   
   return p;
}

template<>
fftw_plan fftw_plan_simple<double, std::complex<double> >( size_t n0, 
                                                           size_t n1, 
                                                           double * in, 
                                                           std::complex<double> * out, 
                                                           int dir )
{
   fftw_plan p = fftw_plan_dft_r2c_2d(n0, n1, in, reinterpret_cast<fftw_complex*>(out),  FFTW_ESTIMATE);
   
   return p;
}

template<>
fftw_plan fftw_plan_simple<std::complex<double>, std::complex<double> >( size_t n0, 
                                                                              size_t n1, 
                                                                              std::complex<double> * in, 
                                                                              std::complex<double> * out, 
                                                                              int dir
                                                                             )
{
   fftw_plan p = fftw_plan_dft_2d(n0, n1, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), dir, FFTW_ESTIMATE);
   
   return p;
}


template<int dim, class Matrix_in, class Matrix_out> 
void fft( Matrix_in & m_in, Matrix_out & m_out, int dir=-1)
{
   
   fftw_plan p;
   
   if(dim == 1)
      p = fftw_plan_simple(m_in.size(), m_in.data(), m_out.data(), dir);
   
   if(dim == 2)
      p = fftw_plan_simple(m_in.rows(), m_in.cols(), m_in.data(), m_out.data(), dir);
   
   fftw_execute(p);
   
   fftw_destroy_plan(p);
}




}//namespace mx

#endif // __mx_vmopfft__

