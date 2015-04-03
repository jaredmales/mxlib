
#ifndef __fft_hpp__
#define __fft_hpp__


#include <fftw3.h>

namespace mx
{

// template<class inT, class outT>
// fftw_plan fftw_plan_simple(size_t N, inT * in, outT* out, int dir)
// {
//    std::cerr << "Must specialize the 1D fft plan\n";
// }
// 
// template<class inT, class outT>
// fftw_plan fftw_plan_simple(size_t N1, size_t N2, inT * in, outT* out, int dir)
// {
//    fftw_plan p;
//    std::cerr << "Must specialize the 2D fft plan\n";
//    return p;
// }
// 
// template<>
// fftw_plan fftw_plan_simple<float, std::complex<float> >( size_t N, 
//                                                            float * in, 
//                                                            std::complex<float> * out, 
//                                                            int dir )
// {
//    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, reinterpret_cast<fftw_complex*>(out),  FFTW_ESTIMATE);
//    
//    return p;
// }
// 
// template<>
// fftw_plan fftw_plan_simple<double, std::complex<double> >( size_t N, 
//                                                            double * in, 
//                                                            std::complex<double> * out, 
//                                                            int dir )
// {
//    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, reinterpret_cast<fftw_complex*>(out),  FFTW_ESTIMATE);
//    
//    return p;
// }
// 
// 
// 
// template<>
// fftw_plan fftw_plan_simple<double, std::complex<double> >( size_t n0, 
//                                                            size_t n1, 
//                                                            double * in, 
//                                                            std::complex<double> * out, 
//                                                            int dir )
// {
//    fftw_plan p = fftw_plan_dft_r2c_2d(n0, n1, in, reinterpret_cast<fftw_complex*>(out),  FFTW_ESTIMATE);
//    
//    return p;
// }
// 
// template<>
// fftw_plan fftw_plan_simple<std::complex<double>, std::complex<double> >( size_t n0, 
//                                                                               size_t n1, 
//                                                                               std::complex<double> * in, 
//                                                                               std::complex<double> * out, 
//                                                                               int dir
//                                                                              )
// {
//    fftw_plan p = fftw_plan_dft_2d(n0, n1, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), dir, FFTW_ESTIMATE);
//    
//    return p;
// }


template<typename arrayOut, typename arrayIn, typename inT = typename arrayIn::Scalar, int cols = arrayIn::ColsAtCompileTime>
struct _fft
{
   void operator()(arrayOut & out, arrayIn & in, int sign = FFTW_FORWARD)
   {
      std::cout << "none\n";
   }
};


template<typename arrayOut, typename arrayIn>
struct _fft<arrayOut, arrayIn, std::complex<float>, -1>
{
   void operator()(arrayOut & out, arrayIn & in, int sign = FFTW_FORWARD)
   {      
      fftwf_plan p = fftwf_plan_dft_2d(in.rows(), in.cols(), reinterpret_cast<fftwf_complex*>(in.data()), reinterpret_cast<fftwf_complex*>(out.data()),  sign, FFTW_ESTIMATE);
      
      fftwf_execute(p);

      fftwf_destroy_plan(p);
   }
};

template<typename arrayOut, typename arrayIn>
struct _fft<arrayOut, arrayIn, std::complex<double>, -1>
{
   void operator()(arrayOut & out, arrayIn & in, int sign = FFTW_FORWARD)
   {      
      fftw_plan p = fftw_plan_dft_2d(in.rows(), in.cols(), reinterpret_cast<fftw_complex*>(in.data()), reinterpret_cast<fftw_complex*>(out.data()),  sign, FFTW_ESTIMATE);
      
      fftw_execute(p);

      fftw_destroy_plan(p);
   }
};


template<typename arrayOut, typename arrayIn>
void fft(arrayOut & out, arrayIn & in, int sign= FFTW_FORWARD)
{
   _fft<arrayOut, arrayIn> dofft;
   dofft(out, in, sign);
}




}//namespace mx

#endif // __fft_hpp__

