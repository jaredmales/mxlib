/** \file templateFFTW.hpp
  * \brief Declares and defines templatized wrappers for the fftw library
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __templateFFTW_hpp__
#define __templateFFTW_hpp__


#include <fftw3.h>

namespace mx
{

template<typename realT>
int fftw_import_system_wisdom();


template<>
inline int fftw_import_system_wisdom<float>()
{
   return ::fftwf_import_system_wisdom();
}

template<>
inline int fftw_import_system_wisdom<double>()
{
   return ::fftw_import_system_wisdom();
}


template<>
inline int fftw_import_system_wisdom<long double>()
{
   return ::fftwl_import_system_wisdom();
}


template<>
inline int fftw_import_system_wisdom<__float128>()
{
   return ::fftwq_import_system_wisdom();
}


template<typename realT>
inline int fftw_import_wisdom_from_filename( const char * filename );


template<>
inline int fftw_import_wisdom_from_filename<float>( const char * filename )
{
   return ::fftwf_import_wisdom_from_filename( filename );
}


template<>
inline int fftw_import_wisdom_from_filename<double>( const char * filename )
{
   return ::fftw_import_wisdom_from_filename( filename );
}


template<>
inline int fftw_import_wisdom_from_filename<long double>( const char * filename )
{
   return ::fftwl_import_wisdom_from_filename( filename );
}


template<>
inline int fftw_import_wisdom_from_filename<__float128>( const char * filename )
{
   return ::fftwq_import_wisdom_from_filename( filename );
}


template<typename realT>
int fftw_export_wisdom_to_filename(const char * filename);


template<>
inline int fftw_export_wisdom_to_filename<float>(const char *filename)
{
   return ::fftwf_export_wisdom_to_filename( filename );
}


template<>
inline int fftw_export_wisdom_to_filename<double>(const char *filename)
{
   return ::fftw_export_wisdom_to_filename( filename );
}


template<>
inline int fftw_export_wisdom_to_filename<long double>(const char *filename)
{
   return ::fftwl_export_wisdom_to_filename( filename );
}


template<>
inline int fftw_export_wisdom_to_filename<__float128>(const char *filename)
{
   return ::fftwq_export_wisdom_to_filename( filename );
}


template<typename realT>
void fftw_make_planner_thread_safe();


template<>
inline void fftw_make_planner_thread_safe<float>()
{
   ::fftwf_make_planner_thread_safe();
}


template<>
inline void fftw_make_planner_thread_safe<double>()
{
   ::fftw_make_planner_thread_safe();
}


template<>
inline void fftw_make_planner_thread_safe<long double>()
{
   ::fftwl_make_planner_thread_safe();
}


template<>
inline void fftw_make_planner_thread_safe<__float128>()
{
   ::fftwl_make_planner_thread_safe();
}


//************ Allocation ************************//

///Call to fftw_malloc, but with type cast.
/**
  * \param n is the number of type realT elements, not the number of bytes.
  */ 
template<typename realT>
realT * fftw_malloc( size_t n);

template<>
inline float * fftw_malloc<float>(size_t n)
{
   return (float *) ::fftwf_malloc(n*sizeof(float));
}

template<>
inline std::complex<float> * fftw_malloc<std::complex<float>>(size_t n)
{
   return (std::complex<float> *) ::fftwf_malloc(n*sizeof(std::complex<float>));
}

template<>
inline double * fftw_malloc<double>(size_t n)
{
   return (double *) ::fftw_malloc(n*sizeof(double));
}

template<>
inline std::complex<double> * fftw_malloc<std::complex<double>>(size_t n)
{
   return (std::complex<double> *) ::fftw_malloc(n*sizeof(std::complex<double>));
}

template<>
inline long double * fftw_malloc<long double>(size_t n)
{
   return (long double *) ::fftwl_malloc(n*sizeof(long double));
}

template<>
inline std::complex<long double> * fftw_malloc<std::complex<long double>>(size_t n)
{
   return (std::complex<long double> *) ::fftwl_malloc(n*sizeof(std::complex<long double>));
}

template<>
inline __float128 * fftw_malloc<__float128>(size_t n)
{
   return (__float128 *) ::fftwq_malloc(n*sizeof(__float128));
}

template<>
inline std::complex<__float128> * fftw_malloc<std::complex<__float128>>(size_t n)
{
   return (std::complex<__float128> *) ::fftwq_malloc(n*sizeof(std::complex<__float128>));
}

//************ De-Allocation ************************//

///Call to fftw_free.
/**
  * \param p is a pointer to the array to de-allocate.
  */ 
template<typename realT>
void fftw_free( realT * p);

template<>
inline void fftw_free<float>( float * p)
{
   ::fftwf_free( p );
}

template<>
inline void fftw_free<std::complex<float>>( std::complex<float> * p)
{
   ::fftwf_free( p );
}

template<>
inline void fftw_free<double>( double * p)
{
   ::fftw_free( p );
}

template<>
inline void fftw_free<std::complex<double>>( std::complex<double> * p)
{
   ::fftw_free( p );
}

template<>
inline void fftw_free<long double>( long double * p)
{
   ::fftwl_free( p );
}

template<>
inline void fftw_free<std::complex<long double>>( std::complex<long double> * p)
{
   ::fftwl_free( p );
}

template<>
inline void fftw_free<__float128>( __float128 * p)
{
   ::fftwq_free( p );
}

template<>
inline void fftw_free<std::complex<__float128>>( std::complex<__float128> * p)
{
   ::fftwq_free( p );
}

//******* Type Specification *********//
template<typename inputDataT, typename ouputDataT>
struct fftwTypeSpec;

template<>
struct fftwTypeSpec<std::complex<float>, std::complex<float>>
{
   typedef float realT;
   typedef std::complex<float> inputDataT;
   typedef std::complex<float> outputDataT;
   
   typedef fftwf_plan planT;
   
};

template<>
struct fftwTypeSpec<std::complex<double>,std::complex<double>>
{
   typedef double realT;
   typedef std::complex<double> inputDataT;
   typedef std::complex<double> outputDataT;
   
   typedef fftw_plan planT;
   
};

template<>
struct fftwTypeSpec<std::complex<long double>,std::complex<long double>>
{
   typedef long double realT;
   typedef std::complex<long double> inputDataT;
   typedef std::complex<long double> outputDataT;
   
   typedef fftwl_plan planT;
   
};

template<>
struct fftwTypeSpec<std::complex<__float128>,std::complex<__float128>>
{
   typedef __float128 realT;
   typedef std::complex<__float128> inputDataT;
   typedef std::complex<__float128> outputDataT;
   
   typedef fftwq_plan planT;
   
};


//****** Plan Destruction *********//

template<typename dataT>
void fftw_destroy_plan( typename fftwTypeSpec<dataT, dataT>::planT plan );

template<>
inline void fftw_destroy_plan<std::complex<float>>( fftwTypeSpec<std::complex<float>,std::complex<float>>::planT plan )
{
   fftwf_destroy_plan(plan);
}

template<>
inline void fftw_destroy_plan<std::complex<double>>( fftwTypeSpec<std::complex<double>, std::complex<double>>::planT plan )
{
   fftw_destroy_plan(plan);
}

template<>
inline void fftw_destroy_plan<std::complex<long double>>( fftwTypeSpec<std::complex<long double>, std::complex<long double>>::planT plan )
{
   fftwl_destroy_plan(plan);
}

template<>
inline void fftw_destroy_plan<std::complex<__float128>>( fftwTypeSpec<std::complex<__float128>, std::complex<__float128>>::planT plan )
{
   fftwq_destroy_plan(plan);
}


//*********  Plan Execution *************//
template<typename inputDataT, typename outputDataT>
void fftw_execute_dft( typename fftwTypeSpec<inputDataT, outputDataT>::planT plan, 
                       inputDataT * in, 
                       outputDataT * out);

template<>
inline void fftw_execute_dft<std::complex<float>,std::complex<float>>( fftwTypeSpec<std::complex<float>, std::complex<float>>::planT plan, 
                                                                      std::complex<float> * in, 
                                                                      std::complex<float> * out )
{
   fftwf_execute_dft( plan, reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
}


template<>
inline void fftw_execute_dft<std::complex<double>,std::complex<double>>( fftwTypeSpec<std::complex<double>, std::complex<double>>::planT plan, 
                                                                         std::complex<double> * in, 
                                                                         std::complex<double> * out )
{
   fftw_execute_dft( plan, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
}

}//namespace mx

#endif // __templateFFTW_hpp__

