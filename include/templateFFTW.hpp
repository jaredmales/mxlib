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

}//namespace mx

#endif // __templateFFTW_hpp__

