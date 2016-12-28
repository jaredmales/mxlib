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


}//namespace mx

#endif // __templateFFTW_hpp__

