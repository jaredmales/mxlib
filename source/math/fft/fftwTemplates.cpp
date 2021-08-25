/** \file fftwTemplates.cpp
  * \brief Definitions of templatized wrappers for the fftw library
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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


#include "math/fft/fftwTemplates.hpp"

namespace mx
{
namespace math
{
namespace fft
{
   
template<>
int fftw_import_system_wisdom<float>()
{
   return ::fftwf_import_system_wisdom();
}

template<>
int fftw_import_system_wisdom<double>()
{
   return ::fftw_import_system_wisdom();
}

template<>
int fftw_import_system_wisdom<long double>()
{
   return ::fftwl_import_system_wisdom();
}

#ifdef HASQUAD
template<>
int fftw_import_system_wisdom<__float128>()
{
   return ::fftwq_import_system_wisdom();
}
#endif

template<>
int fftw_import_wisdom_from_filename<float>( const char * filename )
{
   return ::fftwf_import_wisdom_from_filename( filename );
}

template<>
int fftw_import_wisdom_from_filename<double>( const char * filename )
{
   return ::fftw_import_wisdom_from_filename( filename );
}

template<>
int fftw_import_wisdom_from_filename<long double>( const char * filename )
{
   return ::fftwl_import_wisdom_from_filename( filename );
}

#ifdef HASQUAD
template<>
int fftw_import_wisdom_from_filename<__float128>( const char * filename )
{
   return ::fftwq_import_wisdom_from_filename( filename );
}
#endif

template<>
int fftw_export_wisdom_to_filename<float>(const char *filename)
{
   return ::fftwf_export_wisdom_to_filename( filename );
}

template<>
int fftw_export_wisdom_to_filename<double>(const char *filename)
{
   return ::fftw_export_wisdom_to_filename( filename );
}

template<>
int fftw_export_wisdom_to_filename<long double>(const char *filename)
{
   return ::fftwl_export_wisdom_to_filename( filename );
}

#ifdef HASQUAD
template<>
int fftw_export_wisdom_to_filename<__float128>(const char *filename)
{
   return ::fftwq_export_wisdom_to_filename( filename );
}
#endif

template<>
float * fftw_malloc<float>(size_t n)
{
   return (float *) ::fftwf_malloc(n*sizeof(float));
}

template<>
complexFT * fftw_malloc<complexFT>(size_t n)
{
   return (complexFT *) ::fftwf_malloc(n*sizeof(complexFT));
}

template<>
double * fftw_malloc<double>(size_t n)
{
   return (double *) ::fftw_malloc(n*sizeof(double));
}

template<>
complexDT * fftw_malloc<complexDT>(size_t n)
{
   return (complexDT *) ::fftw_malloc(n*sizeof(complexDT));
}

template<>
long double * fftw_malloc<long double>(size_t n)
{
   return (long double *) ::fftwl_malloc(n*sizeof(long double));
}

template<>
complexLT * fftw_malloc<complexLT>(size_t n)
{
   return (complexLT *) ::fftwl_malloc(n*sizeof(complexLT));
}

#ifdef HASQUAD
template<>
__float128 * fftw_malloc<__float128>(size_t n)
{
   return (__float128 *) ::fftwq_malloc(n*sizeof(__float128));
}

template<>
complexQT * fftw_malloc<complexQT>(size_t n)
{
   return (complexQT *) ::fftwq_malloc(n*sizeof(complexQT));
}
#endif

template<>
void fftw_free<float>( float * p)
{
   ::fftwf_free( p );
}

template<>
void fftw_free<complexFT>( complexFT * p)
{
   ::fftwf_free( p );
}

template<>
void fftw_free<double>( double * p)
{
   ::fftw_free( p );
}

template<>
void fftw_free<complexDT>( complexDT * p)
{
   ::fftw_free( p );
}

template<>
void fftw_free<long double>( long double * p)
{
   ::fftwl_free( p );
}

template<>
void fftw_free<complexLT>( complexLT * p)
{
   ::fftwl_free( p );
}

#ifdef HASQUAD
template<>
void fftw_free<__float128>( __float128 * p)
{
   ::fftwq_free( p );
}

template<>
void fftw_free<complexQT>( complexQT * p)
{
   ::fftwq_free( p );
}
#endif

template<>
void fftw_make_planner_thread_safe<float>()
{
   //::fftwf_make_planner_thread_safe();
}


template<>
void fftw_make_planner_thread_safe<double>()
{
   //::fftw_make_planner_thread_safe();
}


template<>
void fftw_make_planner_thread_safe<long double>()
{
   //::fftwl_make_planner_thread_safe();
}

#ifdef HASQUAD
template<>
void fftw_make_planner_thread_safe<__float128>()
{
   //::fftwl_make_planner_thread_safe();
}
#endif

template<>
void fftw_plan_with_nthreads<float>(int nthreads)
{
   //::fftwf_plan_with_nthreads(nthreads);
}

template<>
void fftw_plan_with_nthreads<double>(int nthreads)
{
   //::fftw_plan_with_nthreads(nthreads);
}

template<>
void fftw_plan_with_nthreads<long double>(int nthreads)
{
   //::fftwl_plan_with_nthreads(nthreads);
}

#ifdef HASQUAD
template<>
void fftw_plan_with_nthreads<__float128>(int nthreads)
{
   //::fftwq_plan_with_nthreads(nthreads);
}
#endif

template<>
fftwTypeSpec<complexFT, complexFT>::planT fftw_plan_dft<complexFT, complexFT>( std::vector<int> n, 
                                                                                      complexFT * in, 
                                                                                      complexFT * out,
                                                                                      int sign,
                                                                                      unsigned flags )
{
   return ::fftwf_plan_dft( n.size(), n.data(), reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out), sign, flags);
}


template<>
fftwTypeSpec<float, complexFT>::planT fftw_plan_dft<float, complexFT>( std::vector<int> n, 
                                                                       float * in, 
                                                                       complexFT * out,
                                                                       int sign,
                                                                       unsigned flags 
                                                                     )
{
   static_cast<void>(sign);
   return ::fftwf_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftwf_complex*>(out), flags);
}

template<>
fftwTypeSpec<complexFT, float>::planT fftw_plan_dft<complexFT, float>( std::vector<int> n, 
                                                                       complexFT * in, 
                                                                       float * out,
                                                                       int sign,
                                                                       unsigned flags 
                                                                     )
{
   static_cast<void>(sign);
   return ::fftwf_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftwf_complex*>(in), out, flags);
}



template<>
fftwTypeSpec<complexDT, complexDT>::planT fftw_plan_dft<complexDT, complexDT>( std::vector<int> n, 
                                                                               complexDT * in, 
                                                                               complexDT * out,
                                                                               int sign,
                                                                               unsigned flags 
                                                                             )
{
   return ::fftw_plan_dft( n.size(), n.data(), reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), sign, flags);
}

template<>
fftwTypeSpec<double, complexDT>::planT fftw_plan_dft<double, complexDT>( std::vector<int> n, 
                                                                         double * in, 
                                                                         complexDT * out,
                                                                         int sign,
                                                                         unsigned flags 
                                                                       )
{
   static_cast<void>(sign);
   return ::fftw_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftw_complex*>(out), flags);
}

template<>
fftwTypeSpec<complexDT, double>::planT fftw_plan_dft<complexDT, double>( std::vector<int> n, 
                                                                         complexDT * in, 
                                                                         double * out,
                                                                         int sign,
                                                                         unsigned flags 
                                                                       )
{
   static_cast<void>(sign);
   return ::fftw_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftw_complex*>(in), out, flags);
}



template<>
fftwTypeSpec<complexLT, complexLT>::planT fftw_plan_dft<complexLT, complexLT>( std::vector<int> n, 
                                                                               complexLT * in, 
                                                                               complexLT * out,
                                                                               int sign,
                                                                               unsigned flags 
                                                                             )
{
   return ::fftwl_plan_dft( n.size(), n.data(), reinterpret_cast<fftwl_complex*>(in), reinterpret_cast<fftwl_complex*>(out), sign, flags);
}

template<>
fftwTypeSpec<long double, complexLT>::planT fftw_plan_dft<long double, complexLT>( std::vector<int> n, 
                                                                                   long double * in, 
                                                                                   complexLT * out,
                                                                                   int sign,
                                                                                   unsigned flags 
                                                                                 )
{
   static_cast<void>(sign);
   return ::fftwl_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftwl_complex*>(out), flags);
}

template<>
fftwTypeSpec<complexLT, long double>::planT fftw_plan_dft<complexLT, long double>( std::vector<int> n, 
                                                                                   complexLT * in, 
                                                                                   long double * out,
                                                                                   int sign,
                                                                                   unsigned flags 
                                                                                 )
{
   static_cast<void>(sign);
   return ::fftwl_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftwl_complex*>(in), out, flags);
}

#ifdef HASQUAD
template<>
fftwTypeSpec<complexQT, complexQT>::planT fftw_plan_dft<complexQT, complexQT>( std::vector<int> n, 
                                                                               complexQT * in, 
                                                                               complexQT * out,
                                                                               int sign,
                                                                               unsigned flags )
{
   return ::fftwq_plan_dft( n.size(), n.data(), reinterpret_cast<fftwq_complex*>(in), reinterpret_cast<fftwq_complex*>(out), sign, flags);
}

template<>
fftwTypeSpec<__float128, complexQT>::planT fftw_plan_dft<__float128, complexQT>( std::vector<int> n, 
                                                                                 __float128 * in, 
                                                                                 complexQT * out,
                                                                                 int sign,
                                                                                 unsigned flags )
{
   return ::fftwq_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftwq_complex*>(out), flags);
}

template<>
fftwTypeSpec<complexQT, __float128>::planT fftw_plan_dft<complexQT, __float128>( std::vector<int> n, 
                                                                                 complexQT * in, 
                                                                                 __float128 * out,
                                                                                 int sign,
                                                                                 unsigned flags )
{
   return ::fftwq_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftwq_complex*>(in), out, flags);
}

#endif

template<>
void fftw_cleanup<float>()
{
   ::fftwf_cleanup();
}

template<>
void fftw_cleanup<double>()
{
   ::fftw_cleanup();
}

template<>
void fftw_cleanup<long double>()
{
   ::fftwl_cleanup();
}

#ifdef HASQUAD
template<>
void fftw_cleanup<__float128>()
{
   ::fftwq_cleanup();
}
#endif

template<>
void fftw_cleanup_threads<float>()
{
   //::fftwf_cleanup_threads();
}

template<>
void fftw_cleanup_threads<double>()
{
   //::fftw_cleanup_threads();
}

template<>
void fftw_cleanup_threads<long double>()
{
   //::fftwl_cleanup_threads();
}

#ifdef HASQUAD
template<>
void fftw_cleanup_threads<__float128>()
{
   //::fftwq_cleanup_threads();
}
#endif

template<>
void fftw_execute_dft<complexFT,complexFT>( fftwTypeSpec<complexFT, complexFT>::planT plan, 
                                            complexFT * in, 
                                            complexFT * out )
{
   ::fftwf_execute_dft( plan, reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
}


template<>
void fftw_execute_dft<float,complexFT>( fftwTypeSpec<float, complexFT>::planT plan, 
                                        float * in, 
                                        complexFT * out )
{
   ::fftwf_execute_dft_r2c( plan, in, reinterpret_cast<fftwf_complex*>(out));
}

template<>
void fftw_execute_dft<complexFT,float>( fftwTypeSpec<complexFT,float>::planT plan, 
                                        complexFT * in,
                                        float * out 
                                      )
{
   ::fftwf_execute_dft_c2r( plan, reinterpret_cast<fftwf_complex*>(in), out );
}

template<>
void fftw_execute_dft<complexDT,complexDT>( fftwTypeSpec<complexDT, complexDT>::planT plan, 
                                            complexDT * in, 
                                            complexDT * out )
{
   ::fftw_execute_dft( plan, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
}

template<>
void fftw_execute_dft<double,complexDT>( fftwTypeSpec<double, complexDT>::planT plan, 
                                         double * in, 
                                         complexDT * out )
{
   ::fftw_execute_dft_r2c( plan, in, reinterpret_cast<fftw_complex*>(out));
}

template<>
void fftw_execute_dft<complexDT, double>( fftwTypeSpec<complexDT,double>::planT plan, 
                                          complexDT * in, 
                                          double * out
                                        )
{
   ::fftw_execute_dft_c2r( plan, reinterpret_cast<fftw_complex*>(in), out);
}


template<>
void fftw_destroy_plan<float>( fftwPlanSpec<float>::planT plan )
{
   ::fftwf_destroy_plan(plan);
}

template<>
void fftw_destroy_plan<double>( fftwPlanSpec<double>::planT plan )
{
   ::fftw_destroy_plan(plan);
}

template<>
void fftw_destroy_plan<long double>( fftwPlanSpec<long double>::planT plan )
{
   ::fftwl_destroy_plan(plan);
}

#ifdef HASQUAD
template<>
void fftw_destroy_plan<__float128>( fftwPlanSpec<__float128>::planT plan )
{
   ::fftwq_destroy_plan(plan);
}

#endif


} //namespace fft
} //namespace math
} //namespace mx


