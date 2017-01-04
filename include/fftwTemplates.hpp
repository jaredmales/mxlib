/** \file fftwTemplates.hpp
  * \brief Declares and defines templatized wrappers for the fftw library
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef __fftwTemplates_hpp__
#define __fftwTemplates_hpp__

#include <complex>
#include <vector>

#include <fftw3.h>

namespace mx
{

typedef std::complex<float>       complexFT;   
typedef std::complex<double>      complexDT;
typedef std::complex<long double> complexLT;
typedef std::complex<__float128>  complexQT;

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
inline complexFT * fftw_malloc<complexFT>(size_t n)
{
   return (complexFT *) ::fftwf_malloc(n*sizeof(complexFT));
}

template<>
inline double * fftw_malloc<double>(size_t n)
{
   return (double *) ::fftw_malloc(n*sizeof(double));
}

template<>
inline complexDT * fftw_malloc<complexDT>(size_t n)
{
   return (complexDT *) ::fftw_malloc(n*sizeof(complexDT));
}

template<>
inline long double * fftw_malloc<long double>(size_t n)
{
   return (long double *) ::fftwl_malloc(n*sizeof(long double));
}

template<>
inline complexLT * fftw_malloc<complexLT>(size_t n)
{
   return (complexLT *) ::fftwl_malloc(n*sizeof(complexLT));
}

template<>
inline __float128 * fftw_malloc<__float128>(size_t n)
{
   return (__float128 *) ::fftwq_malloc(n*sizeof(__float128));
}

template<>
inline complexQT * fftw_malloc<complexQT>(size_t n)
{
   return (complexQT *) ::fftwq_malloc(n*sizeof(complexQT));
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
inline void fftw_free<complexFT>( complexFT * p)
{
   ::fftwf_free( p );
}

template<>
inline void fftw_free<double>( double * p)
{
   ::fftw_free( p );
}

template<>
inline void fftw_free<complexDT>( complexDT * p)
{
   ::fftw_free( p );
}

template<>
inline void fftw_free<long double>( long double * p)
{
   ::fftwl_free( p );
}

template<>
inline void fftw_free<complexLT>( complexLT * p)
{
   ::fftwl_free( p );
}

template<>
inline void fftw_free<__float128>( __float128 * p)
{
   ::fftwq_free( p );
}

template<>
inline void fftw_free<complexQT>( complexQT * p)
{
   ::fftwq_free( p );
}

//****** Plans ******//
template<typename realT>
struct fftwPlanSpec;

template<>
struct fftwPlanSpec<float>
{
   typedef fftwf_plan planT;
};

template<>
struct fftwPlanSpec<double>
{
   typedef fftw_plan planT;
};

template<>
struct fftwPlanSpec<long double>
{
   typedef fftwl_plan planT;
};

template<>
struct fftwPlanSpec<__float128>
{
   typedef fftwq_plan planT;
};


//******* Type Specification *********//
template<typename inputDataT, typename ouputDataT>
struct fftwTypeSpec;

template<>
struct fftwTypeSpec<complexFT, complexFT>
{
   typedef float realT;
   typedef complexFT inputDataT;
   typedef complexFT outputDataT;
   
   typedef fftwPlanSpec<float>::planT planT;
   
};

template<>
struct fftwTypeSpec<float, complexFT>
{
   typedef float realT;
   typedef float inputDataT;
   typedef complexFT outputDataT;
   
   typedef fftwPlanSpec<float>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexFT, float>
{
   typedef float realT;
   typedef complexFT inputDataT;
   typedef float outputDataT;
   
   typedef fftwPlanSpec<float>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexDT,complexDT>
{
   typedef double realT;
   typedef complexDT inputDataT;
   typedef complexDT outputDataT;
   
   typedef fftwPlanSpec<double>::planT planT;
   
};

template<>
struct fftwTypeSpec<double,complexDT>
{
   typedef double realT;
   typedef double inputDataT;
   typedef complexDT outputDataT;
   
   typedef fftwPlanSpec<double>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexDT,double>
{
   typedef double realT;
   typedef complexDT inputDataT;
   typedef double outputDataT;
   
   typedef fftwPlanSpec<double>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexLT,complexLT>
{
   typedef long double realT;
   typedef complexLT inputDataT;
   typedef complexLT outputDataT;
   
   typedef fftwPlanSpec<long double>::planT planT;
   
};

template<>
struct fftwTypeSpec<long double,complexLT>
{
   typedef long double realT;
   typedef long double inputDataT;
   typedef complexLT outputDataT;
   
   typedef fftwPlanSpec<long double>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexLT,long double>
{
   typedef long double realT;
   typedef complexLT inputDataT;
   typedef long double outputDataT;
   
   typedef fftwPlanSpec<long double>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexQT,complexQT>
{
   typedef __float128 realT;
   typedef complexQT inputDataT;
   typedef complexQT outputDataT;
   
   typedef fftwPlanSpec<__float128>::planT planT;
   
};

template<>
struct fftwTypeSpec<__float128,complexQT>
{
   typedef __float128 realT;
   typedef __float128 inputDataT;
   typedef complexQT outputDataT;
   
   typedef fftwPlanSpec<__float128>::planT planT;
   
};

template<>
struct fftwTypeSpec<complexQT,__float128>
{
   typedef __float128 realT;
   typedef complexQT inputDataT;
   typedef __float128 outputDataT;
   
   typedef fftwPlanSpec<__float128>::planT planT;
   
};

//********* Plan Creation **************//
template<typename inputDataT, typename outputDataT>
typename fftwTypeSpec<inputDataT,outputDataT>::planT fftw_plan_dft( std::vector<int> n, 
                                                                    inputDataT * in, 
                                                                    outputDataT * out,
                                                                    int sign,
                                                                    unsigned flags );

template<>
inline fftwTypeSpec<complexFT, complexFT>::planT fftw_plan_dft<complexFT, complexFT>( std::vector<int> n, 
                                                                                      complexFT * in, 
                                                                                      complexFT * out,
                                                                                      int sign,
                                                                                      unsigned flags )
{
   return ::fftwf_plan_dft( n.size(), n.data(), reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out), sign, flags);
}


template<>
inline fftwTypeSpec<float, complexFT>::planT fftw_plan_dft<float, complexFT>( std::vector<int> n, 
                                                                              float * in, 
                                                                              complexFT * out,
                                                                              int sign,
                                                                              unsigned flags )
{
   return ::fftwf_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftwf_complex*>(out), flags);
}

template<>
inline fftwTypeSpec<complexFT, float>::planT fftw_plan_dft<complexFT, float>( std::vector<int> n, 
                                                                              complexFT * in, 
                                                                              float * out,
                                                                              int sign,
                                                                              unsigned flags )
{
   return ::fftwf_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftwf_complex*>(in), out, flags);
}



template<>
inline fftwTypeSpec<complexDT, complexDT>::planT fftw_plan_dft<complexDT, complexDT>( std::vector<int> n, 
                                                                                      complexDT * in, 
                                                                                      complexDT * out,
                                                                                      int sign,
                                                                                      unsigned flags )
{
   return ::fftw_plan_dft( n.size(), n.data(), reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), sign, flags);
}

template<>
inline fftwTypeSpec<double, complexDT>::planT fftw_plan_dft<double, complexDT>( std::vector<int> n, 
                                                                                double * in, 
                                                                                complexDT * out,
                                                                                int sign,
                                                                                unsigned flags )
{
   return ::fftw_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftw_complex*>(out), flags);
}

template<>
inline fftwTypeSpec<complexDT, double>::planT fftw_plan_dft<complexDT, double>( std::vector<int> n, 
                                                                                complexDT * in, 
                                                                                double * out,
                                                                                int sign,
                                                                                unsigned flags )
{
   return ::fftw_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftw_complex*>(in), out, flags);
}



template<>
inline fftwTypeSpec<complexLT, complexLT>::planT fftw_plan_dft<complexLT, complexLT>( std::vector<int> n, 
                                                                                      complexLT * in, 
                                                                                      complexLT * out,
                                                                                      int sign,
                                                                                      unsigned flags )
{
   return ::fftwl_plan_dft( n.size(), n.data(), reinterpret_cast<fftwl_complex*>(in), reinterpret_cast<fftwl_complex*>(out), sign, flags);
}

template<>
inline fftwTypeSpec<long double, complexLT>::planT fftw_plan_dft<long double, complexLT>( std::vector<int> n, 
                                                                                          long double * in, 
                                                                                          complexLT * out,
                                                                                          int sign,
                                                                                          unsigned flags )
{
   return ::fftwl_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftwl_complex*>(out), flags);
}

template<>
inline fftwTypeSpec<complexLT, long double>::planT fftw_plan_dft<complexLT, long double>( std::vector<int> n, 
                                                                                          complexLT * in, 
                                                                                          long double * out,
                                                                                          int sign,
                                                                                          unsigned flags )
{
   return ::fftwl_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftwl_complex*>(in), out, flags);
}

template<>
inline fftwTypeSpec<complexQT, complexQT>::planT fftw_plan_dft<complexQT, complexQT>( std::vector<int> n, 
                                                                                      complexQT * in, 
                                                                                      complexQT * out,
                                                                                      int sign,
                                                                                      unsigned flags )
{
   return ::fftwq_plan_dft( n.size(), n.data(), reinterpret_cast<fftwq_complex*>(in), reinterpret_cast<fftwq_complex*>(out), sign, flags);
}

template<>
inline fftwTypeSpec<__float128, complexQT>::planT fftw_plan_dft<__float128, complexQT>( std::vector<int> n, 
                                                                                        __float128 * in, 
                                                                                        complexQT * out,
                                                                                        int sign,
                                                                                        unsigned flags )
{
   return ::fftwq_plan_dft_r2c( n.size(), n.data(), in, reinterpret_cast<fftwq_complex*>(out), flags);
}

template<>
inline fftwTypeSpec<complexQT, __float128>::planT fftw_plan_dft<complexQT, __float128>( std::vector<int> n, 
                                                                                        complexQT * in, 
                                                                                        __float128 * out,
                                                                                        int sign,
                                                                                        unsigned flags )
{
   return ::fftwq_plan_dft_c2r( n.size(), n.data(), reinterpret_cast<fftwq_complex*>(in), out, flags);
}

//*********  Plan Execution *************//
template<typename inputDataT, typename outputDataT>
void fftw_execute_dft( typename fftwTypeSpec<inputDataT, outputDataT>::planT plan, 
                       inputDataT * in, 
                       outputDataT * out);

template<>
inline void fftw_execute_dft<complexFT,complexFT>( fftwTypeSpec<complexFT, complexFT>::planT plan, 
                                                   complexFT * in, 
                                                   complexFT * out )
{
   ::fftwf_execute_dft( plan, reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out));
}


template<>
inline void fftw_execute_dft<complexDT,complexDT>( fftwTypeSpec<complexDT, complexDT>::planT plan, 
                                                   complexDT * in, 
                                                   complexDT * out )
{
   ::fftw_execute_dft( plan, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
}


//****** Plan Destruction *********//

template<typename realT>
void fftw_destroy_plan( typename fftwPlanSpec<realT>::planT plan );

template<>
inline void fftw_destroy_plan<float>( fftwPlanSpec<float>::planT plan )
{
   ::fftwf_destroy_plan(plan);
}

template<>
inline void fftw_destroy_plan<double>( fftwPlanSpec<double>::planT plan )
{
   ::fftw_destroy_plan(plan);
}

template<>
inline void fftw_destroy_plan<long double>( fftwPlanSpec<long double>::planT plan )
{
   ::fftwl_destroy_plan(plan);
}

template<>
inline void fftw_destroy_plan<__float128>( fftwPlanSpec<__float128>::planT plan )
{
   ::fftwq_destroy_plan(plan);
}



}//namespace mx

#endif // __fftwTemplates_hpp__

