/** \file fftwEnvironment.hpp
  * \brief Declares and defines the fftwEnvironment manager
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

#ifndef fftwEnvironment_hpp
#define fftwEnvironment_hpp

#include <string>

#include "../environment.hpp"


#include "fftwTemplates.hpp"

namespace mx
{

///Return a string corresponding the fftw real floating point type.
/** This is used for wisdom filenames.
  *
  * \tparam realT the real floating point type.
  * 
  * \ingroup fft
  */ 
template<typename realT>
std::string fftw_typename();

template<>
inline std::string fftw_typename<float>()
{
   return "float";
}

template<>
inline std::string fftw_typename<double>()
{
   return "double";
}

template<>
inline std::string fftw_typename<long double>()
{
   return "long_double";
}

#ifdef HASQUAD
template<>
inline std::string fftw_typename<__float128>()
{
   return "quad";
}
#endif

///Create the mxlib standard wisdom filename for the type.
/** Looks for the environment variable MXFFTW_WISDOM.  If not found, then pwd "./" is used.
  *
  * \tparam realT is the real floating point type.
  * 
  * \ingroup fft
  */
template<typename realT>
std::string fftw_wisdom_filename()
{
   std::string path = getEnv("MXFFTW_WISDOM");
   std::string sub = "fftw_wisdom.";
   
   if(path == "") 
   {
      path = ".";
      sub = ".fftw_wisdom.";
   }
   
   std::string filename = path + "/" + sub + fftw_typename<realT>();
   
   return filename;
}

/// Manage the FFTW environment and wisdom using RAII
/** This class manages the FFTW environment.  On construction, wisdom is imported from the 
  * mxlib standard location, and if threads are used threading is initialized and planners are made thread safe.
  *
  * \note the fftw docs recommend against using fftw_make_planner_thread_safe, and says it won't work with omp, but
  * I have not seen any issues using it.
  *
  * On destruction, wisdom is exported and the fftw_cleanup[_threads] function is called.
  *
  * Typically, an object of this type should be created in the main function.  Nothing else needs to be done with it,
  * as it will be destructed on program termination. Example:
   \code
   
   #include <mx/fft/fftwEnvironment.hpp>
   
   int main()
   {
      typedef double realT;
      mx::fftwEnvironment<realT> fftwEnv;
   
      //do stuff . . .
   
      return 0;
   }
   \endcode
  * Note that there is no need to explicitly destroy the object fftwEnv.
  * 
  * \ingroup fft
  */  
template< typename realT, bool threads=true>
struct fftwEnvironment
{
   fftwEnvironment(unsigned nThreads = 1 /**< [in] [optional] the number of threads to use.  This can be changed any time by the program by calling \ref fftw_plan_with_nthreads() */)
   {
      fftw_make_planner_thread_safe<realT>();
      
      if(nThreads = 0) nThreads = 1; //Just to be safe, 1 disables threads but the fftw docs don't say what happens for 0.
      fftw_plan_with_nthreads<realT>(nThreads);

      fftw_import_wisdom_from_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 
   }
   
   ~fftwEnvironment()
   {
      fftw_export_wisdom_to_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 

      fftw_cleanup_threads<realT>();
   }
};


//Partial specialization for no threads.
template< typename realT>
struct fftwEnvironment<realT, false>
{
   fftwEnvironment(unsigned nThreads = 0 )
   {
      fftw_import_wisdom_from_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 
   }
   
   ~fftwEnvironment()
   {
      fftw_export_wisdom_to_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 

      fftw_cleanup<realT>();
   }
};   
   
}//namespace mx

#endif // fftwEnvironment_hpp

