/** \file fftwEnvironment.hpp
  * \brief Declares and defines the fftwEnvironment manager
  * \ingroup fft_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef fftwEnvironment_hpp
#define fftwEnvironment_hpp

#include <string>

#include "../../sys/environment.hpp"


#include "fftwTemplates.hpp"

namespace mx
{
namespace math 
{
namespace fft 
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
std::string fftw_typename<float>();

template<>
std::string fftw_typename<double>();

template<>
std::string fftw_typename<long double>();

#ifdef HASQUAD
template<>
std::string fftw_typename<__float128>();
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
   std::string path = sys::getEnv("MXFFTW_WISDOM");
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
  * \note the fftw docs recommend against using fftw_make_planner_thread_safe.  It does seem to play poorly with omp.
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
template< typename realT, bool threads=false>
struct fftwEnvironment
{
   explicit fftwEnvironment(unsigned nThreads = 0 /**< [in] [optional] the number of threads to use.  This can be changed any time by the program by calling \ref fftw_plan_with_nthreads() */)
   {
      static_cast<void>(nThreads);
      errno = 0;
      int rv = fftw_import_wisdom_from_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 
      
   }
   
   ~fftwEnvironment()
   {
      errno = 0;
      int rv = fftw_export_wisdom_to_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 

      fftw_cleanup<realT>();
   }
};   



//Partial specialization for threads.
template< typename realT>
struct fftwEnvironment<realT, true>
{
   explicit fftwEnvironment(unsigned nThreads = 1)
   {
      fftw_make_planner_thread_safe<realT>();
      
      if(nThreads == 0) nThreads = 1; //Just to be safe, 1 disables threads but the fftw docs don't say what happens for 0.
      fftw_plan_with_nthreads<realT>(nThreads);

      fftw_import_wisdom_from_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 
   }
   
   ~fftwEnvironment()
   {
      fftw_export_wisdom_to_filename<realT>(fftw_wisdom_filename<realT>().c_str()); 

      fftw_cleanup_threads<realT>();
   }
};

}//namespace fft 
}//namespace math
}//namespace mx

#endif // fftwEnvironment_hpp

