/** \file psdFilter.hpp
  * \brief Declares and defines a class for filtering with PSDs
  * \ingroup signal_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018, 2019, 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef psdFilter_hpp
#define psdFilter_hpp

#include <vector>
#include <complex>
#include <Eigen/Dense>

#include "../mxError.hpp"
#include "../math/fft/fft.hpp"
#include "../improc/eigenCube.hpp"

namespace mx
{
namespace sigproc 
{


namespace psdFilterTypes
{

///Types for different ranks in psdFilter
template<typename realT, size_t rank>
struct arrayT;

template<typename realT>
struct arrayT<realT, 1>
{
   typedef std::vector<realT> realArrayT;
   typedef std::vector<std::complex<realT>> complexArrayT;
   
   static void clear( realArrayT & arr)
   {
      arr.clear();
   }
   
   static void clear( complexArrayT & arr)
   {
      arr.clear();
   }
   
};

template<typename realT>
struct arrayT<realT, 2>
{
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> realArrayT;
   typedef Eigen::Array<std::complex<realT>,Eigen::Dynamic, Eigen::Dynamic> complexArrayT;
   
   static void clear( realArrayT & arr)
   {
      arr.resize(0,0);
   }
   
   static void clear( complexArrayT & arr)
   {
      arr.resize(0,0);
   }
};

template<typename realT>
struct arrayT<realT, 3>
{
   typedef improc::eigenCube<realT> realArrayT;
   typedef improc::eigenCube<std::complex<realT>> complexArrayT;
   
   static void clear( realArrayT & arr)
   {
      arr.resize(0,0,0);
   }
   
   static void clear( complexArrayT & arr)
   {
      arr.resize(0,0,0);
   }
};

}

//Forward declaration
template<typename _realT, size_t rank, int cuda = 0>
class psdFilter;

/// A class for filtering noise with PSDs
/** The square-root of the PSD is maintained by this class, either as a pointer to an external array or using internally allocated memory (which will be
  * de-allocated on destruction). 
  * 
  * PSD Requirements: 
  * - the PSD must be in FFT storage order form.  That means including negative frequencies reversed from the end of the end of the array.
  * - the PSD used for this needs to be normalized properly, \ref psds "according to the mxlib standard", to produce filtered noise with the correct statistics.  
  * 
  *
  * Array type varies based on rank. 
  * - For rank==1, the array type is std::vector<realT>
  * - for rank==2, the array type is Eigen::Array<realT, -1, -1>
  * - for rank==3, the array type is mx::improc::eigenCube<realT>
  * and likewise for the complex types.
  * 
  * \tparam _realT real floating type
  * \tparam _rank the rank, or dimension, of the PSD
  *
  * \ingroup psds
  *
  * \todo once fftT has a plan interface with pointers for working memory, use it.
  *
  * \test Scenario: compiling psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
  */
template<typename _realT, size_t _rank>
class psdFilter<_realT, _rank, 0>
{
public:
   
   typedef _realT realT; ///< Real floating point type
   typedef std::complex<_realT> complexT; ///< Complex floating point type.

   static const size_t rank = _rank;

   typedef typename psdFilterTypes::arrayT<realT,rank>::realArrayT realArrayT; ///< std::vector for rank==1, Eigen::Array for rank==2, eigenCube for rank==3.
   typedef typename psdFilterTypes::arrayT<realT,rank>::complexArrayT complexArrayT; ///< std::vector for rank==1, Eigen::Array for rank==2, eigenCube for rank==3.
   
   
protected:

   int m_rows {0}; ///< The number of rows in the filter, and the required number of rows in the noise array.
   int m_cols {0}; ///< The number of columns in the filter, and the required number of columns in the noise array.
   int m_planes {0}; ///< Then number of planes in the filter.
   
   realT m_dFreq1 {1.0}; ///< The frequency scaling of the x-dimension.  Used to scale the output.
   realT m_dFreq2 {1.0}; ///< The frequency scaling of the y-dimension.  Used to scale the output.
   realT m_dFreq3 {1.0}; ///< The frequency scaling of the z-dimension.  Used to scale the output.
   
   realArrayT * m_psdSqrt {nullptr}; ///< Pointer to the real array containing the square root of the PSD.
   bool m_owner {false}; ///< Flag indicates whether or not m_psdSqrt was allocated by this instance, and so must be deallocated.

   mutable complexArrayT m_ftWork;   ///< Working memory for the FFT.  Declared mutable so it can be accessed in the const filter method.
   
   math::fft::fftT< complexT, complexT,rank,0> m_fft_fwd; ///< FFT object for the forward transform.
   math::fft::fftT< complexT, complexT,rank,0> m_fft_back; ///< FFT object for the backward transfsorm.
   
public:
   
   ///C'tor.
   /**
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */ 
   psdFilter();
   
   ///Destructor
   ~psdFilter();
   
protected:

   ///Set the sqaure-root of the PSD to be a pointer to an array containing the square root of the properly normalized PSD.
   /** This does not allocate _npsdSqrt, it merely points to the specified array, which remains your responsibility for deallocation, etc.
     *
     * See the discussion of PSD normalization above.
     * 
     * This private version handles the actual setting of m_psdSqrt, which is rank independent.
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   int psdSqrt( realArrayT * npsdSqrt /**< [in] a pointer to an array containing the square root of the PSD. */ );
   
   /// Set the sqaure-root of the PSD.
   /** This allocates _npsdSqrt and fills it with the values in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This private version handles the actual setting of m_psdSqrt, which is rank independent.
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   int psdSqrt( const realArrayT & npsdSqrt /**< [in] an array containing the square root of the PSD. */ );
   
   ///Set the size of the filter.
   /** Handles allocation of the m_ftWork array and fftw planning.
     *
     * Requires m_psdSqrt to be set first.  This is called by the psdSqrt() and psd() methods.
     * 
     * This version compiles when rank==1
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   template<size_t crank=rank>
   int setSize(typename std::enable_if<crank==1>::type* = 0 );

   ///Set the size of the filter.
   /** Handles allocation of the m_ftWork array and fftw planning.
     *
     * Requires m_psdSqrt to be set first.  This is called by the psdSqrt() and psd() methods.
     * 
     * This version compiles when rank==2
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   template<size_t crank=rank>
   int setSize(typename std::enable_if<crank==2>::type* = 0 );
   
   ///Set the size of the filter.
   /** Handles allocation of the m_ftWork array and fftw planning.
     *
     * Requires m_psdSqrt to be set first.  This is called by the psdSqrt() and psd() methods.
     * 
     * This version compiles when rank==3
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   template<size_t crank=rank>
   int setSize(typename std::enable_if<crank==3>::type* = 0 );
   
public:   

   ///Get the number of rows in the filter
   /**
     * \returns the current value of m_rows.
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */ 
   int rows();
   
   ///Get the number of columns in the filter
   /**
     * \returns the current value of m_cols.
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   int cols();
   
   ///Get the number of planes in the filter
   /**
     * \returns the current value of m_planes.
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   int planes();
   
   ///Set the sqaure-root of the PSD to be a pointer to an array containing the square root of the properly normalized PSD.
   /** This does not allocate _npsdSqrt, it merely points to the specified array, which remains your responsibility for deallocation, etc.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==1
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   template<size_t crank = rank>
   int psdSqrt( realArrayT * npsdSqrt, ///< [in] a pointer to an array containing the square root of the PSD.  
                realT df,              ///< [in] the frequency spacing
                typename std::enable_if<crank==1>::type* = 0
              );
   
   ///Set the square-root of the PSD to be a pointer to an array containing the square root of the properly normalized PSD.
   /** This does not allocate _npsdSqrt, it merely points to the specified array, which remains your responsibility for deallocation, etc.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==2
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   template<size_t crank = rank>
   int psdSqrt( realArrayT * npsdSqrt, ///< [in] a pointer to an array containing the square root of the PSD. 
                realT dk1,             ///< [in] the frequency spacing along dimension 1
                realT dk2,             ///< [in] the frequency spacing along dimension 2
                typename std::enable_if<crank==2>::type* = 0
              );
   
   ///Set the sqaure-root of the PSD to be a pointer to an array containing the square root of the properly normalized PSD.
   /** This does not allocate _npsdSqrt, it merely points to the specified array, which remains your responsibility for deallocation, etc.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==3
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   template<size_t crank = rank>
   int psdSqrt( realArrayT * npsdSqrt, ///< [in] a pointer to an array containing the square root of the PSD. 
                realT dk1,             ///< [in] the frequency spacing along dimension 1
                realT dk2,             ///< [in] the frequency spacing along dimension 2
                realT df,              ///< [in] the frequency spacing along dimension 3
                typename std::enable_if<crank==3>::type* = 0
              );
   
   /// Set the sqaure-root of the PSD.
   /** This allocates _npsdSqrt and fills it with th evalues in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==1
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   template<size_t crank = rank>
   int psdSqrt( const realArrayT & npsdSqrt, ///< [in] an array containing the square root of the PSD.
                realT df,                    ///< [in] the frequency spacing
                typename std::enable_if<crank==1>::type* = 0
              );
   
   ///Set the sqaure-root of the PSD.
   /** This allocates _npsdSqrt and fills it with th evalues in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==2
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   template<size_t crank = rank>
   int psdSqrt( const realArrayT & npsdSqrt, ///< [in] an array containing the square root of the PSD.
                realT dk1,                   ///< [in] the frequency spacing along dimension 1
                realT dk2,                   ///< [in] the frequency spacing along dimension 2
                typename std::enable_if<crank==2>::type* = 0
              );
   
   ///Set the sqaure-root of the PSD.
   /** This allocates _npsdSqrt and fills it with th evalues in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==3
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */  
   template<size_t crank = rank>
   int psdSqrt( const realArrayT & npsdSqrt, ///< [in] an array containing the square root of the PSD.
                realT dk1,                   ///< [in] the frequency spacing along dimension 1
                realT dk2,                   ///< [in] the frequency spacing along dimension 2
                realT df,                    ///< [in] the frequency spacing along dimension 3
                typename std::enable_if<crank==3>::type* = 0
              );
   
   ///Set the sqaure-root of the PSD from the PSD.
   /** This allocates _npsdSqrt and fills it with the square root of the values in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==1
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   template<size_t crank=rank>
   int psd( const realArrayT & npsd, ///< [in] an array containing the PSD
            const realT df,          ///< [in] the frequency spacing
            typename std::enable_if<crank==1>::type* = 0
          );

   ///Set the sqaure-root of the PSD from the PSD.
   /** This allocates _npsdSqrt and fills it with the square root of the values in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==2
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   template<size_t crank=rank>
   int psd( const realArrayT & npsd, ///< [in] an array containing the PSD
            const realT dk1,         ///< [in] the frequency spacing along dimension 1
            const realT dk2,         ///< [in] the frequency spacing along dimension 2
            typename std::enable_if<crank==2>::type* = 0
          );
   
   ///Set the sqaure-root of the PSD from the PSD.
   /** This allocates _npsdSqrt and fills it with the square root of the values in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * This version compiles when rank==3
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */
   template<size_t crank=rank>
   int psd( const realArrayT & npsd, ///< [in] an array containing the PSD
            const realT dk1,         ///< [in] the frequency spacing along dimension 1
            const realT dk2,         ///< [in] the frequency spacing along dimension 2
            const realT df,          ///< [in] the frequency spacing along dimension 3
            typename std::enable_if<crank==3>::type* = 0
          );
   
   ///De-allocate all working memory and reset to initial state.
   /**
     *
     * \test Verify compilation and initialization of the 3 ranks for psdFilter. \ref tests_sigproc_psdFilter_compile "[test doc]" 
     */ 
   void clear();
   
   
   ///Apply the filter.
   /**
     * This version compiles when rank==1
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify filtering and noise normalization. \ref tests_sigproc_psdFilter_filter "[test doc]" 
     */ 
   template<size_t crank=rank>
   int filter( realArrayT & noise,             ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
               realArrayT * noiseIm = nullptr, ///< [out] [optional] an array to fill with the imaginary output of the filter, allowing 2-for-1 calculation.
               typename std::enable_if<crank==1>::type* = 0
             ) const;
   
   ///Apply the filter.
   /**
     * This version compiles when rank==2
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify filtering and noise normalization. \ref tests_sigproc_psdFilter_filter "[test doc]" 
     */ 
   template<size_t crank=rank>
   int filter( realArrayT & noise,             ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
               realArrayT * noiseIm = nullptr, ///< [out] [optional] an array to fill with the imaginary output of the filter, allowing 2-for-1 calculation.
               typename std::enable_if<crank==2>::type* = 0
             ) const;
             
   ///Apply the filter.
   /**
     * This version compiles when rank==3
     * 
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify filtering and noise normalization. \ref tests_sigproc_psdFilter_filter "[test doc]" 
     */ 
   template<size_t crank=rank>
   int filter( realArrayT & noise,             ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
               realArrayT * noiseIm = nullptr, ///< [out] [optional] an array to fill with the imaginary output of the filter, allowing 2-for-1 calculation.
               typename std::enable_if<crank==3>::type* = 0
             ) const;
             
   ///Apply the filter.
   /**
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify filtering and noise normalization. \ref tests_sigproc_psdFilter_filter "[test doc]" 
     */ 
   int operator()( realArrayT & noise /**< [in/out] the noise field of size rows() X cols(), which is filtered in-place. */ ) const;
   
   ///Apply the filter.
   /**
     * \returns 0 on success
     * \returns -1 on error
     * 
     * \test Verify filtering and noise normalization. \ref tests_sigproc_psdFilter_filter "[test doc]" 
     */ 
   int operator()( realArrayT & noise,  ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
                   realArrayT & noiseIm ///< [out] [optional] an array to fill with the imaginary output of the filter, allowing 2-for-1 calculation.
                 ) const; 
};

template<typename realT, size_t rank>
psdFilter<realT,rank>::psdFilter()
{
}

template<typename realT, size_t rank>
psdFilter<realT,rank>::~psdFilter()
{
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
   }
}

template<typename realT, size_t rank>
int psdFilter<realT,rank>::psdSqrt( realArrayT * npsdSqrt )
{
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
   }
   
   m_psdSqrt = npsdSqrt;
   m_owner = false;
   
   setSize();
   
   return 0;
}

template<typename realT, size_t rank>
int psdFilter<realT,rank>::psdSqrt( const realArrayT & npsdSqrt )
{
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
   }
   
   m_psdSqrt = new realArrayT;
   
   (*m_psdSqrt) = npsdSqrt;
   m_owner = true;
   
   setSize();
      
   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::setSize(typename std::enable_if<crank==1>::type* )
{
   if( m_psdSqrt == 0)
   {
      mxError("psdFilter", MXE_PARAMNOTSET, "m_psdSqrt has not been set yet, is still NULL.");
      return -1;
   }
   
   if( m_rows == m_psdSqrt->size() )
   {
      return 0;
   }
   
   m_rows = m_psdSqrt->size();
   m_cols = 1;
   m_planes = 1;
   
   m_ftWork.resize(m_rows);

   m_fft_fwd.plan(m_rows, MXFFT_FORWARD, true);
      
   m_fft_back.plan(m_rows, MXFFT_BACKWARD, true);      

   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::setSize(typename std::enable_if<crank==2>::type* )
{
   if( m_psdSqrt == 0)
   {
      mxError("psdFilter", MXE_PARAMNOTSET, "m_psdSqrt has not been set yet, is still NULL.");
      return -1;
   }
   
   if( m_rows == m_psdSqrt->rows() && m_cols == m_psdSqrt->cols())
   {
      return 0;
   }
   
   m_rows = m_psdSqrt->rows();
   m_cols = m_psdSqrt->cols();
   m_planes = 1;
   
   m_ftWork.resize(m_rows, m_cols);

   m_fft_fwd.plan(m_rows, m_cols, MXFFT_FORWARD, true);
      
   m_fft_back.plan(m_rows, m_cols, MXFFT_BACKWARD, true);      

   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::setSize(typename std::enable_if<crank==3>::type* )
{
   if( m_psdSqrt == 0)
   {
      mxError("psdFilter", MXE_PARAMNOTSET, "m_psdSqrt has not been set yet, is still NULL.");
      return -1;
   }
   
   if( m_rows == m_psdSqrt->rows() && m_cols == m_psdSqrt->cols() && m_planes == m_psdSqrt->planes())
   {
      return 0;
   }
   
   m_rows = m_psdSqrt->rows();
   m_cols = m_psdSqrt->cols();
   m_planes = m_psdSqrt->planes();
   
   m_ftWork.resize(m_rows, m_cols, m_planes);

   m_fft_fwd.plan(m_planes, m_rows, m_cols, MXFFT_FORWARD, true);
      
   m_fft_back.plan(m_planes, m_rows, m_cols, MXFFT_BACKWARD, true);      

   return 0;
}

template<typename realT, size_t rank>
int psdFilter<realT,rank>::rows()
{
   return m_rows;
}
   
template<typename realT, size_t rank>
int psdFilter<realT,rank>::cols()
{
   return m_cols;
}

template<typename realT, size_t rank>
int psdFilter<realT,rank>::planes()
{
   return m_planes;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psdSqrt( realArrayT * npsdSqrt,
                                    realT df,
                                    typename std::enable_if<crank==1>::type*
                                  )
{
   m_dFreq1 = df;
   return psdSqrt(npsdSqrt);
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psdSqrt( realArrayT * npsdSqrt,
                                    realT dk1,
                                    realT dk2,
                                    typename std::enable_if<crank==2>::type*
                                  )
{
   m_dFreq1 = dk1;
   m_dFreq2 = dk2;
   return psdSqrt(npsdSqrt);
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psdSqrt( realArrayT * npsdSqrt,
                                    realT dk1,
                                    realT dk2,
                                    realT df,
                                    typename std::enable_if<crank==3>::type*
                                  )
{
   m_dFreq1 = dk1;
   m_dFreq2 = dk2;
   m_dFreq3 = df;
   return psdSqrt(npsdSqrt);
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psdSqrt( const realArrayT & npsdSqrt,
                                    realT df,
                                    typename std::enable_if<crank==1>::type*
                                 )
{
   m_dFreq1 = df;
   return psdSqrt(npsdSqrt);
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psdSqrt( const realArrayT & npsdSqrt,
                                    realT dk1,
                                    realT dk2,
                                    typename std::enable_if<crank==2>::type*
                                 )
{
   m_dFreq1 = dk1;
   m_dFreq2 = dk2;
   return psdSqrt(npsdSqrt);
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psdSqrt( const realArrayT & npsdSqrt,
                                    realT dk1,
                                    realT dk2,
                                    realT df,
                                    typename std::enable_if<crank==3>::type*
                                 )
{
   m_dFreq1 = dk1;
   m_dFreq2 = dk2;
   m_dFreq3 = df;
   return psdSqrt(npsdSqrt);
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psd( const realArrayT & npsd,
                                const realT df1,
                                typename std::enable_if<crank==1>::type*
                              )
{
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
   }
   
   m_psdSqrt = new realArrayT;
   
   //Vector
   m_psdSqrt->resize(npsd.size());
   for(size_t n=0;n<npsd.size();++n) (*m_psdSqrt)[n] = sqrt(npsd[n]);
   
   m_owner = true;
   
   m_dFreq1 = df1;
   
   setSize();
      
   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psd( const realArrayT & npsd,
                                const realT dk1,
                                const realT dk2,
                                typename std::enable_if<crank==2>::type*
                              )
{
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
   }
   
   m_psdSqrt = new realArrayT;
   
   (*m_psdSqrt) = npsd.sqrt();
   m_owner = true;
   
   m_dFreq1 = dk1;
   m_dFreq2 = dk2;
   
   setSize();
      
   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::psd( const realArrayT & npsd,
                                const realT dk1,
                                const realT dk2,
                                const realT df,
                                typename std::enable_if<crank==3>::type*
                              )
{
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
   }
   
   m_psdSqrt = new realArrayT;
   
   //Cube
   m_psdSqrt->resize(npsd.rows(), npsd.cols(), npsd.planes());
   for(int pp=0;pp <npsd.planes();++pp) 
   {
      for(int cc=0; cc<npsd.cols(); ++cc)
      {
         for(int rr=0; rr<npsd.rows(); ++rr)
         {
            m_psdSqrt->image(pp)(rr,cc) = sqrt(npsd.image(pp)(rr,cc));
         }
      }
   }
   
   m_owner = true;
   
   m_dFreq1 = dk1;
   m_dFreq2 = dk2;
   m_dFreq3 = df;
   
   setSize();
      
   return 0;
}

template<typename realT, size_t rank>
void psdFilter<realT,rank>::clear()
{
   //m_ftWork.resize(0,0);
   psdFilterTypes::arrayT<realT,rank>::clear(m_ftWork);
   
   m_rows = 0;
   m_cols = 0;
   m_planes = 0;
   
   if(m_psdSqrt && m_owner)
   {
      delete m_psdSqrt;
      m_psdSqrt = 0;
   }
}
   
template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::filter( realArrayT & noise, 
                                   realArrayT * noiseIm,
                                   typename std::enable_if<crank==1>::type*
                                 ) const
{
   for(int nn=0; nn< noise.size(); ++nn) m_ftWork[nn] = complexT(noise[nn],0);
   
   //Transform complex noise to Fourier domain.
   m_fft_fwd(m_ftWork.data(), m_ftWork.data() );
   
   //Apply the filter.
   for(int nn=0;nn<m_ftWork.size();++nn) m_ftWork[nn] *= (*m_psdSqrt)[nn];
        
   m_fft_back(m_ftWork.data(), m_ftWork.data());
   
   //Now take the real part, and normalize.
   realT norm = sqrt(noise.size()/m_dFreq1);
   for(int nn=0;nn<m_ftWork.size();++nn) noise[nn] = m_ftWork[nn].real()/norm;
   
   if(noiseIm != nullptr)
   {
      for(int nn=0;nn<m_ftWork.size();++nn) (*noiseIm)[nn] = m_ftWork[nn].imag()/norm;
   }
   
   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::filter( realArrayT & noise, 
                                   realArrayT * noiseIm,
                                   typename std::enable_if<crank==2>::type*
                                 ) const
{
   //Make noise a complex number
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         m_ftWork(ii,jj) = complexT(noise(ii,jj),0);
      }
   }
   
   //Transform complex noise to Fourier domain.
   m_fft_fwd(m_ftWork.data(), m_ftWork.data() );
   
   //Apply the filter.
   m_ftWork *= *m_psdSqrt;
        
   m_fft_back(m_ftWork.data(), m_ftWork.data());
   
   realT norm = sqrt(noise.rows()*noise.cols()/(m_dFreq1*m_dFreq2));
   
   //Now take the real part, and normalize.
   noise = m_ftWork.real()/norm;
   
   if(noiseIm != nullptr)
   {
      *noiseIm = m_ftWork.imag()/norm;
   }
   
   return 0;
}

template<typename realT, size_t rank>
template<size_t crank>
int psdFilter<realT,rank>::filter( realArrayT & noise, 
                                   realArrayT * noiseIm,
                                   typename std::enable_if<crank==3>::type*
                                 ) const
{
   //Make noise a complex number
   for(int pp=0;pp<noise.planes();++pp)
   {
      for(int ii=0;ii<noise.rows();++ii)
      {
         for(int jj=0; jj<noise.cols(); ++jj)
         {
            m_ftWork.image(pp)(ii,jj) = complexT(noise.image(pp)(ii,jj),0);
         }
      }
   }
   
   //Transform complex noise to Fourier domain.
   m_fft_fwd(m_ftWork.data(), m_ftWork.data() );
   
   //Apply the filter.
   for(int pp=0;pp<noise.planes();++pp)  m_ftWork.image(pp) *= m_psdSqrt->image(pp);
        
   m_fft_back(m_ftWork.data(), m_ftWork.data());
   
   //Now take the real part, and normalize.
   
   realT norm = sqrt(m_rows*m_cols*m_planes/(m_dFreq1*m_dFreq2*m_dFreq3));
   for(int pp=0; pp< noise.planes();++pp)  noise.image(pp) = m_ftWork.image(pp).real()/norm;
   
   if(noiseIm != nullptr)
   {
      for(int pp=0; pp< noise.planes();++pp) noiseIm->image(pp) = m_ftWork.image(pp).imag()/norm;
   }
   
   return 0;
}

template<typename realT, size_t rank>
int psdFilter<realT,rank>::operator()( realArrayT & noise ) const
{
   return filter(noise);
}

template<typename realT, size_t rank>
int psdFilter<realT,rank>::operator()( realArrayT & noise,
                                       realArrayT & noiseIm
                                     ) const
{
   return filter(noise, &noiseIm);
}

} //namespace sigproc 
} //namespace mx

#endif //psdFilter_hpp
