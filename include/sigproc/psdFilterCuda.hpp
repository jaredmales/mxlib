/** \file psdFilterCuda.hpp
  * \brief Declares and defines a class for filtering with PSDs on CUDA GPUs
  * \ingroup signal_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef psdFilterCuda_hpp
#define psdFilterCuda_hpp

#include "../cuda/templateCudaPtr.hpp"
#include "../cuda/templateCuda.hpp"

#include "psdFilter.hpp"

#include <cufft.h>
#include <cufftXt.h>
//#include <helper_functions.h>
#include <helper_cuda.h>


namespace mx
{
namespace sigproc 
{


/// A class for filtering noise with PSDs using CUDA cuFFT
/** The square-root of the PSD is maintained by this class, either as a pointer to an external array or using internally allocated memory (which will be
  * de-allocated on destruction. 
  * 
  * PSD Normalization: the PSD used for this needs to be normalized properly to produce filtered noise with the correct statistics.  Given an
  * array of size N which contains the PSD as "power per unit frequency" vs frequency (i.e. the PSD) on some frequency scale with uniform spacing \f$ \Delta f \f$,
  * the normalization is
  \f[
       PSD_{norm} = PSD * N * \Delta f 
  \f]
  * for a 1-dimensional PSD and  
  \f[
       PSD_{norm} = PSD * N_0 * \Delta f_0 * N_1 * \Delta f_1 
  \f]
  * for a 2-dimensional PSD. Remember that these are applied before taking the square root.
  *
  * 
  * \ingroup psd_filter
  *
  * \todo once fftT has a plan interface with pointers for working memory, use it.
  */
template<typename _realT>
class psdFilter<_realT, 1>
{
public:
   
   typedef _realT realT; ///< Real floating point type
   typedef std::complex<_realT> complexT; ///< Complex floating point type.
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> realArrayT; ///< Eigen array type with Scalar==realT 
   typedef Eigen::Array<complexT, Eigen::Dynamic, Eigen::Dynamic> complexArrayT; ///< Eigen array type with Scalar==complexT
   
   typedef realT deviceRealPtrT;
   typedef complexT deviceComplexPtrT;
   
protected:

   int m_rows {0}; ///< The number of rows in the filter, and the required number of rows in the noise array.
   int m_cols {0}; ///< The number of columns in the filter, and the required number of columns in the noise array.

   mx::cuda::cudaPtr<complexT> m_psdSqrt; ///< Pointer to the real array containing the square root of the PSD.

   
   mx::cuda::cudaPtr<complexT> m_scale; ///< the scale factor.
   
   ///Cuda FFT plan.  We only need one since the forward/inverse is part of execution.
   cufftHandle m_fftPlan {0};
   
public:
   
   ///C'tor.
   psdFilter();
   
   ///Destructor
   ~psdFilter();
   
protected:

   ///Set the size of the filter.
   /** Handles allocation of the _ftWork array and fftw planning.
     *
     * Requires _psdSqrt to be set first.  This is called by the psdSqrt() and psd() methods.
     * 
     */
   int setSize();

public:   

   ///Get the number of rows in the filter
   /**
     * \returns the current value of m_rows.
     */ 
   int rows();
   
   ///Get the number of columns in the filter
   /**
     * \returns the current value of m_cols.
     */
   int cols();
   
   ///Set the sqaure-root of the PSD.
   /** This allocates _npsdSqrt and fills it with th evalues in the array.
     *
     * See the discussion of PSD normalization above.
     * 
     * \returns 0 on success
     * \returns -1 on error
     */  
   int psdSqrt( const realArrayT & npsdSqrt /**< [in] an array containing the square root of the PSD.*/ );
   
   int psdSqrt( const cuda::cudaPtr<std::complex<realT>> & npsdSqrt, /**< [in] an array containing the square root of the PSD.*/ 
                size_t rows,
                size_t cols
              );
   
   ///De-allocate all working memory and reset to initial state.
   void clear();
      
   ///Apply the filter.
   /**
     * 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int filter( deviceComplexPtrT * noise           ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
             );
   
   ///Apply the filter.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int operator()( realArrayT & noise /**< [in/out] the noise field of size rows() X cols(), which is filtered in-place. */ );
   
   ///Apply the filter.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int operator()( realArrayT & noise,  ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
                   realArrayT & noiseIm ///< [out] [optional] an array to fill with the imaginary output of the filter, allowing 2-for-1 calculation.
                 ); 
};

template<typename realT>
psdFilter<realT,1>::psdFilter()
{
}

template<typename realT>
psdFilter<realT,1>::~psdFilter()
{
   if( m_fftPlan )
   {
      checkCudaErrors(cufftDestroy(m_fftPlan));
   }
}

template<typename realT>
int psdFilter<realT,1>::setSize()
{
   
   m_scale.resize(1);
   std::complex<realT> scale(1./(m_rows*m_cols),0);
   m_scale.upload(&scale,1);
   
   checkCudaErrors(cufftPlan2d(&m_fftPlan, m_rows, m_cols, CUFFT_C2C));
   
   return 0;
}



template<typename realT>
int psdFilter<realT,1>::rows()
{
   return m_rows;
}
   
template<typename realT>
int psdFilter<realT,1>::cols()
{
   return m_cols;
}
   
// template<typename realT>
// int psdFilter<realT,1>::psdSqrt( realArrayT * npsdSqrt )
// {
//    if(_psdSqrt && _owner)
//    {
//       delete _psdSqrt;
//    }
//    
//    _psdSqrt = npsdSqrt;
//    _owner = false;
//    
//    setSize();
//    
//    return 0;
// }

template<typename realT>
int psdFilter<realT,1>::psdSqrt( const realArrayT & npsdSqrt )
{
   m_rows = npsdSqrt.rows();
   m_cols = npsdSqrt.cols();
   
   
   complexArrayT tmp( m_rows, m_cols);
   tmp.setZero();
   tmp.real() = npsdSqrt;
   
   m_psdSqrt.upload(tmp.data(), m_rows*m_cols);
   
   setSize();
      
   return 0;
}

template<typename realT>
int psdFilter<realT,1>::psdSqrt( const cuda::cudaPtr<std::complex<realT>> & npsdSqrt,
                                 size_t rows,
                                 size_t cols
                               )
{
   m_rows = rows;
   m_cols = cols;
   
   m_psdSqrt.resize(m_rows*m_cols);
      
   ///\todo move this into cudaPtr
   cudaMemcpy( m_psdSqrt.m_devicePtr, npsdSqrt.m_devicePtr, m_rows*m_cols*sizeof(std::complex<realT>), cudaMemcpyDeviceToDevice);
   
   setSize();
      
   return 0;
}

template<typename realT>
void psdFilter<realT,1>::clear()
{
   m_scale.resize(0);
   m_rows = 0;
   m_cols = 0;
      
   m_psdSqrt.resize(0);
}
   
template<typename realT>
int psdFilter<realT,1>::filter( deviceComplexPtrT * noise )
{
   
   //Transform complex noise to Fourier domain.
   
   cufftExecC2C(m_fftPlan, (cuComplex *) noise, (cuComplex *) noise, CUFFT_FORWARD);    
      
   //Apply the filter.
   mx::cuda::pointwiseMul<cuComplex><<<32, 256>>>((cuComplex *) noise, (cufftComplex *) m_psdSqrt.m_devicePtr, m_rows*m_cols);
        
   cufftExecC2C(m_fftPlan, (cuComplex *)noise, (cuComplex *)noise, CUFFT_INVERSE);
   
   
   mx::cuda::scalarMul<cuComplex><<<32, 256>>>((cuComplex *) noise, (cuComplex *) m_scale.m_devicePtr, m_rows*m_cols);
   //Now take the real part, and normalize.
   
   //noise = noise/(noise.rows()*noise.cols());
   
   
   return 0;
}

template<typename realT>
int psdFilter<realT,1>::operator()( realArrayT & noise )
{
   return filter(noise);
}

template<typename realT>
int psdFilter<realT,1>::operator()( realArrayT & noise,
                                  realArrayT & noiseIm
                                )
{
   return filter(noise, &noiseIm);
}

} //namespace sigproc 
} //namespace mx

#endif //psdFilterCuda_hpp
