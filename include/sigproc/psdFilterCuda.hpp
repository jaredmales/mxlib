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

#include "psdfilter.hpp"

namespace mx
{
namespace sigproc 
{

template<typename _realT, int cuda = 0>
class psdFilter;

/// A class for filtering noise with PSDs
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
  * \ingroup psds
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

   
   complexArrayT m_ranWork; ///< Working memory for the noise field on host.
   mx::cuda::cudaPtr<complexT> m_ftWork; ///< Working memory for the FFT on the device.
   
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
   
   ///De-allocate all working memory and reset to initial state.
   void clear();
      
   ///Apply the filter.
   /**
     * 
     * \returns 0 on success
     * \returns -1 on error
     */ 
   int filter( realArrayT & noise,            ///< [in/out] the noise field of size rows() X cols(), which is filtered in-place. 
               realArrayT * noiseIm = nullptr ///< [out] [optional] an array to fill with the imaginary output of the filter, allowing 2-for-1 calculation.
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
psdFilter<realT>::psdFilter()
{
}

template<typename realT>
psdFilter<realT>::~psdFilter()
{
   if( m_fftPlan )
   {
      checkCudaErrors(cufftDestroy(m_fftPlan));
   }
}

template<typename realT>
int psdFilter<realT>::setSize()
{
   if( m_rows * m_cols == m_ftWork.size())
   {
      return 0;
   }
   
   m_ranWork.resize(m_rows, m_cols);
   m_ftWork.resize(m_rows*m_cols);

   checkCudaErrors(cufftPlan2d(&m_fftPlan, m_rows, m_cols, CUFFT_C2C));
   
   return 0;
}



template<typename realT>
int psdFilter<realT>::rows()
{
   return m_rows;
}
   
template<typename realT>
int psdFilter<realT>::cols()
{
   return m_cols;
}
   
// template<typename realT>
// int psdFilter<realT>::psdSqrt( realArrayT * npsdSqrt )
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
int psdFilter<realT>::psdSqrt( const realArrayT & npsdSqrt )
{
   if(m_psdSqrt && m_owner)
   {
      cudaFree(m_psdSqrt);
   }
   
   m_rows = npsdSqrt.rows();
   m_cols = npsdSqrt.cols();
   
   
   m_psdSqrt.upload(npsdSqrt.data(), m_rows*m_cols);
   
   setSize();
      
   return 0;
}

template<typename realT>
void psdFilter<realT>::clear()
{
   m_ranWork.resize(0);
   m_ftWork.resize(0);
   m_rows = 0;
   m_cols = 0;
      
   m_psdSqrt.resize(0);
}
   
template<typename realT>
int psdFilter<realT>::filter( realArrayT & noise, realArrayT * noiseIm )
{
   //Make noise a complex number
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         m_ranWork(ii,jj) = complexT(noise(ii,jj),0);
      }
   }
   
   m_ftWork.upload(m_ranWork, m_rows*m_cols);
   
   
   //Transform complex noise to Fourier domain.
   
   cufftExecC2C(m_fftPlan, (cufftComplex *) complexPupil, (cufftComplex *) complexFocal, CUFFT_FORWARD);    
   
   fft_fwd(_ftWork.data(), _ftWork.data() );
   
   //Apply the filter.
   _ftWork *= *_psdSqrt;
        
   fft_back(_ftWork.data(), _ftWork.data());
   
   //Now take the real part, and normalize.
   noise = _ftWork.real()/(noise.rows()*noise.cols());
   
   if(noiseIm != nullptr)
   {
      *noiseIm = _ftWork.imag()/(noise.rows()*noise.cols());
   }
   
   return 0;
}

template<typename realT>
int psdFilter<realT>::operator()( realArrayT & noise )
{
   return filter(noise);
}

template<typename realT>
int psdFilter<realT>::operator()( realArrayT & noise,
                                  realArrayT & noiseIm
                                )
{
   return filter(noise, &noiseIm);
}

} //namespace sigproc 
} //namespace mx

#endif //psdFilterCuda_hpp
