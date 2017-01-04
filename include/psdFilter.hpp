/** \file psdFilter.hpp
  * \brief Declares and defines a class for filtering with PSDs
  * \ingroup signal_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */


#ifndef __psdFilter_hpp__
#define __psdFilter_hpp__

#include <Eigen/Dense>

#include "mxError.hpp"
#include "fft.hpp"

namespace mx
{
   
/// A class for filtering with PSDs
/** \ingroup psds
  *
  * \todo investigate whether it is really necessary to use the backward transform. 
  */
template<typename _realT>
class psdFilter
{
public:
   
   typedef _realT realT; ///< Real floating point type
   typedef std::complex<_realT> complexT; ///< Complex floating point type.
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> realArrayT; ///< Eigen array type with Scalar==realT 
   typedef Eigen::Array<complexT, Eigen::Dynamic, Eigen::Dynamic> complexArrayT; ///< Eigen array type with Scalar==complexT
   
   
protected:

   int _rows; ///< The number of rows in the filter, and the required number of rows in the noise array.
   int _cols; ///< The number of columns in the filter, and the required number of columns in the noise array.

   realArrayT * _psdSqrt; ///< Pointer to the real array containing the square root of the PSD.
   bool _owner; ///< Flag indicates whether or not _psdSqrt was allocated by this instance, and so must be deallocated.

   complexArrayT _ftWork;   ///< Working memory for the FFT.
   
   fftT< complexT, complexT,2,0> fft_fwd; ///< FFT object for the forward transform.
   fftT< complexT, complexT,2,0> fft_back; ///< FFT object for the backward transfsorm.
   

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
     * \returns the current value of _rows.
     */ 
   int rows();
   
   ///Get the number of columns in the filter
   /**
     * \returns the current value of _cols.
     */
   int cols();
   
   ///Set the sqaure-root of the PSD to be a pointer to an array containing the square root of the PSD.
   /** This does not allocate _npsdSqrt, it merely points to the specified array, which remains your responsibility for deallocation, etc.
     *
     * \param [in] npsdSqrt a pointer to an array containing the square root of the PSD.
     */  
   int psdSqrt( realArrayT * npsdSqrt );
   
   ///Set the sqaure-root of the PSD.
   /** This allocates _npsdSqrt and fills it with th evalues in the array.
     *
     * \param [in] npsdSqrt an array containing the square root of the PSD.
     */  
   int psdSqrt( const realArrayT & npsdSqrt );
   
   ///Set the sqaure-root of the PSD from the PSD.
   /** This allocates _npsdSqrt and fills it with the square root of the values in the array.
     *
     * \param [in] npsdSqrt an array containing the PSD.
     */  
   int psd( const realArrayT & npsd );

   ///De-allocate all working memory and reset to initial state.
   void clear();
   
   
   ///Apply the filter.
   /**
     * \param [in,out] noise is the noise field of the same size as _sqrtPsd, which is filtered.
     */ 
   int filter( realArrayT & noise );
   
   ///Apply the filter.
   /**
     * \param [in,out] noise is the noise field of the same size as _sqrtPsd, which is filtered.
     */ 
   int operator()( realArrayT & noise );
   
};

template<typename realT>
psdFilter<realT>::psdFilter()
{
   _rows = 0;
   _cols = 0;
   
   _psdSqrt = 0;
   _owner = false;
}

template<typename realT>
psdFilter<realT>::~psdFilter()
{
   if(_psdSqrt && _owner)
   {
      delete _psdSqrt;
   }
}

template<typename realT>
int psdFilter<realT>::setSize()
{
   if( _psdSqrt == 0)
   {
      mxError("psdFilter", MXE_PARAMNOTSET, "_psdSqrt has not been set yet, is still NULL.");
      return -1;
   }
   
   if( _rows == _psdSqrt->rows() && _cols == _psdSqrt->cols())
   {
      return 0;
   }
   
   _rows = _psdSqrt->rows();
   _cols = _psdSqrt->cols();
   
   _ftWork.resize(_rows, _cols);

   fft_fwd.plan(_rows, _cols, MXFFT_FORWARD, true);
      
   fft_back.plan(_rows, _cols, MXFFT_BACKWARD, true);      

   return 0;
}



template<typename realT>
int psdFilter<realT>::rows()
{
   return _rows;
}
   
template<typename realT>
int psdFilter<realT>::cols()
{
   return _cols;
}
   
template<typename realT>
int psdFilter<realT>::psdSqrt( realArrayT * npsdSqrt )
{
   if(_psdSqrt && _owner)
   {
      delete _psdSqrt;
   }
   
   _psdSqrt = npsdSqrt;
   _owner = false;
   
   setSize();
   
   return 0;
}

template<typename realT>
int psdFilter<realT>::psdSqrt( const realArrayT & npsdSqrt )
{
   if(_psdSqrt && _owner)
   {
      delete _psdSqrt;
   }
   
   _psdSqrt = new realArrayT;
   
   (*_psdSqrt) = npsdSqrt;
   _owner = true;
   
   setSize();
      
   return 0;
}

template<typename realT>
int psdFilter<realT>::psd( const realArrayT & npsd )
{
   if(_psdSqrt && _owner)
   {
      delete _psdSqrt;
   }
   
   _psdSqrt = new realArrayT;
   
   (*_psdSqrt) = npsd.sqrt();
   _owner = true;
   
   setSize();
      
   return 0;
}

template<typename realT>
void psdFilter<realT>::clear()
{
   _ftWork.resize(0,0);
   _rows = 0;
   _cols = 0;
      
   if(_psdSqrt && _owner)
   {
      delete _psdSqrt;
      _psdSqrt = 0;
   }
}
   
template<typename realT>
int psdFilter<realT>::filter( realArrayT & noise )
{
   //Make noise a complex number
   for(int ii=0;ii<noise.rows();++ii)
   {
      for(int jj=0; jj<noise.cols(); ++jj)
      {
         _ftWork(ii,jj) = complexT(noise(ii,jj),0);
      }
   }
   
   //Transform complex noise to Fourier domain.
   fft_fwd(_ftWork.data(), _ftWork.data() );
   
   //Apply the filter.
   _ftWork *= *_psdSqrt;
        
   fft_back(_ftWork.data(), _ftWork.data());
   
   //Now take the real part, and normalize.
   noise = _ftWork.real()/(noise.rows()*noise.cols());
   
   return 0;
}

template<typename realT>
int psdFilter<realT>::operator()( realArrayT & noise )
{
   return filter(noise);
}

   
}; //namespace mx

#endif //__psdFilter_hpp__
