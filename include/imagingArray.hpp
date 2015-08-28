/** \file imagingArray.hpp
  * \brief Declares and defines a class for managing images 
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */


#ifndef __imagingArray_hpp__
#define __imagingArray_hpp__

#include "templateBLAS.hpp"
#include "fft.hpp"

#include <Eigen/Dense>
namespace mx
{
   
template<typename dataT>
struct fftwAllocator
{
   dataT * alloc(size_t sz)
   {
      return (dataT *) fftw_malloc(sz * sizeof(dataT));
   }

   void free(dataT * ptr)
   {
      fftw_free(ptr);
   }
};

  
///Test whether a type is an imagingArray by testing whether it has a typedef of "is_imagingArray"
/** Used for compile-time determination of type
  * Example usage:
  * \code
  * bool is_eC = is_imagingArray<eigenCube<float> >; //Evaluates to true
  * bool is_not_eC = is_imagingArray<float>; //Evaluates to false
  * \endcode
  */
//This was taken directly from the example at http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
template <typename T>
struct is_imagingArray
{
   // Types "yes" and "no" are guaranteed to have different sizes,
   // specifically sizeof(yes) == 1 and sizeof(no) == 2.
   typedef char yes[1];
   typedef char no[2];
 
   template <typename imageT>
   static yes& test(typename imageT::imagingArrayT*);
 
   template <typename>
   static no& test(...);
 
   // If the "sizeof" of the result of calling test<T>(0) would be equal to sizeof(yes),
   // the first overload worked and T has a nested type named "is_mmatrix".
   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};

template< typename imageT, 
          typename typeT,  
          bool isScalar = is_imagingArray<typeT>::value>
struct imagingArrayInplaceProduct
{
   imagingArrayInplaceProduct()
   {
      static_assert( is_imagingArray<imageT>::value, "imagingArrayInplaceProduct template parameter 1 must be an imagingArray type");
   }
   
   void operator()(imageT & im, typeT & alpha)
   {
      scal<typename imageT::dataT>(im.szX()*im.szY(), (typeT) alpha, im.data(), 1);        
   }
};

template< typename imageT, 
          typename typeT>
struct imagingArrayInplaceProduct<imageT, typeT, true>
{
   imagingArrayInplaceProduct()
   {
      static_assert( is_imagingArray<imageT>::value, "imagingArrayInplaceProduct template parameter 1 must be an imagingArray type");
   }
   
   void operator()(imageT & im, typeT & alpha)
   {
      typedef typename imageT::dataT scalarT;
      //hadp<scalarT>(im.szX()*im.szY(), alpha.data(), im.data()); 
      
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigY(alpha.data(), alpha.szX(), alpha.szY());
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigX(im.data(), im.szX(), im.szY());
      
      eigX *= eigY;
      
   }
};


template< typename imageT, 
          typename typeT,  
          bool isScalar = is_imagingArray<typeT>::value>
struct imagingArrayInplaceDivision
{
   imagingArrayInplaceDivision()
   {
      static_assert( is_imagingArray<imageT>::value, "imagingArrayInplaceDivision template parameter 1 must be an imagingArray type");
   }
   
   void operator()(imageT im, typeT alpha)
   {
      scal<typename imageT::dataT>(im.szX()*im.szY(), ((typeT) (1))/alpha, im.data(), 1);   
   }
};

template< typename imageT, 
          typename typeT>
struct imagingArrayInplaceDivision<imageT, typeT, true>
{
   imagingArrayInplaceDivision()
   {
      static_assert( is_imagingArray<imageT>::value, "imagingArrayInplaceDivision template parameter 1 must be an imagingArray type");
   }
   
   void operator()(imageT im, typeT alpha)
   {
      typedef typename imageT::dataT scalarT;
      hadd<scalarT>(im.szX()*im.szY(), (scalarT) 1.0, alpha.data(),1, im.data(), 1); 
   }
};


template<typename _dataT, class _allocatorT, int cudaGPU>
class imagingArray;

template<typename _dataT, class _allocatorT>
class imagingArray<_dataT, _allocatorT, 0>
{
public:
   
   typedef bool imagingArrayT;
   typedef _dataT dataT;
   typedef _allocatorT allocatorT;
   
protected:
   allocatorT allocator;
   
   dataT * _data;

   int _szX;
   int _szY;

public:

   imagingArray()
   {
      _data = 0;
      _szX = 0;
      _szY = 0;
   }

   ~imagingArray()
   {
      free();
   }
   
   void resize(int szx, int szy)
   {
      if(szx == _szX && szy == _szY) return;
      
      free();

      _szX = szx;
      _szY = szy;

      _data = allocator.alloc(_szX*_szY);
   }

   void free()
   {
      if(_data)
      {
         allocator.free(_data);
         _data = 0;
         _szX = 0;
         _szY = 0;
      }
   }

   int szX()
   {
      return _szX;
   }
      
   int szY()
   {
      return _szY;
   }
      
   dataT * data()
   {
      return _data;
   }
   
   dataT & operator()(int x, int y)
   {
      return _data[x + y*_szX];
   }
   
   template<typename argT>
   void set(argT val)
   {
      int N = _szX*_szY;
      
      dataT dval = (dataT) val;
      for(int i=0;i<N;++i) _data[i] = dval;
      
   }
   
   ///In-place product.  Works for both scalars and imagingArray types.
   template<typename typeT>
   imagingArray & operator*=(typeT & alpha)
   {
      imagingArrayInplaceProduct<imagingArray, typeT> prod;
      prod(*this, alpha);
      return *this;
   }

   ///In-place division.  Works for both scalars and imagingArray types.
   template<typename typeT>
   imagingArray & operator/=(typeT alpha)
   {
      imagingArrayInplaceDivision<imagingArray, typeT> div;
      div(*this, alpha);
      return *this;
   }
   
   dataT sum()
   {
      int N = _szX*_szY;
      dataT _sum = 0;
      
      for(int i=0;i<N;++i) _sum += _data[i];
      
      return _sum;
   }
};






} //namespace mx


#endif //__imagingArray_hpp__
