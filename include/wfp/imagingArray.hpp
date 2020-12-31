/** \file imagingArray.hpp
  * \brief Declares and defines a class for managing images 
  * \ingroup imaging
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */


#ifndef imagingArray_hpp
#define imagingArray_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include "../math/templateBLAS.hpp"
#include "../math/fft/fft.hpp"


namespace mx
{
namespace wfp
{
   
template<typename Scalar>
struct fftwAllocator
{
   Scalar * alloc(size_t sz)
   {
      return (Scalar *) ::fftw_malloc(sz * sizeof(Scalar));
   }

   void free(Scalar * ptr)
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
      math::scal<typename imageT::Scalar>(im.cols()*im.rows(), (typeT) alpha, im.data(), 1);        
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
      typedef typename imageT::Scalar scalarT;
      //hadp<scalarT>(im.cols()*im.rows(), alpha.data(), im.data()); 
      
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigY(alpha.data(), alpha.cols(), alpha.rows());
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigX(im.data(), im.cols(), im.rows());
      
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
      math::scal<typename imageT::Scalar>(im.cols()*im.rows(), ((typeT) (1))/alpha, im.data(), 1);   
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
      typedef typename imageT::Scalar scalarT;
      math::hadd<scalarT>(im.cols()*im.rows(), (scalarT) 1.0, alpha.data(),1, im.data(), 1); 
   }
};



template< typename imageT, 
          typename typeT,  
          bool isScalar = is_imagingArray<typeT>::value>
struct imagingArrayInplaceAdd
{
   imagingArrayInplaceAdd()
   {
      static_assert( is_imagingArray<imageT>::value, "imagingArrayInplaceAdd template parameter 1 must be an imagingArray type");
   }
   
   void operator()(imageT & im, typeT & alpha)
   {
      typedef typename imageT::Scalar scalarT;
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigX(im.data(), im.cols(), im.rows());
      eigX += alpha;
   }
};

template< typename imageT, 
          typename typeT>
struct imagingArrayInplaceAdd<imageT, typeT, true>
{
   imagingArrayInplaceAdd()
   {
      static_assert( is_imagingArray<imageT>::value, "imagingArrayInplaceAdd template parameter 1 must be an imagingArray type");
   }
   
   void operator()(imageT & im, typeT & alpha)
   {
      typedef typename imageT::Scalar scalarT;
      //hadp<scalarT>(im.cols()*im.rows(), alpha.data(), im.data()); 
      
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigY(alpha.data(), alpha.cols(), alpha.rows());
      Eigen::Map<Eigen::Array<scalarT,-1,-1> > eigX(im.data(), im.cols(), im.rows());
      
      eigX += eigY;
      
   }
};



template<typename _Scalar, class _allocatorT, int cudaGPU>
class imagingArray;

template<typename _Scalar, class _allocatorT>
class imagingArray<_Scalar, _allocatorT, 0>
{
public:
   
   typedef bool imagingArrayT;
   typedef _Scalar Scalar;
   typedef _allocatorT allocatorT;
   
protected:
   allocatorT allocator;
   
   Scalar * _data;

   bool _owner;
   
   int _cols;
   int _rows;

public:

   imagingArray()
   {
      _data = 0;
      _owner = false;
      _cols = 0;
      _rows = 0;
   }

   imagingArray(int x, int y)
   {
      _data = 0;
      _owner=false;
      _cols = 0;
      _rows = 0;
      
      resize(x, y);
   }
   
   template<typename eigenT>
   explicit imagingArray(eigenT & im)
   {
      _data = im.data();
      _owner=false;
      _cols = im.cols();
      _rows = im.rows();
   }
   
   ~imagingArray()
   {
      free();
   }
   
   void resize(int szx, int szy)
   {
      if(szx == _cols && szy == _rows && _owner) return;
      
      free();

      _cols = szx;
      _rows = szy;

      _data = allocator.alloc(_cols*_rows);
      _owner = true;
   }

   void free()
   {
      if(_data && _owner)
      {
         allocator.free(_data);
         _data = 0;
         _cols = 0;
         _rows = 0;
      }
   }

   int cols() const
   {
      return _cols;
   }
      
   int rows() const
   {
      return _rows;
   }
      
   Scalar * data()
   {
      return _data;
   }
   
   Scalar * data() const
   {
      return _data;
   }
   
   Scalar operator()(int x, int y) const
   {
      return _data[x + y*_cols];
   }
   
   Scalar & operator()(int x, int y)
   {
      return _data[x + y*_cols];
   }
   
   template<typename argT>
   void set(argT val)
   {
      int N = _cols*_rows;
      
      Scalar dval = (Scalar) val;
      for(int i=0;i<N;++i) _data[i] = dval;
      
   }
   
   void setZero()
   {
      set( (Scalar) 0);
   }
   
   
   imagingArray & operator=(const imagingArray & arr)
   {
      resize(arr.rows(), arr.cols());
      
      for(int i=0; i<_rows; ++i)
      {
         for(int j=0; j<_cols; ++j)
         {
            operator()(i,j) = arr(i,j);
         }
      }
      
      return (*this);
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
   
   ///In-place add.  Works for both scalars and imagingArray types.
   template<typename typeT>
   imagingArray & operator+=(typeT & alpha)
   {
      imagingArrayInplaceAdd<imagingArray, typeT> prod;
      prod(*this, alpha);
      return *this;
   }
   
   Scalar sum()
   {
      int N = _cols*_rows;
      Scalar _sum = 0;
      
      for(int i=0;i<N;++i) _sum += _data[i];
      
      return _sum;
   }
   
   Eigen::Map<Eigen::Array<Scalar,-1,-1> >  eigenMap()
   {
      return Eigen::Map<Eigen::Array<Scalar,-1,-1> >(_data, _rows, _cols);
   }
};


   
template<typename scalarT>
using imagingArrayT = imagingArray<scalarT, fftwAllocator<scalarT>,0>;

// typedef imagingArray<float,mx::fftwAllocator<float>,0> imagingArrayf;
// typedef imagingArray<double,mx::fftwAllocator<double>,0> imagingArrayd;
// typedef imagingArray<std::complex<float>,mx::fftwAllocator<std::complex<float> >,0> imagingArraycf;

} //namespace wfp
} //namespace mx


#endif //imagingArray_hpp
