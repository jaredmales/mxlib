/** \file eigenCube.hpp
  * \brief An image cube with an Eigen API
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup image_processing
  *
  */

#ifndef __eigenCube_hpp__
#define __eigenCube_hpp__

#pragma GCC system_header
#include <Eigen/Dense>


#include "eigenUtils.hpp"
#include "vectorUtils.hpp"

using namespace Eigen;

namespace mx
{
   
/// An image cube with an Eigen API
/** \ingroup image_processing
  */
template<typename dataT>
class eigenCube 
{
public:
   typedef bool is_eigenCube;
   typedef dataT Scalar;

   typedef typename Array<dataT,Dynamic,Dynamic>::Index Index;
   typedef Map<Array<dataT, Dynamic, Dynamic> > imageRef;
   
protected:

   
   Index _rows;
   Index _cols;
   Index _planes;

   dataT * _data;

   bool _owner;

public:

   eigenCube()
   {
      _rows = 0;
      _cols = 0;
      _planes = 0;
      _data = 0;
      _owner = false;
   }

   eigenCube(Index nrows, Index ncols, Index nplanes)
   {
      _rows = nrows;
      _cols = ncols;
      _planes = nplanes;

      _data = new dataT[_rows*_cols*_planes];
      _owner = true;
   }

   eigenCube(dataT * ndata, size_t nrows, size_t ncols, size_t nplanes)
   {
      _rows = nrows;
      _cols = ncols;
      _planes = nplanes;

      _data = ndata;
      _owner = false;
   }

   ~eigenCube()
   {
      if(_owner && _data)
      {
         delete _data;
      }
   }

   void setZero()
   {
      int N = _rows*_cols*_planes;
      
      for(int i=0;i<N;++i) _data[i] = ((dataT) 0);
   }
   
   eigenCube<dataT> & operator=(eigenCube<dataT> & ec)
   {
      resize(ec.rows(), ec.cols(), ec.planes());
      
      int N = _rows*_cols*_planes;
      
      for(int i=0;i<N;++i) _data[i] = ec._data[i];
      
      return *this;
   }
   
   void shallowCopy(eigenCube<dataT> & src, bool takeOwner = false)
   {
      if(_owner && _data)
      {
         delete _data;
         _owner = false;
      }
      
      _rows = src._rows;
      _cols = src._cols;
      _planes = src._planes;
      _data = src._data;
      
      if(takeOwner == true) 
      {
         _owner = true;
         src._owner = false;
      }
      else
      { 
         _owner = false;
      }
   }
   
   void clear()
   {
      if(_owner && _data)
      {
         delete _data;
      }
      
      _rows = 0;
      _cols = 0;
      _planes = 0;
   }
   
   void resize(int r, int c, int p)
   {
      if(_owner && _data)
      {
         delete _data;
      }
      
      _rows = r;
      _cols = c;
      _planes = p;
     
      _data = new dataT[_rows*_cols*_planes];
      _owner = true;
     
   }
   
   void resize(int r, int c)
   {
      resize(r, c, 1);
   }
   
   dataT * data()
   {
      return _data;
   }
   
   const dataT * data() const
   {
      return _data;
   }
   
   const Index rows() const
   {
      return _rows;
   }
   
   const Index cols() const
   {
      return _cols;
   }
   
   const Index planes() const
   {
      return _planes;
   }

   /// Returns a 2D Eigen::Map pointed at the entire cube.
   Map<Array<dataT, Dynamic, Dynamic> > cube()
   {
      return Map<Array<dataT, Dynamic, Dynamic> >(_data, _rows*_cols, _planes);
   }
   
   /// Returns a 2D Eigen::Map pointed at the specified image.
   Map<Array<dataT, Dynamic, Dynamic> > image(Index n)
   {
      return Map<Array<dataT, Dynamic, Dynamic> >(_data + n*_rows*_cols, _rows, _cols);
   }
   
   /// Returns an Eigen::Map-ed vector of the pixels at the given coordinate
   Map<Array<dataT, Dynamic, Dynamic>, Unaligned, Stride<Dynamic, Dynamic> > pixel(Index i, Index j)
   {
      return Map<Array<dataT, Dynamic, Dynamic>, Unaligned, Stride<Dynamic, Dynamic>  >(_data + j*_rows + i, _planes, 1, Stride<Dynamic, Dynamic>(0,_rows*_cols));
   }

   ///Return an Eigen::Map of the cube where each image is a vector
   Map<Array<dataT, Dynamic, Dynamic> > asVectors()
   {
      return Map<Array<dataT, Dynamic, Dynamic> >(_data,  _rows*_cols, _planes);
   }
 
   ///Calculate the covariance matrix of the images in the cube
   void Covar(Matrix<dataT, Dynamic, Dynamic> & cv)
   {
      cv = asVectors().matrix().transpose()*asVectors().matrix();
   }
   
   ///Calculate the mean image of the cube
   template<typename eigenT>
   void mean(eigenT & mim);
      
   ///Calculate the median image of the cube
   template<typename eigenT>
   void median(eigenT & mim);
   
   ///Calculate the sigma clipped mean image of the cube
   template<typename eigenT>
   void sigmaMean(eigenT & mim, Scalar sigma);
   
};


template<typename dataT>
template<typename eigenT>
void eigenCube<dataT>::mean(eigenT & mim)
{
   mim.resize(_rows, _cols);
   
   #pragma omp parallel for schedule(static, 1) num_threads(Eigen::nbThreads())
   for(Index i=0; i < _rows; ++i)
   {
      for(Index j=0;j< _cols; ++j)
      {
         mim(i,j) = pixel(i,j).mean(); 
      }
   }
}

template<typename dataT>
template<typename eigenT>
void eigenCube<dataT>::median(eigenT & mim)
{
   mim.resize(_rows, _cols);
   
   #pragma omp parallel for schedule(static, 10) num_threads(Eigen::nbThreads())
   for(Index i=0; i < _rows; ++i)
   {
      //Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> work;
      std::vector<Scalar> work;
      for(Index j=0;j< _cols; ++j)
      {
         mim(i,j) = eigenMedian(pixel(i,j), &work);             
      }
   }
}

template<typename dataT>
template<typename eigenT>
void eigenCube<dataT>::sigmaMean(eigenT & mim, dataT sigma)
{
   mim.resize(_rows, _cols);

   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<Scalar> work;
      Scalar var;
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.resize(_planes);//work could be smaller after sigmaMean
            int ii=0;
            for(int k =0; k< _planes; ++k) work[k] = (pixel(i,j))(k,0);
            
            mim(i,j) = vectorSigmaMean(work, sigma, var);             
         }
      }
      
      
   }
}











///Test whether a type is an eigenCube by testing whether it has a typedef of "is_eigenCube"
/** Used for compile-time determination of type
  * Example usage:
  * \code
  * bool is_eC = is_eigenCube<eigenCube<float> >; //Evaluates to true
  * bool is_not_eC = is_eigenCube<eigenImagef>; //Evaluates to false
  * \endcode
  */
//This was taken directly from the example at http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
template <typename T>
struct is_eigenCube 
{
   // Types "yes" and "no" are guaranteed to have different sizes,
   // specifically sizeof(yes) == 1 and sizeof(no) == 2.
   typedef char yes[1];
   typedef char no[2];
 
   template <typename imageT>
   static yes& test(typename imageT::is_eigenCube*);
 
   template <typename>
   static no& test(...);
 
   // If the "sizeof" of the result of calling test<T>(0) would be equal to sizeof(yes),
   // the first overload worked and T has a nested type named "is_mmatrix".
   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};

} //namespace mx

#endif //__eigenCube_hpp__
