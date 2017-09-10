/** \file eigenCube.hpp
  * \brief An image cube with an Eigen API
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup image_processing_files
  *
  */

#ifndef eigenCube_hpp
#define eigenCube_hpp

#pragma GCC system_header
#include <Eigen/Dense>


#include "../eigenUtils.hpp"
#include "../math/vectorUtils.hpp"
#include "eigenImage.hpp"



using namespace Eigen;

namespace mx
{
namespace improc 
{
   
/// An image cube with an Eigen-like API
/** \ingroup eigen_image_processing
  */
template<typename dataT>
class eigenCube 
{
public:
   typedef bool is_eigenCube;
   typedef dataT Scalar;

   typedef typename Array<dataT,Dynamic,Dynamic>::Index Index;
   
   typedef Map<Array<dataT, Dynamic, Dynamic> > imageRef;
   
   typedef Array<dataT,Dynamic,Dynamic> imageT;
   
protected:
   
   Index _rows;
   Index _cols;
   Index _planes;

   dataT * _data;

   bool _owner;

public:

   eigenCube();

   eigenCube( Index nrows, 
              Index ncols, 
              Index nplanes
            );

   eigenCube( dataT * ndata, 
              size_t nrows, 
              size_t ncols, 
              size_t nplanes
            );

   ~eigenCube();

   void setZero();
   
   eigenCube<dataT> & operator=(eigenCube<dataT> & ec);
   
   void shallowCopy(eigenCube<dataT> & src, bool takeOwner = false);
   
   void clear();
   
   void resize(int r, int c, int p);
   
   void resize(int r, int c);
   
   dataT * data();
   
   const dataT * data() const;
   
   const Index rows() const;
   
   const Index cols() const;
   
   const Index planes() const;

   /// Returns a 2D Eigen::Map pointed at the entire cube.
   Map<Array<dataT, Dynamic, Dynamic> > cube();
   
   /// Returns a 2D Eigen::Map pointed at the specified image.
   Map<Array<dataT, Dynamic, Dynamic> > image(Index n);
   
   /// Returns an Eigen::Map-ed vector of the pixels at the given coordinate
   Map<Array<dataT, Dynamic, Dynamic>, Unaligned, Stride<Dynamic, Dynamic> > pixel(Index i, Index j);

   ///Return an Eigen::Map of the cube where each image is a vector
   Map<Array<dataT, Dynamic, Dynamic> > asVectors();
 
   ///Calculate the covariance matrix of the images in the cube
   void Covar(Matrix<dataT, Dynamic, Dynamic> & cv);
   
   ///Calculate the mean image of the cube
   /**
     * \tparam eigenT an Eigen-like type.
     */ 
   template<typename eigenT>
   void mean( eigenT & mim /**< [out] the resultant mean image.  Is resized. */);
   
   ///Calculate the mean image of the cube with a mask.
   /**
     * \tparam eigenT an Eigen-like type.
     * \tparam eigenCubeT an eigenCube type.
     */ 
   template<typename eigenT, typename eigenCubeT>
   void mean( eigenT & mim, ///< [out] the resultant mean image.  Is resized. 
              eigenCubeT & mask ///< [in] a mask cube.  Only pixels with value 1 are included in the mean calculation.
            );
   
   ///Calculate the weighted mean image of the cube
   /**
     * \tparam eigenT an Eigen-like type.
     */
   template<typename eigenT>
   void mean( eigenT & mim, ///< [out] the resultant mean image. Is resized. 
              std::vector<dataT> & weights ///< [in] a vector of weights to use for calculating the mean 
            );
   
   ///Calculate the weighted mean image of the cube, with a mask cube.
   /**
     * \tparam eigenT an Eigen-like type.
     * \tparam eigenCubeT an eigenCube type.
     */
   template<typename eigenT, typename eigenCubeT>
   void mean( eigenT & mim, ///< [out] the resultant mean image. Is resized. 
              std::vector<dataT> & weights, ///< [in] a vector of weights to use for calculating the mean 
              eigenCubeT & mask ///< [in] a mask cube.  Only pixels with value 1 are included in the mean calculation.
            );
   
   ///Calculate the median image of the cube
   /**
     * \tparam eigenT an Eigen-like type.
     */
   template<typename eigenT>
   void median( eigenT & mim /**< [out] the resultant median image. Is resized. */);
   
   ///Calculate the sigma clipped mean image of the cube
   /**
     * \tparam eigenT an Eigen-like type.
     */
   template<typename eigenT>
   void sigmaMean( eigenT & mim, ///< [out] the resultant mean image  
                   Scalar sigma ///< [in] the sigma value at which to clip.
                 );
   
   ///Calculate the sigma clipped mean image of the cube, with a mask cube
   /**
     * \tparam eigenT an Eigen-like type.
     * \tparam eigenCubeT an eigenCube type.
     */
   template<typename eigenT, typename eigenCubeT>
   void sigmaMean( eigenT & mim,  ///< [out] the resultant mean image.  Is resized.
                   eigenCubeT & mask, ///< [in] a mask cube.  Only pixels with value 1 are included in the mean calculation.
                   Scalar sigma ///< [in] the sigma value at which to clip.
                 );
   
   ///Calculate the sigma clipped weighted mean image of the cube
   /**
     * \tparam eigenT an Eigen-like type.
     */
   template<typename eigenT>
   void sigmaMean( eigenT & mim, ///< [out] the resultant mean image. Is resized. 
                   std::vector<dataT> & weights, ///< [in] a vector of weights to use for calculating the mean 
                   Scalar sigma ///< [in] the sigma value at which to clip.
                 );
   
   ///Calculate the sigma clipped weighted mean image of the cube, with a mask cube.
   /**
     * \tparam eigenT an Eigen-like type.
     * \tparam eigenCubeT an eigenCube type.
     */
   template<typename eigenT, typename eigenCubeT>
   void sigmaMean( eigenT & mim, ///< [out] the resultant mean image. Is resized. 
                   std::vector<dataT> & weights, ///< [in] a vector of weights to use for calculating the mean 
                   eigenCubeT & mask, ///< [in] a mask cube.  Only pixels with value 1 are included in the mean calculation.
                   Scalar sigma ///< [in] the sigma value at which to clip.
                 );
};

template<typename dataT>
eigenCube<dataT>::eigenCube()
{
   _rows = 0;
   _cols = 0;
   _planes = 0;
   _data = 0;
   _owner = false;
}

template<typename dataT>
eigenCube<dataT>::eigenCube(Index nrows, Index ncols, Index nplanes)
{
   _rows = nrows;
   _cols = ncols;
   _planes = nplanes;

   _data = new dataT[_rows*_cols*_planes];
   _owner = true;
}

template<typename dataT>
eigenCube<dataT>::eigenCube(dataT * ndata, size_t nrows, size_t ncols, size_t nplanes)
{
   _rows = nrows;
   _cols = ncols;
   _planes = nplanes;

   _data = ndata;
   _owner = false;
}

template<typename dataT>
eigenCube<dataT>::~eigenCube()
{
   if(_owner && _data)
   {
      delete _data;
   }
}

template<typename dataT>
void eigenCube<dataT>::setZero()
{
   int N = _rows*_cols*_planes;
   
   for(int i=0;i<N;++i) _data[i] = ((dataT) 0);
}

template<typename dataT>
eigenCube<dataT> & eigenCube<dataT>::operator=(eigenCube<dataT> & ec)
{
   resize(ec.rows(), ec.cols(), ec.planes());
   
   int N = _rows*_cols*_planes;
   
   for(int i=0;i<N;++i) _data[i] = ec._data[i];
   
   return *this;
}

template<typename dataT>
void eigenCube<dataT>::shallowCopy(eigenCube<dataT> & src, bool takeOwner)
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

template<typename dataT>
void eigenCube<dataT>::clear()
{
   if(_owner && _data)
   {
      delete _data;
   }
   
   _rows = 0;
   _cols = 0;
   _planes = 0;
}

template<typename dataT>
void eigenCube<dataT>::resize(int r, int c, int p)
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

template<typename dataT>
void eigenCube<dataT>::resize(int r, int c)
{
   resize(r, c, 1);
}

template<typename dataT>
dataT * eigenCube<dataT>::data()
{
   return _data;
}

template<typename dataT>
const dataT * eigenCube<dataT>::data() const
{
   return _data;
}

template<typename dataT>
const typename eigenCube<dataT>::Index eigenCube<dataT>::rows() const
{
   return _rows;
}

template<typename dataT>
const typename eigenCube<dataT>::Index eigenCube<dataT>::cols() const
{
   return _cols;
}

template<typename dataT>
const typename eigenCube<dataT>::Index eigenCube<dataT>::planes() const
{
   return _planes;
}

template<typename dataT>
Map<Array<dataT, Dynamic, Dynamic> > eigenCube<dataT>::cube()
{
   return Map<Array<dataT, Dynamic, Dynamic> >(_data, _rows*_cols, _planes);
}

template<typename dataT>
Map<Array<dataT, Dynamic, Dynamic> > eigenCube<dataT>::image(Index n)
{
   return Map<Array<dataT, Dynamic, Dynamic> >(_data + n*_rows*_cols, _rows, _cols);
}

template<typename dataT>
Map<Array<dataT, Dynamic, Dynamic>, Unaligned, Stride<Dynamic, Dynamic> > eigenCube<dataT>::pixel(Index i, Index j)
{
   return Map<Array<dataT, Dynamic, Dynamic>, Unaligned, Stride<Dynamic, Dynamic>  >(_data + j*_rows + i, _planes, 1, Stride<Dynamic, Dynamic>(0,_rows*_cols));
}

template<typename dataT>
Map<Array<dataT, Dynamic, Dynamic> > eigenCube<dataT>::asVectors()
{
   return Map<Array<dataT, Dynamic, Dynamic> >(_data,  _rows*_cols, _planes);
}

template<typename dataT>
void eigenCube<dataT>::Covar(Matrix<dataT, Dynamic, Dynamic> & cv)
{
   cv = asVectors().matrix().transpose()*asVectors().matrix();
}
   
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
template<typename eigenT, typename eigenCubeT>
void eigenCube<dataT>::mean( eigenT & mim, 
                             eigenCubeT & mask
                           )
{
   mim.resize(_rows, _cols);
   
   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<dataT> work;
      
      #pragma omp for schedule(static, 1) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            
            work.clear();
            
            for(Index k=0; k<_planes; ++k)
            {
            
               if( (mask.pixel(i,j))(k,0) == 1)
               {
                  work.push_back( (pixel(i,j))(k,0) );
               }
            }
            if(work.size() > 0.75*_planes)
            {
               mim(i,j) = math::vectorMean(work); 
            }
            else mim(i,j) = std::numeric_limits<dataT>::quiet_NaN();   
            
         }
      }
   }
}

template<typename dataT>
template<typename eigenT>
void eigenCube<dataT>::mean( eigenT & mim, 
                             std::vector<dataT> & weights
                           )
{
   mim.resize(_rows, _cols);

   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<Scalar> work;
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.resize(_planes);//work could be smaller after sigmaMean
            int ii=0;
            for(int k =0; k< _planes; ++k) work[k] = (pixel(i,j))(k,0);
            
            mim(i,j) = math::vectorMean(work, weights);             
         }
      }
   }
}

template<typename dataT>
template<typename eigenT, typename eigenCubeT>
void eigenCube<dataT>::mean( eigenT & mim, 
                             std::vector<dataT> & weights, 
                             eigenCubeT & mask 
                           )
{
   mim.resize(_rows, _cols);

   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<Scalar> work, wwork;
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.clear();
            wwork.clear();
            
            int ii=0;
            
            for(Index k =0; k< _planes; ++k) 
            {
               if( (mask.pixel(i,j))(k,0) == 1)
               {
                  work.push_back( (pixel(i,j))(k,0) );
                  wwork.push_back( weights[k] );
               }
            }
            if(work.size() > 0.75*_planes)
            {
               mim(i,j) = math::vectorMean(work, wwork); 
            }
            else mim(i,j) = std::numeric_limits<dataT>::quiet_NaN();
         }
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
         mim(i,j) = imageMedian(pixel(i,j), &work);             
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
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.resize(_planes);//work could be smaller after sigmaMean
            int ii=0;
            for(int k =0; k< _planes; ++k) work[k] = (pixel(i,j))(k,0);
            
            mim(i,j) = math::vectorSigmaMean(work, sigma);             
         }
      }
      
      
   }
}

template<typename dataT>
template<typename eigenT, typename eigenCubeT>
void eigenCube<dataT>::sigmaMean(eigenT & mim, eigenCubeT & mask, dataT sigma)
{
   mim.resize(_rows, _cols);

   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<Scalar> work;
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.clear();
            int ii=0;
            
            for(Index k =0; k< _planes; ++k) 
            {
               if( (mask.pixel(i,j))(k,0) == 1)
               {
                  work.push_back( (pixel(i,j))(k,0) );
               }
            }
            if(work.size() > 0.75*_planes)
            {
               mim(i,j) = math::vectorSigmaMean(work, sigma); 
            }
            else mim(i,j) = std::numeric_limits<dataT>::quiet_NaN();
            
         }
      }
      
      
   }
}

template<typename dataT>
template<typename eigenT>
void eigenCube<dataT>::sigmaMean(eigenT & mim, std::vector<dataT> & weights, dataT sigma)
{
   mim.resize(_rows, _cols);

   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<Scalar> work;
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.resize(_planes);//work could be smaller after sigmaMean
            int ii=0;
            for(int k =0; k< _planes; ++k) work[k] = (pixel(i,j))(k,0);
            
            mim(i,j) = math::vectorSigmaMean(work, weights, sigma);             
         }
      }
      
      
   }
}

template<typename dataT>
template<typename eigenT, typename eigenCubeT>
void eigenCube<dataT>::sigmaMean(eigenT & mim, std::vector<dataT> & weights, eigenCubeT & mask, dataT sigma)
{
   mim.resize(_rows, _cols);

   #pragma omp parallel num_threads(Eigen::nbThreads())
   {
      std::vector<Scalar> work, wwork;
      
      #pragma omp for schedule(static, 10) 
      for(Index i=0; i < _rows; ++i)
      {
         for(Index j=0;j< _cols; ++j)
         {
            work.clear();
            wwork.clear();
            
            int ii=0;
            
            for(Index k =0; k< _planes; ++k) 
            {
               if( (mask.pixel(i,j))(k,0) == 1)
               {
                  work.push_back( (pixel(i,j))(k,0) );
                  wwork.push_back( weights[k] );
               }
            }
            if(work.size() > 0.75*_planes)
            {
               mim(i,j) = math::vectorSigmaMean(work, wwork, sigma); 
            }
            else mim(i,j) = std::numeric_limits<dataT>::quiet_NaN();
         }
      }
   }
}







} //namespace improc 

} //namespace mx

#endif //eigenCube_hpp
