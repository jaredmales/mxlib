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


#include "../math/vectorUtils.hpp"
#include "eigenImage.hpp"
#include "imageUtils.hpp"

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

   typedef typename Eigen::Array<dataT,Eigen::Dynamic,Eigen::Dynamic>::Index Index;

   typedef Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > imageRef;

   typedef Eigen::Array<dataT,Eigen::Dynamic,Eigen::Dynamic> imageT;

protected:

   Index _rows;
   Index _cols;
   Index _planes;

   dataT * m_data {nullptr};

   bool _owner;

public:

   eigenCube();

   /// C'tor which will allocate space..
   /** 
     */
   eigenCube( Index nrows,  ///< [in] Number of rows in the cube
              Index ncols,  ///< [in] Number of columns the cube
              Index nplanes ///< [in] Number of planes in the cube
            );

   /// C'tor taking an existing array as an argument.
   /** The existing array is used as is, as in a map, and ownership is not taken.  You are
     * responsible for memory management, e.g. free-ing this array.
     */
   eigenCube( dataT * ndata, ///< [in] Allocated array with a cube of nrows x ncols x nplanes
              size_t nrows,  ///< [in] Number of rows in the cube
              size_t ncols,  ///< [in] Number of columns the cube
              size_t nplanes ///< [in] Number of planes in the cube
            );

   ~eigenCube();

   void setZero();

   eigenCube<dataT> & operator=(const eigenCube<dataT> & ec);

   void shallowCopy(eigenCube<dataT> & src, bool takeOwner = false);

   void clear();

   void resize(int r, int c, int p);

   void resize(int r, int c);

   dataT * data();

   const dataT * data() const;

   const Index rows() const;

   const Index cols() const;

   const Index planes() const;

   /// Returns a 2D Eigen::Eigen::Map pointed at the entire cube.
   Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > cube();

   /// Returns a 2D Eigen::Eigen::Map pointed at the specified image.
   Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > image(Index n /**< [in] the image number */);

   /// Returns a 2D Eigen::Eigen::Map pointed at the specified image.
   const Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > image(Index n /**< [in] the image number */) const;
   
   /// Returns an Eigen::Eigen::Map-ed vector of the pixels at the given coordinate
   Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > pixel(Index i, Index j);

   ///Return an Eigen::Eigen::Map of the cube where each image is a vector
   Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > asVectors();

   ///Calculate the covariance matrix of the images in the cube
   void Covar(Eigen::Matrix<dataT, Eigen::Dynamic, Eigen::Dynamic> & cv);

   ///Calculate the sum image of the cube
   /**
     * \tparam eigenT an Eigen-like type.
     */
   template<typename eigenT>
   void sum( eigenT & mim /**< [out] the resultant sum image.  Is resized. */);
   
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
              eigenCubeT & mask, ///< [in] a mask cube.  Only pixels with value 1 are included in the mean calculation.
              double minGoodFract = 0.0 ///< [in] [optional] the minimum fraction of good pixels, if not met then the pixel is NaN-ed.
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
              eigenCubeT & mask, ///< [in] a mask cube.  Only pixels with value 1 are included in the mean calculation.
              double minGoodFract = 0.0 ///< [in] [optional] the minimum fraction of good pixels, if not met then the pixel is NaN-ed.
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
                   Scalar sigma, ///< [in] the sigma value at which to clip.
                   double minGoodFract = 0.0 ///< [in] [optional] the minimum fraction of good pixels, if not met then the pixel is NaN-ed.
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
                   Scalar sigma, ///< [in] the sigma value at which to clip.
                   double minGoodFract = 0.0 ///< [in] [optional] the minimum fraction of good pixels, if not met then the pixel is NaN-ed.
                 );
   
   
};


template<typename dataT>
eigenCube<dataT>::eigenCube()
{
   _rows = 0;
   _cols = 0;
   _planes = 0;
   _owner = false;
}

template<typename dataT>
eigenCube<dataT>::eigenCube(Index nrows, Index ncols, Index nplanes)
{
   _rows = nrows;
   _cols = ncols;
   _planes = nplanes;

   m_data = new dataT[_rows*_cols*_planes];
   _owner = true;
}

template<typename dataT>
eigenCube<dataT>::eigenCube(dataT * ndata, size_t nrows, size_t ncols, size_t nplanes)
{
   _rows = nrows;
   _cols = ncols;
   _planes = nplanes;

   m_data = ndata;
   _owner = false;
}

template<typename dataT>
eigenCube<dataT>::~eigenCube()
{
   if(_owner && m_data)
   {
      delete[] m_data;
   }
}

template<typename dataT>
void eigenCube<dataT>::setZero()
{
   int N = _rows*_cols*_planes;

   for(int i=0;i<N;++i) m_data[i] = ((dataT) 0);
}

template<typename dataT>
eigenCube<dataT> & eigenCube<dataT>::operator=(const eigenCube<dataT> & ec)
{
   resize(ec.rows(), ec.cols(), ec.planes());

   int N = _rows*_cols*_planes;

   for(int i=0;i<N;++i) m_data[i] = ec.m_data[i];

   return *this;
}

template<typename dataT>
void eigenCube<dataT>::shallowCopy(eigenCube<dataT> & src, bool takeOwner)
{
   if(_owner && m_data)
   {
      delete m_data;
      _owner = false;
   }

   _rows = src._rows;
   _cols = src._cols;
   _planes = src._planes;
   m_data = src.m_data;

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
   if(_owner && m_data)
   {
      delete[] m_data;
   }

   _rows = 0;
   _cols = 0;
   _planes = 0;
}

template<typename dataT>
void eigenCube<dataT>::resize(int r, int c, int p)
{
   if(_owner && m_data)
   {
      delete[] m_data;
   }

   _rows = r;
   _cols = c;
   _planes = p;

   m_data = new dataT[_rows*_cols*_planes];
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
   return m_data;
}

template<typename dataT>
const dataT * eigenCube<dataT>::data() const
{
   return m_data;
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
Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > eigenCube<dataT>::cube()
{
   return Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> >(m_data, _rows*_cols, _planes);
}

template<typename dataT>
Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > eigenCube<dataT>::image(Index n)
{
   return Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> >(m_data + n*_rows*_cols, _rows, _cols);
}

template<typename dataT>
const Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > eigenCube<dataT>::image(Index n) const
{
   return Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> >(m_data + n*_rows*_cols, _rows, _cols);
}

template<typename dataT>
Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > eigenCube<dataT>::pixel(Index i, Index j)
{
   return Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>  >(m_data + j*_rows + i, _planes, 1, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(0,_rows*_cols));
}

template<typename dataT>
Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> > eigenCube<dataT>::asVectors()
{
   return Eigen::Map<Eigen::Array<dataT, Eigen::Dynamic, Eigen::Dynamic> >(m_data,  _rows*_cols, _planes);
}

template<typename dataT>
void eigenCube<dataT>::Covar(Eigen::Matrix<dataT, Eigen::Dynamic, Eigen::Dynamic> & cv)
{
   cv = asVectors().matrix().transpose()*asVectors().matrix();
}

template<typename dataT>
template<typename eigenT>
void eigenCube<dataT>::sum(eigenT & mim)
{
   mim.resize(_rows, _cols);

   #pragma omp parallel for schedule(static, 1) num_threads(Eigen::nbThreads())
   for(Index i=0; i < _rows; ++i)
   {
      for(Index j=0;j< _cols; ++j)
      {
         mim(i,j) = pixel(i,j).sum();
      }
   }
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
                             eigenCubeT & mask,
                             double minGoodFract
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
            if(work.size() > minGoodFract*_planes)
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
            //int ii=0;
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
                             eigenCubeT & mask,
                             double minGoodFract
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

            //int ii=0;

            for(Index k =0; k< _planes; ++k)
            {
               if( (mask.pixel(i,j))(k,0) == 1)
               {
                  work.push_back( (pixel(i,j))(k,0) );
                  wwork.push_back( weights[k] );
               }
            }
            if(work.size() > minGoodFract*_planes)
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
      //Eigen::Eigen::Array<Scalar, Eigen::Eigen::Dynamic, Eigen::Eigen::Dynamic> work;
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
            for(int k =0; k< _planes; ++k) work[k] = (pixel(i,j))(k,0);

            mim(i,j) = math::vectorSigmaMean(work, sigma);
         }
      }


   }
}

template<typename dataT>
template<typename eigenT, typename eigenCubeT>
void eigenCube<dataT>::sigmaMean( eigenT & mim, 
                                  eigenCubeT & mask, 
                                  dataT sigma,
                                  double minGoodFract
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
            work.clear();
            int ii=0;

            for(Index k =0; k< _planes; ++k)
            {
               if( (mask.pixel(i,j))(k,0) == 1)
               {
                  work.push_back( (pixel(i,j))(k,0) );
               }
            }
            if(work.size() > minGoodFract*_planes)
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
void eigenCube<dataT>::sigmaMean( eigenT & mim, 
                                  std::vector<dataT> & weights, 
                                  eigenCubeT & mask, 
                                  dataT sigma,
                                  double minGoodFract
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
            if(work.size() > minGoodFract*_planes)
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
