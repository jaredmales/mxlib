#ifndef __eigenUtils_hpp__
#define __eigenUtils_hpp__

namespace mx
{

#include <Eigen/Dense>
#include <cmath>
#include <sofa.h>
#include "geo.h"
extern "C"
{
#include <cblas.h>
}
  
#include <templateLapack.hpp>

/// Fills in the cells of an Eigen 2D Array with their radius from the center
/** \ingroup image_processing
  *
  * \param m is the allocated Eigen Array
  * \param xc is the x center
  * \param yc is the y center
  * \param scale [optional] is a scaling to apply to each value (default = 1)
  */  
template<class eigenT, typename arithT = typename eigenT::Scalar> 
void radiusImage( eigenT & m, 
                arithT xc, 
                arithT yc, 
                arithT scale=1
              )
{
   arithT f_x, f_y;

   size_t dim1 = m.rows();
   size_t dim2 = m.cols();
   
   for(size_t i=0; i < dim1; i++)
   {
      f_x = (i-xc)*(i-xc);
      
      for(size_t j=0; j < dim2; j++)
      {
         f_y = (j-yc)*(j-yc);
         
         m(i,j) = std::sqrt( f_x + f_y)*scale;
      }
   }
   
}

/// Fills in the cells of Eigen 2D Array with their radius from the canonical center
/** \ingroup image_processing
  *
  * The center is @f$ (x_c, y_c) = (0.5*(dim_1-1), 0.5*(dim_2 -1)) @f$.
  *
  * \param m is the allocated Eigen Array
  * \param scale [optional] is a scaling to apply to each value (default = 1)
  */  
template<class eigenT, typename arithT = typename eigenT::Scalar> 
void radiusImage(eigenT & m, arithT scale=1)
{
   arithT xc, yc;
   
   xc = 0.5*(m.rows()-1);
   yc = 0.5*(m.cols()-1);
   
   radiusImage(m, xc, yc, scale);
}




/// Fills in the cells of an Eigen 2D Array with their angle relative to the center
/** \ingroup image_processing
  *
  * \param m is the allocated Eigen Array
  * \param xc is the x center
  * \param yc is the y center
  * \param scale [optional] is a scaling to apply to each value (default is \ref DR2D)
  */  
template<class eigenT, typename arithT = typename eigenT::Scalar> 
void angleImage( eigenT & m, 
                arithT xc, 
                arithT yc, 
                arithT scale= DR2D
               )
{
   arithT f_x, f_y;

   size_t dim1 = m.rows();
   size_t dim2 = m.cols();
   
   for(size_t i=0; i < dim1; i++)
   {
      f_x = (i-xc);
      
      for(size_t j=0; j < dim2; j++)
      {
         f_y = (j-yc);
         
         m(i,j) = fmod(atan2(f_y, f_x) + D2PI, D2PI)  *scale;
      }
   }
   
}

/// Fills in the cells of Eigen 2D Array with their angle relative the canonical center
/** \ingroup image_processing
  *
  * The center is @f$ (x_c, y_c) = (0.5*(dim_1-1), 0.5*(dim_2 -1)) @f$.
  *
  * \param m is the allocated Eigen Array
  * \param scale [optional] is a scaling to apply to each value (default = \ref DR2D)
  */  
template<class eigenT, typename arithT = typename eigenT::Scalar> 
void angleImage(eigenT & m, arithT scale= DR2D)
{
   arithT xc, yc;
   
   xc = 0.5*(m.rows()-1);
   yc = 0.5*(m.cols()-1);
   
   angleImage(m, xc, yc, scale);
  
}


/// Fills in the cells of an Eigen 2D Array with their radius amd angle relative to the center
/** \ingroup image_processing
  *
  * \param m is the allocated Eigen Array
  * \param xc is the x center
  * \param yc is the y center
  * \param scale [optional] is a scaling to apply to each value (default is \ref DR2D)
  */  
template<class eigenT, typename arithT = typename eigenT::Scalar> 
void radAngImage( eigenT & rIm,
                  eigenT & qIm,
                  arithT xc, 
                  arithT yc,
                  arithT rscale = 1,
                  arithT qscale= DR2D
                 )
{
   arithT f_x, f_y;

   size_t dim1 = rIm.rows();
   size_t dim2 = rIm.cols();
   
   for(size_t i=0; i < dim1; ++i)
   {
      f_x = ( ((arithT)i)-xc);
      
      for(size_t j=0; j < dim2; ++j)
      {
         f_y = (((arithT)j)-yc);
         rIm(i,j) = std::sqrt( f_x*f_x + f_y*f_y)*rscale;
         qIm(i,j) = fmod(atan2(f_y, f_x) +D2PI, D2PI) *qscale;
      }
   }
}

template<typename eigenT>
void imageRegionIndices( vector<size_t> & idx, 
                         eigenT &rIm, 
                         eigenT &qIm,
                         typename eigenT::Scalar xcen,
                         typename eigenT::Scalar ycen,
                         typename eigenT::Scalar min_r, 
                         typename eigenT::Scalar max_r,
                         typename eigenT::Scalar min_q, 
                         typename eigenT::Scalar max_q)
{

   int min_x = -max_r, max_x = max_r, min_y = -max_r, max_y = max_r;

   if(max_q == 0) max_q = 360.;
   
   size_t msize = ((DPI*(max_r*max_r - min_r*min_r)) * (max_q-min_q)/360.) *1.01 + 1;
   
   //This was tested, this is slightly faster than resize with an erase.
   idx.reserve(msize);
   
   size_t x0 = xcen+min_x;
   size_t x1 = xcen+max_x;
   size_t y0 = ycen+min_y;
   size_t y1 = ycen+max_y;
   
   for(size_t i = x0; i< x1; ++i)
   {
      for(size_t j = y0; j< y1; ++j)
      { 
         if(rIm(i,j) >= min_r && rIm(i,j) <= max_r && qIm(i,j) >= min_q && qIm(i,j) <= max_q) 
         {
            idx.push_back(j*rIm.rows() + i);
         }
      }
   }
}

template<typename dataT>
static int eigenMedian_compare (const void * a, const void * b)
{
  if( *(dataT*)a < *(dataT*)b) return -1;
  if( *(dataT*)a > *(dataT*)b) return 1;
  
  return 0;
   
}

#if 1
template<typename eigenT>
typename eigenT::Scalar eigenMedian(const eigenT & mat, Eigen::Array<typename eigenT::Scalar, Eigen::Dynamic, Eigen::Dynamic> * work =0)
{
   typename eigenT::Scalar med;
   
   bool localWork = false;
   if(work == 0) 
   {
      work = new Eigen::Array<typename eigenT::Scalar, Eigen::Dynamic, Eigen::Dynamic>;
      localWork = true;
   }
   
   *work = mat;
   
   qsort(work->data(), mat.size(), sizeof(typename eigenT::Scalar), &eigenMedian_compare<typename eigenT::Scalar>);
   
   if(mat.size()%2 == 0)
   {
      med = 0.5*( (*work)((int) floor(0.5*mat.size()) - 1) + (*work)((int) floor(0.5*mat.size())));
   }
   else
   {
      med = (*work)((int) floor(0.5*mat.size()));
   }
      
   if(localWork) delete work;
   
   return med;
}
 
#else

template<typename eigenT>
typename eigenT::Scalar eigenMedian(const eigenT & mat, std::vector<typename eigenT::Scalar> * work =0)
{
   typename eigenT::Scalar med;
   
   bool localWork = false;
   if(work == 0) 
   {
      work = new std::vector<typename eigenT::Scalar>;
      localWork = true;
   }
   
   work->resize(mat.size());
   
   int ii = 0;
   for(int i=0;i<mat.rows();++i)
   {
      for(int j=0; j<mat.cols();++j)
      {
         (*work)[ii] = mat(i,j);
         ++i;
      }
   }
   
   qsort(work->data(), mat.size(), sizeof(typename eigenT::Scalar), &eigenMedian_compare<typename eigenT::Scalar>);
   
   if(mat.size()%2 == 0)
   {
      med = 0.5*( (*work)[(int) floor(0.5*mat.size()) - 1] + (*work)[(int) floor(0.5*mat.size())]);
   }
   else
   {
      med = (*work)[(int) floor(0.5*mat.size())];
   }
      
   if(localWork) delete work;
   
   return med;
} 
#endif   


/// Calculates the lower triangular part of the covariance matrix of ims.
/** Uses cblas_ssyrk.  cv is resized to ims.rows() X ims.rows().
  * Calculates \f$ cv = AA^T \f$.
  * 
  * \param ims is the eigen matrix/array to calculate the covariance of
  * \param cv is the eigen matrix/array where to store the result
  *
  * \tparam eigenT1 is the eigen matrix/array type of cv.
  * \tparam eigenT2 is the eigen matrix/array type of ims
  */ 
template<typename eigenT1, typename eigenT2>
void eigen_covar_ssyrk(eigenT1 &cv, eigenT2 &ims)
{
   cv.resize(ims.rows(), ims.rows());
   
   cblas_ssyrk(/*const enum CBLAS_ORDER Order*/ CblasColMajor, /*const enum CBLAS_UPLO Uplo*/ CblasLower,
                 /*const enum CBLAS_TRANSPOSE Trans*/ CblasNoTrans, /*const int N*/ims.rows(), /*const int K*/ ims.cols(),
                 /*const float alpha*/ 1.0, /*const float *A*/ims.data(), /*const int lda*/ ims.rows(),
                 /*const float beta*/ 0., /*float *C*/ cv.data(), /*const int ldc*/ cv.rows());
   
}   

template<typename eigenT1, typename eigenT2>
void eigen_covar_dsyrk(eigenT1 &cv, eigenT2 &ims)
{
   cv.resize(ims.rows(), ims.rows());
   
   cblas_dsyrk(/*const enum CBLAS_ORDER Order*/ CblasColMajor, /*const enum CBLAS_UPLO Uplo*/ CblasLower,
                 /*const enum CBLAS_TRANSPOSE Trans*/ CblasNoTrans, /*const int N*/ims.rows(), /*const int K*/ ims.cols(),
                 /*const float alpha*/ 1.0, /*const float *A*/ims.data(), /*const int lda*/ ims.rows(),
                 /*const float beta*/ 0., /*float *C*/ cv.data(), /*const int ldc*/ cv.rows());
   
}   

template<typename eigenT>
int eigenSYEVR(eigenT &X, eigenT &eigvec, eigenT &eigval) 
{
        /*
                This function calculates the eigenvalues and eigenvectors of 
                the n*n symmetric matrix X. 
                The matrices have to be in Fortran vector format.
                The eigenvectors will be put columnwise in the n*n matrix eigvec,
                where the corresponding eigenvalues will be put in the vector 
                eigval (length n of course). Only the lower triangle of the matrix
                X is used. The content of X is not changed.
                
                This function first queries the Lapack routines for optimal workspace 
                sizes. These memoryblocks are then allocated and the decomposition is 
                calculated using the Lapack function "dsyevr". The allocated memory 
                is then freed. 
        */

   typedef typename eigenT::Scalar dataT;
   
   dataT *WORK;
   int *ISUPPZ, *IWORK;
   int  numeig, info, sizeWORK, sizeIWORK;
           
   int n = X.rows();
   eigvec.resize(n,n);
   eigval.resize(n,1);
   
   /*  Use a copy of X so we don't need to change its value or use its memoryblock */
   eigenT Xc = X;//     Xc=malloc(n*n*sizeof(float));
                
   /*  The support of the eigenvectors. We will not use this but the routine needs it  */
   ISUPPZ = (int *) malloc (2*n*sizeof(dataT));
        
   /*  Allocate temporarily minimally allowed size for workspace arrays */
   WORK = (dataT *) malloc (26*n*sizeof(dataT));
   IWORK = (int *) malloc (10*n*sizeof(int));
                
   /*  Check for NULL-pointers.  */
   if ((ISUPPZ==NULL)||(WORK==NULL)||(IWORK==NULL)) 
   {
      printf("malloc failed in eigen_decomposition\n"); 
      return 2;
   }
        
   /*  Query the Lapack routine for optimal sizes for workspace arrays  */
   info=syevr<dataT>('V', 'A', 'L', n, Xc.data(), n, 0, 0, 0, 0, lamch<dataT>('S'), &numeig, eigval.data(), eigvec.data(), n, ISUPPZ, WORK, -1, IWORK, -1);
   sizeWORK = (int)WORK[0]; 
   sizeIWORK = IWORK[0]; 
        
   /*  Free previous allocation and reallocate preferable workspaces, Check result  */
   free(WORK);
   free(IWORK);
   WORK = (dataT *) malloc (sizeWORK*sizeof(dataT));
   IWORK = (int *) malloc (sizeIWORK*sizeof(int));
   if ((WORK==NULL)||(IWORK==NULL)) 
   {
      printf("malloc failed in eigen_decomposition\n"); 
      return 2;
   }
   printf("starting\n");
   fflush(stdout);
        
   /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
   info=syevr<dataT>('V', 'A', 'L', n, Xc.data(), n, 0, 0, 0, 0, lamch<dataT>('S'), &numeig, eigval.data(), eigvec.data(), n, ISUPPZ, WORK, sizeWORK, IWORK, sizeIWORK);
        
    /*  Cleanup and exit  */
   free(WORK); free(IWORK); free(ISUPPZ);
   return info;
}       

}//namespace mx

#endif //__eigenUtils_hpp__
