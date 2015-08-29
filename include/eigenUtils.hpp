#ifndef __eigenUtils_hpp__
#define __eigenUtils_hpp__

#include <Eigen/Dense>
#include <cmath>
#include <sofa.h>
#include "geo.h"

#include "templateBLAS.hpp"
#include "templateLapack.hpp"


namespace mx
{


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

///Get the vector indices of an annular region in an image
/** \ingroup image_processing
  * 
  * \param rIm is a radius image of the type produced by \ref radiusImage
  * \param qIm is an angle image of the type produce by \ref angleImage
  * \param xcen is the x center of the image
  * \param ycen is the y center of the image
  * \param min_r is the minimum radius of the region
  * \param max_r is the maximum radius of the region
  * \param min_q is the minimum angle of the region
  * \param max_q is the maximum angle of the region
  * 
  * \returns a vector containing the 1D indices of the region defined by the input parameters
  */
template<typename eigenT>
std::vector<size_t> imageRegionIndices( eigenT &rIm, 
                         eigenT &qIm,
                         typename eigenT::Scalar xcen,
                         typename eigenT::Scalar ycen,
                         typename eigenT::Scalar min_r, 
                         typename eigenT::Scalar max_r,
                         typename eigenT::Scalar min_q, 
                         typename eigenT::Scalar max_q)
{

   std::vector<size_t> idx;
   
   int min_x = -max_r, max_x = max_r, min_y = -max_r, max_y = max_r;

   if(max_q == 0) max_q = 360.;
   
   size_t msize = ((DPI*(max_r*max_r - min_r*min_r)) * (max_q-min_q)/360.) *1.01 + 1;
   
   //This was tested, this is slightly faster than resize with an erase.
   idx.reserve(msize);
   
   int x0 = xcen+min_x;
   if(x0 < 0) x0 = 0;
   int x1 = xcen+max_x+1;
   if(x1 > rIm.rows()) x1 = rIm.rows();
   int y0 = ycen+min_y;
   if(y0 < 0) y0 = 0;
   int y1 = ycen+max_y+1;
   if(y1 > rIm.cols()) y1 = rIm.cols();
   
   for(size_t i = x0; i< x1; ++i)
   {
      for(size_t j = y0; j< y1; ++j)
      { 
         if(rIm(i,j) >= min_r && rIm(i,j) < max_r && qIm(i,j) >= min_q && qIm(i,j) < max_q) 
         {
            idx.push_back(i*rIm.cols() + j);
         }
      }
   }
   
   
   return idx;
}





template<typename imageTout, typename imageTin, typename coeffT>
void cutImageRegion(imageTout & imout, const imageTin & imin,  coeffT & coeffs, bool resize = true)
{
   if(resize)
   {
      imout.resize(coeffs.size(),1);
   }
   
   #pragma omp parallel for schedule(static, 1)
   for(int i=0;i<coeffs.size();++i)
   {
      imout(i) = imin(coeffs[i]);
   }
   
}
 
template<typename imageTout, typename imageTin, typename coeffT>
void insertImageRegion(imageTout imout, const imageTin & imin,  coeffT & coeffs)
{
   #pragma omp parallel for schedule(static, 1)
   for(int i=0;i<coeffs.size();++i)
   {
      imout(coeffs[i]) = imin(i);
   }
   
} 


template<typename eigenT, typename eigenTin>
void removeRowsAndCols(eigenT & out, const eigenTin & in, int st, int w)
{
   
   out.resize(in.rows() - w, in.cols() - w);
   
   out.topLeftCorner(st,st) = in.topLeftCorner(st,st);
   
   out.bottomLeftCorner(in.rows()-(st+w), st) = in.bottomLeftCorner(in.rows()-(st+w), st);
   
   out.topRightCorner(st, in.cols()-(st+w))  = in.topRightCorner(st, in.cols()-(st+w));
   
   out.bottomRightCorner(in.rows()-(st+w),in.cols()-(st+w)) = in.bottomRightCorner(in.rows()-(st+w),in.cols()-(st+w));
}







template<typename eigenT, typename eigenTin>
void removeRows(eigenT & out,  const eigenTin & in, int st, int w)
{
   
   out.resize(in.rows() - w, in.cols());
   
   out.topLeftCorner(st,in.cols()) = in.topLeftCorner(st,in.cols());
   
   out.bottomLeftCorner(in.rows()-(st+w), in.cols()) = in.bottomLeftCorner(in.rows()-(st+w), in.cols());
   
}

template<typename eigenT, typename eigenTin>
void removeCols(eigenT & out,  const eigenTin & in, int st, int w)
{
   
   out.resize(in.rows(), in.cols() - w);
   
   out.topLeftCorner(in.rows(), st) = in.topLeftCorner(in.rows(), st);
   
   out.topRightCorner(in.rows(),in.cols()-(st+w)) = in.topRightCorner(in.rows(),in.cols()-(st+w));
   
}   
   
template<typename dataT>
static int eigenMedian_compare (const void * a, const void * b)
{
  if( *(dataT*)a < *(dataT*)b) return -1;
  if( *(dataT*)a > *(dataT*)b) return 1;
  
  return 0;
   
}



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
         ++ii;
      }
   }

   int n = 0.5*mat.size();
   
   nth_element(work->begin(), work->begin()+n, work->end());
   
   med = (*work)[n];
   
   if(mat.size()%2 == 0)
   {
      //nth_element(work->begin(), work->begin()+n-1, work->end());
      med = 0.5*(med + *std::max_element(work->begin(), work->begin()+n)); //(*work)[n-1]);
   }
         
   if(localWork) delete work;
   
   return med;
} 


/// Calculates the lower triangular part of the covariance matrix of ims.
/** Uses cblas_ssyrk.  cv is resized to ims.cols() X ims.cols().
  * Calculates \f$ cv = A^T*A\f$.
  * 
  * \param ims is the eigen matrix/array (images as columns) to calculate the covariance of
  * \param cv is the eigen matrix/array where to store the result
  *
  * \tparam eigenT1 is the eigen matrix/array type of cv.
  * \tparam eigenT2 is the eigen matrix/array type of ims
  */ 
template<typename eigenT1, typename eigenT2>
void eigenSYRK(eigenT1 &cv, const eigenT2 &ims)
{
   cv.resize(ims.cols(), ims.cols());
   
   syrk<typename eigenT1::Scalar>(/*const enum CBLAS_ORDER Order*/ CblasColMajor, /*const enum CBLAS_UPLO Uplo*/ CblasLower,
                 /*const enum CBLAS_TRANSPOSE Trans*/ CblasTrans, /*const int N*/ims.cols(), /*const int K*/ ims.rows(),
                 /*const float alpha*/ 1.0, /*const float *A*/ims.data(), /*const int lda*/ ims.rows(),
                 /*const float beta*/ 0., /*float *C*/ cv.data(), /*const int ldc*/ cv.rows());
   
}   



/// A struct to hold the working memory for eigenSYEVR and maintain it between calls if desired.
template<typename sizeT, typename intT, typename floatT>
struct syevrMem
{
   sizeT sizeISuppZ;
   sizeT sizeWork;
   sizeT sizeIWork;
   
   intT *iSuppZ;
   
   floatT *work;
   
   intT *iWork;
   
   syevrMem()
   {
      sizeISuppZ = 0;
      sizeWork = 0;
      sizeIWork = 0;
      
      iSuppZ = 0;
      work = 0;
      iWork = 0;
   }
   
   ~syevrMem()
   {
      if(iSuppZ) ::free(iSuppZ);
      if(work) ::free(work);
      if(iWork) ::free(iWork);

   }
        
   void free()
   {
      if(iSuppZ) ::free(iSuppZ);
      sizeISuppZ = 0;
      
      if(work) ::free(work);
      sizeWork = 0;
      
      if(iWork) ::free(iWork);
      sizeIWork = 0;
   }
   
};

/// Calculate select eigenvalues and eigenvectors of an Eigen Array
/** Uses the templateLapack wrapper for syevr.
  * 
  * \tparam cvT is the scalar type of X (a.k.a. the covariance matrix)
  * \tparam calcT is the type in which to calculate the eigenvectors/eigenvalues
  * 
  * \param [in] X is a square matrix which is either upper or lower (default) triangular
  * \param [out] eigvec will contain the eigenvectors as columns
  * \param [out] eigval will contain the eigenvalues
  * \param [in] ev0 is the first desired eigenvalue default 0
  * \param [in] ev1 if >= ev0 thenthis is the last desired eigenvalue.  If -1 all eigenvalues are returned.
  * \param [in] UPLO specifies whether X is upper ('U') or lower ('L') triangular.  Default is ('L').
  * \param [in] mem holds the working memory arrays, can be re-passed to avoid unnecessary re-allocations
  * 
  * \retval the return code from syevr.
  */
template<typename cvT, typename calcT> //, typename eigenT>
int eigenSYEVR( Eigen::Array<cvT, Eigen::Dynamic, Eigen::Dynamic> &X, 
                Eigen::Array<calcT, Eigen::Dynamic, Eigen::Dynamic> &eigvec, 
                Eigen::Array<calcT, Eigen::Dynamic, Eigen::Dynamic> &eigval, 
                int ev0=0, 
                int ev1=-1, 
                char UPLO = 'L',
                syevrMem<int, int, calcT> * mem = 0
              ) 
{     
   //calcT *WORK;
   //int *ISUPPZ, *IWORK;
   int  numeig, info;//, sizeWORK, sizeIWORK;
   char RANGE = 'A';
   
   int localMem = 0;
   
   if(mem == 0)
   {
      mem = new syevrMem<int, int, calcT>;
      localMem = 1;
   }
   
   
   int n = X.rows();
   
   int IL = 1;
   int IU = n;
   if(ev0 >= 0 && ev1 >= ev0)
   {
      RANGE = 'I';
      IL = ev0+1; //This is FORTRAN, after all
      IU = ev1;
   }
   
   
   eigvec.resize(n,IU-IL+1);
   eigval.resize(n, 1); 
      
   //Copy X, casting to calcT
   Eigen::Array<calcT, Eigen::Dynamic, Eigen::Dynamic> Xc = X.template cast<calcT>();
   
   if( mem->sizeISuppZ < 2*n)
   {
      if(mem->iSuppZ) free(mem->iSuppZ);
         
      mem->sizeISuppZ = 2*n;
      mem->iSuppZ = (int *) malloc (mem->sizeISuppZ*sizeof(int));
   
      if ( mem->iSuppZ==NULL ) 
      {
         printf("malloc failed in eigenSYEVR\n"); 
         return 2;
      }
   }
   
   
   //  Allocate minimum allowed sizes for workspace
   int sizeWork = 26*n;
   calcT * work = (calcT *) malloc (sizeWork*sizeof(calcT));
   
   
   int sizeIWork = 10*n;
   int * iWork = (int *) malloc (sizeIWork*sizeof(int));

   //  Query for optimum sizes for workspace 
   info=syevr<calcT>('V', RANGE, UPLO, n, Xc.data(), n, 0, 0, IL, IU, lamch<float>('S'), &numeig, eigval.data(), eigvec.data(), n, mem->iSuppZ, work, -1, iWork, -1);

   // Now allocate optimum sizes
   /* -- tested increasing by x10, didn't improve performance at all 
    */
   if( mem->sizeWork <  ((int) work[0])*(1))
   {
      if(mem->work) free(mem->work);
  
      mem->sizeWork = ((int) work[0])*1;
      mem->work = (calcT *) malloc ((mem->sizeWork)*sizeof(calcT));
   }
   free(work);
   
   
   if(mem->sizeIWork < iWork[0]*1)
   {
      if(mem->iWork) free(mem->iWork);

      mem->sizeIWork = iWork[0]*1; 
      mem->iWork = (int *) malloc ((mem->sizeIWork)*sizeof(int));
   }   
   free(iWork);
   

   if ((mem->work==NULL)||(mem->iWork==NULL)) 
   {
      printf("malloc failed in eigenSYEVR\n"); 
      return 2;
   }
                
   // Now actually do the calculationg
   info=syevr<calcT>('V', RANGE, UPLO, n, Xc.data(), n, 0, 0, IL, IU, lamch<float>('S'), &numeig, eigval.data(), eigvec.data(), n, mem->iSuppZ, mem->work, mem->sizeWork, mem->iWork, mem->sizeIWork);     
   
    /*  Cleanup and exit  */
   //free(WORK); 
   //free(IWORK); 
   //free(ISUPPZ);
      
   if(localMem) delete mem;
   
   return info;
}       


template<typename imOutT, typename imInT>
void padImage(imOutT & imOut, imInT & imIn, int nrows, int ncols)
{
   imOut = imOutT::Zero(nrows, ncols);
   
   int stRow = (nrows - imIn.rows())/2;
   int stCol = (ncols - imIn.cols())/2;
   
   imOut.block(stRow, stCol, imIn.rows(),imIn.cols()) = imIn;
}

template<typename imOutT, typename imInT>
void padImageRef(imOutT imOut, imInT imIn, int nrows, int ncols)
{
   imOut = imOutT::Zero(nrows, ncols);
   
   int stRow = (nrows - imIn.rows())/2;
   int stCol = (ncols - imIn.cols())/2;
   
   imOut.block(stRow, stCol, imIn.rows(),imIn.cols()) = imIn;
}


template<typename imOutT, typename imInT>
void cutImage(imOutT & imOut, imInT & imIn, int nrows, int ncols)
{
   imOut = imOutT::Zero(nrows, ncols);
   
   int stRow = (imIn.rows()-nrows)/2;
   int stCol = (imIn.cols()-ncols)/2;
   
   imOut = imIn.block(stRow, stCol, nrows, ncols);
}


}//namespace mx

#endif //__eigenUtils_hpp__
