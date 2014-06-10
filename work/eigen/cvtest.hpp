

#include <Eigen/Dense>
#include "timeUtils.hpp"
#include "eigenImage.hpp"
#include "eigenUtils.hpp"
#include "pout.hpp"

extern "C"
{
#include "/usr/local/atlas/include/cblas.h"
}

#define NIMS (1000)
#define NPIXS (512*512)

using namespace Eigen;
using namespace mx;

void fullprodcv();



template<typename eigenT>
void diagprodcv(eigenT &ims, eigenImagef &cv)
{
   double t0, t1;

//    eigenImagef im1(NIMS, NPIXS);
// 
//    im1.setRandom();
// 
//    eigenImagef cv(NIMS, NIMS);//, cvt(NIMS, NIMS);
   
//    t0 = get_curr_time();
//    cvt = im1.matrix()*im1.matrix().transpose();
//    t1 = get_curr_time();
//    pout(t1-t0);

   //This is not multi-threaded
    t0 = get_curr_time();
    cv.setZero();
    cv.matrix().selfadjointView<Eigen::Lower>().rankUpdate(ims.matrix(), 1);
    t1 = get_curr_time();

//    t0 = get_curr_time();
//    
//    
//    for(int i=0;i<cv.rows(); ++i)
//    {
//       #pragma omp parallel for schedule(static, 1) num_threads(Eigen::nbThreads())
//       for(int j=0;j<cv.cols(); ++j)
//       {
//          if(j > i) break;
//          
//          cv(i,j) = im1.matrix().row(i).dot(im1.matrix().row(j));
//       }
//    }
//    t1 = get_curr_time();
         
   
   pout(t1-t0);
   
   //pout(cv); //, "\n\n\n");

//    pout("\n");
//    pout(cvt);
   
   
}

template<typename eigenT>
void diagprodcv_ssyrk(eigenT &ims, eigenT &cv)
{
   
   double t0, t1;

//    eigenImagef im1(NIMS, NPIXS);
// 
//    im1.setRandom();
// 
//    eigenImagef cv(NIMS, NIMS);
//    cv.setZero();
   
   t0 = get_curr_time();

   eigen_covar_ssyrk(cv, ims);
   
//    cblas_ssyrk(/*const enum CBLAS_ORDER Order*/ CblasColMajor, /*const enum CBLAS_UPLO Uplo*/ CblasLower,
//                  /*const enum CBLAS_TRANSPOSE Trans*/ CblasNoTrans, /*const int N*/ims.rows(), /*const int K*/ ims.cols(),
//                  /*const float alpha*/ 1.0, /*const float *A*/ims.data(), /*const int lda*/ ims.rows(),
//                  /*const float beta*/ 0., /*float *C*/ cv.data(), /*const int ldc*/ cv.rows());
   
   t1 = get_curr_time();
  
   pout(t1-t0);
  
  //pout(cv);
   
}   

