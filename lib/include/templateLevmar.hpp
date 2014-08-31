
#ifndef __templateLevmar_hpp__
#define __templateLevmar_hpp__

extern "C"
{
#include "levmar.h"
}

template<typename floatT>
int levmar_dif( void (*func)(floatT *p, floatT *hx, int m, int n, void *adata),
                floatT *p, 
                floatT *x, 
                int m, 
                int n, 
                int itmax, 
                floatT *opts,
                floatT *info, 
                floatT *work, 
                floatT *covar, 
                void *adata)
{
   return -1;
}

template<>
int levmar_dif<double>( void (*func)(double *p, double *hx, int m, int n, void *adata),
                        double *p, 
                        double *x, 
                        int m, 
                        int n, 
                        int itmax, 
                        double *opts,
                        double *info, 
                        double *work, 
                        double *covar, 
                        void *adata)
{
   return dlevmar_dif(func,p,x,m,n,itmax,opts,info,work,covar,adata);
}

template<>
int levmar_dif<float>( void (*func)(float *p, float *hx, int m, int n, void *adata),
                        float *p, 
                        float *x, 
                        int m, 
                        int n, 
                        int itmax, 
                        float *opts,
                        float *info, 
                        float *work, 
                        float *covar, 
                        void *adata)
{
   return slevmar_dif(func,p,x,m,n,itmax,opts,info,work,covar,adata);
}


template<typename floatT>
int levmar_der( void (*func)(floatT *p, floatT *hx, int m, int n, void *adata),
                void (*jacf)(floatT *p, floatT *j, int m, int n, void *adata),
                floatT *p, 
                floatT *x, 
                int m, 
                int n, 
                int itmax, 
                floatT *opts,
                floatT *info, 
                floatT *work, 
                floatT *covar, 
                void *adata)
{
   return -1;
}

template<>
int levmar_der<double>( void (*func)(double *p, double *hx, int m, int n, void *adata),
                        void (*jacf)(double *p, double *j, int m, int n, void *adata),
                        double *p, 
                        double *x, 
                        int m, 
                        int n, 
                        int itmax, 
                        double *opts,
                        double *info, 
                        double *work, 
                        double *covar, 
                        void *adata)
{
   return dlevmar_der(func,jacf,p,x,m,n,itmax,opts,info,work,covar,adata);
}


template<>
int levmar_der<float>( void (*func)(float *p, float *hx, int m, int n, void *adata),
                       void (*jacf)(float *p, float *j, int m, int n, void *adata),
                       float *p, 
                       float *x, 
                       int m, 
                       int n, 
                       int itmax, 
                       float *opts,
                       float *info, 
                       float *work, 
                       float *covar, 
                       void *adata)
{
   return slevmar_der(func,jacf,p,x,m,n,itmax,opts,info,work,covar,adata);
}

#endif // __templateLevmar_hpp__

