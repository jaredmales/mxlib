/** \file templateLevmar.cpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Templatized wrappers to the levmar minimization routines (definitions).
 * \ingroup fitting_files
 *
 */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#include "math/fit/templateLevmar.hpp"

extern "C"
{
#include "../../vendor/levmar-2.6/levmar.h"
}

namespace mx
{
namespace math
{
namespace fit 
{
   
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

} //namespace mx
} //namespace math
} //namespace fit 
