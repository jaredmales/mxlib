/** \file templateLevmar.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Templatized wrappers to the levmar minimization routines..
 * \ingroup fitting_files
 *
 */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __templateLevmar_hpp__
#define __templateLevmar_hpp__


//*************************************************************************************************//
// copied from lm.h:

/* work arrays size for ?levmar_der and ?levmar_dif functions.
 * should be multiplied by sizeof(double) or sizeof(float) to be converted to bytes
 */
#define LM_DER_WORKSZ(npar, nmeas) (2*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))
#define LM_DIF_WORKSZ(npar, nmeas) (4*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))

/* work arrays size for ?levmar_bc_der and ?levmar_bc_dif functions.
 * should be multiplied by sizeof(double) or sizeof(float) to be converted to bytes
 */
#define LM_BC_DER_WORKSZ(npar, nmeas) (2*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))
#define LM_BC_DIF_WORKSZ(npar, nmeas) LM_BC_DER_WORKSZ((npar), (nmeas)) /* LEVMAR_BC_DIF currently implemented using LEVMAR_BC_DER()! */

/* work arrays size for ?levmar_lec_der and ?levmar_lec_dif functions.
 * should be multiplied by sizeof(double) or sizeof(float) to be converted to bytes
 */
#define LM_LEC_DER_WORKSZ(npar, nmeas, nconstr) LM_DER_WORKSZ((npar)-(nconstr), (nmeas))
#define LM_LEC_DIF_WORKSZ(npar, nmeas, nconstr) LM_DIF_WORKSZ((npar)-(nconstr), (nmeas))

/* work arrays size for ?levmar_blec_der and ?levmar_blec_dif functions.
 * should be multiplied by sizeof(double) or sizeof(float) to be converted to bytes
 */
#define LM_BLEC_DER_WORKSZ(npar, nmeas, nconstr) LM_LEC_DER_WORKSZ((npar), (nmeas)+(npar), (nconstr))
#define LM_BLEC_DIF_WORKSZ(npar, nmeas, nconstr) LM_LEC_DIF_WORKSZ((npar), (nmeas)+(npar), (nconstr))

/* work arrays size for ?levmar_bleic_der and ?levmar_bleic_dif functions.
 * should be multiplied by sizeof(double) or sizeof(float) to be converted to bytes
 */
#define LM_BLEIC_DER_WORKSZ(npar, nmeas, nconstr1, nconstr2) LM_BLEC_DER_WORKSZ((npar)+(nconstr2), (nmeas)+(nconstr2), (nconstr1)+(nconstr2))
#define LM_BLEIC_DIF_WORKSZ(npar, nmeas, nconstr1, nconstr2) LM_BLEC_DIF_WORKSZ((npar)+(nconstr2), (nmeas)+(nconstr2), (nconstr1)+(nconstr2))

#define LM_OPTS_SZ    	 5 /* max(4, 5) */
#define LM_INFO_SZ    	 10
#define LM_ERROR         -1
#define LM_INIT_MU    	 1E-03
#define LM_STOP_THRESH	 1E-17
#define LM_DIFF_DELTA    1E-06
#define LM_VERSION       "2.6 (November 2011)"

namespace mx
{
namespace math
{
namespace fit 
{
   
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
                void *adata
              );

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
                        void *adata
                      );

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
                        void *adata
                     );

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
                void *adata
              );

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
                        void *adata
                      );

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
                       void *adata
                     );

} //namespace mx
} //namespace math
} //namespace fit 

#endif // __templateLevmar_hpp__

