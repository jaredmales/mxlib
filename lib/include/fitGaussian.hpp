/** \file fitGaussian.hpp
 * \author Jared R. Males
 * \brief Tools for fitting Gaussians to data.
 * \ingroup peak_fit
 *
 */

#ifndef __fitGaussian_hpp__
#define __fitGaussian_hpp__

#include "levmarInterface.hpp"
#include "gaussian.hpp"

namespace mx
{

/** \addtogroup peak_fit
  * @{
  */

///\ref levmarInterface fitter structure for the symmetric Gaussian.
template<typename _floatT>
struct gaussian2D_sym_fitter
{
   typedef _floatT floatT;
   
   static const int nparams = 5;
   
   static void func(floatT *p, floatT *hx, int m, int n, void *adata)
   {
      array2Fit<floatT> * arr = (array2Fit<floatT> *) adata;
   
      size_t idx_mat, idx_dat;

      idx_dat = 0;
   
      for(int i=0;i<arr->nx; i++)
      {
         for(int j=0;j<arr->ny;j++)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = gaussian2D<floatT>(i,j,p[0],p[1], p[2], p[3], p[4]) - arr->data[idx_mat];
            
            idx_dat++;
         }
      }
   }
};

///\ref levmarInterface fitter structure for the general elliptical Gaussian.
template<typename _floatT>
struct gaussian2D_gen_fitter
{
   typedef _floatT floatT;
   
   static const int nparams = 7;
   
   //The finite difference method appears to be ~x2 faster
   //typedef bool hasJacobian;
   
   static void func(floatT *p, floatT *hx, int m, int n, void *adata)
   {
      array2Fit<floatT> * arr = (array2Fit<floatT> *) adata;
   
      size_t idx_mat, idx_dat;

      //Check for positive-definiteness of {{a b}{b c}}
      if( p[4]*p[6] - p[5]*p[5] <= 0)
      {
         idx_dat = 0;
         //If it's not positive-definite, then we just fill in with the value of the image itself.
         for(int i=0;i<arr->nx; ++i)
         {
            for(int j=0;j<arr->ny; ++j)
            { 
               idx_mat = i+j*arr->nx;
   
               hx[idx_dat] = arr->data[idx_mat];
                        
               ++idx_dat;
            }
         }         
         return;
      }
      
      //If positive-definite, now actually calculate
      idx_dat = 0;
   
      for(int i=0;i<arr->nx; ++i)
      {
         for(int j=0;j<arr->ny; ++j)
         { 
            idx_mat = i+j*arr->nx;
   
            hx[idx_dat] = gaussian2D<floatT>(i,j, p[0], p[1], p[2], p[3], p[4], p[5], p[6]) - arr->data[idx_mat];
                        
            ++idx_dat;
         }
      }
   }
   
   static void jacf(floatT *p, floatT *jacob, int m, int n, void *adata)
   {
      array2Fit<floatT> * arr = (array2Fit<floatT> *) adata;
   
      size_t idx_mat, idx_dat;

      floatT j_tmp[7];
      
      //Check for positive-definiteness of {{a b}{b c}}
      if( p[4]*p[6] - p[5]*p[5] <= 0)
      {
         idx_dat = 0;
         //If it's not positive-definite, then we just fill in with 0s.
         for(int i=0;i<arr->nx; ++i)
         {
            for(int j=0;j<arr->ny; ++j)
            {    
               for(int k=0; k< 7; ++k) jacob[idx_dat++] = 0;
            }
         }         
         return;
      }
      
      //If positive-definite, now actually calculate
      idx_dat = 0;
   
      for(int i=0;i<arr->nx; ++i)
      {
         for(int j=0;j<arr->ny; ++j)
         { 
            idx_mat = i+j*arr->nx;
   
            gaussian2D_jacobian<floatT>(j_tmp, i,j, p[0], p[1], p[2], p[3], p[4], p[5], p[6]); 
            for(int k=0;k<7;++k) jacob[idx_dat++] = j_tmp[k];
         }
      }
   }
   
   void paramNormalizer(floatT * p, int dir)
   {
      //Prepare for fit
      if(dir == 1)
      {
         floatT a, b, c;
         gaussian2D_rot2gen(a,b,c,p[4], p[5], p[6]);
         p[4] = a;
         p[5] = b;
         p[6] = c;
         return;
      }
      
      //Convert after fit
      if(dir == -1)
      {
         floatT sigma_x, sigma_y, theta;
         gaussian2D_gen2rot(sigma_x, sigma_y, theta, p[4], p[5], p[6]);
         p[4] = sigma_x;
         p[5] = sigma_y;
         p[6] = theta;
         return;
      }
   }
};


///Class to manage fitting a 2D Gaussian to data via the \ref levmarInterface
/** In addition to the requirements on fitterT specified by \ref levmarInterface
  * this class also requires this definition in fitterT
  * \code
  * static const int nparams = 7;
  * \endcode
  * where the number 7 is replaced by the number of parameters that fitterT expects to fit. 
  *
  * \tparam fitterT a type meeting the above requirements.
  */
template<typename fitterT>
class fitGaussian2D : public levmarInterface<fitterT>
{
   
public:
   
   typedef typename fitterT::floatT floatT;

   static const int nparams = fitterT::nparams;

   array2Fit<float> arr;
   
   void initialize()
   {
      this->allocate_params(nparams);
      this->adata = &arr;      
   }
   
   fitGaussian2D()
   {
      initialize();
   }
      
   ~fitGaussian2D()
   {
   }
   
   void setGuess(floatT G0, floatT A, floatT x0, floatT y0, floatT sigma_x, floatT sigma_y, floatT theta)
   {
      if(nparams != 7) return;
      
      this->p[0] = G0;
      this->p[1] = A;
      this->p[2] = x0;
      this->p[3] = y0;      
      this->p[4] = sigma_x;
      this->p[5] = sigma_y;
      this->p[6] = theta;
   }
   
   void setArray(floatT *data, int nx, int ny)
   {
      arr.data = data;
      arr.nx = nx;
      arr.ny = ny;
      
      this->n = nx*ny;
      
   }
   
   int fit()
   {
      fitterT fitter;
      fitter.paramNormalizer(this->p, 1);
      
      levmarInterface<fitterT>::fit();
      
      fitter.paramNormalizer(this->p, -1);
   }
      
   
};


///@}

} //namespace mx

#endif //__fitGaussian_hpp__

