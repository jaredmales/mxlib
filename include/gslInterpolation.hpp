/** \file gslInterpolation.hpp
  * \brief Wrappers for using the GNU Scientific Library 1-D interpolation functions
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __gslInterpolation_hpp__
#define __gslInterpolation_hpp__

#include <gsl/gsl_interp.h>

namespace mx
{

struct gslInterpolator
{
   gsl_interp * interp;
   gsl_interp_accel * acc;
   
   double *_xin;
   double *_yin;
   
   void setup(const gsl_interp_type * interpT, double *xin, double *yin, size_t Nin)
   {
      interp =  gsl_interp_alloc(interpT, Nin);
      acc = gsl_interp_accel_alloc ();

      gsl_interp_init(interp, xin, yin, Nin);
   
      gsl_interp_accel_reset(acc);
      
      _xin = xin;
      _yin = yin;
   }
   
   void setup(const gsl_interp_type * interpT, std::vector<double> & xin, std::vector<double> & yin)
   {
      setup(interpT, xin.data(), yin.data(), xin.size());
   }
   
   gslInterpolator()
   {
   }
   
   gslInterpolator(const gsl_interp_type * interpT, double *xin, double *yin, size_t Nin)
   {
      setup(interpT, xin,yin,Nin);
   }
   
   gslInterpolator(const gsl_interp_type * interpT, std::vector<double> & xin, std::vector<double> & yin)
   {
      setup(interpT, xin.data(), yin.data(), xin.size());
   }
   
   ~gslInterpolator()
   {
      gsl_interp_free(interp);
      gsl_interp_accel_free (acc);
   }
   
   double interpolate(double & x)
   {
      double y;
      gsl_interp_eval_e (interp, _xin, _yin, x, acc, &y);
      return y;
   }
   
};
///Interpolate a 1-D data X vs Y discrete function onto a new X axis
/**
  * \param interpT of the <a href="https://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html#Interpolation-Types">gsl interpolation types</a>.
  * \param [in] xin the input x-axis
  * \param [in] yin the input y-values
  * \param [in] Nin the size of the input x-y data series
  * \param [in] xout the desired x-axis
  * \param [out] yout the output interpolated y-values
  * \param [in] Nout the size of the output x-y axis
  * 
  * \retval
  */
int gsl_interpolate( const gsl_interp_type * interpT, 
                     double *xin, 
                     double *yin, 
                     size_t Nin, 
                     double *xout, 
                     double *yout, 
                     size_t Nout )
{
   gsl_interp * interp =  gsl_interp_alloc(interpT, Nin);
   gsl_interp_accel * acc = gsl_interp_accel_alloc ();

   gsl_interp_init(interp, xin, yin, Nin);
   
   gsl_interp_accel_reset(acc);
   
   for(int i=0;i<Nout; ++i)
   {
      int e = gsl_interp_eval_e (interp, xin, yin, xout[i], acc, &yout[i]);
      //std::cout << xout[i] << " " << yin[i] << "\n";
   }
   
   gsl_interp_free(interp);
   gsl_interp_accel_free (acc);
}

///Interpolate a 1-D data X vs Y discrete function onto a new X axis (vector version)
int gsl_interpolate( const gsl_interp_type * interpT,
                     std::vector<double> & xin,
                     std::vector<double> & yin,
                     std::vector<double> & xout,
                     std::vector<double> & yout )
{
   yout.resize(xout.size());
   
   return gsl_interpolate(interpT, xin.data(), yin.data(), xin.size(), xout.data(), yout.data(), xout.size());
}

} //namespace mx

#endif //__gslInterpolation_hpp__
