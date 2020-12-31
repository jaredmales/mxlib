/** \file levmarInterface.hpp
 * \author Jared R. Males
 * \brief A c++ interface to the templatized levmar minimization routines..
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

#ifndef levmarInterface_hpp
#define levmarInterface_hpp

#include <iostream>

#include "../../mxlib.hpp"

#include "templateLevmar.hpp"
#include "../../sys/timeUtils.hpp"


namespace mx
{
namespace math
{
namespace fit
{

//Forwards
template <typename T>
struct hasJacobian;

///A templatized interface to the levmar package
/** Requires a fitter class, which conforms to one of the following minimum specifications.
  * To use the finite difference jacobian calculation: 
  * \code
  * //This will cause the levmar_dif routine to be used
  * template<typename _realT>
  * struct dif_fitter
  * {
  *    typedef _realT realT; //required
  * 
  *    static void func(realT *p, realT *hx, int m, int n, void *adata)
  *     {
  *        //do stuff here . . .
  *     }
  * };
  * \endcode
  * If you wish to provide your own jacobian:
  * \code 
  * //This will cause the levmar_der routine to be used
  * template<typename _realT>
  * struct der_fitter
  * {
  *    typedef _realT realT; //required
  * 
  *    typdef bool hasJacobian; //this signals that jacf exists and should be used.
  * 
  *    static void func(realT *p, realT *hx, int m, int n, void *adata)
  *    {
  *        //do stuff here . . .
  *    }
  * 
  *    // This is the jacobian.
  *    static void jacf(realT *p, realT *j, int m, int n, void *adata)
  *    {
  *        //do stuff here . . .
  *    }
  * 
  * };
  * 
  * \endcode
  *
  * Note that if your fitter has a jacobian function (jacf), you must define
  * \code
  * typedef bool hasJacobian;
  * \endcode
  * for it to be used.
  * 
  * \tparam fitterT a class with at least the minimum interface described above.
  * 
  * \ingroup fitting
  */
template<class fitterT>
class levmarInterface
{

public:
   
   typedef typename fitterT::realT realT;
   
   
protected:
   realT *p; ///<  Parameter array.  On input is the initial estimates. On output has the estimated solution. 
   realT *init_p; ///< Parameter array on input, saved for comparison.
   
   int m;     ///<  Parameter vector dimension (i.e. number of unknowns)
   bool own_p; ///< Flag indicating whether the p array is owned by this object (for de-allocation).

   int itmax; ///< Maximum number of iterations, default is 100
      
public: 
   realT *x; ///<  I: measurement vector. NULL implies a zero vector 
   
   int n;     ///< I: measurement vector dimension 
   


   /// Options passed to the minimization routines.  See \ref set_opts for details.
   realT opts[LM_OPTS_SZ];    
   
   /// Information regarding the minimization. 
   /** See the levmar source code for documentation.  These fields are accessed by
     * \ref get_initial_norm, \ref get_final_norm, \ref get_iterations, \ref get_reason_code,
     * \ref get_reason_string, \ref get_fevals, \ref get_jevals, \ref get_nlinsys
     */ 
   realT info[LM_INFO_SZ];
   
   /// Elapsed time of the fitting procedure
   double deltaT;
                 
   /// Working memory passed to the levmar routines.
   /** From the levmar documentation: at least LM_DER/DIF_WORKSZ() reals large, allocated if NULL
     * Here this is always allocated by a call to \ref allocate_work.
     */ 
   realT *work;
   
   ///The current size of the work array
   int work_sz;  
   
   ///Covariance matrix corresponding to LS solution; mxm.
   /** Here this is allocated to size m-x-m by allocate_work, but if you don't want it 
     * allcoated and calculated 
     * set \ref getCovar to false and NULL will be passed.
     */ 
   realT *covar; 
   
   ///The current size of the covar array
   int covar_sz;  
   
   ///Controls whether the covar array is allocated. 
   bool getCovar;
   
   ///Pointer to possibly additional data, passed uninterpreted to func & jacf. 
   /** Set to NULL if not needed.
    */
   void *adata;     
   
private:
   ///Initialization common to all constructors
   void initialize();
   
public:
   
   /// Default constructor.
   /** With this constructor, you must set the parameter, data, and adata before calling fit().
     */
   levmarInterface();
   
   /// Setup constructor
   /** with this constructor fit() will work immediately.
     *
     */
   levmarInterface( realT *i_p, ///< [in] pointer to the initial parameter guess
                    realT *i_x, ///< [in] pointer to the data (can be NULL)
                    int i_m, ///< [in] the size of i_p
                    int i_n, ///< [in] the size of i_x
                    void *i_adata ///< [in] pointer to auxiliary data (can be NULL)
                  );

   ///Destructor
   /** Frees the work and covar matrices.
     */
   ~levmarInterface();
   
   ///Set number of parameters, but don't allocate
   void nParams(int i_m /**< [in] the number of parameters */);
   
   ///Get the current number of parameters
   /** \returns the current number of parameters (m)
     */
   int nParams();
   
   ///Allocate parameters array based on previous call to \ref nParams
   void allocate_params();
   
   ///Set number of parameters and allocate
   void allocate_params(int i_m);
   
   ///Point the parameter pointer at an externally allocated array
   void point_params(realT * i_p);
   
   ///Point the parameter pointer at an externally allocated array
   void point_params(realT * i_p, int i_m);
   
   ///Copy parameters to the parameter array
   /** This assumes that either the array was allocated (i.e. with \ref allocate_params) or
     * that \ref point_params has been called
     */ 
   void set_params(realT * i_p);
   
   ///Get current pointer array address
   realT * get_params();
   
   
   ///Set the maximum number of iterations
   /** Sets itmax.  Initialization default is itmax = 100
     * 
     * \param i_itmax the new value of itmax to set 
     */
   void set_itmax(int i_itmax);
   
   ///Get the maximum number of iterations
   int get_itmax();
   
   ///Allocate the work and covar matrices
   /** Uses a function object specialized for whether or not there is a jacobian
     * to determine the size of work. 
     */
   void allocate_work();
   
   ///Set one of the minimization options to val
   /** The options correspond to:
     * 
     * 0: the scale factor of the initial \f$ \mu \f$
     * 
     * 1: \f$ \epsilon_1 \f$ Stopping threshold for ||J^T e||_inf
     * 
     * 2: \f$ \epsilon_1 \f$ Stopping threshold for ||Dp||_2
     * 
     * 3: \f$ \epsilon_1 \f$ Stopping threshold for ||e||_2
     * 
     * 4: \f$ \Delta_{diff}\f$ Stepsize for finite differences
     * 
     * \param n is the option number
     * \param val is the value to set
     */
   void set_opts(int n, realT val);

   ///Set one or all of the minimization options to the default.
   /** See \ref set_opts for discription of the options.
     *
     * \param n the option number.  Pass -1 to set all options to defaults.
     */
   void set_opts_default(int n = -1);

   ///Get the current value of an option.
   /** See \ref set_opts for a description of the options
     *
     * \param n the option number
     */ 
   realT get_opts(int n);
   
   ///Perform the fit
   /** This calls \ref allocate_work, and then dispatches the levmar routine appropriate for fitterT.
     */
   int fit();
      
   ///Returns the L2-norm before minimization occurs
   realT get_initial_norm();
   
   ///Returns the L2-norm at the end of the minimization
   realT get_final_norm();
   
   ///Get the number of iterations taken during the minimization
   int get_iterations();
   
   ///Get a code specifying the reason minimization terminated.
   int get_reason_code();
   
   ///Get the descriptive string describing the reason minimization terminated.
   std::string get_reason_string();

   ///Get the number of function evaluations during the minimization
   int get_fevals();
   
   ///Get the number of jacobian evaluations during the minimization
   int get_jevals();

   ///Get the number of linear system solutions during the minimization
   int get_nlinsys();
   
   ///Get the elapsed time of the fit
   double get_deltaT()
   {
      return deltaT;
   }
   
   //Status reports
   
   ///Output current parameters to a stream
   /** Prints a formatted list of all current fit parameters.
     *
     * \tparam iosT is a std::ostream-like type.
     * \tparam comment is a comment character to start each line.  Can be '\0'.
     */ 
   template<typename iosT, char comment='#'>
   iosT & dumpParameters( iosT & ios /**< [in] a std::ostream-like stream. */);
   
   ///Dump the parameter vector to stdout.
   /** 
     * \tparam comment is a comment character to start each line.  Can be '\0'.
     */
   template<char comment='#'>
   std::ostream & dumpParameters();
   
   ///Output current parameters to a stream
   /** Prints a formatted list of all current fit parameters.
     *
     * \tparam iosT is a std::ostream-like type.
     * \tparam comment is a comment character to start each line.  Can be '\0'.
     */
   template<typename iosT, char comment='#'>
   iosT & dumpReport( iosT & ios, ///< [in] a std::ostream-like stream.
                      bool dumpParams = true ///< [in] [optional] whether or not to dump the parameters.
                    );
      
   ///Dump a status report to stdout
   /** 
     * \tparam comment is a comment character to start each line.  Can be '\0'.
     */ 
   template<char comment='#'>
   std::ostream & dumpReport( bool dumpParams = true /**< [in] [optional] whether or not to dump the parameters.*/);

};

template<class fitterT>
void levmarInterface<fitterT>::initialize()
{
   p=0;
   own_p = false;
   
   init_p = 0;
   
   m=0;
      
   x=0;
   n=0;
   set_opts_default(-1);//set all opts to defaults.
   work=0;
   work_sz = 0;
   covar=0;
   covar_sz = 0;
   getCovar = true;
   adata=0;
   
   for(int i=0;i<LM_INFO_SZ; ++i) info[i] = 0;
   deltaT = 0;
   
   itmax = 100;
}

template<class fitterT>
levmarInterface<fitterT>::levmarInterface()
{
   initialize();
}

template<class fitterT>
levmarInterface<fitterT>::levmarInterface(typename fitterT::realT *i_p,
                                              typename fitterT::realT *i_x,
                                              int i_m,
                                              int i_n,
                                              void * i_adata)
{
   initialize();
   
   p = i_p;
   x = i_x;
   m = i_m;
   n = i_n;
   adata = i_adata;
}

template<class fitterT>
levmarInterface<fitterT>::~levmarInterface()
{
   if(p && own_p) free(p);
   
   if(init_p) free(init_p);
   
   if(work) free(work);
   
   if(covar) free(covar);
}
   
template<class fitterT>
void levmarInterface<fitterT>::nParams(int i_m)
{
   //If we own and have allocated p, then de-alloc
   if(p && own_p) 
   {
      free(p);
      p = 0;
   }
   
   m = i_m;
   
   //Also allocate the init_p storage for initial guess
   if( init_p) 
   {
      free(init_p);
   }
   
   init_p =  (typename fitterT::realT *) malloc(sizeof(typename fitterT::realT) * m);

}

template<class fitterT>
int levmarInterface<fitterT>::nParams()
{
   return m;
}

   
template<class fitterT>   
void levmarInterface<fitterT>::allocate_params()
{
   if(p && own_p)
   {
      free(p);
   }
   
   p = (typename fitterT::realT *) malloc(sizeof(typename fitterT::realT) * m);

   own_p = true;
}

template<class fitterT>   
void levmarInterface<fitterT>::allocate_params(int i_m)
{
   nParams(i_m);
   
   allocate_params();
}
   
template<class fitterT>
void levmarInterface<fitterT>::point_params(realT * i_p)
{
   if(p && own_p)
   {
      free(p);
   }
   
   p = i_p;
}


template<class fitterT>
void levmarInterface<fitterT>::point_params(realT * i_p, int i_m)
{
   nParams(i_m);
   point_params(i_p);
}

template<class fitterT>
void levmarInterface<fitterT>::set_itmax(int i_itmax)
{
   itmax = i_itmax;
}
   
template<class fitterT>
int levmarInterface<fitterT>::get_itmax()
{
   return itmax;
}
   
template<class fitterT>
typename fitterT::realT * levmarInterface<fitterT>::get_params()
{
   return p;
}

template<class fitterT>
void levmarInterface<fitterT>::set_params(realT * i_p)
{
   for(int i=0;i<m;i++) p[i] = i_p[i];
}

   
//Functor which  is used by allocate() if hasJacobian is false 
template<class fitterT, bool jacf = hasJacobian<fitterT>::value>
struct levmar_allocate_size
{
   typedef typename fitterT::realT realT;
   
   int operator()(int m, int n)
   {
      return LM_DIF_WORKSZ(m,n);
   }
};

//Functor which  is used by allocate() if hasJacobian is true 
template<class fitterT>
struct levmar_allocate_size<fitterT, true>
{
   typedef typename fitterT::realT realT;
   
   int operator()(int m, int n)
   {
      return LM_DER_WORKSZ(m,n);
   }
};

template<class fitterT>
void levmarInterface<fitterT>::allocate_work()
{
   //Create function object to get allocation size, which depends on whether there is Jacobian.
   levmar_allocate_size<fitterT> alloc_sz;
   
   if(work_sz < alloc_sz(m,n) || !work)
   {
      if(work) free(work);
      
      work_sz = alloc_sz(m,n);
      work = (realT *) malloc( work_sz * sizeof(realT));
   }
   
   //Allocate if covar is desired and unallocated.
   if(getCovar)
   {
      if(covar_sz < m*m || !covar)
      {
         if(covar) free(covar);
      
         covar_sz = m*m;
         covar = (realT *) malloc(covar_sz * sizeof(realT));
      }
   }
   else
   {
      //If covar is not desired, de-allocate if allocated before.
      if(covar) free(covar);
      covar = 0;
   }
}


template<class fitterT>
void levmarInterface<fitterT>::set_opts(int n, realT val)
{
   opts[n] = val;
}


template<class fitterT>
void levmarInterface<fitterT>::set_opts_default(int n)
{
   //See the source in lm_core.c
   
   if(n == 0 || n == -1)
   {
      opts[0] = LM_INIT_MU;
   }
   
   if(n > 0 && n < 4)
   {
      opts[n] = LM_STOP_THRESH;
   }
   else if(n == -1)
   {
      opts[1] = LM_STOP_THRESH;
      opts[2] = LM_STOP_THRESH;
      opts[3] = LM_STOP_THRESH;
   }
   
   if(n == 4  || n == -1)
   {
      opts[4] = LM_DIFF_DELTA;
   }
}

template<class fitterT>
typename fitterT::realT levmarInterface<fitterT>::get_opts(int n)
{
   return opts[n];
}
   
//Functor which  is used by fit() if hasJacobian is false 
template<class fitterT, bool jacf = hasJacobian<fitterT>::value>
struct do_levmar
{
   typedef typename fitterT::realT realT;
   
   int operator()(realT *p, 
                  realT *x, 
                  int m, 
                  int n, 
                  int itmax, 
                  realT *opts,
                  realT *info, 
                  realT *work, 
                  realT *covar, 
                  void *adata)
   {
      return levmar_dif<realT>( &fitterT::func, p, x, m, n, itmax, opts, info, work, covar, adata);
   }
};

//Functor which  is used by fit() if hasJacobian is true 
template<class fitterT>
struct do_levmar<fitterT, true>
{
   typedef typename fitterT::realT realT;
      
   int operator()(realT *p, 
                  realT *x, 
                  int m, 
                  int n, 
                  int itmax, 
                  realT *opts,
                  realT *info, 
                  realT *work, 
                  realT *covar, 
                  void *adata)
   {
      return levmar_der<realT>( &fitterT::func, &fitterT::jacf, p, x, m, n, itmax, opts, info, work, covar, adata);
   }
};


template<class fitterT>
int levmarInterface<fitterT>::fit()
{
   realT * _opts;
   
   allocate_work();
  
   if(opts[0] == 0) _opts= 0;
   else _opts = opts;
  
   //These may not be updated by the levmar library
   info[8] = 0;
   info[9] = 0;
   
   
   if( !init_p) 
   {
      init_p =  (typename fitterT::realT *) malloc(sizeof(typename fitterT::realT) * m);
   }
   
   for(int i = 0; i < m; ++i) init_p[i] = p[i];
   
   double t0 = sys::get_curr_time();
   
   //Create one of the above functors, which depends on whether fitterT has a Jacobian.
   do_levmar<fitterT> fitter;
   
   fitter(p,x,m,n,itmax,_opts,info,work,covar,adata);
   
   deltaT = sys::get_curr_time() - t0;
   
   return 0;
}



template<class fitterT>
typename fitterT::realT levmarInterface<fitterT>::get_initial_norm()
{
   return info[0];
}
   
template<class fitterT>
typename fitterT::realT levmarInterface<fitterT>::get_final_norm()
{
   return info[1];
}

template<class fitterT>
int levmarInterface<fitterT>::get_iterations()
{
   return (int) info[5];
}
   
template<class fitterT>
int levmarInterface<fitterT>::get_reason_code()
{
   return (int) info[6];
}
  
template<class fitterT>
int levmarInterface<fitterT>::get_fevals()
{
   return (int) info[7];
}
   
template<class fitterT>
int levmarInterface<fitterT>::get_jevals()
{
   return (int) info[8];
}

template<class fitterT>
int levmarInterface<fitterT>::get_nlinsys()
{
   return (int) info[9];
}

template<class fitterT>
std::string levmarInterface<fitterT>::get_reason_string()
{
   switch(get_reason_code())
   {
      case 1:
         return "stopped by small gradient J^T e";
         break;
      case 2:
         return "stopped by small Dp";
         break;
      case 3:
         return "stopped by itmax";
         break;
      case 4:
         return "singular matrix. Restart from current p with increased mu";
         break;
      case 5:
         return "no further error reduction is possible. Restart with increased mu";
         break;
      case 6:
         return "stopped by small ||e||_2";
         break;
      case 7:
         return "stopped by invalid (i.e. NaN or Inf) \"func\" values. This is a user error";
         break;
      default:
         return "unknown reason code";
   }
     
}

template<class fitterT>
template<typename iosT, char comment>
iosT & levmarInterface<fitterT>::dumpParameters( iosT & ios )
{
   //This causes the stream to not output a '\0'
   char c[] = {comment, '\0'};
   
   ios << c << "Current parameters (initial):\n";
   for(int i=0;i<m;i++) 
   {
      ios << c;
      ios <<  "p[" << i << "] = " << p[i] << " (" << init_p[i] << ")\n";
   }
   
   return ios;
}

template<class fitterT>
template<char comment>
std::ostream & levmarInterface<fitterT>::dumpParameters()
{
   return dumpParameters<std::ostream,comment>(std::cout);
}

template<class fitterT>
template<typename iosT, char comment>
iosT & levmarInterface<fitterT>::dumpReport( iosT & ios,
                                             bool dumpParams
                                           )
{
   char c[] = {comment, '\0'}; //So a '\0' won't be written to stream.
   
   ios << c << "--------------------------------------\n";
   ios << c << "mx::math::fit::levmarInterface Results \n";
   ios << c << "--------------------------------------\n";
   if(dumpParams) dumpParameters<iosT,comment>(ios);
   ios << c << "Reason for termination: " << get_reason_string() << "\n";
   ios << c << "Initial norm: " << get_initial_norm() << "\n";
   ios << c << "Final norm: " << get_final_norm() << "\n";
   ios << c << "Number of iterations: " << get_iterations() << "\n";
   ios << c << "Function evals: " << get_fevals() << "\n";
   ios << c << "Jacobian evals: " << get_jevals() << "\n";
   ios << c << "Elapsed time: " << get_deltaT() << " secs\n";
   dumpGitStatus<iosT,comment>(ios);
   return ios;
   
}
   
template<class fitterT>
template<char comment>
std::ostream & levmarInterface<fitterT>::dumpReport( bool dumpParams )
{
   return dumpReport<std::ostream, comment>(std::cout);
   
}
   
   
///Test whether a function type has a Jacobian function by testing whether it has a typedef of "hasJacobian"
/** Used for compile-time determination of whether the fitter has a Jacobian.
  * 
  * \ingroup fitting
  */
template <typename T>
struct hasJacobian //This was taken directly from the example at http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
{
   // Types "yes" and "no" are guaranteed to have different sizes,
   // specifically sizeof(yes) == 1 and sizeof(no) == 2.
   typedef char yes[1];
   typedef char no[2];
 
   template <typename fitterT>
   static yes& test(typename fitterT::hasJacobian*);
 
   template <typename>
   static no& test(...);
 
   /// If hasJacobian<fitterT>::value == true, then fitterT has a Jacobian and the appropriate levmar routines are used.
   /// If ::value == false, then the numerical derivatives are calculated.
   // If the "sizeof" of the result of calling test<T>(0) would be equal to sizeof(yes),
   // the first overload worked and T has a nested type named "hasJacobian".
   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};

///Wrapper for a native array to pass to \ref levmarInterface
/**
  * \ingroup fitting
  */ 
template<typename realT>
struct array2Fit
{
   realT * data {0}; ///Pointer to the array
   size_t nx {0}; ///X dimension of the array
   size_t ny {0}; ///Y dimension of the array
};


} //namespace fit 
} //namespace math
} //namespace mx

#endif //levmarInterface_hpp

