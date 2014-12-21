/** \file levmarInterface.hpp
 * \author Jared R. Males
 * \brief A templatized interface to the levmar minimization routines..
 * \ingroup fitting
 *
 */

#ifndef __levmarInterface_hpp__
#define __levmarInterface_hpp__


#include "templateLevmar.hpp"

namespace mx
{

//Test whether a function type has a Jacobian function by testing whether it has a typedef of "hasJacobian"
/* Used for compile-time determination of type
  */
//This was taken directly from the example at http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
template <typename T>
struct hasJacobian
{
   // Types "yes" and "no" are guaranteed to have different sizes,
   // specifically sizeof(yes) == 1 and sizeof(no) == 2.
   typedef char yes[1];
   typedef char no[2];
 
   template <typename fitterT>
   static yes& test(typename fitterT::hasJacobian*);
 
   template <typename>
   static no& test(...);
 
   // If the "sizeof" of the result of calling test<T>(0) would be equal to sizeof(yes),
   // the first overload worked and T has a nested type named "is_mmatrix".
   static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};

///Wrapper for a native array to pass to \ref levmarInterface
template<typename floatT>
struct array2Fit
{
   floatT * data;
   size_t nx;
   size_t ny;
};

/** \addtogroup fitting
  * @{
  */

///An interface to the levmar package
/** Requires a fitter class, which conforms to one of the following minimum specifications.
  * To use the finite difference jacobian calculation: 
  * \code
  * //This will cause the levmar_dif routine to be used
  * template<typename _floatT>
  * struct dif_fitter
  * {
  *    typedef _floatT floatT; //required
  * 
  *    static void func(floatT *p, floatT *hx, int m, int n, void *adata)
  *     {
  *        //do stuff here . . .
  *     }
  * };
  * \endcode
  * If you wish to provide your own jacobian:
  * \code 
  * //This will cause the levmar_der routine to be used
  * template<typename _floatT>
  * struct der_fitter
  * {
  *    typedef _floatT floatT; //required
  * 
  *    typdef bool hasJacobian; //this signals that jacf exists and should be used.
  * 
  *    static void func(floatT *p, floatT *hx, int m, int n, void *adata)
  *    {
  *        //do stuff here . . .
  *    }
  * 
  *    static void jacf(floatT *p, floatT *j, int m, int n, void *adata)
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
  */
template<class fitterT>
class levmarInterface
{

public:
   
   typedef typename fitterT::floatT floatT;
   
protected:
   floatT *p; ///<  Parameter array.  On input is the initial estimates. On output has the estimated solution. 
   int m;     ///<  Parameter vector dimension (i.e. #unknowns)
   bool own_p; ///< Flag indicating whether the p array is owned by this object (for de-allocation).

   int itmax; ///< Maximum number of iterations, default is 100
      
public: 
   floatT *x; ///<  I: measurement vector. NULL implies a zero vector 
   
   int n;     ///< I: measurement vector dimension 
   


   /// Options passed to the minimization routines.  See \ref set_opts for details.
   floatT opts[LM_OPTS_SZ];    
   
   /// Information regarding the minimization. 
   /** See the levmar source code for documentation.  These fields are accessed by
     * \ref get_initial_norm, \ref get_final_norm, \ref get_iterations, \ref get_reason_code,
     * \ref get_reason_string, \ref get_fevals, \ref get_jevals, \ref get_nlinsys
     */ 
   floatT info[LM_INFO_SZ];
   
   /// Elapsed time of the fitting procedure
   double deltaT;
                 
   /// Working memory passed to the levmar routines.
   /** From the levmar documentation: at least LM_DER/DIF_WORKSZ() reals large, allocated if NULL
     * Here this is always allocated by a call to \ref allocate_work.
     */ 
   floatT *work;
   
   ///The current size of the work array
   int work_sz;  
   
   ///Covariance matrix corresponding to LS solution; mxm.
   /** Here this is allocated to size m-x-m by allocate_work, but if you don't want it 
     * allcoated and calculated 
     * set \ref getCovar to false and NULL will be passed.
     */ 
   floatT *covar; 
   
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
   levmarInterface();
   
   /// Setup constructor
   /** After construction ::fit() will work.
     *
     * \param i_p pointer to the initial parameter guess
     * \param i_x pointer to the data (can be NULL)
     * \param i_m the size of i_p
     * \param i_n the size of i_x
     * \param i_adata pointer to auxiliary data (can be NULL)
     */
   levmarInterface( floatT *i_p,
                     floatT *i_x,
                     int i_m,
                     int i_n,
                     void *i_adata);

   ///Destructor
   /** Frees the work and covar matrices.
     */
   ~levmarInterface();
   
   ///Set number of parameters, but don't allocate
   void set_nparams(int i_m);
   
   ///Get the current number of parameters
   int get_nparams();
   
   ///Allocate parameters array based on previous call to \ref set_nparams
   void allocate_params();
   
   ///Set number of parameters and allocate
   void allocate_params(int i_m);
   
   ///Point the parameter pointer at an externally allocated array
   void point_params(floatT * i_p);
   
   ///Point the parameter pointer at an externally allocated array
   void point_params(floatT * i_p, int i_m);
   
   ///Copy parameters to the parameter array
   /** This assumes that either the array was allocated (i.e. with \ref allocate_params) or
     * that \ref point_params has been called
     */ 
   void set_params(floatT * i_p);
   
   ///Get current pointer array address
   floatT * get_params();
   
   
   ///Set the maximum number of iterations
   /** Sets itmax.  Initizlization default is itmax = 100
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
   /** The options correspond too:
     * 
     * 0: the scale factor the initial \f$ \mu \f$
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
   void set_opts(int n, floatT val);

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
   floatT get_opts(int n);
   
   ///Perform the fit
   /** This calls \ref allocate_work, and then dispatches the levmar routine appropriate for fitterT.
     */
   int fit();
      
   ///Returns the L2-norm before minimization occurs
   floatT get_initial_norm();
   
   ///Returns the L2-norm at the end of the minimization
   floatT get_final_norm();
   
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
   
   ///Dump the parameter vector to stdout.
   void dump_params()
   {
      std::cout << "Current parameters:\n";
      for(int i=0;i<m;i++) std::cout << "p[" << i << "] = " << p[i] << "\n";
   }
   
   ///Dump a status report to stdout
   void dump_report()
   {
      std::cout << "-------------------------------------------\n";
      dump_params();
      std::cout << "Reason for termination: " << get_reason_string() << "\n";
      std::cout << "Initial norm: " << get_initial_norm() << "\n";
      std::cout << "Final norm: " << get_final_norm() << "\n";
      std::cout << "Number of iterations: " << get_iterations() << "\n";
      std::cout << "Function evals: " << get_fevals() << "\n";
      std::cout << "Jacobian evals: " << get_jevals() << "\n";
      std::cout << "Elapsed time: " << get_deltaT() << " secs\n";
   }
      
};

template<class fitterT>
void levmarInterface<fitterT>::initialize()
{
   p=0;
   own_p = false;
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
levmarInterface<fitterT>::levmarInterface(typename fitterT::floatT *i_p,
                                              typename fitterT::floatT *i_x,
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
   if(work) free(work);
   if(covar) free(covar);
}
   
template<class fitterT>
void levmarInterface<fitterT>::set_nparams(int i_m)
{
   //If we own and have allocated p, then de-alloc
   if(p && own_p) 
   {
      free(p);
      p = 0;
   }
   
   m = i_m;
}

template<class fitterT>
int levmarInterface<fitterT>::get_nparams()
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
   
   p = (typename fitterT::floatT *) malloc(sizeof(typename fitterT::floatT) * m);

   own_p = true;
}

template<class fitterT>   
void levmarInterface<fitterT>::allocate_params(int i_m)
{
   set_nparams(i_m);
   
   allocate_params();
}
   
template<class fitterT>
void levmarInterface<fitterT>::point_params(floatT * i_p)
{
   if(p && own_p)
   {
      free(p);
   }
   
   p = i_p;
}


template<class fitterT>
void levmarInterface<fitterT>::point_params(floatT * i_p, int i_m)
{
   set_nparams(i_m);
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
typename fitterT::floatT * levmarInterface<fitterT>::get_params()
{
   return p;
}

template<class fitterT>
void levmarInterface<fitterT>::set_params(floatT * i_p)
{
   for(int i=0;i<m;i++) p[i] = i_p[i];
}

   
//Functor which  is used by allocate() if hasJacobian is false 
template<class fitterT, bool jacf = hasJacobian<fitterT>::value>
struct levmar_allocate_size
{
   typedef typename fitterT::floatT floatT;
   
   int operator()(int m, int n)
   {
      return LM_DIF_WORKSZ(m,n);
   }
};

//Functor which  is used by allocate() if hasJacobian is true 
template<class fitterT>
struct levmar_allocate_size<fitterT, true>
{
   typedef typename fitterT::floatT floatT;
   
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
      work = (floatT *) malloc( LM_DIF_WORKSZ(m,n) * sizeof(floatT));
   }
   
   //Allocate if covar is desired and unallocated.
   if(getCovar)
   {
      if(covar_sz < m*m || !covar)
      {
         if(covar) free(covar);
      
         covar_sz = m*m;
         covar = (floatT *) malloc(covar_sz * sizeof(floatT));
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
void levmarInterface<fitterT>::set_opts(int n, floatT val)
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
typename fitterT::floatT levmarInterface<fitterT>::get_opts(int n)
{
   return opts[n];
}
   
//Functor which  is used by fit() if hasJacobian is false 
template<class fitterT, bool jacf = hasJacobian<fitterT>::value>
struct do_levmar
{
   typedef typename fitterT::floatT floatT;
   
   int operator()(floatT *p, 
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
      return ::levmar_dif<floatT>( &fitterT::func, p, x, m, n, itmax, opts, info, work, covar, adata);
   }
};

//Functor which  is used by fit() if hasJacobian is true 
template<class fitterT>
struct do_levmar<fitterT, true>
{
   typedef typename fitterT::floatT floatT;
      
   int operator()(floatT *p, 
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
      return ::levmar_der<floatT>( &fitterT::func, &fitterT::jacf, p, x, m, n, itmax, opts, info, work, covar, adata);
   }
};


template<class fitterT>
int levmarInterface<fitterT>::fit()
{
   floatT * _opts;
   
   allocate_work();
  
   if(opts[0] == 0) _opts= 0;
   else _opts = opts;
  
   //These may not be updated by the levmar library
   info[8] = 0;
   info[9] = 0;
   
   
   double t0 = get_curr_time();
   
   //Create one of the above functors, which depends on whether fitterT has a Jacobian.
   do_levmar<fitterT> fitter;
   
   
   
   fitter(p,x,m,n,itmax,_opts,info,work,covar,adata);
   
   deltaT = get_curr_time() - t0;
}



template<class fitterT>
typename fitterT::floatT levmarInterface<fitterT>::get_initial_norm()
{
   return info[0];
}
   
template<class fitterT>
typename fitterT::floatT levmarInterface<fitterT>::get_final_norm()
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

///@}

} //namespace mx

#endif //__levmarInterface_hpp__
