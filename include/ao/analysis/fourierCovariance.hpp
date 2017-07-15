/** \file fourierCovariance.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculation of the modal covariance in the Fourier basis.
  * \ingroup mxAO_files
  * 
  */

#ifndef __fourierCovariance_hpp__
#define __fourierCovariance_hpp__


#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;


#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


#include "../../math/func/jinc.hpp"
#include "../../fourierModes.hpp"
#include "../../improc/fitsFile.hpp"
#include "../../improc/eigenCube.hpp"
#include "../../mxlib_uncomp_version.h"
#include "../../mxlib.h"
#include "../../ompLoopWatcher.hpp"
#include "../../timeUtils.hpp"

#include "aoAtmosphere.hpp"
#include "aoPSDs.hpp"
#include "aoSystem.hpp"
#include "mxaoa_version.h"


namespace mx
{
namespace AO
{
   
#ifndef WSZ

/** \def WSZ
  * \brief The size of the gsl_integration workspaces.
  */ 
#define WSZ 100000

#endif


//Forward Declaration
template<typename floatT, typename aosysT>
struct fourierCovariance;

///Worker for the azimuthal integral (in phi) for the basic Fourier mode covariance.
/** 
  * \param phi the angle at which to evaluate the integrand
  * \param params a pointer to a object of type fourierCovariance<floatT, aosyT>
  *
  * \tparam floatT a floating point type used for all calculations
  * \tparam aosysT the type of the AO system structure 
  * 
  * \ingroup mxAOAnalytic
  */
template<typename floatT, typename aosysT>
floatT phiInt_basic (floatT phi, void * params) 
{
   
   fourierCovariance<floatT, aosysT> * Pp = (fourierCovariance<floatT, aosysT> *) params;
   
   floatT L0 = Pp->aosys->atm.L_0(); 
   floatT D = Pp->aosys->D();

   int p = Pp->p;
   floatT m = Pp->m;
   floatT n = Pp->n;

   int pp = Pp->pp;
   floatT mp = Pp->mp;
   floatT np = Pp->np;

   floatT k = Pp->k;

   floatT cosp = cos(phi);
   floatT sinp = sin(phi);

   /*** no prime ***/
   floatT kmn_p, kmn_m;
   floatT Ji_mn_p, Ji_mn_m;
   floatT Q_mn;//, Qs_mn;

   kmn_p = sqrt( pow(k*cosp + m/D, 2) + pow(k*sinp + n/D, 2));
   kmn_m = sqrt( pow(k*cosp - m/D, 2) + pow(k*sinp - n/D, 2));

   Ji_mn_p = math::func::jinc(pi<floatT>()*D*kmn_p);
   Ji_mn_m = math::func::jinc(pi<floatT>()*D*kmn_m);

   
   floatT N = 1./sqrt(0.5 + p*math::func::jinc(2*pi<floatT>()*sqrt(m*m+n*n)));
   
   Q_mn = N*(Ji_mn_p + p*Ji_mn_m);

   /*** primed ***/
   floatT kmpnp_p, kmpnp_m;
   floatT Ji_mpnp_p, Ji_mpnp_m;
   floatT Q_mpnp;//, Qs_mpnp;

   kmpnp_p = sqrt( pow(k*cosp + mp/D, 2) + pow(k*sinp + np/D, 2));
   kmpnp_m = sqrt( pow(k*cosp - mp/D, 2) + pow(k*sinp - np/D, 2));

   
   Ji_mpnp_p = math::func::jinc(pi<floatT>()*D*kmpnp_p);
   Ji_mpnp_m = math::func::jinc(pi<floatT>()*D*kmpnp_m);

   floatT Np = 1./sqrt(0.5 + pp*math::func::jinc(2*pi<floatT>()*sqrt(mp*mp+np*np)));
   
   Q_mpnp = Np*(Ji_mpnp_p + pp*Ji_mpnp_m);
  
   floatT P = Pp->aosys->psd(Pp->aosys->atm, k); //vonKarmanPSD(k, D, L0, Pp->subPiston, Pp->subTipTilt);
   
   return P*k*(Q_mn*Q_mpnp);
}

///Worker for the azimuthal integral (in phi) for the modified Fourier mode covariance.
/** 
  * \param phi the angle at which to evaluate the integrand
  * \param params a pointer to a object of type fourierCovariance<floatT, aosyT>
  *
  * \tparam floatT a floating point type used for all calculations
  * \tparam aosysT the type of the AO system structure 
  */
template<typename floatT, typename aosysT>
floatT phiInt_mod (floatT phi, void * params) 
{
   
   fourierCovariance<floatT, aosysT> * Pp = (fourierCovariance<floatT, aosysT> *) params;
   
   floatT L0 = Pp->aosys->atm.L_0(); 
   floatT D = Pp->aosys->D();

   int p = Pp->p;
   floatT m = Pp->m;
   floatT n = Pp->n;

   int pp = Pp->pp;
   floatT mp = Pp->mp;
   floatT np = Pp->np;

   floatT k = Pp->k;

   floatT cosp = cos(phi);
   floatT sinp = sin(phi);

   /*** no prime ***/
   floatT kmn_p, kmn_m;
   floatT Ji_mn_p, Ji_mn_m;
   floatT Qc_mn, Qs_mn;

   kmn_p = sqrt( pow(k*cosp + m/D, 2) + pow(k*sinp + n/D, 2));
   kmn_m = sqrt( pow(k*cosp - m/D, 2) + pow(k*sinp - n/D, 2));

   Ji_mn_p = math::func::jinc(pi<floatT>()*D*kmn_p);
   Ji_mn_m = math::func::jinc(pi<floatT>()*D*kmn_m);


   /*** primed ***/
   floatT kmpnp_p, kmpnp_m;
   floatT Ji_mpnp_p, Ji_mpnp_m;
   floatT Qc_mpnp, Qs_mpnp;

   kmpnp_p = sqrt( pow(k*cosp + mp/D, 2) + pow(k*sinp + np/D, 2));
   kmpnp_m = sqrt( pow(k*cosp - mp/D, 2) + pow(k*sinp - np/D, 2));

   
   Ji_mpnp_p = math::func::jinc(pi<floatT>()*D*kmpnp_p);
   Ji_mpnp_m = math::func::jinc(pi<floatT>()*D*kmpnp_m);
   
   floatT QQ;
   
   if(p == pp)
   {
      QQ = 2.0*( Ji_mn_p*Ji_mpnp_p + Ji_mn_m*Ji_mpnp_m);
   }
   else
   {
      QQ = 2.0*( Ji_mn_p*Ji_mpnp_m + Ji_mn_m*Ji_mpnp_p);
   }
  
   floatT P = Pp->aosys->psd(Pp->aosys->atm, k); //vonKarmanPSD(k, D, L0, Pp->subPiston, Pp->subTipTilt);
   return P*k*QQ;
}



///Worker function for the radial integral in the covariance calculation
/** 
  * \param k the spatial frequency at which to evaluate the integrand
  * \param params a pointer to a object of type fourierCovariance<floatT, aosyT>
  *
  * \tparam floatT a floating point type used for all calculations.  As of Nov 2016 must be double due to gsl_integration.
  * \tparam aosysT the type of the AO system structure 
  */
template<typename floatT, typename aosysT>
floatT kInt (floatT k, void * params) 
{
   fourierCovariance<floatT, aosysT> * Pp = (fourierCovariance<floatT, aosysT> *) params;
     
   floatT result, error;

   gsl_function func;
   
   //Here choose between basic and modified Fourier modes.
   if(Pp->useBasic)
   {
      func.function = &phiInt_basic<floatT, aosysT>;
   }
   else
   {
      func.function = &phiInt_mod<floatT, aosysT>;
   }
   func.params = Pp;
   
   Pp->k = k;
   
   //Tolerances of 1e-4 seem to be all we can ask for
   gsl_integration_qag(&func, 0., 2*pi<double>(), 1e-10, 1e-4, WSZ, GSL_INTEG_GAUSS21, Pp->phi_w, &result, &error);

   
   return result;
}

///Structure to manage the Fourier mode covariance calculation
/** 
  * \tparam floatT a floating point type used for all calculations.  As of Nov 2016 must be double due to gsl_integration.
  * \tparam aosysT the type of the AO system structure 
  */
template<typename floatT, typename aosysT>
struct fourierCovariance
{
   ///Pointer to an AO system, which contains the relevant spatial PSD of turbulence.
   aosysT * aosys;

   ///Flag controlling use of basic or modified Fourier modes.  If true, the basic sin/cos modes are used.  If false (default), the modified modes are u sed.
   bool useBasic;
      
   ///p-index of the unprimed mode.  +/-1 for modified modes.  If basic, then +1==>cosine, -1==>sine. 
   int p;
   
   ///The m-index of the unprimed mode, corresponding to the \f$ k_u = m/D \f$ component of spatial frequency.
   floatT m;
   
   ///The n-indexof the unprimed mode, corresponding to the \f$ k_v = n/D \f$ component of spatial frequency.
   floatT n;

   ///p-index of the primed mode.  +/-1 for modified modes.  If basic, then +1==>cosine, -1==>sine.    
   int pp;
   
   ///The m-index of the primed mode, corresponding to the \f$ k_u = m/D \f$ component of spatial frequency.
   floatT mp;
   
   ///The n-indexof the primed mode, corresponding to the \f$ k_v = n/D \f$ component of spatial frequency.
   floatT np;
   
   ///Spatiall frequency being calculated, passed for use in the integrand worker functions.
   floatT k;

   ///Absolute tolerance for the radial integral.  Default is 1e-7.
   floatT absTol;
   
   ///Relative tolerance for the radial integral. Default is 1e-7.
   floatT relTol;
   
   ///Working memory for the azimuthal integral.
   gsl_integration_workspace * phi_w;
   
   ///Working memory for the radial integral.
   gsl_integration_workspace * k_w;

   
   ///Constructor
   fourierCovariance()
   {
      useBasic = false;
   
      absTol = 1e-7;
      relTol = 1e-7;
      
      phi_w = gsl_integration_workspace_alloc (WSZ);
      k_w = gsl_integration_workspace_alloc (WSZ);
   }
   
   ///Destructor
   ~fourierCovariance()
   {
      gsl_integration_workspace_free(phi_w);
      gsl_integration_workspace_free(k_w);
   }

   ///Calculate the covariance between the two modes.
   /** \todo document me
     * \todo handle gsl errors
     */ 
   floatT getVariance(floatT & error)
   {
      floatT result;
  
      gsl_function func;
      func.function = &kInt<floatT, aosysT>;
      func.params = this;
   
      gsl_set_error_handler_off();
      
      int ec = gsl_integration_qagiu (&func, 0, absTol, relTol, WSZ, k_w, &result, &error);
      
      return result;
   }
   
};


template<typename floatT>
int fourierVarVec( const std::string & fname,
                   int N,
                   floatT D,
                   floatT L_0,
                   bool subPist,
                   bool subTilt,
                   floatT absTol,
                   floatT relTol,
                   bool modifed=true)
{
   aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > aosys;
   
   aosys.load_MagAO_model();

   //This is just a normalization parameter in this context.
   aosys.atm.r_0(1.0, 0.5e-6);
   aosys.lam_wfs(0.5e-6); //So that output just needs to be divided by r0^(5/3)

   aosys.D( D );
   aosys.atm.L_0( L_0 );
   aosys.psd.subPiston( subPist );
   aosys.psd.subTipTilt( subTilt );

   
   std::vector<floatT> var(N,0);
   
   mx::ompLoopWatcher<> watcher(N, std::cout);
   
   #pragma omp parallel
   {   
      fourierCovariance<floatT, aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > > Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      floatT result, error;
   
      #pragma omp for schedule(static,5)
      for(int i=0; i< N; ++i)
      { 
         Pp.p = +1;
         Pp.m = i+1;
         Pp.n = 0;
         
         Pp.pp = +1;
         Pp.mp = i+1;
         Pp.np = 0;
            
         result = Pp.getVariance(error);
   
         var[i] = result;
         
         watcher.incrementAndOutputStatus();
      }
   } 

   std::ofstream fout;
   fout.open(fname);
   
   for(int i=0; i<N; ++i)
   {
      fout << i+1 << " " << var[i] << "\n";
   }
   
   fout.close();
   
}

template<typename floatT>
int fourierVarMap( const std::string & fname,
                   int N,
                   floatT D,
                   floatT L_0,
                   bool subPist,
                   bool subTilt,
                   floatT absTol,
                   floatT relTol,
                   bool modifed=true)
{
   std::vector<mx::fourierModeDef> ml;
   mx::makeFourierModeFreqs_Rect(ml, N);
   
      
   aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > aosys;
   
   aosys.load_MagAO_model();

   //This is just a normalization parameter in this context.
   aosys.atm.r_0(1.0, 0.5e-6);
   aosys.lam_wfs(0.5e-6); //So that output just needs to be divided by r0^(5/3)

   aosys.D( D );
   aosys.atm.L_0( L_0 );
   aosys.psd.subPiston( subPist );
   aosys.psd.subTipTilt( subTilt );

   int psz = ml.size();
   
   Eigen::Array<floatT, -1,-1> var;
   var.resize(N+1, N+1);
   var.setZero();
   
   mx::ompLoopWatcher<> watcher(0.5*psz, std::cout);
   
   #pragma omp parallel
   {   
      fourierCovariance<floatT, aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > > Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      floatT result, error;
   
      #pragma omp for schedule(static,5)
      for(int i=0; i< psz; ++i)
      { 
         if(ml[i].p == -1) continue; 
            
         Pp.p = ml[i].p;
         Pp.m = ml[i].m;
         Pp.n = ml[i].n;
         
         Pp.pp = ml[i].p;
         Pp.mp = ml[i].m;
         Pp.np = ml[i].n;
            
         result = Pp.getVariance(error);
   
         var( 0.5*var.rows() + ml[i].m, 0.5*var.cols() + ml[i].n) = result;
         var( 0.5*var.rows() - ml[i].m, 0.5*var.cols() - ml[i].n) = result;
      
         watcher.incrementAndOutputStatus();
      }
   } 
   
   improc::fitsHeader head;
   head.append("DIAMETER", aosys.D(), "Diameter in meters");
   head.append("L0", aosys.atm.L_0(), "Outer scale (L_0) in meters");
   head.append("SUBPIST", aosys.psd.subPiston(), "Piston subtractioon true/false flag");
   head.append("SUBTILT", aosys.psd.subTipTilt(), "Tip/Tilt subtractioon true/false flag");
   head.append("ABSTOL", absTol, "Absolute tolerance in qagiu");
   head.append("RELTOL", relTol, "Relative tolerance in qagiu");
   
   fitsHeaderGitStatus(head, "mxlib_comp",  mxlib_compiled_git_sha1(), mxlib_compiled_git_repo_modified());
   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   fitsHeaderGitStatus(head, "mxaoanalytic",  MXAOANALYTIC_CURRENT_SHA1, MXAOANALYTIC_REPO_MODIFIED);
   
   improc::fitsFile<floatT> ff;
   ff.write(fname + ".fits", var, head);
   
}

template<typename floatT>
int fourierCovarMap( const std::string & fname,
                     int N,
                     floatT D,
                     floatT L_0,
                     bool subPist,
                     bool subTilt,
                     floatT absTol,
                     floatT relTol,
                     bool modifed=true)
{
   std::vector<mx::fourierModeDef> ml;
   mx::makeFourierModeFreqs_Rect(ml, N);
   
      
   aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > aosys;
   
   aosys.load_MagAO_model();

   //This is just a normalization parameter in this context.
   aosys.atm.r_0(1.0, 0.5e-6);

   aosys.D( D );
   aosys.atm.L_0( L_0 );
   aosys.psd.subPiston( subPist );
   aosys.psd.subTipTilt( subTilt );

   int psz = ml.size();
   
   Eigen::Array<floatT,-1,-1> covar( psz, psz);
   covar.setZero();
   
   //int ncalc = 0.5*( psz*psz - psz);
   
   mx::ompLoopWatcher<> watcher(psz, std::cout);
   
   #pragma omp parallel
   {   
      fourierCovariance<floatT, aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > > Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      floatT result, error;
   
      #pragma omp for schedule(static,5)
      for(int i=0; i< psz; ++i)
      { 
         for(int j=i; j< psz; ++j)
         {
            Pp.p = ml[i].p;
            Pp.m = ml[i].m;
            Pp.n = ml[i].n;
         
            Pp.pp = ml[j].p;
            Pp.mp = ml[j].m;
            Pp.np = ml[j].n;
            result = Pp.getVariance(error);
   
            covar(i,j) = result;
            
         }
         watcher.incrementAndOutputStatus();
      }
   } 
   
   improc::fitsHeader head;
   head.append("DIAMETER", aosys.D(), "Diameter in meters");
   head.append("L0", aosys.atm.L_0(), "Outer scale (L_0) in meters");
   head.append("SUBPIST", aosys.psd.subPiston(), "Piston subtractioon true/false flag");
   head.append("SUBTILT", aosys.psd.subTipTilt(), "Tip/Tilt subtractioon true/false flag");
   head.append("ABSTOL", absTol, "Absolute tolerance in qagiu");
   head.append("RELTOL", relTol, "Relative tolerance in qagiu");
   
   fitsHeaderGitStatus(head, "mxlib_comp",  mxlib_compiled_git_sha1(), mxlib_compiled_git_repo_modified());
   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   fitsHeaderGitStatus(head, "mxaoanalytic",  MXAOANALYTIC_CURRENT_SHA1, MXAOANALYTIC_REPO_MODIFIED);
   
   improc::fitsFile<floatT> ff;
   ff.write(fname + ".fits", covar, head);
   
}




template<typename floatT>
int fourierCovarMapSeparated( const std::string & fname,
                              int N,
                              floatT D,
                              floatT L_0,
                              bool subPist,
                              bool subTilt,
                              floatT absTol,
                              floatT relTol,
                              bool modifed=true)
{
   std::vector<mx::fourierModeDef> ml;
   mx::makeFourierModeFreqs_Rect(ml, N);
      
   aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > aosys;
   
   aosys.load_MagAO_model();

   //This is just a normalization parameter in this context.
   aosys.atm.r_0(1.0, 0.5e-6);

   aosys.D( D );
   aosys.atm.L_0( L_0 );
   aosys.psd.subPiston( subPist );
   aosys.psd.subTipTilt( subTilt );

   int psz = 0.5*ml.size();
   
   Eigen::Array<floatT,-1,-1> covar_pp( psz, (int)(.5*psz)), covar_ppp(psz, (int)(0.5*psz));
   covar_pp.setZero();
   covar_ppp.setZero();
   
   #pragma omp parallel
   {   

      fourierCovariance<floatT, aoSystem<floatT, vonKarmanSpectrum<floatT>, pywfsUnmod<floatT> > > Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      floatT result, error;
   
      #pragma omp for
      for(int i=0; i< psz; ++i)
      {  
         for(int j=0; j<= 0.5*i; ++j)
         {
            for(int k=0; k< 2; ++k)
            {
               Pp.p = ml[i].p;
               Pp.m = ml[i].m;
               Pp.n = ml[i].n;
         
               Pp.pp = ml[2*j + k].p;
               Pp.mp = ml[2*j + k].m;
               Pp.np = ml[2*j + k].n;
               result = Pp.getVariance(error);
   
               if( Pp.p == Pp.pp)
               {
                  covar_pp(i, j) = result;
               }
               else
               {
                  covar_ppp(i, j) = result;
               }
            }
         }
      }
   } 
   
   improc::fitsHeader head;
   head.append("DIAMETER", aosys.D(), "Diameter in meters");
   head.append("L0", aosys.atm.L_0(), "Outer scale (L_0) in meters");
   head.append("SUBPIST", aosys.psd.subPiston(), "Piston subtractioon true/false flag");
   head.append("SUBTILT", aosys.psd.subTipTilt(), "Tip/Tilt subtractioon true/false flag");
   head.append("ABSTOL", absTol, "Absolute tolerance in qagiu");
   head.append("RELTOL", relTol, "Relative tolerance in qagiu");
   
   fitsHeaderGitStatus(head, "mxlib_comp",  mxlib_compiled_git_sha1(), mxlib_compiled_git_repo_modified());
   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   fitsHeaderGitStatus(head, "mxaoanalytic",  MXAOANALYTIC_CURRENT_SHA1, MXAOANALYTIC_REPO_MODIFIED);
   
   improc::fitsFile<floatT> ff;
   ff.write(fname + "_pp.fits", covar_pp, head);
   ff.write(fname + "_ppp.fits", covar_ppp, head);
   
}

template<typename realT>
void calcKLCoeffs( const std::string & outFile,
                   const std::string & cvFile )
{
   improc::fitsFile<realT> ff;
   
   Eigen::Array<realT,-1,-1> cvT, cv, evecs, evals;
   
   ff.read(cv, cvFile);

   //cvT = cv.block(0,0, 1000,1000);//.transpose();
   
   std::cerr << cvT.rows() << " " << cvT.cols() << "\n";
   std::cerr << "1\n";
   syevrMem<int, int, double> mem;
   
   double t0 = get_curr_time<double>();
   int info = eigenSYEVR<double,double>(evecs, evals, 0, cv, 0, -1, 'U', &mem);
   double t1 = get_curr_time<double>();

   std::cerr << "2\n";
   
   if(info !=0 ) 
   {
      std::cerr << "info =" << info << "\n";
      exit(0);
   }

   std::cerr << "Time = " << t1-t0 << " secs\n";
   
   //Normalize the eigenvectors
   for(int i=0;i< evecs.cols(); ++i)
   {
      evecs.col(i) = evecs.col(i)/sqrt(fabs(evals(i)));
   }
   
   ff.write(outFile, evecs);
}

template<typename eigenArrT1, typename eigenArrT2, typename eigenArrT3>
void makeKL( eigenArrT1 & kl,
             eigenArrT2 & evecs,
             eigenArrT3 && rvecs )
{
   
   int tNims = evecs.rows();
   int tNpix = rvecs.rows();

   int n_modes = tNims;
   
   //Now calculate KL images
   /*
    *  KL = E^T * R  ==> C = A^T * B
    */
   math::gemm<typename eigenArrT1::Scalar>(CblasColMajor, CblasTrans, CblasTrans, n_modes, tNpix,
                              tNims, 1., evecs.data(), evecs.rows(), rvecs.data(), rvecs.rows(),
                                 0., kl.data(), kl.rows());
   
}

template<typename realT>
void makeFKL( const std::string & outFile,
              const std::string & coeffs,
              int N,
              int pupSize )
{
   improc::fitsFile<realT> ff;
   Eigen::Array<realT, -1, -1> evecs;
   
   ff.read(evecs, coeffs);
   
   improc::eigenCube<realT> Rims;
   makeFourierBasis_Rect(Rims, pupSize, N, MX_FOURIER_MODIFIED);
   
   
   std::cout << Rims.planes() << " " << evecs.cols() << "\n";
   
   
   
   
   Eigen::Array<realT,-1,-1> kl;
   kl.resize( Rims.planes(), Rims.rows()*Rims.cols());
   
   std::cerr  << 1 << "\n";
   makeKL( kl, evecs, Rims.cube());
   std::cerr  << 2 << "\n";
   
   Eigen::Array<realT,-1,-1> klT = kl.transpose();
   //kl.resize(0,0);
   //Rims.resize(0,0);
   //evecs.resize(0,0);
      
   improc::eigenCube<realT> klims(klT.data(), Rims.rows(), Rims.cols(), Rims.planes());
   
   improc::eigenCube<realT> klimsR;
   klimsR.resize( klims.rows(), klims.cols(), klims.planes());
   
   for(int i=0; i< klims.planes(); ++i)
   {
      klimsR.image(i) = klims.image(klims.planes()-1-i);
   }
   
   ff.write(outFile, klimsR);
   std::cerr  << 3 << "\n";
}

} //namespace AO
} //namespace mx

#endif //__fourierCovariance_hpp__

