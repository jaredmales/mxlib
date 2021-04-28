/** \file fourierCovariance.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculation of the modal covariance in the Fourier basis.
  * \ingroup mxAO_files
  *
  */

#ifndef fourierCovariance_hpp
#define fourierCovariance_hpp


#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"
#include "../../sigproc/fourierModes.hpp"
#include "../../ioutils/fits/fitsFile.hpp"
#include "../../improc/eigenImage.hpp"
#include "../../improc/eigenCube.hpp"
#include "../../mxlib_uncomp_version.h"
#include "../../ipc/ompLoopWatcher.hpp"
#include "../../sys/timeUtils.hpp"
#include "../../math/eigenLapack.hpp"

#include "../../math/func/airyPattern.hpp"



#include "aoAtmosphere.hpp"
#include "aoPSDs.hpp"
#include "aoSystem.hpp"
#include "varmapToImage.hpp"


namespace mx
{
namespace AO
{
namespace analysis
{

#ifndef WSZ

/** \def WSZ
  * \brief The size of the gsl_integration workspaces.
  */
#define WSZ 100000

#endif


//Forward Declaration
template<typename realT, typename aosysT>
struct fourierCovariance;

///Worker for the azimuthal integral (in phi) for the basic Fourier mode covariance.
/**
  * \param phi the angle at which to evaluate the integrand
  * \param params a pointer to a object of type fourierCovariance<realT, aosyT>
  *
  * \tparam realT a floating point type used for all calculations
  * \tparam aosysT the type of the AO system structure
  *
  * \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT phiInt_basic (realT phi, void * params)
{

   fourierCovariance<realT, aosysT> * Pp = (fourierCovariance<realT, aosysT> *) params;

   realT D = Pp->aosys->D();

   int p = Pp->p;
   realT m = Pp->m;
   realT n = Pp->n;

   int pp = Pp->pp;
   realT mp = Pp->mp;
   realT np = Pp->np;

   realT k = Pp->k;

   realT cosp = cos(phi);
   realT sinp = sin(phi);

   /*** no prime ***/
   realT kmn_p, kmn_m;
   realT Ji_mn_p, Ji_mn_m;
   realT Q_mn;//, Qs_mn;

   kmn_p = sqrt( pow(k*cosp + m/D, 2) + pow(k*sinp + n/D, 2));
   kmn_m = sqrt( pow(k*cosp - m/D, 2) + pow(k*sinp - n/D, 2));

   Ji_mn_p = math::func::jinc(math::pi<realT>()*D*kmn_p);
   Ji_mn_m = math::func::jinc(math::pi<realT>()*D*kmn_m);


   realT N = 1./sqrt(0.5 + p*math::func::jinc(math::two_pi<realT>()*sqrt(m*m+n*n)));

   Q_mn = N*(Ji_mn_p + p*Ji_mn_m);

   /*** primed ***/
   realT kmpnp_p, kmpnp_m;
   realT Ji_mpnp_p, Ji_mpnp_m;
   realT Q_mpnp;//, Qs_mpnp;

   kmpnp_p = sqrt( pow(k*cosp + mp/D, 2) + pow(k*sinp + np/D, 2));
   kmpnp_m = sqrt( pow(k*cosp - mp/D, 2) + pow(k*sinp - np/D, 2));


   Ji_mpnp_p = math::func::jinc(math::pi<realT>()*D*kmpnp_p);
   Ji_mpnp_m = math::func::jinc(math::pi<realT>()*D*kmpnp_m);

   realT Np = 1./sqrt(0.5 + pp*math::func::jinc(math::two_pi<realT>()*sqrt(mp*mp+np*np)));

   Q_mpnp = Np*(Ji_mpnp_p + pp*Ji_mpnp_m);

   realT P = Pp->aosys->psd(Pp->aosys->atm, k, 1.0); //vonKarmanPSD(k, D, L0, Pp->subPiston, Pp->subTipTilt);

   return P*k*(Q_mn*Q_mpnp);
}

///Worker for the azimuthal integral (in phi) for the modified Fourier mode covariance.
/**
  * \param phi the angle at which to evaluate the integrand
  * \param params a pointer to a object of type fourierCovariance<realT, aosyT>
  *
  * \tparam realT a floating point type used for all calculations
  * \tparam aosysT the type of the AO system structure
  */
template<typename realT, typename aosysT>
realT phiInt_mod (realT phi, void * params)
{

   fourierCovariance<realT, aosysT> * Pp = (fourierCovariance<realT, aosysT> *) params;

   //realT L0 = Pp->aosys->atm.L_0();
   realT D = Pp->aosys->D();

   int p = Pp->p;
   realT m = Pp->m;
   realT n = Pp->n;

   int pp = Pp->pp;
   realT mp = Pp->mp;
   realT np = Pp->np;

   realT k = Pp->k;

   realT mnCon = Pp->mnCon;

   realT cosp = cos(phi);
   realT sinp = sin(phi);

   /*** no prime ***/
   realT kmn_p, kmn_m;
   realT Ji_mn_p, Ji_mn_m;
   
   kmn_p = sqrt( pow(k*cosp + m/D, 2) + pow(k*sinp + n/D, 2));
   kmn_m = sqrt( pow(k*cosp - m/D, 2) + pow(k*sinp - n/D, 2));

   Ji_mn_p = math::func::jinc(math::pi<realT>()*D*kmn_p);
   Ji_mn_m = math::func::jinc(math::pi<realT>()*D*kmn_m);


   /*** primed ***/
   realT kmpnp_p, kmpnp_m;
   realT Ji_mpnp_p, Ji_mpnp_m;
   
   kmpnp_p = sqrt( pow(k*cosp + mp/D, 2) + pow(k*sinp + np/D, 2));
   kmpnp_m = sqrt( pow(k*cosp - mp/D, 2) + pow(k*sinp - np/D, 2));


   Ji_mpnp_p = math::func::jinc(math::pi<realT>()*D*kmpnp_p);
   Ji_mpnp_m = math::func::jinc(math::pi<realT>()*D*kmpnp_m);

   realT QQ;

   if(p == pp)
   {
      QQ = 2.0*( Ji_mn_p*Ji_mpnp_p + Ji_mn_m*Ji_mpnp_m);
   }
   else
   {
      QQ = 2.0*( Ji_mn_p*Ji_mpnp_m + Ji_mn_m*Ji_mpnp_p);
   }

   realT P = Pp->aosys->psd(Pp->aosys->atm, k, 1.0); //psd(Pp->aosys->atm, k, Pp->aosys->lam_sci(), Pp->aosys->lam_wfs(), 1.0);

   if(mnCon > 0 )
   {
      if( k*D <= mnCon )
      {
         P *= pow(math::two_pi<realT>()*Pp->aosys->atm.v_wind()* k * (Pp->aosys->minTauWFS()+Pp->aosys->deltaTau()),2);
      }
   }


   return P*k*QQ;
}



///Worker function for the radial integral in the covariance calculation
/**
  * \param k the spatial frequency at which to evaluate the integrand
  * \param params a pointer to a object of type fourierCovariance<realT, aosyT>
  *
  * \tparam realT a floating point type used for all calculations.  As of Nov 2016 must be double due to gsl_integration.
  * \tparam aosysT the type of the AO system structure
  */
template<typename realT, typename aosysT>
realT kInt (realT k, void * params)
{
   fourierCovariance<realT, aosysT> * Pp = (fourierCovariance<realT, aosysT> *) params;

   realT result, error;

   gsl_function func;

   //Here choose between basic and modified Fourier modes.
   if(Pp->useBasic)
   {
      func.function = &phiInt_basic<realT, aosysT>;
   }
   else
   {
      func.function = &phiInt_mod<realT, aosysT>;
   }
   func.params = Pp;

   Pp->k = k;

   //Tolerances of 1e-4 seem to be all we can ask for
   gsl_integration_qag(&func, 0., math::two_pi<double>(), 1e-10, 1e-4, WSZ, GSL_INTEG_GAUSS21, Pp->phi_w, &result, &error);


   return result;
}

/// Structure to manage the Fourier mode covariance calculation, passed to integration functions
/**
  * \tparam realT a floating point type used for all calculations.  As of Nov 2016 must be double due to gsl_integration.
  * \tparam aosysT the type of the AO system structure
  */
template<typename realT, typename aosysT>
struct fourierCovariance
{
   ///Pointer to an AO system, which contains the relevant spatial PSD of turbulence.
   aosysT * aosys {nullptr};

   ///Flag controlling use of basic or modified Fourier modes.  If true, the basic sin/cos modes are used.  If false (default), the modified modes are u sed.
   bool useBasic {false};

   ///p-index of the unprimed mode.  +/-1 for modified modes.  If basic, then +1==>cosine, -1==>sine.
   int p;

   ///The m-index of the unprimed mode, corresponding to the \f$ k_u = m/D \f$ component of spatial frequency.
   realT m;

   ///The n-indexof the unprimed mode, corresponding to the \f$ k_v = n/D \f$ component of spatial frequency.
   realT n;

   ///p-index of the primed mode.  +/-1 for modified modes.  If basic, then +1==>cosine, -1==>sine.
   int pp;

   ///The m-index of the primed mode, corresponding to the \f$ k_u = m/D \f$ component of spatial frequency.
   realT mp;

   ///The n-indexof the primed mode, corresponding to the \f$ k_v = n/D \f$ component of spatial frequency.
   realT np;

   ///Spatial frequency being calculated, passed for use in the integrand worker functions.
   realT k;

   ///The maximum controlled value of spatial frequency (k*D).  If < 0 then not controlled.
   realT mnCon {0};

   ///Absolute tolerance for the radial integral.  Default is 1e-7.
   realT absTol {1e-7};

   ///Relative tolerance for the radial integral. Default is 1e-7.
   realT relTol {1e-7};

   ///Working memory for the azimuthal integral.
   gsl_integration_workspace * phi_w;

   ///Working memory for the radial integral.
   gsl_integration_workspace * k_w;

   ///Constructor
   fourierCovariance()
   {
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
   realT getVariance(realT & error)
   {
      realT result;

      gsl_function func;
      func.function = &kInt<realT, aosysT>;
      func.params = this;

      gsl_set_error_handler_off();

      int ec = gsl_integration_qagiu (&func, 0, absTol, relTol, WSZ, k_w, &result, &error);

      //Finally, convert to square sampling.
      return result;
   }

};

///Calculate a vector of Fourier mode variances.
template<typename realT, typename aosysT>
int fourierVarVec( const std::string & fname,
                   int N,
                   aosysT & aosys,
                   realT absTol,
                   realT relTol,
                   bool modifed=true
                 )
{

   std::vector<realT> var(N,0);

   ipc::ompLoopWatcher<> watcher(N, std::cerr);

   realT mnCon = 0;
   if( aosys.d_min() > 0)
   {
      mnCon = floor( aosys.D()/aosys.d_min()/2.0);
   }

   #pragma omp parallel
   {
      fourierCovariance<realT, aosysT> Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      Pp.mnCon = mnCon;

      realT result, error;

      #pragma omp for schedule(dynamic,5)
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
      realT D = aosys.D();
      realT k = (i+1)/D;
      realT P = aosys.psd(aosys.atm, k, aosys.lam_sci(), 0, aosys.lam_wfs(), 1.0);

      if(mnCon > 0 )
      {
         if( k*D < mnCon )
         {
            P *= pow(math::two_pi<realT>()*aosys.atm.v_wind()* k * (aosys.minTauWFS()+aosys.deltaTau()),2);
         }
      }


      fout << i+1 << " " << var[i] << " " << P/pow(D,2) << "\n";
   }

   fout.close();

   return 0;
}

///Calculate a map of Fourier variances by convolution with the PSD
/** Uses the Airy pattern for the circularly unobstructed aperture.
  *
  * \returns 0 on success
  * \returns -1 on error
  */
template<typename realT, typename aosysT>
int fourierPSDMap( improc::eigenImage<realT> & var, ///< [out] The variance estimated by convolution with the PSD
                   improc::eigenImage<realT> & psd, ///< [out] the PSD map
                   int N, ///< [in] the number of components to analyze
                   int overSample,
                   aosysT & aosys ///< [in[ the AO system defining the PSD characteristics.
                 )
{
   N *= overSample;
   psd.resize(2*N + 1, 2*N+1);

   realT mnCon = 0;
   if( aosys.d_min() > 0)
   {
      mnCon = floor( aosys.D()/aosys.d_min()/2.0);
   }

   for(int i=0; i<=N; ++i)
   {
      for(int j=-N; j<=N; ++j)
      {

         realT D = aosys.D();
         realT k = sqrt( pow(i,2) + pow(j,2))/D/overSample;

         realT P = aosys.psd(aosys.atm, k, aosys.lam_sci(), 0, aosys.lam_wfs(), 1.0);

         if(mnCon > 0 )
         {
            if( k*D < mnCon )
            {
               P *= pow(math::two_pi<realT>()*aosys.atm.v_wind()* k * (aosys.minTauWFS()+aosys.deltaTau()),2);
            }
         }

         psd(N+i, N + j) = P/pow(D*overSample,2);
         psd(N-i, N - j) = P/pow(D*overSample,2);
      }
   }

   //Create Airy PSF for convolution with PSD psd.
   Eigen::Array<realT, -1,-1> psf;
   psf.resize(2*N+3,2*N+3);
   for(int i=0;i<psf.rows();++i)
   {
      for(int j=0;j<psf.cols();++j)
      {
         psf(i,j) = mx::math::func::airyPattern(sqrt( pow( i-floor(.5*psf.rows()),2) + pow(j-floor(.5*psf.cols()),2))/overSample);
      }
   }

   mx::AO::analysis::varmapToImage(var, psd, psf);

   return 0;
}

template<typename realT>
int fourierCovarMap( const std::string & fname, ///< [out] the path where the output FITS file will be written
                     int N, ///< [in] the linear number of Fourier modes across the aperture.  The Nyquist frequency is set by N/2.
                     realT D,
                     realT L_0,
                     bool subPist,
                     bool subTilt,
                     realT absTol,
                     realT relTol,
                     bool modified=true
                   )
{
   std::vector<mx::sigproc::fourierModeDef> ml;
   mx::sigproc::makeFourierModeFreqs_Rect(ml, N);


   aoSystem<realT, vonKarmanSpectrum<realT>, pywfsUnmod<realT> > aosys;

   aosys.loadMagAOX();

   //This is just a normalization parameter in this context.
   aosys.atm.r_0(1.0, 0.5e-6);

   aosys.D( D );
   aosys.atm.L_0( L_0 );
   aosys.psd.subPiston( subPist );
   aosys.psd.subTipTilt( subTilt );

   int psz = ml.size();

   Eigen::Array<realT,-1,-1> covar( psz, psz);
   covar.setZero();

   //int ncalc = 0.5*( psz*psz - psz);

   ipc::ompLoopWatcher<> watcher(psz, std::cout);

   std::cerr << "Starting . . .\n";
   #pragma omp parallel
   {
      fourierCovariance<realT, aoSystem<realT, vonKarmanSpectrum<realT>, pywfsUnmod<realT> > > Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      if(!modified) Pp.useBasic = true;

      realT result, error;

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

   fits::fitsHeader head;
   head.append("DIAMETER", aosys.D(), "Diameter in meters");
   head.append("NSUBAP", N, "Linear number of s.f. sampled");
   head.append("L0", aosys.atm.L_0(0), "Outer scale (L_0) in meters");
   head.append("SUBPIST", aosys.psd.subPiston(), "Piston subtractioon true/false flag");
   head.append("SUBTILT", aosys.psd.subTipTilt(), "Tip/Tilt subtractioon true/false flag");
   head.append("ABSTOL", absTol, "Absolute tolerance in qagiu");
   head.append("RELTOL", relTol, "Relative tolerance in qagiu");

   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   

   fits::fitsFile<realT> ff;
   ff.write(fname + ".fits", covar, head);

   return 0;
}




template<typename realT, typename aosysT>
int fourierCovarMapSeparated( const std::string & fname,
                              int N,
                              aosysT & aosys,
                              realT absTol,
                              realT relTol,
                              bool modified=true)
{
   std::vector<mx::sigproc::fourierModeDef> ml;
   mx::sigproc::makeFourierModeFreqs_Rect(ml, N);

   int psz = 0.5*ml.size();

   Eigen::Array<realT,-1,-1> covar_pp( (int) (0.5*psz), (int)(.5*psz)), covar_ppp( (int) (0.5*psz), (int)(0.5*psz));
   covar_pp.setZero();
   covar_ppp.setZero();

   realT mnCon = 0;
   if( aosys.d_min() > 0)
   {
      mnCon = floor( aosys.D()/aosys.d_min()/2.0);
   }

   ipc::ompLoopWatcher<> watcher((psz+1)*0.125*(psz+1)*2, std::cout);

   #pragma omp parallel
   {

      fourierCovariance<realT, aosysT > Pp;
      Pp.absTol = absTol;
      Pp.relTol = relTol;
      Pp.aosys = &aosys;

      if(!modified) Pp.useBasic = true;

      Pp.mnCon = mnCon;

      realT result, error;

      #pragma omp for schedule(dynamic,5)
      for(int i=0; i< psz; i+=2)
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
                  covar_pp(i/2, j) = result;
               }
               else
               {
                  covar_ppp(i/2, j) = result;
               }
               watcher.incrementAndOutputStatus();
            }
         }
      }
   }

   fits::fitsHeader head;
   head.append("DIAMETER", aosys.D(), "Diameter in meters");
   head.append("L0", aosys.atm.L_0(), "Outer scale (L_0) in meters");
   head.append("SUBPIST", aosys.psd.subPiston(), "Piston subtractioon true/false flag");
   head.append("SUBTILT", aosys.psd.subTipTilt(), "Tip/Tilt subtractioon true/false flag");
   head.append("ABSTOL", absTol, "Absolute tolerance in qagiu");
   head.append("RELTOL", relTol, "Relative tolerance in qagiu");

   
   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   
   fits::fitsFile<realT> ff;
   ff.write(fname + "_pp.fits", covar_pp, head);
   ff.write(fname + "_ppp.fits", covar_ppp, head);

   return 0;
}

template<typename realT>
void calcKLCoeffs( const std::string & outFile,
                   const std::string & cvFile )
{
   fits::fitsFile<realT> ff;

   Eigen::Array<realT,-1,-1> cvT, cv, evecs, evals;

   ff.read(cv, cvFile);

   //cvT = cv.block(0,0, 1000,1000);//.transpose();

   std::cerr << cvT.rows() << " " << cvT.cols() << "\n";
   std::cerr << "1\n";
   math::syevrMem<double> mem;

   double t0 = sys::get_curr_time();

   int info = math::eigenSYEVR<double,double>(evecs, evals, cv, 0, -1, 'U', &mem);

   double t1 = sys::get_curr_time();

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
   fits::fitsFile<realT> ff;
   Eigen::Array<realT, -1, -1> evecs;

   ff.read(evecs, coeffs);

   improc::eigenCube<realT> Rims;
   sigproc::makeFourierBasis_Rect(Rims, pupSize, N, MX_FOURIER_MODIFIED);

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

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //fourierCovariance_hpp
