/** \file zernikeCovariance.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculation of the modal covariance in the zernike basis.
  * \ingroup mxAO_files
  *
  */

#ifndef zernikeCovariance_hpp
#define zernikeCovariance_hpp


#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"
#include "../../sigproc/zernike.hpp"
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

/// Structure to manage the zernike mode covariance calculation, passed to integration functions
/**
  * \tparam realT a floating point type used for all calculations.  As of Dec 2023 must be double due to gsl_integration.
  * \tparam aosysT the type of the AO system structure
  */
template<typename realT, typename aosysT>
struct zernikeCovariance
{
   ///Pointer to an AO system, which contains the relevant spatial PSD of turbulence.
   aosysT * m_aosys {nullptr};

   ///The n-index of the unprimed mode.
   realT m_n;

   ///The m-index of the unprimed mode.
   realT m_m;

   ///The n-indexof the primed mode, corresponding to the \f$ k_v = n/D \f$ component of spatial frequency.
   realT m_np;

   ///The m-index of the primed mode, corresponding to the \f$ k_u = m/D \f$ component of spatial frequency.
   realT m_mp;

   ///Spatial frequency being calculated, passed for use in the integrand worker functions.
   realT m_k;

   ///Absolute tolerance for the radial integral.  Default is 1e-10. 
   realT m_kIntEpsAbs {1e-10};

   ///Relative tolerance for the radial integral. Default is 0, meaning absolute is used.
   realT m_kIntEpsRel {0};

    ///Absolute tolerance for the azimuthal integral.  Default is 1e-10.
   realT m_phiIntEpsAbs {1e-10};

   ///Relative tolerance for the azimuthal integral. Default is 0, meaning absolute is used.
   realT m_phiIntEpsRel {0};

protected:
   size_t m_workspaceSize {1000000};

   ///Working memory for the azimuthal integral.
   gsl_integration_workspace * m_workspacePhi {nullptr};

   ///Working memory for the radial integral.
   gsl_integration_workspace * m_workspaceK {nullptr};

public:

     ///Constructor
     zernikeCovariance()
     {
         m_workspacePhi = gsl_integration_workspace_alloc (m_workspaceSize);
         m_workspaceK = gsl_integration_workspace_alloc (m_workspaceSize);
     }
 
     ///Destructor
     ~zernikeCovariance()
     {
         gsl_integration_workspace_free(m_workspacePhi);
         gsl_integration_workspace_free(m_workspaceK);
     }
 
     void workspaceSize( size_t wsz )
     {
         gsl_integration_workspace_free(m_workspacePhi);
         gsl_integration_workspace_free(m_workspaceK);
 
         m_workspacePhi = gsl_integration_workspace_alloc (m_workspaceSize);
         m_workspaceK = gsl_integration_workspace_alloc (m_workspaceSize);
     }
 
     size_t workspaceSize()
     {
         return m_workspaceSize;
     }
 
     gsl_integration_workspace * workspacePhi()
     {
         return m_workspacePhi;
     }

     ///Calculate the covariance between the two modes.
     /** \todo document me
       * \todo handle gsl errors
       */
     realT getCovariance(realT & error)
     {
         realT result;

         if(m_aosys == nullptr)
         {
             mxThrowException(err::paramnotset, "zernikeCovariance::getCovariance", "AO system not setup (aosys is nullptr)");
         }

       gsl_function func;
       func.function = &kInt;
       func.params = this;

       gsl_set_error_handler_off();

       int ec = gsl_integration_qagiu (&func, 0, m_kIntEpsAbs, m_kIntEpsRel, m_workspaceSize, m_workspaceK, &result, &error);

       return result;
    }

    ///Worker function for the radial integral in the covariance calculation
    /**
      * \returns the value of the integrand at spatial frequency \p k.
      */
    static realT kInt( realT k,      ///< [in] the spatial frequency at which to evaluate the integrand, meters
                       void * params ///< [in] a pointer to a object of type zernikeCovariance<realT, aosyT>
                     );

    /// Worker for the azimuthal integral (in phi) for the zernike mode covariance.
    /**
      * \return the value of the integrand at the angle \p phi
      */
    static realT phiInt( realT phi,    ///< [in] the angle at which to evaluate the integrand, radians
                         void * params ///< [in] a pointer to a object of type zernikeCovariance<realT, aosyT>
                       );
};


template<typename realT, typename aosysT>
realT zernikeCovariance<realT, aosysT>::kInt( realT k, 
                                              void * params
                                            )
{
   zernikeCovariance<realT, aosysT> * Pp = (zernikeCovariance<realT, aosysT> *) params;

   realT result, error;

   gsl_function func;

   func.function = &phiInt;
   
   func.params = Pp;

   Pp->m_k = k;

   gsl_integration_qag(&func, 0., math::two_pi<double>(), Pp->m_phiIntEpsAbs, Pp->m_phiIntEpsRel, Pp->workspaceSize(), GSL_INTEG_GAUSS61, Pp->workspacePhi(), &result, &error);


   return result;
}

template<typename realT, typename aosysT>
realT zernikeCovariance<realT, aosysT>::phiInt( realT phi, 
                                                void * params
                                              )
{
    zernikeCovariance<realT, aosysT> * Pp = (zernikeCovariance<realT, aosysT> *) params;

    realT D = Pp->m_aosys->D();

    realT k = Pp->m_k;
   
    /*** no prime ***/
    realT m = Pp->m_m;
    realT n = Pp->m_n;

    std::complex<realT> Q_mn = zernikeQ(k*D/2.0, phi, n, m);

    /*** primed ***/
    realT mp = Pp->m_mp;
    realT np = Pp->m_np;

    std::complex<realT> Q_mpnp = zernikeQ(k*D/2.0, phi, np, mp);

    /*** The product ***/
    realT QQ = real(conj(Q_mn)*Q_mpnp);

    realT P = Pp->m_aosys->psd(Pp->m_aosys->atm, 0, k, Pp->m_aosys->lam_sci(), Pp->m_aosys->lam_wfs(), 1.0);

    return P*k*QQ;
}


#if 0

///Calculate a vector of zernike mode variances.
template<typename realT, typename aosysT>
int zernikeVarVec( const std::string & fname,
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
   if( aosys.d_min(0) > 0)
   {
      mnCon = floor( aosys.D()/aosys.d_min(0)/2.0);
   }

   #pragma omp parallel
   {
      zernikeCovariance<realT, aosysT> Pp;
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

   realT D = aosys.D();
   realT r_0 = aosys.atm.r_0();
   realT tot = 1.0299*pow(D/r_0, 5./3.);

   realT sum = 0;

   std::cout << 0 << " " << 0 << " " << 0 << " " <<  tot << "\n";
   for(int i=0; i<N; ++i)
   {
      realT k = (i+1)/D;

      realT P = aosys.psd(aosys.atm, 0, k, 1.0);// aosys.lam_sci(), 0, aosys.lam_wfs(), 1.0);

      if(mnCon > 0 )
      {
         if( k*D < mnCon )
         {
            P *= pow(math::two_pi<realT>()*aosys.atm.v_wind()* k * (aosys.minTauWFS(0)+aosys.deltaTau()),2);
         }
      }

      sum += var[i];
      fout << i+1 << " " << var[i] << " " << P/pow(D,2) << " " << tot-sum << "\n";
   }

   fout.close();

   return 0;
}

///Calculate a map of zernike variances by convolution with the PSD
/** Uses the Airy pattern for the circularly unobstructed aperture.
  *
  * \returns 0 on success
  * \returns -1 on error
  */
template<typename realT, typename aosysT>
int zernikePSDMap( improc::eigenImage<realT> & var, ///< [out] The variance estimated by convolution with the PSD
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
int zernikeCovarMap( const std::string & fname, ///< [out] the path where the output FITS file will be written
                     int N, ///< [in] the linear number of zernike modes across the aperture.  The Nyquist frequency is set by N/2.
                     realT D,
                     realT L_0,
                     bool subPist,
                     bool subTilt,
                     realT absTol,
                     realT relTol,
                     bool modified=true
                   )
{
   std::vector<mx::sigproc::zernikeModeDef> ml;
   mx::sigproc::makezernikeModeFreqs_Rect(ml, N);


   aoSystem<realT, vonKarmanSpectrum<realT>> aosys;

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
      zernikeCovariance<realT, aoSystem<realT, vonKarmanSpectrum<realT> > > Pp;
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
int zernikeCovarMapSeparated( const std::string & fname,
                              int N,
                              aosysT & aosys,
                              realT absTol,
                              realT relTol,
                              bool modified=true)
{
   std::vector<mx::sigproc::zernikeModeDef> ml;
   mx::sigproc::makezernikeModeFreqs_Rect(ml, N);

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

      zernikeCovariance<realT, aosysT > Pp;
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
   sigproc::makezernikeBasis_Rect(Rims, pupSize, N, MX_zernike_MODIFIED);

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

#endif

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //zernikeCovariance_hpp
