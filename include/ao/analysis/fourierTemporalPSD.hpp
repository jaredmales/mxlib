/** \file fourierTemporalPSD.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculation of the temporal PSD of Fourier modes.
  * \ingroup mxAOm_files
  *
  */

//***********************************************************************//
// Copyright 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef fourierTemporalPSD_hpp
#define fourierTemporalPSD_hpp

#include <iostream>
#include <fstream>

#include <sys/stat.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Eigen/Dense>

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"
#include "../../math/func/airyPattern.hpp"
#include "../../math/vectorUtils.hpp"
#include "../../ioutils/fits/fitsFile.hpp"
#include "../../sigproc/fourierModes.hpp"
#include "../../sigproc/psdVarMean.hpp"
#include "../../ioutils/stringUtils.hpp"
#include "../../ioutils/readColumns.hpp"
#include "../../ioutils/binVector.hpp"
#include "../../ioutils/fileUtils.hpp"

#include "../../ipc/ompLoopWatcher.hpp"
#include "../../mxError.hpp"
#include "../../math/gslInterpolation.hpp"

#include "aoSystem.hpp"
#include "aoPSDs.hpp"
#include "wfsNoisePSD.hpp"
#include "clAOLinearPredictor.hpp"
#include "clGainOpt.hpp"
#include "varmapToImage.hpp"
#include "speckleAmpPSD.hpp"

#include "aoConstants.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

#ifndef WSZ

/** \def WFZ
  * \brief Size of the GSL integration workspace
  */
#define WSZ 100000

#endif

enum basis : unsigned int { basic, ///< The basic sine and cosine Fourier modes
                            modified, ///< The modified Fourier basis from \cite males_guyon_2017
                            projected_basic,
                            projectedm_modified
                          };


//Forward declaration
template<typename realT, typename aosysT>
realT F_basic (realT kv, void * params) ;

//Forward declaration
template<typename realT, typename aosysT>
realT F_mod (realT kv, void * params);

//Forward declaration
template<typename realT, typename aosysT>
realT Fm_projMod (realT kv, void * params);

///Class to manage the calculation of temporal PSDs of the Fourier modes in atmospheric turbulence.
/** Works with both basic (sines/cosines) and modified Fourier modes.
  *
  * \tparam realT is a real floating point type for calculations.  Currently must be double due to gsl_integration.
  * \tparam aosysT is an AO system type, usually of type ao_system.
  *
  * \todo Split off the integration parameters in a separate structure.
  * \todo once integration parameters are in a separate structure, make this a class with protected members.
  * \todo GSL error handling
  *
  * \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
struct fourierTemporalPSD
{
   ///Pointer to an AO system structure.
   aosysT * m_aosys {nullptr};

   realT m_f {0}; ///< the current temporal frequency0
   realT m_m {0}; ///< the spatial frequency m index
   realT m_n {0}; ///< the spatial frequency n index
   realT m_cq {0}; ///< The cosine of the wind direction
   realT m_sq {0}; ///< The sine of the wind direction
   realT m_spatialFilter  {false}; ///< Flag indicating if a spatial filter is applied
   
   int m_p {1}; ///< The parity of the mode, +/- 1.  If _useBasis==MXAO_FTPSD_BASIS_BASIC then +1 indicates cosine, -1 indicates sine.
   int _layer_i; ///< The index of the current layer.

   int _useBasis; ///< Set to  MXAO_FTPSD_BASIS_BASIC/MODIFIED/PROJECTED_* to use the basic sin/cos modes, the modified Fourier modes, or a projection of them.

   ///Worskspace for the gsl integrators, allocated to WSZ if constructed as worker (with allocate == true).
   gsl_integration_workspace * _w;

   realT _absTol; ///< The absolute tolerance to use in the GSL integrator
   realT _relTol; ///< The relative tolerance to use in the GSL integrator

  int m_mode_i; ///< Projected basis mode index

  Eigen::Array<realT, -1, -1> m_modeCoeffs; ///< Coeeficients of the projection onto the Fourier modes
  realT m_minCoeffVal;


  std::vector<realT> Jps;
  std::vector<realT> Jms;
  std::vector<int> ps;
  std::vector<realT> ms;
  std::vector<realT> ns;

   void initProjection()
   {
      Jps.resize(m_modeCoeffs.cols());
      Jms.resize(m_modeCoeffs.cols());
      ps.resize(m_modeCoeffs.cols());
      ms.resize(m_modeCoeffs.cols());
      ns.resize(m_modeCoeffs.cols());

      for(int i=0; i< m_modeCoeffs.cols(); ++i)
      {
         int m, n, p;
         sigproc::fourierModeCoordinates( m, n, p, i);
         ps[i] = p;
         ms[i] = m;
         ns[i] = n;
      }
  }

public:
   ///Default c'tor
   fourierTemporalPSD();

   ///Constructor with workspace allocation
   /**
     * \param allocate if true, then the workspace for GSL integrators is allocated.
     */
   explicit fourierTemporalPSD(bool allocate);

   ///Destructor
   /** Frees GSL workspace if it was allocated.
     */
   ~fourierTemporalPSD();

protected:
   ///Initialize parameters to default values.
   void initialize();

public:

   /** \name GSL Integration Tolerances
     * For good results it seems that absolute tolerance (absTol) needs to be 1e-10.  Lower tolerances cause some frequencies to drop out, etc.
     * Relative tolerance (relTol) seems to be less sensitive, and 1e-4 works on cases tested as of 1 Jan, 2017.
     *
     * See the documentation for the GSL Library integrators at (https://www.gnu.org/software/gsl/manual/htmlm_node/QAGI-adaptive-integration-on-infinite-intervals.html)
     * @{
     */

   ///Set absolute tolerance
   /**
     * \param at is the new absolute tolerance.
     */
   void absTol(realT at);

   ///Get the current absolute tolerance
   /**
     * \returns _absTol
     */
   realT absTol();

   ///Set relative tolerance
   /**
     * \param rt is the new relative tolerance.
     */
   void relTol(realT rt);

   ///Get the current relative tolerance
   /**
     * \returns _relTol
     */
   realT relTol();

   ///@}

   ///Determine the frequency of the highest V-dot-k peak
   /**
     * \param m the spatial frequency u index
     * \param n the spatial frequency v index
     *
     * \return the frequency of the fastest peak
     */
   realT fastestPeak( int m,
                       int n );


   ///Calculate the temporal PSD for a Fourier mode for a single layer.
   /**
     *
     * \todo implement error checking.
     * \todo need a way to track convergence failures in integral without throwing an error.
     * \todo need better handling of averaging for the -17/3 extension.
     *
     */
   int singleLayerPSD( std::vector<realT> &PSD,  ///< [out] the calculated PSD
                       std::vector<realT> &freq, ///< [in] the populated temporal frequency grid defining the frequencies at which the PSD is calculated
                       realT m,                  ///< [in] the first index of the spatial frequency
                       realT n,                  ///< [in] the second index of the spatial frequency
                       int layer_i,              ///< [in] the index of the layer, for accessing the atmosphere parameters
                       int p,                    ///< [in] sets which mode is calculated (if basic modes, p = -1 for sine, p = +1 for cosine)
                       realT fmax = 0            ///< [in] [optional] set the maximum temporal frequency for the calculation.  The PSD is filled in with a -17/3 power law past this frequency.
                     );      

   ///\cond multilayerm_parallel
   //Conditional to exclude from Doxygen.

protected:
   //Type to allow overloading of the multiLayerPSD workers based on whether they are parallelized or not.
   template<bool m_parallel>
   struct isParallel{};

   //Parallelized version of multiLayerPSD, with OMP directives.
   int m_multiLayerPSD( std::vector<realT> & PSD,
                       std::vector<realT> & freq,
                       realT m,
                       realT n,
                       int p,
                       realT fmax,
                       isParallel<true> parallel );

   //Non-Parallelized version of multiLayerPSD, without OMP directives.
   int m_multiLayerPSD( std::vector<realT> & PSD,
                       std::vector<realT> & freq,
                       realT m,
                       realT n,
                       int p,
                       realT fmax,
                       isParallel<false> parallel );

   ///\endcond

public:
   ///Calculate the temporal PSD for a Fourier mode in a multi-layer model.
   /**
     *
     * \tparam parallel controls whether layers are calculated in parallel.  Default is true.  Set to false if this is called inside a parallelized loop, as in \ref makePSDGrid.
     *
     * \todo implement error checking
     * \todo handle reports of convergence failures form singleLayerPSD when implemented.
     *
     */
   template<bool parallel=true>
   int multiLayerPSD( std::vector<realT> & PSD,  ///< [out] the calculated PSD
                      std::vector<realT> & freq, ///< [in] the populated temporal frequency grid defining at which frequencies the PSD is calculated
                      realT m,                   ///< [in] the first index of the spatial frequency
                      realT n,                   ///< [in] the second index of the spatial frequency
                      int p,                     ///< [in] sets which mode is calculated (if basic modes, p = -1 for sine, p = +1 for cosine)
                      realT fmax = 0             ///< [in] [optional] set the maximum temporal frequency for the calculation. The PSD is filled in
                                                 /// with a -17/3 power law past this frequency.  If 0, then it is taken to be 150 Hz + 2*fastestPeak(m,n).
                    );

   ///Calculate PSDs over a grid of spatial frequencies.
   /** The grid of spatial frequencies is square, set by the maximum value of m and n.
     *
     * The PSDs are written as mx::binVector binary files to a directory.  We do not use FITS since
     * this adds overhead and cfitisio handles parallelization poorly due to the limitation on number of created file pointers.
     *
     */
   void makePSDGrid( const std::string & dir, ///< [in] the directory for output of the PSDs
                     int mnMax,               ///< [in] the maximum value of m and n in the grid.
                     realT dFreq,             ///< [in] the temporal frequency spacing.
                     realT maxFreq,           ///< [in] the maximum temporal frequency to calculate
                     realT fmax = 0           ///< [in] [optional] set the maximum temporal frequency for the calculation. The PSD is filled in with a -17/3 power law past
                                              /// this frequency.  If 0, then it is taken to be 150 Hz + 2*fastestPeak(m,n).
                   );

   ///Analyze a PSD grid under closed-loop control.
   /** This always analyzes the simple integrator, and can also analyze the linear preditor controller.
     */
   int analyzePSDGrid( const std::string & subDir,         ///< [out] the sub-directory of psdDir where to write the results.
                       const std::string & psdDir,         ///< [in]  the directory containing the grid of PSDs.
                       int mnMax,                          ///< [in]  the maximum value of m and n in the grid.
                       int mnCon,                          ///< [in]  the maximum value of m and n which can be controlled.
                       int lpNc,                           ///< [in]  the number of linear predictor coefficients to analyze.  If 0 then LP is not analyzed.
                       std::vector<realT> & mags,          ///< [in]  the guide star magnitudes to analyze for.
                       int lifetimeTrials = 0,             ///< [in]  [optional] number of trials used for calculating speckle lifetimes.  If 0, lifetimes are not calculated. 
                       bool uncontrolledLifetimes = false, ///< [in] [optional] flag controlling whether lifetimes are calculate for uncontrolled modes.
                       bool writePSDs = false              ///< [in]  [optional] flag controlling if resultant PSDs are saved
                     );



   /** \name Disk Storage
     * These methods handle writing to and reading from disk.  The calculated PSDs are store in the mx::BinVector binary format.
     *
     * A grid of PSDs is specified by its directory name.  The directory contains one frequency file (freq.binv), and a set of
     * PSD files, named according to psd_\<m\>_\<n\>_.binv.
     *
     *
     * @{
     */
   ///Get the frequency scale for a PSD grid.
   /**
     */
   void getGridFreq( std::vector<realT> & freq, ///< [out] the vector to populate with the frequency scale.
                     const std::string & dir  ///< [in] specifies the directory containing the grid.
                   );

   ///Get a single PSD from a PSD grid.
   /**
     */
   void getGridPSD( std::vector<realT> & psd, ///< [out] the vector to populate with the PSD.
                    const std::string & dir, ///< [in] specifies the directory containing the grid.
                    int m,  ///< [in] specifies the u component of spatial frequency.
                    int n  ///< [in] specifies the v component of spatial frequency.
                  );

   ///Get both the frequency scale and a single PSD from a PSD grid.
   /**
     */
   void getGridPSD( std::vector<realT> & freq, ///< [out] the vector to populate with the frequency scale.
                    std::vector<realT> & psd, ///< [out] the vector to populate with the PSD.
                    const std::string & dir, ///< [in] specifies the directory containing the grid.
                    int m, ///< [in] specifies the u component of spatial frequency.
                    int n ///< [in] specifies the v component of spatial frequency.
                  );

   ///@}
};

template<typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::fourierTemporalPSD()
{
   m_aosys = nullptr;
   initialize();
}

template<typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::fourierTemporalPSD(bool allocate)
{
   m_aosys = nullptr;
   initialize();

   if(allocate)
   {
      _w = gsl_integration_workspace_alloc (WSZ);
   }
}

template<typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::~fourierTemporalPSD()
{
   if(_w)
   {
      gsl_integration_workspace_free (_w);
   }
}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::initialize()
{
   _useBasis = basis::modified;
   _w = 0;

   _absTol = 1e-10;
   _relTol = 1e-4;
}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::absTol(realT at)
{
   _absTol = at;
}

template<typename realT, typename aosysT>
realT fourierTemporalPSD<realT, aosysT>::absTol()
{
   return _absTol;
}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::relTol(realT rt)
{
   _relTol = rt;
}

template<typename realT, typename aosysT>
realT fourierTemporalPSD<realT, aosysT>::relTol()
{
   return _relTol;
}

template<typename realT, typename aosysT>
realT fourierTemporalPSD<realT, aosysT>::fastestPeak( int m,
                                                        int n )
{
   realT ku, kv, vu,vv;

   ku = ( (realT) m / m_aosys->D());
   kv = ( (realT) n / m_aosys->D());

   realT f, fmax = 0;

   for(size_t i=0; i< m_aosys->atm.n_layers(); ++i)
   {
      vu = m_aosys->atm.layer_v_wind(i) * cos(m_aosys->atm.layer_dir(i));
      vv = m_aosys->atm.layer_v_wind(i) * sin(m_aosys->atm.layer_dir(i));

      f = fabs(ku*vu + kv*vv);
      if(f > fmax) fmax = f;
   }

   return fmax;
}

template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::singleLayerPSD( std::vector<realT> &PSD,
                                                        std::vector<realT> &freq,
                                                        realT m,
                                                        realT n,
                                                        int layer_i,
                                                        int p,
                                                        realT fmax )
{
   if(fmax == 0) fmax = freq[freq.size()-1];

   realT v_wind = m_aosys->atm.layer_v_wind(layer_i);
   realT q_wind = m_aosys->atm.layer_dir(layer_i);


   //Rotate the basis
   realT cq = cos(q_wind);
   realT sq = sin(q_wind);


   realT scale = 2*(1/v_wind); //Factor of 2 for negative frequencies.

   //We'll get the occasional failure to reach tolerance error, just ignore them all for now.
   gsl_set_error_handler_off();

   //Create a local instance so that we're reentrant
   fourierTemporalPSD<realT, aosysT> params(true);

   params.m_aosys = m_aosys;
   params._layer_i = layer_i;
   params.m_m = m*cq + n*sq;
   params.m_n = -m*sq + n*cq;
   params.m_cq = cq; //for de-rotating ku and kv for spatial filtering
   params.m_sq = sq; //for de-rotation ku and kv for spatial filtering
   if(m_aosys->spatialFilter_ku() < std::numeric_limits<realT>::max() || m_aosys->spatialFilter_kv() < std::numeric_limits<realT>::max()) params.m_spatialFilter = true; 
   
   params.m_p = p;


   params.m_mode_i = m_mode_i;
   params.m_modeCoeffs = m_modeCoeffs;
   params.m_minCoeffVal = m_minCoeffVal;

   realT result, error;

   //Setup the GSL calculation
   gsl_function func;
   switch(_useBasis)
   {
      case basis::basic: //MXAO_FTPSD_BASIS_BASIC:
         func.function = &F_basic<realT, aosysT>;
         break;
      case basis::modified: //MXAO_FTPSD_BASIS_MODIFIED:
         func.function = &F_mod<realT, aosysT>;
         break;
      case basis::projected_basic: //MXAO_FTPSD_BASIS_PROJECTED_BASIC:
         mxError("fourierTemporalPSD::singleLayerPSD", MXE_NOTIMPL, "Projected basic-basis modes are not implemented.");
         break;
      case basis::projectedm_modified:// MXAO_FTPSD_BASIS_PROJECTED_MODIFIED:
         params.Jps = Jps;
         params.Jms = Jms;
         params.ms = ms;
         params.ns = ns;
         params.ps = ps;


         func.function = &Fm_projMod<realT, aosysT>;

         break;
      default:
         mxError("fourierTemporalPSD::singleLayerPSD", MXE_INVALIDARG, "value of _useBasis isn not valid.");
         return -1;
   }

   func.params = &params;


   //Here we only calculate up to fmax.
   size_t i=0;
   while( freq[i] <= fmax )
   {
      params.m_f = freq[i];

      int ec = gsl_integration_qagi (&func, _absTol, _relTol, WSZ, params._w, &result, &error);

      if(ec == GSL_EDIVERGE)
      {
         std::cerr << "GSL_EDIVERGE:" << p << " " << freq[i] << " " << v_wind << " " << m << " " << n << " " << m_m << " " << m_n << "\n";
         std::cerr << "ignoring . . .\n";
      }

      PSD[i] = scale*result;

      ++i;
      if(i >= freq.size()) break;
   }

   //Now fill in from fmax to the actual max frequency with a -(alpha+2) power law.
   size_t j=i;
   
   if(j == freq.size()) return 0;
   
   //First average result for last 50.
   PSD[j] = PSD[i-50] * pow( freq[i-50]/freq[j], m_aosys->atm.alpha()+2);//seventeen_thirds<realT>());
   for(size_t k=49; k> 0; --k)
   {
      PSD[j] +=  PSD[i-k] * pow( freq[i-k]/freq[j], m_aosys->atm.alpha()+2); //seventeen_thirds<realT>());
   }
   PSD[j] /= 50.0;
   ++j;
   ++i;
   if(j == freq.size()) return 0;
   while(j < freq.size())
   {
      PSD[j] = PSD[i-1] * pow( freq[i-1]/freq[j], m_aosys->atm.alpha()+2); //seventeen_thirds<realT>());
      ++j;
   }


   return 0;
}


template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::m_multiLayerPSD( std::vector<realT> & PSD,
                                                       std::vector<realT> & freq,
                                                       realT m,
                                                       realT n,
                                                       int p,
                                                       realT fmax,
                                                       isParallel<true> parallel )
{
   static_cast<void>(parallel);
   
   #pragma omp parallel
   {
      //Records each layer PSD
      std::vector<realT> single_PSD(freq.size());

      #pragma omp for
      for(size_t i=0; i< m_aosys->atm.n_layers(); ++i)
      {
         singleLayerPSD(single_PSD, freq, m, n, i, p, fmax);

         //Now add the single layer PSD to the overall PSD, weighted by Cn2
         #pragma omp critical
         for(size_t j=0; j<freq.size(); ++j)
         {
           PSD[j] += m_aosys->atm.layer_Cn2(i)*single_PSD[j];
         }
      }
   }

   return 0;
}

template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::m_multiLayerPSD( std::vector<realT> & PSD,
                                                        std::vector<realT> & freq,
                                                        realT m,
                                                        realT n,
                                                        int p,
                                                        realT fmax,
                                                        isParallel<false>  parallel )
{
   static_cast<void>(parallel);
   
   //Records each layer PSD
   std::vector<realT> single_PSD(freq.size());

   for(size_t i=0; i< m_aosys->atm.n_layers(); ++i)
   {
      singleLayerPSD(single_PSD, freq, m, n, i, p, fmax);

      //Now add the single layer PSD to the overall PSD, weighted by Cn2
      for(size_t j=0; j<freq.size(); ++j)
      {
         PSD[j] += m_aosys->atm.layer_Cn2(i)*single_PSD[j];
      }
   }

   return 0;
}


template<typename realT, typename aosysT>
template<bool parallel>
int fourierTemporalPSD<realT, aosysT>::multiLayerPSD( std::vector<realT> & PSD,
                                                       std::vector<realT> & freq,
                                                       realT m,
                                                       realT n,
                                                       int p,
                                                       realT fmax )
{
   //PSD is zeroed every time to make sure we don't accumulate on repeated calls
   for(size_t j=0;j<PSD.size(); ++j) PSD[j] = 0;

   if(fmax == 0)
   {
      fmax = 150 + 2*fastestPeak(m, n);
   }
  
   return m_multiLayerPSD( PSD, freq, m, n, p, fmax, isParallel<parallel>());


}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::makePSDGrid( const std::string & dir,
                                                     int mnMax,
                                                     realT dFreq,
                                                     realT maxFreq,
                                                     realT fmax
                                                   )
{
   std::vector<realT> freq;

   std::vector<sigproc::fourierModeDef> spf;

   std::string fn;

   sigproc::makeFourierModeFreqs_Rect(spf, 2*mnMax);

   //Calculate number of samples, and make sure we get to at least maxFreq
   int N = (int) maxFreq/dFreq;
   if( N * dFreq < maxFreq) N += 1;


   /*** Dump Params to file ***/
   mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

   std::ofstream fout;
   fn = dir + '/' + "params.txt";
   fout.open(fn);

   fout << "#---------------------------\n";
   m_aosys->dumpAOSystem(fout);
   fout << "#---------------------------\n";
   fout << "# PSD Grid Parameters\n";
   fout << "#    absTol " << _absTol << '\n';
   fout << "#    relTol " << _relTol << '\n';
   fout << "#    useBasis " << _useBasis << '\n';
   fout << "#    makePSDGrid call:\n";
   fout << "#       mnMax = " << mnMax << '\n';
   fout << "#       dFreq = " << dFreq << '\n';
   fout << "#       maxFreq = " << maxFreq << '\n';
   fout << "#       fmax = " << fmax << '\n';
   fout << "#---------------------------\n";


   fout.close();

   //Create frequency scale.
   math::vectorScale(freq, N, dFreq, dFreq); //offset from 0 by dFreq, so f=0 not included

   fn = dir + '/' + "freq.binv";

   ioutils::writeBinVector( fn, freq);

   size_t nLoops = 0.5*spf.size();

   ipc::ompLoopWatcher<> watcher(nLoops, std::cout);

   #pragma omp parallel
   {
      std::vector<realT> PSD;
      PSD.resize(N);
      std::string fname;

      int m, n;

      #pragma omp for
      for(size_t i=0; i<nLoops; ++i)
      {
         m = spf[i*2].m;
         n = spf[i*2].n;

         if(fabs((realT)m/m_aosys->D()) >= m_aosys->spatialFilter_ku() || fabs((realT)n/m_aosys->D()) >= m_aosys->spatialFilter_kv()) 
         {
            watcher.incrementAndOutputStatus();
            continue;
         }
         
         multiLayerPSD<false>( PSD, freq, m, n, 1, fmax);

         fname = dir + '/' + "psd_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";

         ioutils::writeBinVector( fname, PSD);

         watcher.incrementAndOutputStatus();
      }
   }
}

template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::analyzePSDGrid( const std::string & subDir,
                                                       const std::string & psdDir,
                                                       int mnMax,
                                                       int mnCon,
                                                       int lpNc,
                                                       std::vector<realT> & mags,
                                                       int lifetimeTrials,
                                                       bool uncontrolledLifetimes,
                                                       bool writePSDs
                                                     )
{

   std::string dir = psdDir + "/" + subDir;

   /*** Dump Params to file ***/
   mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

   std::ofstream fout;
   std::string fn = dir + '/' + "params.txt";
   fout.open(fn);

   fout << "#---------------------------\n";
   m_aosys->dumpAOSystem(fout);
   fout << "#---------------------------\n";
   fout << "# Analysis Parameters\n";
   fout << "#    mnMax    = " << mnMax << '\n';
   fout << "#    mnCon    = " << mnCon << '\n';
   fout << "#    lpNc     = " << lpNc << '\n';
   fout << "#    mags     = ";
   for(size_t i=0; i<mags.size()-1; ++i) fout << mags[i] << ",";
   fout << mags[mags.size()-1] << '\n';
   fout << "#    lifetimeTrials = " << lifetimeTrials << '\n';
   fout << "#    uncontrolledLifetimes = " << uncontrolledLifetimes << '\n';
   fout << "#    writePSDs = " << std::boolalpha << writePSDs << '\n';
   //fout << "#    intTimes = ";
   //for(int i=0; i<intTimes.size()-1; ++i) fout << intTimes[i] << ",";
   //fout << intTimes[intTimes.size()-1] << '\n';

   fout.close();

   //**** Calculating A Variance Map ****//

   realT fs = 1.0/m_aosys->tauWFS();
   realT tauWFS = m_aosys->tauWFS();
   realT deltaTau = m_aosys->deltaTau();
   
   std::vector<sigproc::fourierModeDef> fms;

   sigproc::makeFourierModeFreqs_Rect(fms, 2*mnMax);
   size_t nModes = 0.5*fms.size();

   Eigen::Array<realT, -1, -1> gains, vars, speckleLifetimes, gains_lp, vars_lp, speckleLifetimes_lp;

   gains.resize(2*mnMax+1, 2*mnMax+1);
   vars.resize(2*mnMax+1, 2*mnMax+1);
   speckleLifetimes.resize(2*mnMax+1, 2*mnMax+1);
   
   gains(mnMax, mnMax) = 0;
   vars(mnMax, mnMax) = 0;
   speckleLifetimes(mnMax, mnMax) = 0;
   
   gains_lp.resize(2*mnMax+1, 2*mnMax+1);
   vars_lp.resize(2*mnMax+1, 2*mnMax+1);
   speckleLifetimes_lp.resize(2*mnMax+1, 2*mnMax+1);
   
   gains_lp(mnMax, mnMax) = 0;
   vars_lp(mnMax, mnMax) = 0;
   speckleLifetimes_lp(mnMax, mnMax) = 0;

   bool doLP = false;
   if(lpNc > 1) doLP = true;
   Eigen::Array<realT, -1, -1> lpC;

   if(doLP)
   {
      lpC.resize(nModes, lpNc);
      lpC.setZero();
   }

   std::vector<realT> S_si, S_lp;

   if(writePSDs)
   {
      for(size_t s=0; s< mags.size(); ++s)
      {
         std::string psdOutDir = dir + "/" + "outputPSDs_" + ioutils::convertToString(mags[s]) + "_si";
         mkdir(psdOutDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

         if(doLP)
         {
            psdOutDir = dir + "/" + "outputPSDs_" + ioutils::convertToString(mags[s]) + "_lp";
            mkdir(psdOutDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         }
      }
   }
   
   
   ipc::ompLoopWatcher<> watcher(nModes*mags.size(), std::cout);

   for(size_t s = 0; s < mags.size(); ++s)
   {
      #pragma omp parallel
      {
         realT localMag = mags[s];

         realT var0;

         realT gopt, var;

         realT gopt_lp, var_lp;

         std::vector<realT> tfreq; //The frequency scale of the PSDs
         std::vector<realT> tPSDp; //The open-loop turbulence PSD for a Fourier mode
         
         //**< Get the frequency grid, and nyquist limit it to f_s/2
         getGridPSD( tfreq, tPSDp, psdDir, 0, 1 ); //To get the freq grid

         size_t imax = 0;
         while( tfreq[imax] <= 0.5*fs ) 
         {
            ++imax;
            if(imax > tfreq.size()-1) break;
         }
         
         if(imax < tfreq.size()-1 && tfreq[imax] <= 0.5*fs*(1.0 + 1e-7)) ++imax;
         
         tfreq.erase(tfreq.begin() + imax, tfreq.end());
         //**>
         
         std::vector<realT> tPSDn; //The open-loop WFS noise PSD         
         tPSDn.resize(tfreq.size()); 

         //**< Setup the controllers
         mx::AO::analysis::clAOLinearPredictor<realT> tflp;
         mx::AO::analysis::clGainOpt<realT> go_si(tauWFS, deltaTau);
         mx::AO::analysis::clGainOpt<realT> go_lp(tauWFS, deltaTau);

         go_si.f(tfreq);
         go_lp.f(tfreq);

         realT gmax = 0;
         realT gmax_lp = 0;
         //**>
         
         int m, n;
         
         
         //**< For use in lifetime calculations
         sigproc::psdVarMean<sigproc::psdVarMeanParams<realT>> pvm;
         std::vector<std::complex<realT>> ETFxn;
         std::vector<std::complex<realT>> NTFxn;
         
         if(lifetimeTrials > 0)
         {
            ETFxn.resize(tfreq.size());
            NTFxn.resize(tfreq.size());
    
            std::string tfOutFile = dir + "/" + "outputTF_" + ioutils::convertToString(mags[s]) + "_si/";
            std::cerr << "Creating " << tfOutFile << "\n";
            ioutils::createDirectories(tfOutFile);
            
            if(doLP)
            {
               tfOutFile = dir + "/" + "outputTF_" + ioutils::convertToString(mags[s]) + "_lp/";
               std::cerr << "Creating " << tfOutFile << "\n";
               ioutils::createDirectories(tfOutFile);
            }
   
         }
         
         //**>
         
         //want to schedule dynamic with small chunks so maximal processor usage,
         //otherwise we can end up with a small number of cores being used at the end
         #pragma omp for schedule(dynamic, 5)
         for(size_t i=0; i<nModes; ++i)
         {
            //Determine the spatial frequency at this step
            m = fms[2*i].m;
            n = fms[2*i].n;
            
            if(fabs((realT)m/m_aosys->D()) >= m_aosys->spatialFilter_ku() || fabs((realT)n/m_aosys->D()) >= m_aosys->spatialFilter_kv())
            {
               gains( mnMax + m, mnMax + n ) = 0;
               gains( mnMax - m, mnMax - n ) = 0;
               
               gains_lp( mnMax + m, mnMax + n ) = 0;
               gains_lp( mnMax - m, mnMax - n ) = 0;
               
               vars( mnMax + m, mnMax + n) = 0;
               vars( mnMax - m, mnMax - n ) = 0;
               
               vars_lp( mnMax + m, mnMax + n) = 0;
               vars_lp( mnMax - m, mnMax - n ) = 0;
               speckleLifetimes( mnMax + m, mnMax + n ) = 0;
               speckleLifetimes( mnMax - m, mnMax - n ) = 0;
               speckleLifetimes_lp( mnMax + m, mnMax + n ) = 0;
               speckleLifetimes_lp( mnMax - m, mnMax - n ) = 0;
            }
            else
            {
         
               realT k = sqrt(m*m + n*n)/m_aosys->D();
               
               //Get the WFS noise PSD (which is already resized to match tfreq)
               wfsNoisePSD<realT>( tPSDn, m_aosys->beta_p(m,n), m_aosys->Fg(localMag), tauWFS, m_aosys->npix_wfs((size_t) 0), m_aosys->Fbg((size_t) 0), m_aosys->ron_wfs((size_t) 0));
               
               //**< Get the open-loop turb. PSD
               getGridPSD( tPSDp, psdDir, m, n );
               
               //Get integral of entire open-loop PSD
               var0 = sigproc::psdVar( tfreq, tPSDp);
               
               //erase points above Nyquist limit
               tPSDp.erase(tPSDp.begin() + imax, tPSDp.end());
               
               //And now determine the variance which has been erased.
               //limVar is the out-of-band variance, which we add back in for completeness
               realT limVar = var0 - sigproc::psdVar( tfreq, tPSDp);
               //**>
               
               //Get the variance of the theoretical PSD
               //var0 = m_aosys->psd(m_aosys->atm, k, m_aosys->lam_sci(), m_aosys->lam_wfs(), m_aosys->secZeta() ) / pow(m_aosys->D(),2);
               
               
               //**< Determine if we're inside the hardwarecontrol limit
               bool inside = false;
               
               if( m_aosys->circularLimit() )
               {
                  if( m*m + n*n <= mnCon*mnCon) inside = true;
               }
               else
               {
                  if(fabs(m) <= mnCon && fabs(n) <= mnCon) inside = true;
               }
               //**>
               
               
               if(inside)
               {
                  gmax = 0;
                  gopt = go_si.optGainOpenLoop(var, tPSDp, tPSDn, gmax);
               
                  var += limVar;
                  
                  if(doLP)
                  {
                     tflp.regularizeCoefficients( gmax_lp, gopt_lp, var_lp, go_lp, tPSDp, tPSDn, lpNc);
                     for(int n=0; n< lpNc; ++n) lpC(i,n) = go_lp.a(n);
                     
                     var_lp += limVar;
                  }
                  else
                  {
                     gopt_lp = 0;
                  }
               }
               else
               {
                  gopt = 0;
                  var = var0;
                  var_lp = var0;
                  gopt_lp = 0;
                  go_lp.a(std::vector<realT>({1}));
                  go_lp.b(std::vector<realT>({1}));
               }
               
               //**< Determine if closed-loop is making a difference:
               
               //mult used to be 1.3 when using the vonKarman PSD instead of open-loop integral
               // 1.3 was a hack since we haven't yet done the convolution...
               
               realT mult = 1;
               if(gopt > 0 && var > mult*var0)
               {
                  gopt = 0;
                  var = var0;
               }
               
               if(gopt_lp > 0 && var_lp > mult*var0)
               {
                  gopt_lp = 0;
                  var_lp = var0;
               }
               //**>
               
               //**< Fill in the gain and variance maps
               gains( mnMax + m, mnMax + n ) = gopt;
               gains( mnMax - m, mnMax - n ) = gopt;
               
               gains_lp( mnMax + m, mnMax + n ) = gopt_lp;
               gains_lp( mnMax - m, mnMax - n ) = gopt_lp;
               
               vars( mnMax + m, mnMax + n) = var;
               vars( mnMax - m, mnMax - n ) = var;
               
               vars_lp( mnMax + m, mnMax + n) = var_lp;
               vars_lp( mnMax - m, mnMax - n ) = var_lp;
               //**>
                
               //**< Calulcate Speckle Lifetimes
               if(lifetimeTrials > 0 && ( uncontrolledLifetimes || inside ))
               {
                  std::vector<realT> spfreq, sppsd;
                                    
                  if(gopt > 0)
                  {
                     for(size_t i=0;i<tfreq.size();++i)
                     {
                        ETFxn[i] = go_si.clETF(i, gopt);
                        NTFxn[i] = go_si.clNTF(i, gopt);
                     }
                  }
                  else
                  {
                     for(size_t i=0;i<tfreq.size();++i)
                     {
                        ETFxn[i] = 1;
                        NTFxn[i] = 0;
                     }
                  }
                  
                  
                  std::string tfOutFile = dir + "/" + "outputTF_" + ioutils::convertToString(mags[s]) + "_si/";
                  
                  std::string etfOutFile = tfOutFile +  "etf_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
                  ioutils::writeBinVector( etfOutFile, ETFxn);
               
                  std::string ntfOutFile = tfOutFile +  "ntf_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
                  ioutils::writeBinVector( ntfOutFile, NTFxn);
                  
                  if(i==0) //Write freq on the first one
                  {
                     std::string fOutFile = tfOutFile + "freq.binv";
                     ioutils::writeBinVector(fOutFile, tfreq);
                  }
#if 0
                  speckleAmpPSD( spfreq, sppsd, tfreq, tPSDp, ETFxn, tPSDn, NTFxn, lifetimeTrials);
                  realT spvar = sigproc::psdVar(spfreq, sppsd);
                  
                  realT splifeT = 100.0;
                  realT error;
                  //realT tau = sppsd[0]/(2*spvar);// bins[1]/(2.*tfreq.back())*vars[1]/vars[0];
                  realT tau = pvm(error, spfreq, sppsd, splifeT) * (splifeT)/spvar;
                     
                  speckleLifetimes( mnMax + m, mnMax + n ) = tau;
                  speckleLifetimes( mnMax - m, mnMax - n ) = tau;
#endif
                  
                  if(doLP)
                  {
                     if(gopt_lp > 0)
                     {
                        for(size_t i=0;i<tfreq.size();++i)
                        {
                           ETFxn[i] = go_lp.clETF(i, gopt_lp);
                           NTFxn[i] = go_lp.clNTF(i, gopt_lp);
                        }
                     }
                     else
                     {
                        for(size_t i=0;i<tfreq.size();++i)
                        {
                           ETFxn[i] = 1;
                           NTFxn[i] = 0;
                        }
                     }   
                     
                     std::string tfOutFile = dir + "/" + "outputTF_" + ioutils::convertToString(mags[s]) + "_lp/";
                     
                     std::string etfOutFile = tfOutFile +  "etf_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
                     ioutils::writeBinVector( etfOutFile, ETFxn);
                     
                     std::string ntfOutFile = tfOutFile +  "ntf_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
                     ioutils::writeBinVector( ntfOutFile, NTFxn);
                     
                     if(i==0) //Write freq on the first one
                     {
                        std::string fOutFile = tfOutFile + "freq.binv";
                        ioutils::writeBinVector(fOutFile, tfreq);
                     }
                  
#if 0
                     speckleAmpPSD( spfreq, sppsd, tfreq, tPSDp, ETFxn, tPSDn, NTFxn, lifetimeTrials);
                     realT spvar = sigproc::psdVar(spfreq, sppsd);
                  
                     //realT tau = sppsd[0]/(2*spvar);// bins[1]/(2.*tfreq.back())*vars[1]/vars[0];
                     realT tau = pvm(error, spfreq, sppsd, splifeT) * (splifeT)/spvar;
                  
                     speckleLifetimes_lp( mnMax + m, mnMax + n ) = tau;
                     speckleLifetimes_lp( mnMax - m, mnMax - n ) = tau;
#endif
                  }
                  
               }
               
               //**>
               
               //Calculate the controlled PSDs and output
               if(writePSDs)
               {
                  std::string psdOutFile = dir + "/" + "outputPSDs_" + ioutils::convertToString(mags[s]) + "_si/";
                  psdOutFile +=  "psd_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
               
                  std::vector<realT> psdOut(tPSDp.size());
               
                  //Calculate the output PSD if gains are applied
                  if(gopt > 0)
                  {
                     realT ETF, NTF;
               
                     for(size_t i=0; i< tfreq.size(); ++i)
                     {
                        go_si.clTF2( ETF, NTF, i,gopt);
                        psdOut[i] = tPSDp[i]*ETF + tPSDn[i]*NTF;
                     }
                  }
                  else //otherwise just copy
                  {
                     psdOut = tPSDp;
                  }
               
                  ioutils::writeBinVector( psdOutFile, psdOut);
               
                  if(i==0) //Write freq on the first one
                  {
                     psdOutFile = dir + "/" + "outputPSDs_" + ioutils::convertToString(mags[s]) + "_si/freq.binv";
                     ioutils::writeBinVector(psdOutFile, tfreq);
                  }
                  
                  
                  
                  if(doLP)
                  {
                     psdOutFile = dir + "/" + "outputPSDs_" + ioutils::convertToString(mags[s]) + "_lp/";
                     psdOutFile +=  "psd_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
               
                     //Calculate the output PSD if gains are applied
                     if(gopt_lp > 0)
                     {
                        realT ETF, NTF;
               
                        for(size_t i=0; i < tfreq.size(); ++i)
                        {
                           go_lp.clTF2( ETF, NTF, i, gopt_lp);
                           psdOut[i] = tPSDp[i]*ETF + tPSDn[i]*NTF;
                        }
                     }
                     else //otherwise just copy
                     {
                        psdOut = tPSDp;
                     }
               
                     ioutils::writeBinVector( psdOutFile, psdOut);
               
                     if(i==0)
                     {
                        psdOutFile = dir + "/" + "outputPSDs_" + ioutils::convertToString(mags[s]) + "_lp/freq.binv";
                        ioutils::writeBinVector(psdOutFile, tfreq);
                     }
                  }
               }
            }
            watcher.incrementAndOutputStatus();
            
         } //omp for i..nModes
      }//omp Parallel

      Eigen::Array<realT, -1,-1> cim, psf;

      //Create Airy PSF for convolution with variance map.
      psf.resize(77,77);
      for(int i=0;i<psf.rows();++i)
      {
         for(int j=0;j<psf.cols();++j)
         {
            psf(i,j) = mx::math::func::airyPattern(sqrt( pow( i-floor(.5*psf.rows()),2) + pow(j-floor(.5*psf.cols()),2)));
         }
      }

      fits::fitsFile<realT> ff;
      std::string fn = dir + "/gainmap_" + ioutils::convertToString(mags[s]) + "_si.fits";
      ff.write( fn, gains);

      fn = dir + "/varmap_" + ioutils::convertToString(mags[s]) + "_si.fits";
      ff.write( fn, vars);


      //Perform convolution for uncontrolled modes
      /*mx::AO::analysis::varmapToImage(cim, vars, psf);

      //Now switch to uncovolved variance for controlled modes
      for(int m=0; m<= mnMax; ++m)
      {
         for(int n=-mnMax; n<= mnMax; ++n)
         {
            if(gains(mnMax + m, mnMax + n) > 0)
            {
               cim( mnMax + m, mnMax + n) = vars( mnMax + m, mnMax + n);
               cim( mnMax - m, mnMax - n ) = vars( mnMax + m, mnMax + n);
            }
         }
      }*/
      cim = vars;

      realT S = exp(-1*cim.sum());
      S_si.push_back(S);
      cim /= S;

      fn = dir + "/contrast_" + ioutils::convertToString(mags[s]) + "_si.fits";
      ff.write( fn, cim);

      
      if(lifetimeTrials > 0)
      {
         fn = dir + "/speckleLifetimes_" + ioutils::convertToString(mags[s]) + "_si.fits";
         ff.write( fn, speckleLifetimes);
      }
      
      
      if(doLP)
      {
         fn = dir + "/gainmap_" + ioutils::convertToString(mags[s]) + "_lp.fits";
         ff.write( fn, gains_lp);

         fn = dir + "/lpcmap_" + ioutils::convertToString(mags[s]) + "_lp.fits";
         ff.write( fn, lpC);

         fn = dir + "/varmap_" + ioutils::convertToString(mags[s]) + "_lp.fits";
         ff.write( fn, vars_lp);

         /*mx::AO::analysis::varmapToImage(cim, vars_lp, psf);

         //Now switch to unconvolved variance for controlled modes
         for(int m=0; m<= mnMax; ++m)
         {
            for(int n=-mnMax; n<= mnMax; ++n)
            {
               if(gains_lp(mnMax + m, mnMax + n) > 0)
               {
                  cim( mnMax + m, mnMax + n) = vars_lp( mnMax + m, mnMax + n);
                  cim( mnMax - m, mnMax - n ) = vars_lp( mnMax + m, mnMax + n);
               }
            }
         }*/
         cim = vars_lp;

         realT S = exp(-1*cim.sum());
         S_lp.push_back(S);
         cim /= S;

         fn = dir + "/contrast_" + ioutils::convertToString(mags[s]) + "_lp.fits";
         ff.write( fn, cim);
         
         if(lifetimeTrials > 0)
         {
            fn = dir + "/speckleLifetimes_" + ioutils::convertToString(mags[s]) + "_lp.fits";
            ff.write( fn, speckleLifetimes_lp);
         }
      }

   }//s (mag)

   fn = dir + "/strehl_si.txt";
   fout.open(fn);
   for(size_t i=0;i<mags.size(); ++i)
   {
      fout << mags[i] << " " << S_si[i] << "\n";
   }

   fout.close();

   if(doLP)
   {
      fn = dir + "/strehl_lp.txt";
      fout.open(fn);
      for(size_t i=0;i<mags.size(); ++i)
      {
         fout << mags[i] << " " << S_lp[i] << "\n";
      }

      fout.close();
   }


   return 0;
}



template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::getGridFreq( std::vector<realT> & freq,
                                                     const std::string & dir )
{
   std::string fn;
   fn = dir + '/' + "freq.binv";
   ioutils::readBinVector(freq, fn);
}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::getGridPSD( std::vector<realT> & psd,
                                                     const std::string & dir,
                                                     int m,
                                                     int n )
{
   std::string fn;
   fn = dir + '/' + "psd_" + ioutils::convertToString(m) + '_' + ioutils::convertToString(n) + ".binv";
   ioutils::readBinVector(psd, fn);
}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::getGridPSD( std::vector<realT> & freq,
                                                     std::vector<realT> & psd,
                                                     const std::string & dir,
                                                     int m,
                                                     int n )
{
   getGridFreq( freq, dir );
   getGridPSD(psd, dir, m, n);
}

///Worker function for GSL Integration for the basic sin/cos Fourier modes.
/** \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT F_basic (realT kv, void * params)
{
   fourierTemporalPSD<realT, aosysT> * Fp = (fourierTemporalPSD<realT, aosysT> *) params;

   realT f = Fp->m_f;
   realT v_wind = Fp->m_aosys->atm.layer_v_wind(Fp->_layer_i);

   realT D = Fp->m_aosys->D();
   realT m = Fp->m_m;
   realT n = Fp->m_n;
   int p = Fp->m_p;

   realT ku = f/v_wind;

   realT kp = sqrt( pow(ku + m/D,2) + pow(kv + n/D,2) );
   realT kpp = sqrt( pow(ku - m/D,2) + pow(kv - n/D,2) );

   realT Q1 = math::func::jinc(math::pi<realT>()*D*kp);

   realT Q2 = math::func::jinc(math::pi<realT>()*D*kpp);

   realT Q = (Q1 + p*Q2);

   realT P =  Fp->m_aosys->psd(Fp->m_aosys->atm, sqrt( pow(ku,2) + pow(kv,2)),  Fp->m_aosys->secZeta());

   return P*Q*Q ;
}

///Worker function for GSL Integration for the modified Fourier modes.
/** \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT F_mod (realT kv, void * params)
{
   fourierTemporalPSD<realT, aosysT> * Fp = (fourierTemporalPSD<realT, aosysT> *) params;
   

   realT f = Fp->m_f;
   realT v_wind = Fp->m_aosys->atm.layer_v_wind(Fp->_layer_i);

   realT ku = f/v_wind;

   if(Fp->m_spatialFilter)
   {
      //de-rotate the spatial frequency vector back to pupil coordinates
      realT dku = ku*Fp->m_cq - kv*Fp->m_sq;
      realT dkv = ku*Fp->m_sq + kv*Fp->m_cq;
      //Return if spatially filtered
      if(fabs(dku) >= Fp->m_aosys->spatialFilter_ku()) return 0;
   
      if(fabs(dkv) >= Fp->m_aosys->spatialFilter_kv()) return 0;
   }
   
   realT D = Fp->m_aosys->D();
   realT m = Fp->m_m;
   realT n = Fp->m_n;
   
   realT kp = sqrt( pow(ku + m/D,2) + pow(kv + n/D,2) );
   realT kpp = sqrt( pow(ku - m/D,2) + pow(kv - n/D,2) );

   realT Jp = math::func::jinc(math::pi<realT>()*D*kp);

   realT Jm = math::func::jinc(math::pi<realT>()*D*kpp);

   realT QQ = 2*(Jp*Jp + Jm*Jm);

   realT P =  Fp->m_aosys->psd(Fp->m_aosys->atm, sqrt( pow(ku,2) + pow(kv,2)), Fp->m_aosys->lam_sci(), Fp->m_aosys->lam_wfs(), Fp->m_aosys->secZeta() );

   return P*QQ ;
}



///Worker function for GSL Integration for a basis projected onto the modified Fourier modes.
/** \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT Fm_projMod (realT kv, void * params)
{
   fourierTemporalPSD<realT, aosysT> * Fp = (fourierTemporalPSD<realT, aosysT> *) params;

   realT f = Fp->m_f;
   realT v_wind = Fp->m_aosys->atm.layer_v_wind(Fp->_layer_i);

   realT q_wind = Fp->m_aosys->atm.layer_dir(Fp->_layer_i);

   //For rotating the basis
   realT cq = cos(q_wind);
   realT sq = sin(q_wind);

   realT D = Fp->m_aosys->D();

   int mode_i = Fp->m_mode_i;

   realT m;
   realT n;
   int p, pp;// = Fp->m_p;

   realT ku = f/v_wind;

   realT kp, km, Jp, Jpp, Jm, Jmp, QQ;

   QQ = 0;

   for(int i=0; i< Fp->m_modeCoeffs.cols(); ++i)
   {
      if(Fp->m_modeCoeffs(i, Fp->m_modeCoeffs.cols()-1 - mode_i) == 0) continue;

      m = Fp->ms[i]*cq + Fp->ns[i]*sq;
      n = -Fp->ms[i]*sq + Fp->ns[i]*cq;

      kp = sqrt( pow(ku + m/D,2) + pow(kv + n/D,2) );
      km = sqrt( pow(ku - m/D,2) + pow(kv - n/D,2) );

      Fp->Jps[i] = math::func::jinc(math::pi<realT>()*D*kp);
      Fp->Jms[i] = math::func::jinc(math::pi<realT>()*D*km);

   }

   int N = Fp->m_modeCoeffs.cols();

   realT sc;

   for(int i=0; i< N ; ++i)
   {
      if( fabs(Fp->m_modeCoeffs(i, Fp->m_modeCoeffs.cols()-1 - mode_i)) < Fp->m_minCoeffVal) continue;

      Jp = Fp->Jps[i];
      Jm = Fp->Jms[i];
      p = Fp->ps[i];

      QQ += 2*( Jp*Jp + Jm*Jm)*Fp->m_modeCoeffs(i,Fp->m_modeCoeffs.cols()-1 - mode_i)*Fp->m_modeCoeffs(mode_i,i);

      for(int j=(i+1); j < N; ++j)
      {
         if( fabs(Fp->m_modeCoeffs(j,Fp->m_modeCoeffs.cols()-1 - mode_i)) < Fp->m_minCoeffVal) continue;
         sc =  Fp->m_modeCoeffs(i,Fp->m_modeCoeffs.cols()-1 - mode_i)*Fp->m_modeCoeffs(j,Fp->m_modeCoeffs.cols()-1 - mode_i);

         //if(fabs(sc) < m_minCoeffProduct) continue;

         //if( sc*sc < 1e-2) continue;
         Jpp = Fp->Jps[j];

         Jmp = Fp->Jms[j];

         pp = Fp->ps[j];

         if(p == pp)
         {
            QQ += 2*2*( Jp*Jpp + Jm*Jmp)*sc;
         }
         else
         {
            QQ += 2*2*( Jp*Jmp + Jm*Jpp)*sc;
         }

      }
   }


   //realT QQ = 2*(Jp*Jp + Jm*Jm);


   realT P =  Fp->m_aosys->psd(Fp->m_aosys->atm, sqrt( pow(ku,2) + pow(kv,2)), Fp->m_aosys->lam_sci(), Fp->m_aosys->lam_wfs(), Fp->m_aosys->secZeta() );

   return P*QQ ;
}

/*
extern template
struct fourierTemporalPSD<float, aoSystem<float, vonKarmanSpectrum<float>, std::ostream>>;
*/
extern template
struct fourierTemporalPSD<double, aoSystem<double, vonKarmanSpectrum<double>, std::ostream>>;
/*
extern template
struct fourierTemporalPSD<long double, aoSystem<long double, vonKarmanSpectrum<long double>, std::ostream>>;

#ifdef HASQUAD
extern template
struct fourierTemporalPSD<__float128, aoSystem<__float128, vonKarmanSpectrum<__float128>, std::ostream>>;
#endif
*/

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //fourierTemporalPSD_hpp
