/** \file fourierTemporalPSD.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculation of the temporal PSD of Fourier modes.
  * \ingroup mxAO_files
  * 
  */

#ifndef fourierTemporalPSD_hpp
#define fourierTemporalPSD_hpp

#include <iostream>
#include <fstream>

#include <sys/stat.h>

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Eigen/Dense>

#include "../../math/func/jinc.hpp"
#include "../../vectorUtils.hpp"
#include "../../stringUtils.hpp"
#include "../../sigproc/fourierModes.hpp"
#include "../../readColumns.hpp"
#include "../../binVector.hpp"

#include "../../ompLoopWatcher.hpp"
#include "../../mxError.hpp"
#include "../../gslInterpolation.hpp"

#include "wfsNoisePSD.hpp"
#include "clAOLinearPredictor.hpp"
#include "clGainOpt.hpp"
#include "varmapToImage.hpp"

#include "aoConstants.hpp"
using namespace mx::AO::constants;

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
                            projected_modified
                          };
                          

//Forward declaration
template<typename realT, typename aosysT>
realT F_basic (realT kv, void * params) ;

//Forward declaration
template<typename realT, typename aosysT>
realT F_mod (realT kv, void * params); 

//Forward declaration
template<typename realT, typename aosysT>
realT F_projMod (realT kv, void * params); 

///Class to manage the calculation of temporal PSDs of the Fourier modes in atmospheric turbulence.
/** Works with both basic (sines/cosines) and modified Fourier modes.
  *
  * \tparam realT is a real floating point type for calculations.  Currently must be double due to gsl_integration.
  * \tparam aosysT is an AO system type, usually of type ao_system.
  * 
  * \todo Split off the integration parameters in a separate structure.
  * \todo once integration parameters are in a separate structure, make this a class with protected members.
  * 
  * \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
struct fourierTemporalPSD
{
   ///Pointer to an AO system structure.
   aosysT * _aosys;
  
   realT _f; ///< the current temporal frequency 
   realT _m; ///< the spatial frequency m index
   realT _n; ///< the spatial frequency n index
   int _p; ///< The parity of the mode, +/- 1.  If _useBasis==MXAO_FTPSD_BASIS_BASIC then +1 indicates cosine, -1 indicates sine.
   int _layer_i; ///< The index of the current layer.
   
   int _useBasis; ///< Set to  MXAO_FTPSD_BASIS_BASIC/MODIFIED/PROJECTED_* to use the basic sin/cos modes, the modified Fourier modes, or a projection of them.
   
   ///Worskspace for the gsl integrators, allocated to WSZ if constructed as worker (with allocate == true).
   gsl_integration_workspace * _w;
   
   realT _absTol; ///< The absolute tolerance to use in the GSL integrator 
   realT _relTol; ///< The relative tolerance to use in the GSL integrator
  
  int _mode_i; ///< Projected basis mode index

  Eigen::Array<realT, -1, -1> _modeCoeffs; ///< Coeeficients of the projection onto the Fourier modes 

  
  
  std::vector<realT> Jps;
  std::vector<realT> Jms;
  std::vector<int> ps;
  std::vector<realT> ms;
  std::vector<realT> ns;
   
   void initProjection()
   {
      Jps.resize(_modeCoeffs.cols());
      Jms.resize(_modeCoeffs.cols());
      ps.resize(_modeCoeffs.cols());
      ms.resize(_modeCoeffs.cols());
      ns.resize(_modeCoeffs.cols());
         
      for(int i=0; i< _modeCoeffs.cols(); ++i)
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
   fourierTemporalPSD(bool allocate);
   
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
     * See the documentation for the GSL Library integrators at (https://www.gnu.org/software/gsl/manual/html_node/QAGI-adaptive-integration-on-infinite-intervals.html)
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
     * \param[out] PSD is the calculated PSD
     * \param[in] freq is the populated temporal frequency grid defining the frequencies at which the PSD is calculated
     * \param[in] m is the first index of the spatial frequency
     * \param[in] n is the second index of the spatial frequency
     * \param[in] layer_i the index of the layer, for accessing the atmosphere parameters
     * \param[in] p sets which mode is calculated (if basic modes, p = -1 for sine, p = +1 for cosine)
     * \param[in] fmax [optional] set the maximum temporal frequency for the calculation.  The PSD is filled in 
     *                             with a -17/3 power law past this frequency.
     * 
     * \todo implement error checking.
     * \todo need a way to track convergence failures in integral without throwing an error. 
     * \todo need better handling of averaging for the -17/3 extension.
     *  
     */ 
   int singleLayerPSD( std::vector<realT> &PSD, 
                       std::vector<realT> &freq, 
                       realT m, 
                       realT n, 
                       int layer_i, 
                       int p, 
                       realT fmax = 0);

   ///\cond multilayer_parallel
   //Conditional to exclude from Doxygen.
   
protected:
   //Type to allow overloading of the multiLayerPSD workers based on whether they are parallelized or not.
   template<bool _parallel>
   struct isParallel{};

   //Parallelized version of multiLayerPSD, with OMP directives.
   int _multiLayerPSD( std::vector<realT> & PSD,
                       std::vector<realT> & freq, 
                       realT m,
                       realT n,
                       int p,
                       realT fmax,
                       isParallel<true> parallel );

   //Non-Parallelized version of multiLayerPSD, without OMP directives.   
   int _multiLayerPSD( std::vector<realT> & PSD,
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
   int multiLayerPSD( std::vector<realT> & PSD, ///< [out] the calculated PSD
                      std::vector<realT> & freq, ///< [in] the populated temporal frequency grid defining at which frequencies the PSD is calculated
                      realT m, ///< [in] the first index of the spatial frequency
                      realT n, ///< [in] the second index of the spatial frequency
                      int p, ///< [in] sets which mode is calculated (if basic modes, p = -1 for sine, p = +1 for cosine)
                      realT fmax = 0 ///< [in] [optional] set the maximum temporal frequency for the calculation. The PSD is filled in 
                                     /// with a -17/3 power law past this frequency.  If 0, then it is taken to be 150 Hz + 2*fastestPeak(m,n).
                    );

   ///Calculate PSDs over a grid of spatial frequencies.
   /** The grid of spatial frequencies is square, set by the maximum value of m and n.
     * 
     * The PSDs are written as mx::binVector binary files to a directory.  We do not use FITS since
     * this adds overhead and cfitisio handles parallelization poorly due to the limitation on number of created file pointers.
     * 
     */ 
   void makePSDGrid( std::string dir, ///< [in] the directory for output of the PSDs
                     int mnMax, ///< [in] the maximum value of m and n in the grid. 
                     realT dFreq, ///< [in] the temporal frequency spacing.
                     realT maxFreq, ///< [in] the maximum temporal frequency to calculate
                     realT fmax = 0 ///< [in] [optional] set the maximum temporal frequency for the calculation. The PSD is filled in with a -17/3 power law past 
                                    /// this frequency.  If 0, then it is taken to be 150 Hz + 2*fastestPeak(m,n).
                   ); 
   
   ///Analyze a PSD grid under closed-loop control.
   /** This always analyzes the simple integrator, and can also analyze the linear preditor controller.
     */
   int analyzePSDGrid( std::string subDir, ///< [in] the sub-directory of psdDir where to write the results.
                       std::string psdDir,  ///< [in] the directory containing the grid of PSDs. 
                       int mnMax, ///< [in] the maximum value of m and n in the grid.
                       int mnCon, ///< [in] the maximum value of m and n which can be controlled.
                       int lpNc, ///< [in] the number of linear predictor coefficients to analyze.  If 0 then LP is not analyzed.
                       std::vector<realT> & mags, ///< [in] the guide star magnitudes to analyze for.
                       std::vector<int> & intTimes ///< [in] the integration times, in units of aosysT::minTauWFS.
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
     * \param [out] freq is the vector to populate with the frequency scale.
     * \param [in] dir specifies the directory containing the grid.
     */ 
   void getGridFreq( std::vector<realT> & freq,
                     const std::string & dir );

   ///Get a single PSD from a PSD grid.
   /** 
     * \param [out] psd is the vector to populate with the PSD.
     * \param [in] dir specifies the directory containing the grid.
     * \param [in] m specifies the u component of spatial frequency.
     * \param [in] n specifies the v component of spatial frequency.
     */ 
   void getGridPSD( std::vector<realT> & psd,
                    const std::string & dir,
                    int m,
                    int n );
   
   ///Get both the frequency scale and a single PSD from a PSD grid.
   /** 
     * \param [out] freq is the vector to populate with the frequency scale.
     * \param [out] psd is the vector to populate with the PSD.
     * \param [in] dir specifies the directory containing the grid.
     * \param [in] m specifies the u component of spatial frequency.
     * \param [in] n specifies the v component of spatial frequency.
     */ 
   void getGridPSD( std::vector<realT> & freq,
                    std::vector<realT> & psd,
                    const std::string & dir,
                    int m,
                    int n );
       
   ///@}
};

template<typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::fourierTemporalPSD()
{
   initialize();
}

template<typename realT, typename aosysT>
fourierTemporalPSD<realT, aosysT>::fourierTemporalPSD(bool allocate)
{
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
   _useBasis = basis::modified;//  MXAO_FTPSD_BASIS_MODIFIED;
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
   
   ku = ( (realT) m / _aosys->D());
   kv = ( (realT) n / _aosys->D());
   
   realT f, fmax = 0;
   
   for(int i=0; i< _aosys->atm.n_layers(); ++i)
   {
      vu = _aosys->atm.layer_v_wind(i) * cos(_aosys->atm.layer_dir(i));
      vv = _aosys->atm.layer_v_wind(i) * sin(_aosys->atm.layer_dir(i));
      
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

   realT v_wind = _aosys->atm.layer_v_wind(layer_i);
   realT q_wind = _aosys->atm.layer_dir(layer_i);
   realT r_0 = _aosys->atm.r_0();
   realT D = _aosys->D();

   //Rotate the basis
   realT cq = cos(q_wind);
   realT sq = sin(q_wind);

   
   realT scale = 2*(1/v_wind); //Factor of 2 for negative frequencies.

   //We'll get the occasional failure to reach tolerance error, just ignore them all for now.
   gsl_set_error_handler_off();
   
   //Create a local instance so that we're reentrant
   fourierTemporalPSD<realT, aosysT> params(true);
   
   params._aosys = _aosys;
   params._layer_i = layer_i;
   params._m = m*cq + n*sq;
   params._n = -m*sq + n*cq;
   params._p = p;

   
   params._mode_i = _mode_i;
   params._modeCoeffs = _modeCoeffs;

   
   realT result, error, result2;

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
      case basis::projected_modified:// MXAO_FTPSD_BASIS_PROJECTED_MODIFIED:
         params.Jps = Jps;
         params.Jms = Jms;
         params.ms = ms;
         params.ns = ns;
         params.ps = ps;

   
         func.function = &F_projMod<realT, aosysT>;
         
         break; 
      default:
         mxError("fourierTemporalPSD::singleLayerPSD", MXE_INVALIDARG, "value of _useBasis isn not valid.");
         return -1;
   }
   
   func.params = &params;


   //Here we only calculate up to fmax.
   int i=0;      
   while(freq[i] <= fmax && i < freq.size())
   {
      params._f = freq[i];
   
      int ec = gsl_integration_qagi (&func, _absTol, _relTol, WSZ, params._w, &result, &error);
      
      if(ec == GSL_EDIVERGE)
      {
         std::cerr << "GSL_EDIVERGE:" << p << " " << freq[i] << " " << v_wind << " " << m << " " << n << " " << _m << " " << _n << "\n";
         std::cerr << "ignoring . . .\n";
      }
   
      PSD[i] = scale*result;
   
      ++i;
   }
   
   //Now fill in from fmax to the actual max frequency with a -17/3 power law.
   int j=i;
   
   //First average result for last 50.
   PSD[j] = PSD[i-50] * pow( freq[i-50]/freq[j], seventeen_thirds<realT>());
   for(int k=49; k> 0; --k)
   {
      PSD[j] +=  PSD[i-k] * pow( freq[i-k]/freq[j], seventeen_thirds<realT>());
   }
   PSD[j] /= 50.0;
   ++j;
   ++i;
   while(j < freq.size()) 
   {
      PSD[j] = PSD[i-1] * pow( freq[i-1]/freq[j], seventeen_thirds<realT>());
      ++j;
   }
   
   
   
   return 0;
}
 

template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::_multiLayerPSD( std::vector<realT> & PSD,
                                                       std::vector<realT> & freq, 
                                                       realT m,
                                                       realT n,
                                                       int p,
                                                       realT fmax,
                                                       isParallel<true> parallel )
{
   #pragma omp parallel
   {
      //Records each layer PSD
      std::vector<realT> single_PSD(freq.size());

      #pragma omp for
      for(int i=0; i< _aosys->atm.n_layers(); ++i)
      {
         singleLayerPSD(single_PSD, freq, m, n, i, p, fmax);

         //Now add the single layer PSD to the overall PSD, weighted by Cn2
         #pragma omp critical
         for(int j=0; j<freq.size(); ++j)
         {
           PSD[j] += _aosys->atm.layer_Cn2(i)*single_PSD[j]; 
         }
      } 
   }
   
   return 0;
}

template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::_multiLayerPSD( std::vector<realT> & PSD,
                                                        std::vector<realT> & freq, 
                                                        realT m,
                                                        realT n,
                                                        int p,
                                                        realT fmax,
                                                        isParallel<false>  parallel )
{
   //Records each layer PSD
   std::vector<realT> single_PSD(freq.size());

   for(int i=0; i< _aosys->atm.n_layers(); ++i)
   {
      singleLayerPSD(single_PSD, freq, m, n, i, p, fmax);

      //Now add the single layer PSD to the overall PSD, weighted by Cn2
      for(int j=0; j<freq.size(); ++j)
      {
         PSD[j] += _aosys->atm.layer_Cn2(i)*single_PSD[j]; 
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
   for(int j=0;j<PSD.size(); ++j) PSD[j] = 0;

   if(fmax == 0) 
   {
      fmax = 150 + 2*fastestPeak(m, n);
   }
   
   return _multiLayerPSD( PSD, freq, m, n, p, fmax, isParallel<parallel>());
   

}
      
template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::makePSDGrid( std::string dir,
                                                      int mnMax,
                                                      realT dFreq,
                                                      realT maxFreq,
                                                      realT fmax )
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
   _aosys->dumpAOSystem(fout);
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
   mx::vectorScale(freq, N, dFreq, dFreq); //offset from 0 by dFreq, so f=0 not included
   
   fn = dir + '/' + "freq.binv";

   mx::writeBinVector( fn, freq);
   
   int nLoops = 0.5*spf.size();
   
   mx::ompLoopWatcher<> watcher(nLoops, std::cout);

   #pragma omp parallel
   {
      std::vector<realT> PSD;
      PSD.resize(N);
      std::string fname;
      
      int m, n;
      
      #pragma omp for 
      for(int i=0; i<nLoops; ++i)
      {
         m = spf[i*2].m;
         n = spf[i*2].n;
         
         multiLayerPSD<false>( PSD, freq, m, n, 1, fmax);
         
         fname = dir + '/' + "psd_" + mx::convertToString(m) + '_' + mx::convertToString(n) + ".binv";
      
         mx::writeBinVector( fname, PSD);
         
         watcher.incrementAndOutputStatus();
      }
   }
}
  
template<typename realT, typename aosysT>
int fourierTemporalPSD<realT, aosysT>::analyzePSDGrid( std::string subDir,
                                                       std::string psdDir,
                                                       int mnMax,
                                                       int mnCon,
                                                       int lpNc,
                                                       std::vector<realT> & mags,
                                                       std::vector<int> & intTimes
                                                     )
{

   std::string dir = psdDir + "/" + subDir;
   
   /*** Dump Params to file ***/
   mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   
   std::ofstream fout;
   std::string fn = dir + '/' + "params.txt";
   fout.open(fn);
   
   fout << "#---------------------------\n";
   _aosys->dumpAOSystem(fout);
   fout << "#---------------------------\n";
   fout << "# Analysis Parameters\n";
   fout << "#    mnMax    = " << mnMax << '\n';
   fout << "#    mnCon    = " << mnCon << '\n';
   fout << "#    lpNc     = " << lpNc << '\n';
   fout << "#    mags     = "; 
   for(int i=0; i<mags.size()-1; ++i) fout << mags[i] << ",";
   fout << mags[mags.size()-1] << '\n';
   fout << "#    intTimes = "; 
   for(int i=0; i<intTimes.size()-1; ++i) fout << intTimes[i] << ",";
   fout << intTimes[intTimes.size()-1] << '\n';
   
   fout.close();
   
   //**** Calculating A Variance Map ****//
   
   realT fs = 1.0/_aosys->minTauWFS();
   
   std::vector<sigproc::fourierModeDef> fms;
   
   sigproc::makeFourierModeFreqs_Rect(fms, 2*mnMax);
   int nModes = 0.5*fms.size();

   Eigen::Array<realT, -1, -1> gains, vars, gains_lp, vars_lp;
   
   gains.resize(2*mnMax+1, 2*mnMax+1);
   vars.resize(2*mnMax+1, 2*mnMax+1);
   
   gains(mnMax, mnMax) = 0;
   vars(mnMax, mnMax) = 0;

   gains_lp.resize(2*mnMax+1, 2*mnMax+1);
   vars_lp.resize(2*mnMax+1, 2*mnMax+1);
   
   gains_lp(mnMax, mnMax) = 0;
   vars_lp(mnMax, mnMax) = 0;
      
   bool doLP = false;
   if(lpNc > 1) doLP = true;
   Eigen::Array<realT, -1, -1> lpC;  
   
   if(doLP) 
   {
      lpC.resize(nModes, lpNc);
      lpC.setZero();
   }
   
   std::vector<realT> S_si, S_lp;
   
   mx::ompLoopWatcher<> watcher(nModes*mags.size()*intTimes.size(), std::cout);

   for(int j=0; j<intTimes.size(); ++j)
   {      
      if(intTimes[j] == 0) continue;
      
      S_si.clear();
      S_lp.clear();
      
      for(int s = 0; s < mags.size(); ++s)
      {
         //_aosys->starMag(mags[s]);
   
        #pragma omp parallel 
        { 
            realT localMag = mags[s];
            realT localIntTime = intTimes[j];
            
            realT var0;
            
            realT gopt, var;
      
            realT gopt_lp, var_lp;

            std::vector<realT> tfreq, tPSDp;
            std::vector<realT> tPSDn;      
      
            getGridPSD( tfreq, tPSDp, psdDir, 0, 1 ); //To get the freq grid
            
            int imax = 0;
            while( tfreq[imax] <= 0.5*fs/localIntTime && imax < tfreq.size()) ++imax;
            tfreq.erase(tfreq.begin() + imax, tfreq.end());
   
            tPSDn.resize(tfreq.size());
   
            mx::AO::analysis::clAOLinearPredictor<realT> tflp;
            
            mx::AO::analysis::clGainOpt<realT> go_si(localIntTime/fs, 1.5/fs);
            mx::AO::analysis::clGainOpt<realT> go_lp(localIntTime/fs, 1.5/fs);
            
            go_si.f(tfreq);
            go_lp.f(tfreq);
            
            realT gmax = 0;
            realT gmax_lp = 0;
            double t0, t1, t00, t11;
         
            int m, n;
            
            #pragma omp for schedule(dynamic, 5) //want to schedule dynamic with small chunks so maximal processor usage, otherwise we can end up with a small number of cores being used at the end
            for(int i=0; i<nModes; ++i)
            {
              
               //if( fms[i].p == -1 ) continue;
               m = fms[2*i].m;
               n = fms[2*i].n;
              
               wfsNoisePSD( tPSDn, (realT) _aosys->beta_p(m,n), _aosys->Fg(localMag), (realT) (localIntTime/fs), (realT) _aosys->npix_wfs(), (realT) _aosys->Fbg(), (realT) _aosys->ron_wfs());
               
               realT k = sqrt(m*m + n*n)/_aosys->D();
               
               getGridPSD( tPSDp, psdDir, m, n );
               tPSDp.erase(tPSDp.begin() + imax, tPSDp.end());
               var0 = _aosys->psd(_aosys->atm, k,0,1.0)*pow(_aosys->atm.lam_0()/_aosys->lam_wfs(),2) / pow(_aosys->D(),2); //
               
               if(fabs(m) <= mnCon && fabs(n) <= mnCon)
               {
                  gmax = 0;
                  gopt = go_si.optGainOpenLoop(var, tPSDp, tPSDn, gmax);
            
                  if(doLP)
                  {
                     tflp.regularizeCoefficients( gmax_lp, gopt_lp, var_lp, go_lp, tPSDp, tPSDn, lpNc);
                     for(int n=0; n< lpNc; ++n) lpC(i,n) = go_lp.a(n);
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
                  
               if(gopt > 0 && var > var0)
               {
                  gopt = 0;
                  var = var0;
               }

               if(gopt_lp > 0 && var_lp > var0)
               {
                  gopt_lp = 0;
                  var_lp = var0;
               }
         
               gains( mnMax + m, mnMax + n ) = gopt;
               gains( mnMax - m, mnMax - n ) = gopt;
      
               gains_lp( mnMax + m, mnMax + n ) = gopt_lp;
               gains_lp( mnMax - m, mnMax - n ) = gopt_lp;
         
               vars( mnMax + m, mnMax + n) = var;
               vars( mnMax - m, mnMax - n ) = var;
      
               vars_lp( mnMax + m, mnMax + n) = var_lp;
               vars_lp( mnMax - m, mnMax - n ) = var_lp;
               
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
   
         mx::improc::fitsFile<realT> ff;
         std::string fn = dir + "/gainmap_" + mx::convertToString<int>(mags[s]) + "_si_" + mx::convertToString<int>(intTimes[j]) + ".fits"; 
         ff.write( fn, gains);
   
         fn = dir + "/varmap_" + mx::convertToString<int>(mags[s]) + "_si_" + mx::convertToString<int>(intTimes[j]) + ".fits";
         ff.write( fn, vars);
   
         mx::AO::analysis::varmapToImage(cim, vars, psf);
         realT S = exp( -1*vars.sum());
         S_si.push_back(S);
         cim /= S;
         
         fn = dir + "/contrast_" + mx::convertToString<int>(mags[s]) + "_si_" + mx::convertToString<int>(intTimes[j]) + ".fits";
         ff.write( fn, cim);
         
         if(doLP)
         {
            fn = dir + "/gainmap_" + mx::convertToString<int>(mags[s]) + "_lp_" + mx::convertToString<int>(intTimes[j]) + ".fits"; 
            ff.write( fn, gains_lp);
   
            fn = dir + "/lpcmap_" + mx::convertToString<int>(mags[s]) + "_lp_" + mx::convertToString<int>(intTimes[j]) + ".fits"; 
            ff.write( fn, lpC);
            
            fn = dir + "/varmap_" + mx::convertToString<int>(mags[s]) + "_lp_" + mx::convertToString<int>(intTimes[j]) + ".fits";
            ff.write( fn, vars_lp);
         
            mx::AO::analysis::varmapToImage(cim, vars_lp, psf);
            S = exp( -1*vars_lp.sum());
            S_lp.push_back(S);
            
            cim /= S;
            fn = dir + "/contrast_" + mx::convertToString<int>(mags[s]) + "_lp_" + mx::convertToString<int>(intTimes[j]) + ".fits";
            ff.write( fn, cim);
         }
         
      }//s (mag)
      
      std::string fn = dir + "/strehl_si_" + mx::convertToString<int>(intTimes[j]) + ".txt";
      fout.open(fn);
      for(int i=0;i<mags.size(); ++i)
      {
         fout << mags[i] << " " << S_si[i] << "\n";
      }
      
      fout.close();
      
      if(doLP)
      {
         fn = dir + "/strehl_lp_" + mx::convertToString<int>(intTimes[j]) + ".txt";
         fout.open(fn);
         for(int i=0;i<mags.size(); ++i)
         {
            fout << mags[i] << " " << S_lp[i] << "\n";
         }
      
         fout.close();
      }
      
      
   }//j (intTimes)
   
   return 0;
}



template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::getGridFreq( std::vector<realT> & freq,
                                                      const std::string & dir )
{
   std::string fn;      
   fn = dir + '/' + "freq.binv";
   mx::readBinVector(freq, fn);
}

template<typename realT, typename aosysT>
void fourierTemporalPSD<realT, aosysT>::getGridPSD( std::vector<realT> & psd,
                                                     const std::string & dir,
                                                     int m,
                                                     int n )
{
   std::string fn;
   fn = dir + '/' + "psd_" + mx::convertToString(m) + '_' + mx::convertToString(n) + ".binv";
   mx::readBinVector(psd, fn);
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
   
   realT f = Fp->_f;
   realT v_wind = Fp->_aosys->atm.layer_v_wind(Fp->_layer_i); 
   
   //realT r_0 = Fp->_aosys->atm.r_0(); 
   //realT L_0 = Fp->_aosys->atm.L_0(); 
   realT D = Fp->_aosys->D(); 
   realT m = Fp->_m;
   realT n = Fp->_n; 
   int p = Fp->_p;

   realT ku = f/v_wind;
   
   realT kp = sqrt( pow(ku + m/D,2) + pow(kv + n/D,2) );
   realT kpp = sqrt( pow(ku - m/D,2) + pow(kv - n/D,2) );
   
   realT Q1 = math::func::jinc(pi<realT>()*D*kp);
   
   realT Q2 = math::func::jinc(pi<realT>()*D*kpp);
   
   realT Q = (Q1 + p*Q2);
   
   realT P =  Fp->_aosys->psd(Fp->_aosys->atm, sqrt( pow(ku,2) + pow(kv,2)), 0, Fp->_aosys->lam_wfs());
   
   return P*Q*Q ;
}

///Worker function for GSL Integration for the modified Fourier modes.
/** \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT F_mod (realT kv, void * params) 
{
   fourierTemporalPSD<realT, aosysT> * Fp = (fourierTemporalPSD<realT, aosysT> *) params;
   
   realT f = Fp->_f;
   realT v_wind = Fp->_aosys->atm.layer_v_wind(Fp->_layer_i); 
   
   //realT r_0 = Fp->_aosys->atm.r_0(); 
   //realT L_0 = Fp->_aosys->atm.L_0(); 
   realT D = Fp->_aosys->D(); 
   realT m = Fp->_m;
   realT n = Fp->_n; 
   int p = Fp->_p;

   realT ku = f/v_wind;
   
   realT kp = sqrt( pow(ku + m/D,2) + pow(kv + n/D,2) );
   realT kpp = sqrt( pow(ku - m/D,2) + pow(kv - n/D,2) );
   
   realT Jp = math::func::jinc(pi<realT>()*D*kp);
   
   realT Jm = math::func::jinc(pi<realT>()*D*kpp);
   
   realT QQ = 2*(Jp*Jp + Jm*Jm);
   
   realT P =  Fp->_aosys->psd(Fp->_aosys->atm, sqrt( pow(ku,2) + pow(kv,2)), Fp->_aosys->lam_sci(), 0, Fp->_aosys->lam_wfs(), Fp->_aosys->secZeta() );
   
   return P*QQ ;
}
 


///Worker function for GSL Integration for a basis projected onto the modified Fourier modes.
/** \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT F_projMod (realT kv, void * params) 
{
   fourierTemporalPSD<realT, aosysT> * Fp = (fourierTemporalPSD<realT, aosysT> *) params;
   
   realT f = Fp->_f;
   realT v_wind = Fp->_aosys->atm.layer_v_wind(Fp->_layer_i); 
   
   realT q_wind = Fp->_aosys->atm.layer_dir(Fp->_layer_i);
   
   //For rotating the basis
   realT cq = cos(q_wind);
   realT sq = sin(q_wind);

   
   //realT r_0 = Fp->_aosys->atm.r_0(); 
   //realT L_0 = Fp->_aosys->atm.L_0(); 
   realT D = Fp->_aosys->D(); 
   
   int mode_i = Fp->_mode_i;
   
   realT m; 
   realT n; 
   int p, pp;// = Fp->_p;

   realT ku = f/v_wind;
   
   realT kp, kpp, km, kmp, Jp, Jpp, Jm, Jmp, QQ;
   
   QQ = 0;
   
   for(int i=0; i< Fp->_modeCoeffs.cols(); ++i)
   {
      if(Fp->_modeCoeffs(mode_i,i) == 0) continue;
      
      m = Fp->ms[i]*cq + Fp->ns[i]*sq;
      n = -Fp->ms[i]*sq + Fp->ns[i]*cq;
   
      kp = sqrt( pow(ku + m/D,2) + pow(kv + n/D,2) );
      km = sqrt( pow(ku - m/D,2) + pow(kv - n/D,2) );
   
      Fp->Jps[i] = math::func::jinc(pi<realT>()*D*kp);
      Fp->Jms[i] = math::func::jinc(pi<realT>()*D*km);
      
   }
   
   int N = Fp->_modeCoeffs.cols();
   
   realT sc;
   
   for(int i=0; i< N ; ++i)
   {
      if( abs(Fp->_modeCoeffs(mode_i,i)) == 0) break; //continue;
      
      Jp = Fp->Jps[i];
      Jm = Fp->Jms[i];
      p = Fp->ps[i];
      
      QQ += 2*( Jp*Jp + Jm*Jm)*Fp->_modeCoeffs(mode_i,i)*Fp->_modeCoeffs(mode_i,i);
       
      for(int j=(i+1); j < N; ++j)
      {
         if( abs(Fp->_modeCoeffs(mode_i,j)) == 0) break; //continue;
         sc =  Fp->_modeCoeffs(mode_i,i)*Fp->_modeCoeffs(mode_i, j);
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
   
   
   realT P =  Fp->_aosys->psd(Fp->_aosys->atm, sqrt( pow(ku,2) + pow(kv,2)), Fp->_aosys->lam_sci(), 0, Fp->_aosys->lam_wfs(), Fp->_aosys->secZeta() );
   
   return P*QQ ;
}

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //fourierTemporalPSD_hpp
