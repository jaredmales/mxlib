/** \file zernikeTemporalPSD.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Calculation of the temporal PSD of Zernike modes.
  * \ingroup mxAOm_files
  *
  */

//***********************************************************************//
// Copyright 2022 Jared R. Males (jaredmales@gmail.com)
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

#ifndef zernikeTemporalPSD_hpp
#define zernikeTemporalPSD_hpp

#include <iostream>
#include <fstream>
#include <cmath>

#include <sys/stat.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Eigen/Dense>

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"
#include "../../math/func/airyPattern.hpp"
#include "../../math/vectorUtils.hpp"
#include "../../ioutils/fits/fitsFile.hpp"
#include "../../sigproc/zernike.hpp"
#include "../../ioutils/stringUtils.hpp"
#include "../../ioutils/readColumns.hpp"
#include "../../ioutils/binVector.hpp"
#include "../../ioutils/fileUtils.hpp"

#include "../../ipc/ompLoopWatcher.hpp"
#include "../../mxError.hpp"

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



//Forward declaration
template<typename realT, typename aosysT>
realT F_zernike (realT kv, void * params);


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
template<typename _realT, typename aosysT>
struct zernikeTemporalPSD
{
   typedef _realT realT;
   typedef std::complex<realT> complexT;
   
   ///Pointer to an AO system structure.
   aosysT * m_aosys {nullptr};

protected:
    realT m_apertureRatio {1}; /// Aperture ratio, < 1, used for segments.
    realT m_apertureRatio4 {1};   /// Aperture ratio ^4 used for calculations

public:
    void apertureRatio( realT ar )
    {
        m_apertureRatio = ar;
        m_apertureRatio4 = pow(ar,4);
    }

    realT apertureRatio4()
    {
        return m_apertureRatio4;
    }

   realT m_f {0}; ///< the current temporal frequency
   realT m_zern_j {0}; ///< the current mode number
   int m_zern_m {0}; ///< The current mode m
   int m_zern_n {0}; ///< The current mode n

   realT m_cq {0}; ///< The cosine of the wind direction
   realT m_sq {0}; ///< The sine of the wind direction
   realT m_spatialFilter  {false}; ///< Flag indicating if a spatial filter is applied
   
   int _layer_i; ///< The index of the current layer.

   ///Worskspace for the gsl integrators, allocated to WSZ if constructed as worker (with allocate == true).
   gsl_integration_workspace * _w;

   realT _absTol; ///< The absolute tolerance to use in the GSL integrator
   realT _relTol; ///< The relative tolerance to use in the GSL integrator



  std::vector<realT> Jps;
  std::vector<realT> Jms;
  std::vector<int> ps;
  std::vector<realT> ms;
  std::vector<realT> ns;


public:
   ///Default c'tor
   zernikeTemporalPSD();

   ///Constructor with workspace allocation
   /**
     * \param allocate if true, then the workspace for GSL integrators is allocated.
     */
   explicit zernikeTemporalPSD(bool allocate);

   ///Destructor
   /** Frees GSL workspace if it was allocated.
     */
   ~zernikeTemporalPSD();

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


   ///Calculate the temporal PSD for a Zernike mode for a single layer.
   /**
     *
     * \todo implement error checking.
     * \todo need a way to track convergence failures in integral without throwing an error.
     * \todo need better handling of averaging for the -17/3 extension.
     *
     */
   int singleLayerPSD( std::vector<realT> &PSD,  ///< [out] the calculated PSD
                       std::vector<realT> &freq, ///< [in] the populated temporal frequency grid defining the frequencies at which the PSD is calculated
                       int zern_j,                  ///< [in] 
                       int layer_i,              ///< [in] the index of the layer, for accessing the atmosphere parameters
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
                       int zern_j,
                       realT fmax,
                       isParallel<true> parallel );

   //Non-Parallelized version of multiLayerPSD, without OMP directives.
   int m_multiLayerPSD( std::vector<realT> & PSD,
                       std::vector<realT> & freq,
                       int zern_j,
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
                      int zern_j,
                      realT fmax = 0             ///< [in] [optional] set the maximum temporal frequency for the calculation. The PSD is filled in
                                                 /// with a -17/3 power law past this frequency.  If 0, then it is taken to be 150 Hz + 2*fastestPeak(m,n).
                    );


};

template<typename realT, typename aosysT>
zernikeTemporalPSD<realT, aosysT>::zernikeTemporalPSD()
{
   m_aosys = nullptr;
   initialize();
}

template<typename realT, typename aosysT>
zernikeTemporalPSD<realT, aosysT>::zernikeTemporalPSD(bool allocate)
{
   m_aosys = nullptr;
   initialize();

   if(allocate)
   {
      _w = gsl_integration_workspace_alloc (WSZ);
   }
}

template<typename realT, typename aosysT>
zernikeTemporalPSD<realT, aosysT>::~zernikeTemporalPSD()
{
   if(_w)
   {
      gsl_integration_workspace_free (_w);
   }
}

template<typename realT, typename aosysT>
void zernikeTemporalPSD<realT, aosysT>::initialize()
{
   _w = 0;

   _absTol = 1e-10;
   _relTol = 1e-4;
}

template<typename realT, typename aosysT>
void zernikeTemporalPSD<realT, aosysT>::absTol(realT at)
{
   _absTol = at;
}

template<typename realT, typename aosysT>
realT zernikeTemporalPSD<realT, aosysT>::absTol()
{
   return _absTol;
}

template<typename realT, typename aosysT>
void zernikeTemporalPSD<realT, aosysT>::relTol(realT rt)
{
   _relTol = rt;
}

template<typename realT, typename aosysT>
realT zernikeTemporalPSD<realT, aosysT>::relTol()
{
   return _relTol;
}


template<typename realT, typename aosysT>
int zernikeTemporalPSD<realT, aosysT>::singleLayerPSD( std::vector<realT> &PSD,
                                                       std::vector<realT> &freq,
                                                       int zern_j,
                                                       int layer_i,
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
   zernikeTemporalPSD<realT, aosysT> params(true);

   params.m_aosys = m_aosys;
   params.m_apertureRatio = m_apertureRatio;
   params.m_apertureRatio4 = m_apertureRatio4;
   params._layer_i = layer_i;
   params.m_zern_j = zern_j;
   sigproc::noll_nm(params.m_zern_n, params.m_zern_m, params.m_zern_j);
   params.m_cq = cq; //for de-rotating ku and kv for spatial filtering
   params.m_sq = sq; //for de-rotation ku and kv for spatial filtering
   if(m_aosys->spatialFilter_ku() < std::numeric_limits<realT>::max() || m_aosys->spatialFilter_kv() < std::numeric_limits<realT>::max()) params.m_spatialFilter = true; 
   


   realT result, error;

   //Setup the GSL calculation
   gsl_function func;
   func.function = &F_zernike<realT, aosysT>;


   func.params = &params;


   //Here we only calculate up to fmax.
   size_t i=0;
   while( freq[i] <= fmax )
   {
      params.m_f = freq[i];

      int ec = gsl_integration_qagi (&func, _absTol, _relTol, WSZ, params._w, &result, &error);

      if(ec == GSL_EDIVERGE)
      {
         std::cerr << "GSL_EDIVERGE:" << " " << freq[i] << " " << v_wind << " " << zern_j << "\n";
         std::cerr << "ignoring . . .\n";
      }

      PSD[i] = scale*result;

      ++i;
      if(i >= freq.size()) break;
   }
/*
   //Now fill in from fmax to the actual max frequency with a -(alpha+2) power law.
   size_t j=i;
   
   if(j == freq.size()) return 0;
   
   //First average result for last 50.
   PSD[j] = PSD[i-50] * pow( freq[i-50]/freq[j], m_aosys->atm.alpha(layer_i)+2);//seventeen_thirds<realT>());
   for(size_t k=49; k> 0; --k)
   {
      PSD[j] +=  PSD[i-k] * pow( freq[i-k]/freq[j], m_aosys->atm.alpha(layer_i)+2); //seventeen_thirds<realT>());
   }
   PSD[j] /= 50.0;
   ++j;
   ++i;
   if(j == freq.size()) return 0;
   while(j < freq.size())
   {
      PSD[j] = PSD[i-1] * pow( freq[i-1]/freq[j], m_aosys->atm.alpha(layer_i)+2); //seventeen_thirds<realT>());
      ++j;
   }
*/

   return 0;
}


template<typename realT, typename aosysT>
int zernikeTemporalPSD<realT, aosysT>::m_multiLayerPSD( std::vector<realT> & PSD,
                                                       std::vector<realT> & freq,
                                                       int zern_j,
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
         singleLayerPSD(single_PSD, freq, zern_j, i, fmax);

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
int zernikeTemporalPSD<realT, aosysT>::m_multiLayerPSD( std::vector<realT> & PSD,
                                                        std::vector<realT> & freq,
                                                        int zern_j,
                                                        realT fmax,
                                                        isParallel<false>  parallel )
{
   static_cast<void>(parallel);
   
   //Records each layer PSD
   std::vector<realT> single_PSD(freq.size());

   for(size_t i=0; i< m_aosys->atm.n_layers(); ++i)
   {
      singleLayerPSD(single_PSD, freq, zern_j, i, fmax);

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
int zernikeTemporalPSD<realT, aosysT>::multiLayerPSD( std::vector<realT> & PSD,
                                                       std::vector<realT> & freq,
                                                       int zern_j,
                                                       realT fmax )
{
   //PSD is zeroed every time to make sure we don't accumulate on repeated calls
   for(size_t j=0;j<PSD.size(); ++j) PSD[j] = 0;

   /*if(fmax == 0)
   {
      fmax = 150 + 2*fastestPeak(m, n);
   }*/
  
   return m_multiLayerPSD( PSD, freq, zern_j, fmax, isParallel<parallel>());


}


                             


///Worker function for GSL Integration for the Zernike modes.
/** \ingroup mxAOAnalytic
  */
template<typename realT, typename aosysT>
realT F_zernike (realT kv, void * params)
{
    zernikeTemporalPSD<realT, aosysT> * Fp = (zernikeTemporalPSD<realT, aosysT> *) params;

    realT f = Fp->m_f;
    realT v_wind = Fp->m_aosys->atm.layer_v_wind(Fp->_layer_i);
    realT D = Fp->m_aosys->D();

    realT apertureRatio4 = Fp->apertureRatio4();

    realT ku = f/v_wind;

    realT phi = atan(kv/ku);

    realT k = sqrt( pow(ku,2) + pow(kv,2) );
  
    realT Q2norm = apertureRatio4 * sigproc::zernikeQNorm(k*D/2.0, phi, Fp->m_zern_n, Fp->m_zern_m);

    realT P =  Fp->m_aosys->psd(Fp->m_aosys->atm, Fp->_layer_i, k, Fp->m_aosys->lam_sci(), Fp->m_aosys->lam_wfs(), Fp->m_aosys->secZeta() );


    return P*Q2norm ;
}




/*extern template
struct zernikeTemporalPSD<float, aoSystem<float, vonKarmanSpectrum<float>, std::ostream>>;*/


//extern template
//struct zernikeTemporalPSD<double, aoSystem<double, vonKarmanSpectrum<double>, std::ostream>>;

/*
extern template
struct zernikeTemporalPSD<long double, aoSystem<long double, vonKarmanSpectrum<long double>, std::ostream>>;

#ifdef HASQUAD
extern template
struct zernikeTemporalPSD<__float128, aoSystem<__float128, vonKarmanSpectrum<__float128>, std::ostream>>;
#endif
*/

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //zernikeTemporalPSD_hpp
