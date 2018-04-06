/** \file filterNoise1D.hpp
  * \brief Filtering noise in 1 dimension
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
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

#ifndef filterNoise1D_hpp
#define filterNoise1D_hpp

#include "../randomT.hpp"
#include "psdFilter.hpp"

namespace mx
{
namespace sigproc 
{

///Class to manage filtering of 1-D noise time-series with a PSD
template<typename _realT>
struct filterNoise1D
{
   typedef _realT realT;
   
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> realArrayT;
   
   mx::normDistT<realT> m_normVar;

   bool m_seeded {false};
   
   mx::sigproc::psdFilter<realT> m_psdFilt;
   
   size_t m_nsamps {0};
   
   int m_oversamp {0};
   
   realT m_scaleFactor{1};
   
   realArrayT m_noise;
   
   filterNoise1D();
   
   filterNoise1D( int oversamp,
                  realArrayT * psdSqrt
                );
   
   int initialize( int oversamp,
                   realArrayT * psdSqrt
                 );
   
   static int makeSqrtVonKarman( realArrayT & sqrtPD,
                                 size_t nsamps,
                                 int oversamp,
                                 realT dt,
                                 realT alpha,
                                 realT T0 = 0,
                                 realT t0 = 0
                               );
   
   
   int initNoise();
   
   int genNoise( realArrayT & cnoise );
   
   realT scaleFactor();
   
   int scaleFactor( realT sf );
   
   realT measureScaleFactor( realT expectVar,
                             int nTrials
                           );
};

template<typename realT>
filterNoise1D<realT>::filterNoise1D()
{
   return;
}

template<typename realT>
filterNoise1D<realT>::filterNoise1D( int oversamp,
                                         realArrayT * psdSqrt
                                       )
{
   initialize(oversamp, psdSqrt);
}

template<typename realT>
int filterNoise1D<realT>::initialize( int oversamp,
                                        realArrayT * psdSqrt
                                      )
{
   m_oversamp = oversamp;
 
   m_nsamps = psdSqrt->rows() / oversamp;
   
   m_psdFilt.psdSqrt( psdSqrt );
   
   if(!m_seeded)
   {
      m_normVar.seed();
      m_seeded = true;
   }
   
   return 0;
}

template<typename realT>
int filterNoise1D<realT>::makeSqrtVonKarman( realArrayT & sqrtPSD,
                                             size_t nsamps,
                                             int oversamp,
                                             realT dt,
                                             realT alpha,
                                             realT T0,
                                             realT t0
                                           )
{
   realArrayT freq;
   realArrayT PSD;
   
   freq.resize( nsamps*oversamp, 1);

   mx::sigproc::frequency_grid1D( freq, dt );
   
   sqrtPSD.resize(freq.rows(), freq.cols());
   
   mx::sigproc::vonKarman_psd( sqrtPSD, freq, alpha,  T0, t0, 1);   
   
   realT df = freq(1,0) - freq(0,0);
   
   for(int r=0; r < sqrtPSD.rows(); ++r) sqrtPSD(r,0) = sqrt( sqrtPSD(r,0) * sqrtPSD.rows()*df);
   
   return 0;
}

template<typename realT>
int filterNoise1D<realT>::initNoise()
{
   m_noise.resize( m_psdFilt.rows(), 1); //a no-op if not a change.
   
   for(int ii = 0; ii < m_noise.rows(); ++ii) m_noise(ii,0) = m_normVar;
   
   return 0;
}

template<typename realT>
int filterNoise1D<realT>::genNoise( realArrayT & cnoise )
{
   initNoise();
   
   m_psdFilt( m_noise );
      
   size_t start = 0.5*(m_psdFilt.rows()-1) - 0.5*m_nsamps;
   
   cnoise = m_noise.block( start, 0, m_nsamps, 1)*m_scaleFactor;
   
   return 0;
}

template<typename realT>
realT filterNoise1D<realT>::scaleFactor()
{
   return m_scaleFactor;
}
   
template<typename realT>
int filterNoise1D<realT>::scaleFactor( realT sf )
{
   m_scaleFactor = sf;
   return 0;
}
   
template<typename realT>
realT filterNoise1D<realT>::measureScaleFactor( realT expectVar,
                                                int nTrials
                                              )
{
   realArrayT cnoise;
   
   realT dv = 0;
   realT n = 0;
   for(int i=0; i< nTrials; ++i)
   {
      genNoise(cnoise);
   
      dv += sqrt(cnoise.square().sum()/cnoise.rows());
      ++n;
   }
   
   m_scaleFactor = expectVar / (dv/n);
   return m_scaleFactor;
}

} //namespace sigproc 
} //namespace mx
#endif //filterNoise1D_hpp
