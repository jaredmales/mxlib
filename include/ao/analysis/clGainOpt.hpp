/** \file clGainOpt.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to manage closed loop gain optimization.
  * \ingroup mxAO_files
  *
  */

//***********************************************************************//
// Copyright 2016-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef clGainOpt_hpp
#define clGainOpt_hpp

#ifdef MX_INCLUDE_BOOST
#include <boost/math/tools/minima.hpp>
#endif

#include <type_traits>

#include <Eigen/Dense>

#include "../../sys/timeUtils.hpp"

#include "../../math/constants.hpp"

//#define ALLOW_F_ZERO

namespace mx
{
namespace AO
{
namespace analysis
{

//forward declaration of worker functor
template<typename realT>
struct clGainOptOptGain_OL;


/// A class to manage optimizing closed-loop gains
/**
  * \tparam _realT the real floating point type in which to do all arithmetic.  Is used to define the complex type as well.
  *
  * \ingroup mxAOAnalytic
  */
template<typename _realT>
struct clGainOpt
{
   typedef _realT realT; ///< The real data type
   typedef std::complex<_realT> complexT; ///< The complex data type

protected:

   int _N; ///< Number of integrations in the (optional) moving average.  Default is 1.
   realT _Ti; ///< The loop sampling interval
   realT _tau; ///< The loop delay

   std::vector<realT> _b; ///< Vector of FIR coefficients
   std::vector<realT> _a; ///< Vector of IIR coefficients

   std::vector<realT> _f; ///< Vector of frequencies

   bool _fChanged; ///< True if frequency or max size of _a and _b changes

   bool _changed; ///< True if any of the members which make up the basic transfer functions are changed

   Eigen::Array<realT, -1, -1> _cs;
   Eigen::Array<realT, -1, -1> _ss;

   std::vector<std::complex<realT> > _H_dm, _H_wfs, _H_ma, _H_del, _H_con;

public:
   /** Parameters for stability analysis
     * @{
     */

   realT _maxFindMin; ///< The Minimum value for the maximum stable gain finding algorithm.

   ///@}

   /** Parameters for minimization finding
     * @{
     */

   realT _minFindMin; ///< The Minimum value for the minimum finding algorithm.
   realT _minFindMaxFact; ///< The maximum value, as a multiplicative factor of maximum gain
   int _minFindBits; ///< The bits of precision to use for minimum finding. Defaults to std::numeric_limits<realT>::digits.
   uintmax_t _minFindMaxIter; ///< The maximum iterations allowed for minimization.

   ///@}


   /// Default c'tor.
   clGainOpt();

   /// C'tor setting the loop timings.
   /**
     */
   clGainOpt( realT Ti,  ///< [in] the desired loop sampling interval.
              realT tau ///< [in] the desired loop delay.
            );

   /// Initialize this instance.
   void init();

   /// Get the number of integrations in the (optional) moving average
   /**
     * \returns the current value of _N.
     */
   int N();

   /// Set the number of integrations in the moving average
   /**
     */
   void N( int newN /**< [in] the value of _N. */);

   /// Get the loop sampling interval
   /**
     * \returns the current value of _Ti.
     */
   realT Ti();

   /// Set the loop sampling interval
   /**
     */
   void Ti( realT newTi /**< [in]  the new value of _Ti. */);

   /// Get the loop delay
   /**
     * \returns the current value of _tau.
     */
   realT tau();

   /// Set the loop delay
   /**
     */
   void tau( realT newTau /**< [in] the new value of _tau.*/);

   /// Set the vector of FIR coefficients
   /**
     */
   void b( const std::vector<realT> & newB  /**< [in] a vector of coefficients, which is copied to _b.*/);

   /// Set the vector of FIR coefficients
   /**
     */
   void b( const Eigen::Array<realT, -1, -1> & newB  /**< [in] a column-vector Eigen::Array of coefficients, which is copied to _b.*/);

   /// Set the vector of IIR coefficients
   /**
     */
   void a ( const std::vector<realT> & newA  /**< [in] a vector of coefficients, which is copied to _a. */);

   /// Set the vector of IIR coefficients
   /**
     */
   void a( const Eigen::Array<realT, -1, -1> & newA  /**< [in] a column-vector Eigen::Array of coefficients, which is copied to _a.*/);

   realT a( size_t i) { return _a[i];}

   /// Set the FIR and IIR coefficients so that the control law is a leaky integrator.
   /** Set remember to 1.0 for a pure integrator control law.
     *
     */
   void setLeakyIntegrator( realT remember  /**< [in] a number usually close to 1 setting the amount "remembered" from previous iterations.*/);

   /// Set the vector of frequencies
   /**
     */
   void f ( realT * newF, ///< [in] a pointer to an array containing the new frequencies
            size_t nF     ///< [in] the number of elements of size(realT) in newF.
          );

   /// Set the vector of frequencies
   /**
     */
   void f ( std::vector<realT> & newF  /**< [in] a vector containing the new frequencies */);

   ///Get the size of the frequency vector
   /**
     * \returns _f.size()
     */
   size_t f_size()
   {
      return _f.size();
   }

   /// Get the i-th value of frequency.
   /** No range checks are conducted.
     *
     * \returns the value of _f[i]
     *
     */
   realT f( size_t i /**< [in] the index of the frequency to return*/ );

   /// Calculate the open-loop transfer function
   /**
     * \return the complex value of the open-loop transfer function at f.
     */
   complexT olXfer( int fi,           ///< [in] the index of the frequency
                    complexT & H_dm,  ///< [out] the transfer function of the DM
                    complexT & H_del, ///< [out] the delay transfer function
                    complexT & H_con  ///< [out] the controller transfer function.
                  );

   /// Calculate the open-loop transfer function
   /**
     * \overload
     *
     * \returns the complex value of the open-loop transfer function at f[fi].
     */
   complexT olXfer( int fi  /**< [in] the index of the frequency */);


   /// Return the closed loop error transfer function (ETF) at frequency f for gain g.
   /**
     * \returns the closed loop ETF at f and g.
     */
   complexT clETF( int fi,  ///< [in] the index of the frequency at which to calculate the ETF
                   realT g  ///< [in] the loop gain.
                 );

   /// Return the closed loop error transfer function (ETF) phase at frequency f for gain g.
   /**
     * \returns the phase of the closed loop ETF at f and g.
     */
   realT clETFPhase( int fi,  ///< [in] the index of the frequency at which to calculate the ETF
                     realT g  ///< [in] the loop gain.
                   );
   
   /// Return the norm of the closed loop error transfer function (ETF) at frequency f for gain g.
   /**
     * \returns the norm of the closed loop ETF at f and g.
     */
   realT clETF2( int fi, ///< [in] the index of the frequency at which to calculate the ETF.
                 realT g ///< [in] the loop gain.
               );

   /// Return the closed loop noise transfer function (NTF) at frequency f for gain g.
   /**
     * \returns the closed loop NTF at f and g.
     */
   complexT clNTF( int fi,  ///< [in] the index of the frequency at which to calculate the NTF
                   realT g  ///< [in] the loop gain.
                 );

   /// Return the norm of the closed loop noise transfer function (NTF) at frequency f for gain g.
   /**
     * \returns the value of the closed loop NTF at f and g.
     */
   realT clNTF2( int fi, ///< [in] the index of the frequency at which to calculate the NTF
                 realT g ///< [in] the loop gain.
               );

   /// Return the norm of the closed loop transfer functions at frequency f for gain g.
   /** Calculates both the error transfer function (ETF) and the noise transfer function (NTF).
     * This minimizes the various complex number operations, compared to calling both clETF2 and clNTF2.
     *
     */
   void clTF2( realT & ETF, ///< [out] is set to the ETF at f and g
               realT & NTF, ///< [out] is set to the NTF at f and g
               int fi,      ///< [in] is the index of the frequency at which to calculate the ETF
               realT g      ///< [in] is the loop gain.
             );


   /// Calculate the closed loop variance given open-loop PSDs and gain
   /** Calculates the following quantities.
     \f[
     \sigma_{err}^2 = \sum_i \left| ETF(f_i) \right|^2 PSD_{err}(fi) \Delta f\\
     \sigma_{noise}^2 = \sum_i \left| NTF(f_i) \right|^2 PSD_{noise}(fi) \Delta f\\
     \sigma^2 = \sigma_{err}^2 + \sigma_{noise}^2
     \f]
     * \f$ \sigma^2 \f$ is returned, and \f$ \sigma_{err}^2 \f$ and \f$ \sigma_{noise}^2 \f$ are available as the optional
     * arguments varErr and varNoise.
     *
     * \returns the total variance (error + noise) in closed loop
     */
   realT clVariance( realT & varErr,                ///< [out] the variance in the residual process error.
                     realT & varNoise,              ///< [out] the variance in the residual measurement noise.
                     std::vector<realT> & PSDerr,   ///< [in] the open-loop process error PSD.
                     std::vector<realT> & PSDnoise, ///< [in] the open-loop measurement noise PSD.
                     realT g                        ///< [in] the gain.
                   );

   /// Calculate the closed loop variance given open-loop PSDs and gain
   /** Overload of clVariance without the varErr and varNoise output parameters.
     *
     * \overload
     *
     * \returns the total variance (error + noise) in closed loop
     */
   realT clVariance( std::vector<realT> & PSDerr,   ///< [in] the open-loop process error PSD.
                     std::vector<realT> & PSDnoise, ///< [in] the open-loop measurement noise PSD.
                     realT g                        ///< [in] the gain.
                   );

   /// Find the maximum stable gain for the loop parameters
   /**
     *
     * \returns the maximum stable gain for the loop parameters
     */


   /// Find the maximum stable gain for the loop parameters
   /** Conducts a search along the Nyquist contour of the open-loop transfer function to find
     * the most-negative crossing of the real axis.
     *
     * \returns the maximum stable gain for the loop parameters
     */
   realT maxStableGain( realT & ll, ///< [in.out] the lower limit used for the search
                        realT & ul ///< [in.out] the upper limit used for hte search
                      );

   /// Find the maximum stable gain for the loop parameters
   /** Conducts a search along the Nyquist contour of the open-loop transfer function to find
     * the most-negative crossing of the real axis.
     *
     * This version allows constant arguments.
     * \overload
     *
     * \returns the maximum stable gain for the loop parameters
     */
   realT maxStableGain( const realT & ll, ///< [in] the lower limit used for the search
                        const realT & ul ///< [in] the upper limit used for hte search
                      );

   /// Find the maximum stable gain for the loop parameters
   /** Conducts a search along the Nyquist contour of the open-loop transfer function to find
     * the most-negative crossing of the real axis.
     *
     * This version uses _maxFindMin for the lower limit and no upper limit.
     *
     * \overload
     *
     * \returns the maximum stable gain for the loop parameters
     */
   realT maxStableGain();

   ///Return the optimum closed loop gain given an open loop PSD
   /** Uses _gmax for the upper limit.
     * \returns the optimum gain
     */
   realT optGainOpenLoop( realT & var,                    ///< [out] the variance at the optimum gain
                          std::vector<realT> & PSDerr,    ///< [in] open loop error PSD
                          std::vector<realT> & PSDnoise   ///< [in] open loop measurement noise PSD
                        );

   ///Return the optimum closed loop gain given an open loop PSD
   /**
     * \returns the optimum gain.
     */
   realT optGainOpenLoop( realT & var,                   ///< [out] the variance at the optimum gain
                          std::vector<realT> & PSDerr,   ///< [in] open loop error PSD
                          std::vector<realT> & PSDnoise, ///< [in] open loop measurement noise PSD
                          realT & gmax                   ///< [in] maximum gain to consider.  If 0, then _gmax is used.
                        );
   
   ///Calculate the pseudo open-loop PSD given a closed loop PSD
   /**
     * \returns 0 on success
     */
   int pseudoOpenLoop( std::vector<realT> & PSD, ///< [in.out] input closed loop PSD, on output contains the pseudo open loop error PSD ,
                       realT g                   ///< [in] the loop gain when PSD was measured.
                     );

   int nyquist( std::vector<realT> & re,
                std::vector<realT> & im,
                realT g
              );
};

template<typename realT>
clGainOpt<realT>::clGainOpt()
{
   init();
}

template<typename realT>
clGainOpt<realT>::clGainOpt(realT Ti, realT tau)
{
   init();

   _Ti = Ti;
   _tau = tau;
}

template<typename realT>
void clGainOpt<realT>::init()
{
   _N = 1;

   setLeakyIntegrator(1.0);

   _Ti = 1./1000.;
   _tau = 2.5*_Ti;

   _maxFindMin = 0.0;

   _minFindMin = 1e-9;
   _minFindMaxFact = 0.999;
   _minFindBits = std::numeric_limits<realT>::digits;
   _minFindMaxIter = 1000;

   _fChanged = true;
   _changed = true;

}

template<typename realT>
int clGainOpt<realT>::N()
{
   return _N;
}

template<typename realT>
void clGainOpt<realT>::N( int newN )
{
   _N = newN;
   _changed = true;
}

template<typename realT>
realT clGainOpt<realT>::Ti()
{
   return _Ti;
}

template<typename realT>
void clGainOpt<realT>::Ti( realT newTi)
{
   _Ti = newTi;
   _changed = true;
}

template<typename realT>
realT clGainOpt<realT>::tau()
{
   return _tau;
}

template<typename realT>
void clGainOpt<realT>::tau(realT newTau)
{
   _tau = newTau;
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::b( const std::vector<realT> & newB)
{
   if( newB.size() > (size_t) _cs.cols() )
   {
      _fChanged = true;
   }

   _b = newB;
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::b( const Eigen::Array<realT, -1, -1>  & newB)
{
   if( newB.cols() > _cs.cols() )
   {
      _fChanged = true;
   }

   _b.resize( newB.cols());

   for(size_t i=0; i< _b.size(); ++i)
   {
      _b[i] = newB(0,i);
   }

   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::a ( const std::vector<realT> & newA)
{
   if( newA.size()+1 > (size_t) _cs.cols())
   {
      _fChanged = true;
   }

   _a = newA;
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::a( const Eigen::Array<realT, -1, -1>  & newA)
{
   if( newA.cols() +1 > _cs.cols())
   {
      _fChanged = true;
   }

   _a.resize( newA.cols());

   for(size_t i=0; i< _a.size(); ++i)
   {
      _a[i] = newA(0,i);
   }

   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::setLeakyIntegrator(realT remember)
{
   _b.resize(1);
   _b[0] = 1.0;

   _a.resize(1);
   _a[0] = remember;

   _fChanged = true;
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::f(realT * newF, size_t nF)
{
   _f.resize(nF);
   for(int i=0; i< nF; ++i) _f[i] = newF[i];

   _fChanged = true;
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::f(std::vector<realT> & newF)
{
   _f = newF;
   _fChanged = true;

   _changed = true;
}

template<typename realT>
realT clGainOpt<realT>::f(size_t i)
{

   return _f[i];
}

template<typename realT>
std::complex<realT> clGainOpt<realT>::olXfer(int fi)
{
   complexT H_dm;
   complexT H_del;
   complexT H_con;

   return olXfer(fi, H_dm, H_del, H_con);
}

//If PRECALC_TRIG is defined, then the cosine and sine tables are pre-calculated and used instead of repeated exp(-i) calls.
//This is much much faster, though uses more memory.  In general, only undefine this for testing or debugging.
#define PRECALC_TRIG

template<typename realT>
std::complex<realT> clGainOpt<realT>::olXfer(int fi, complexT & H_dm, complexT & H_del, complexT & H_con)
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0)
   {
   #else
   if(_f[fi] < 0)
   {
   #endif

      H_dm = 0;
      H_del = 0;
      H_con = 0;
      return 0;
   }

   #ifdef PRECALC_TRIG
   if(_fChanged)
   {
      size_t jmax = std::max(_a.size()+1, _b.size());

      _cs.resize(_f.size(), jmax);
      _ss.resize(_f.size(), jmax);

      for(size_t i=0; i< _f.size(); ++i)
      {
         _cs(i,0) = 1.0;
         _ss(i,0) = 0.0;

         for(size_t j = 1; j<jmax; ++j)
         {
            _cs(i,j) = cos(math::two_pi<realT>()*_f[i]*_Ti*realT(j));
            _ss(i,j) = sin(math::two_pi<realT>()*_f[i]*_Ti*realT(j));
         }
      }

      _fChanged = false;
   }
   #endif

   if(_changed)
   {
      _H_dm.resize(_f.size(), 0);
      _H_wfs.resize(_f.size(), 0);
      _H_ma.resize(_f.size(), 0);
      _H_del.resize(_f.size(), 0);
      _H_con.resize(_f.size(), 0);

      size_t jmax = std::min(_a.size(), _b.size());

      //#pragma omp parallel for
      for(size_t i=0; i<_f.size(); ++i)
      {
         #ifndef ALLOW_F_ZERO
         if(_f[i] <= 0) continue;
         #else
         if(_f[i] < 0) continue;
         #endif

         complexT s = complexT(0.0, math::two_pi<realT>()*_f[i]);

         complexT expsT = exp(-s*_Ti);

         if(_f[i] == 0)
         {
            _H_dm[i] = std::complex<realT>(1,0);
         }
         else
         {
            _H_dm[i] = (realT(1) - expsT)/(s*_Ti);
         }
         
         _H_wfs[i] = _H_dm[i];

         _H_ma[i] = 1;  //realT(1./_N)*(realT(1) - pow(expsT,_N))/(realT(1) - expsT);

         _H_del[i] = exp(-s*_tau);

         complexT FIR = complexT(_b[0],0);

         complexT IIR = complexT(0.0, 0.0);
         for(size_t j = 1; j < jmax; ++j)
         {
            #ifdef PRECALC_TRIG
            realT cs = _cs(i,j);
            realT ss = _ss(i,j);
            FIR += _b[j]*complexT(cs, -ss);
            IIR += _a[j-1]*complexT(cs, -ss);
            #else
            complexT expZ = exp( -s*_Ti*realT(j));
            FIR += _b[j]*expZ;
            IIR += _a[j-1]*expZ;
            #endif
         }

         for(size_t jj=jmax; jj< _a.size()+1; ++jj)
         {
            #ifdef PRECALC_TRIG
            realT cs = _cs(i,jj);
            realT ss = _ss(i,jj);
            IIR += _a[jj-1]*complexT(cs, -ss);
            #else
            complexT expZ = exp( -s*_Ti*realT(jj));
            IIR += _a[jj-1]*expZ;
            #endif
         }

         for(size_t jj=jmax; jj<_b.size(); ++jj)
         {
            #ifdef PRECALC_TRIG
            realT cs = _cs(i,jj);
            realT ss = _ss(i,jj);
            FIR += _b[jj]*complexT(cs, -ss);
            #else
            complexT expZ = exp( -s*_Ti*realT(jj));
            FIR += _b[jj]*expZ;
            #endif
         }

         _H_con[i] = FIR/( realT(1.0) - IIR);

         /*if( i == 0 || i == 1)
         {
            std::cerr << i << " " << _f[fi] << " " << s << " " << expsT << " " << _H_wfs[i] << " " << _H_dm[i] << " " << _H_con[i] << " " << _H_del[i] << "\n";
            //exit(0);
         }/**/

      }

      _changed = false;
   }

   H_dm = _H_dm[fi];
   H_del = _H_del[fi];//*_H_ma[fi];
   H_con = _H_con[fi];

   return ( _H_dm[fi]*_H_wfs[fi]*_H_del[fi]*_H_con[fi]);

}

template<typename realT>
std::complex<realT> clGainOpt<realT>::clETF(int fi, realT g)
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0) return 0;
   #else
   if(_f[fi] < 0) return 0;
   #endif

   return (realT(1)/( realT(1) + g*olXfer(fi)));
}

template<typename realT>
realT clGainOpt<realT>::clETFPhase(int fi, realT g)
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0) return 0;
   #else
   if(_f[fi] < 0) return 0;
   #endif

   return std::arg((realT(1)/( realT(1) + g*olXfer(fi))));
}

template<typename realT>
realT clGainOpt<realT>::clETF2(int fi, realT g)
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0) return 0;
   #else
   if(_f[fi] < 0) return 0;
   #endif

   return norm(realT(1)/( realT(1) + g*olXfer(fi)));
}

template<typename realT>
std::complex<realT> clGainOpt<realT>::clNTF( int fi, 
                                             realT g
                                           )
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0) return 0;
   #else
   if(_f[fi] < 0) return 0;
   #endif

   complexT H_dm, H_del, H_con;

   complexT olX =  olXfer(fi, H_dm, H_del, H_con); //H_dm*H_wfs*H_ma*H_del*H_con;

   return -(H_dm*H_del*g*H_con)/(realT(1) + g*olX);

}

template<typename realT>
realT clGainOpt<realT>::clNTF2( int fi, 
                                realT g
                              )
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0) return 0;
   #else
   if(_f[fi] < 0) return 0;
   #endif

   complexT H_dm, H_del, H_con;

   complexT olX =  olXfer(fi, H_dm, H_del, H_con); //H_dm*H_wfs*H_ma*H_del*H_con;

   complexT NTF = -(H_dm*H_del*g*H_con)/(realT(1) + g*olX);

   return norm(NTF);
}

template<typename realT>
void clGainOpt<realT>::clTF2( realT & ETF, 
                              realT & NTF, 
                              int fi, 
                              realT g
                            )
{
   #ifndef ALLOW_F_ZERO
   if(_f[fi] <= 0)
   #else
   if(_f[fi] < 0)
   #endif
   {
      ETF = 0;
      NTF = 0;
      return;
   }

   complexT H_dm, H_del, H_con;

   complexT olX =  olXfer(fi, H_dm, H_del, H_con); //H_dm*H_wfs*H_ma*H_del*H_con;

   if(_f[fi] == 0)
   {
   }
   
   ETF = norm(realT(1)/( realT(1) + g*olX));
   NTF = norm(-(H_dm*H_del*g*H_con) /( realT(1) + g*olX) );

   /*if(_f[fi] == 0)
   {
      std::cerr << "ETF: " << ETF << " NTF: " << NTF << "\n";
   }*/
}

template<typename realT>
realT clGainOpt<realT>::clVariance( realT & varErr,
                                    realT & varNoise,
                                    std::vector<realT> & PSDerr,
                                    std::vector<realT> & PSDnoise,
                                    realT g
                                  )
{
   if(_f.size() != PSDerr.size() || _f.size() != PSDnoise.size())
   {
      std::cerr << "clVariance: Frequency grid and PSDs must be same size." << std::endl;
      return -1;
   }

   realT ETF, NTF, df;

   varErr = 0;
   varNoise = 0;

   df = _f[1] - _f[0];

   for(size_t i=0; i< PSDerr.size(); ++i)
   {
      if(g == 0)
      {
         ETF = 1;
         NTF = 0;
      }
      else
      {
         clTF2(ETF, NTF, i, g);
      }
      varErr += ETF*PSDerr[i]*df;
      varNoise += NTF*PSDnoise[i]*df;
   }

   return varErr + varNoise;
}

template<typename realT>
realT clGainOpt<realT>::clVariance( std::vector<realT> & PSDerr,
                                    std::vector<realT> & PSDnoise,
                                    realT g
                                  )
{
   realT varErr;
   realT varNoise;

   return clVariance(varErr, varNoise, PSDerr, PSDnoise, g );
}


template<typename realT>
realT clGainOpt<realT>::maxStableGain( realT & ll, realT & ul)
{
   static_cast<void>(ul);
   
   std::vector<realT> re, im;

   if(ll == 0) ll = _maxFindMin;

   nyquist( re, im, 1.0);

   int gi_c = re.size()-1;

   for(int gi=re.size()-2; gi >= 0; --gi)
   {
      if( -1.0/re[gi] < ll) continue;

      if( ( re[gi] < 0) && ( im[gi+1] >= 0 && im[gi] < 0) )
      {
         //Check for loop back in Nyquist diagram
         if( re[gi] <= re[gi_c] ) gi_c = gi;
      }
   }

   return -1.0/ re[gi_c];

}

template<typename realT>
realT maxStableGain(const realT & ll, const realT & ul)
{
   realT rll = ll;
   realT rul = ul;

   maxStableGain(rll, rul);
}

template<typename realT>
realT clGainOpt<realT>::maxStableGain()
{
   realT ul = 0;
   realT ll = _maxFindMin;

   return maxStableGain(ll, ul);
}

//Implement the minimization, allowing pre-compiled specializations
namespace impl 
{
   
template<typename realT>
realT optGainOpenLoop( clGainOptOptGain_OL<realT> & olgo,
                       realT & var,
                       realT & gmax,
                       realT & minFindMin,
                       realT & minFindMaxFact,
                       int minFindBits,
                       uintmax_t minFindMaxIter
                     )
{
   #ifdef MX_INCLUDE_BOOST
   realT gopt;

   try
   {
      std::pair<realT,realT> brack;
      brack = boost::math::tools::brent_find_minima<clGainOptOptGain_OL<realT>, realT>(olgo, minFindMin, minFindMaxFact*gmax, minFindBits, minFindMaxIter);
      gopt = brack.first;
      var = brack.second;
   }
   catch(...)
   {
      std::cerr << "optGainOpenLoop: No root found\n";
      gopt = minFindMaxFact*gmax;
      var = 0;
   }
   
   return gopt;
   #else
   static_assert(std::is_fundamental<realT>::value || !std::is_fundamental<realT>::value, "impl::optGainOpenLoop<realT> is not specialized for type realT, and MX_INCLUDE_BOOST is not defined, so I can't just use boost.");
   return 0;
   #endif
}

template<>
float optGainOpenLoop<float>( clGainOptOptGain_OL<float> & olgo,
                              float & var,
                              float & gmax,
                              float & minFindMin,
                              float & minFindMaxFact,
                              int minFindBits,
                              uintmax_t minFindMaxIter
                            );

template<>
double optGainOpenLoop<double>( clGainOptOptGain_OL<double> & olgo,
                                double & var,
                                double & gmax,
                                double & minFindMin,
                                double & minFindMaxFact,
                                int minFindBits,
                                uintmax_t minFindMaxIter
                              );

template<>
long double optGainOpenLoop<long double>( clGainOptOptGain_OL<long double> & olgo,
                                          long double & var,
                                          long double & gmax,
                                          long double & minFindMin,
                                          long double & minFindMaxFact,
                                          int minFindBits,
                                          uintmax_t minFindMaxIter
                                        );

#ifdef HASQUAD
template<>
__float128 optGainOpenLoop<__float128>( clGainOptOptGain_OL<__float128> & olgo,
                                        __float128 & var,
                                        __float128 & gmax,
                                        __float128 & minFindMin,
                                        __float128 & minFindMaxFact,
                                        int minFindBits,
                                        uintmax_t minFindMaxIter
                                      );
#endif

} //namespace impl

template<typename realT>
realT clGainOpt<realT>::optGainOpenLoop( realT & var,
                                         std::vector<realT> & PSDerr,
                                         std::vector<realT> & PSDnoise)
{
   realT gmax = 0;
   return optGainOpenLoop(var, PSDerr, PSDnoise, gmax);
}


template<typename realT>
realT clGainOpt<realT>::optGainOpenLoop( realT & var,
                                         std::vector<realT> & PSDerr,
                                         std::vector<realT> & PSDnoise,
                                         realT & gmax )
{
   clGainOptOptGain_OL<realT> olgo;
   olgo.go = this;
   olgo.PSDerr = &PSDerr;
   olgo.PSDnoise = &PSDnoise;

   if(gmax <= 0) gmax = maxStableGain();
   
   return impl::optGainOpenLoop( olgo, var, gmax, _minFindMin, _minFindMaxFact, _minFindBits, _minFindMaxIter);
}

template<typename realT>
int clGainOpt<realT>::pseudoOpenLoop( std::vector<realT> & PSD, realT g)
{
   realT e;
   for(int f=0; f< _f.size(); ++f)
   {
      e = clETF2(f, g);

      if(e > 0) PSD[f] = PSD[f]/ e;
   }

   return 0;
}

template<typename realT>
int clGainOpt<realT>::nyquist( std::vector<realT> & re,
                               std::vector<realT> & im,
                               realT g
                             )
{
   re.resize(_f.size());
   im.resize(_f.size());

   complexT etf;

   for(size_t f=0; f< _f.size(); ++f)
   {
      etf = g*olXfer(f);// clETF(f, g);
      re[f] = real(etf);
      im[f] = imag(etf);
   }
   
   return 0;
}

//------------ Workers ---------------------

///Bisection worker struct for finding optimum closed loop gain from open loop PSDs
template<typename realT>
struct clGainOptOptGain_OL
{
   clGainOpt<realT> * go;
   std::vector<realT> * PSDerr;
   std::vector<realT> * PSDnoise;

   realT operator()(const realT & g)
   {
      return go->clVariance(*PSDerr, *PSDnoise, g);
   }
};


//Explicit Instantiation
extern template
class clGainOpt<float>;

extern template
class clGainOpt<double>;

extern template
class clGainOpt<long double>;

#ifdef HASQUAD
extern template
class clGainOpt<__float128>;
#endif

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //clGainOpt_hpp
