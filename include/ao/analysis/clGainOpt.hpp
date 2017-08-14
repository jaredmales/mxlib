/** \file clGainOpt.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to manage closed loop gain optimization.
  * \ingroup mxAO_files
  * 
  */

#ifndef clGainOpt_hpp
#define clGainOpt_hpp




//#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/constants/constants.hpp>

#include <Eigen/Dense>

#include <mx/timeUtils.hpp>

using namespace boost::math::constants;

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
   /** Paramaters for stability analysis
     * @{
     */
   
   realT _maxFindMin; ///< The Minimum value for the maximum stable gain finding algorithm.
   //int _maxFindBits; ///< The bits of precision to use for finding maximum stable gain. Defaults to 8.
   //boost::uintmax_t _maxFindMaxIter; ///< The maximum iterations allowed for finding maximum stable gain.
   //realT _dg; ///< Gain stepsize for finite difference derivative calculation. Default = 0.01;
   //int _stabULMaxIt; ///< Maximum number of iterations for establishing upper limit.  Default = 10;
   
   ///@}
   
   /** Parameters for minimization finding
     * @{
     */ 
   
   realT _minFindMin; ///< The Minimum value for the minimum finding algorithm.
   realT _minFindMaxFact; ///< The maximum value, as a multiplicative factor of maximum gain
   int _minFindBits; ///< The bits of precision to use for minimum finding. Defaults to std::numeric_limits<realT>::digits.
   boost::uintmax_t _minFindMaxIter; ///< The maximum iterations allowed for minimization.
   
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
   
   /// Return the norm of the closed loop error transfer function (ETF) at frequency f for gain g.
   /**
     * \returns the norm of the closed loop ETF at f and g.
     */
   realT clETF2( int fi, ///< [in] the index of the frequency at which to calculate the ETF. 
                 realT g ///< [in] the loop gain. 
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
   realT maxStableGain( realT & ll, ///< [in/out] the lower limit used for the search 
                        realT & ul ///< [in/out] the upper limit used for hte search
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
   realT optGainOpenLoop( std::vector<realT> & PSDerr,    ///< [in] open loop error PSD
                          std::vector<realT> & PSDnoise   ///< [in] open loop measurement noise PSD
                        );
      
   ///Return the optimum closed loop gain given an open loop PSD
   /**
     * \returns the optimum gain.
     */ 
   realT optGainOpenLoop( std::vector<realT> & PSDerr,   ///< [in] open loop error PSD
                          std::vector<realT> & PSDnoise, ///< [in] open loop measurement noise PSD
                          realT & gmax                   ///< [in] maximum gain to consider.  If 0, then _gmax is used.
                        );

   ///Calculate the pseudo open-loop PSD given a closed loop PSD
   /**
     * \returns 0 on success
     */ 
   int pseudoOpenLoop( std::vector<realT> & PSD, ///< [in/out] input closed loop PSD, on output contains the pseudo open loop error PSD ,
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
   //_maxFindBits = 8;//std::numeric_limits<realT>::digits;
   //_maxFindMaxIter = 1000;
   //_dg = 0.01;
   //_stabULMaxIt = 10;
   
   
   _minFindMin = 1e-9;
   _minFindMaxFact = 0.98;
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
   if( newB.size() > _cs.cols() && newB.size() > _cs.cols())
   {
      _fChanged = true;
   }
   
   _b = newB;
   _changed = true;
}
   
template<typename realT>
void clGainOpt<realT>::b( const Eigen::Array<realT, -1, -1>  & newB)
{
   if( newB.cols() > _cs.cols() && newB.cols() > _cs.cols())
   {
      _fChanged = true;
   }
   
   _b.resize( newB.cols());
   
   for(int i=0; i< _b.size(); ++i)
   {
      _b[i] = newB(0,i);
   }
   
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::a ( const std::vector<realT> & newA)
{
   if( newA.size()+1 > _cs.cols())
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
   
   for(int i=0; i< _a.size(); ++i)
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

template<typename realT>
std::complex<realT> clGainOpt<realT>::olXfer(int fi, complexT & H_dm, complexT & H_del, complexT & H_con)
{
   if(_f[fi] <= 0) 
   {
      H_dm = 0;
      H_del = 0;
      H_con = 0;
      return 0;
   }
   
   if(_fChanged)
   {
      int jmax = std::max(_a.size()+1, _b.size());
      
      //std::cerr << "_fChanged  " << jmax << std::endl;
      
      _cs.resize(_f.size(), jmax);
      _ss.resize(_f.size(), jmax);
      
      //#pragma omp parallel for
      for(int i=0; i< _f.size(); ++i)
      {
         _cs(i,0) = 1.0;
         _ss(i,0) = 0.0;
         
         for(int j = 1; j<jmax; ++j)
         {
            _cs(i,j) = cos(two_pi<realT>()*_f[i]*_Ti*realT(j));
            _ss(i,j) = sin(two_pi<realT>()*_f[i]*_Ti*realT(j));
         }
      }
      
      _fChanged = false;
   }
   
   if(_changed)
   {

      _H_dm.resize(_f.size());
      _H_wfs.resize(_f.size());
      _H_ma.resize(_f.size());
      _H_del.resize(_f.size());
      _H_con.resize(_f.size());
               
      int jmax = std::min(_a.size(), _b.size());

      //#pragma omp parallel for
      for(int i=0; i<_f.size(); ++i)
      {
         if(_f[i] <=0 ) continue;
         
         complexT s = complexT(0.0, two_pi<realT>()*_f[i]);
   
         complexT expsT = exp(-s*_Ti);
   
         _H_dm[i] = (realT(1) - expsT)/(s*_Ti);
   
         _H_wfs[i] = _H_dm[i];//(realT(1) - expsT)/(s*_Ti);
   
         _H_ma[i] = 1;//realT(1./_N)*(realT(1) - pow(expsT,_N))/(realT(1) - expsT);
         
         _H_del[i] = exp(-s*_tau); //complexT(_tau,0));
         
         complexT FIR = complexT(_b[0],0);
         complexT IIR = complexT(0.0, 0.0); 
        
         for(int j = 1; j < jmax; ++j)
         {
            //realT cs = _cs(i,j);//cos(2.*pi<realT>()*_f[i]*_Ti*realT(j));
            //realT ss = _ss(i,j);//sin(2.*pi<realT>()*_f[i]*_Ti*realT(j));
         
            complexT expZ = exp( -s*_Ti*realT(j));
            
            FIR += _b[j]*expZ; //complexT(cs, -ss);
            IIR += _a[j-1]*expZ;//complexT(cs, -ss);
         }
           
         for(int jj=jmax; jj< _a.size()+1; ++jj)
         {
            //realT cs = _cs(i,jj);//cos(2.*pi<realT>()*_f[i]*_Ti*realT(jj));
            //realT ss = _ss(i,jj);//sin(2.*pi<realT>()*_f[i]*_Ti*realT(jj));

            complexT expZ = exp( -s*_Ti*realT(jj));
            
            IIR += _a[jj-1]*expZ; //complexT(cs, -ss);
         }
        
         for(int jj=jmax; jj<_b.size(); ++jj)
         {
            //realT cs = _cs(i,jj);//cos(2.*pi<realT>()*_f[i]*_Ti*realT(jj));
            //realT ss = _ss(i,jj);//sin(2.*pi<realT>()*_f[i]*_Ti*realT(jj));

            complexT expZ = exp( -s*_Ti*realT(jj));
            
            FIR += _b[jj]*expZ; //complexT(cs, -ss);
         }
         
         _H_con[i] = FIR/( realT(1.0) - IIR);

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
   if(_f[fi] <= 0) return 0;
      
   return (realT(1)/( realT(1) + g*olXfer(fi)));
}

template<typename realT>
realT clGainOpt<realT>::clETF2(int fi, realT g)
{
   if(_f[fi] <= 0) return 0;
      
   return norm(realT(1)/( realT(1) + g*olXfer(fi)));
}
   
template<typename realT>   
realT clGainOpt<realT>::clNTF2(int fi, realT g)
{
   if(_f[f] <= 0) return 0;
   
   complexT H_dm, H_del, H_con;
   
   complexT olX =  olXfer(fi, H_dm, H_del, H_con); //H_dm*H_wfs*H_ma*H_del*H_con;      
   
   complexT NTF = -(H_dm*H_del*g*H_con)/(realT(1) + g*olX);
   
   return norm(NTF);
}

template<typename realT>
void clGainOpt<realT>::clTF2(realT & ETF, realT & NTF, int fi, realT g)
{      
   if(_f[fi] <= 0) 
   {
      ETF = 0;
      NTF = 0;
      return;
   }
         
   complexT H_dm, H_del, H_con;
   
   complexT olX =  olXfer(fi, H_dm, H_del, H_con); //H_dm*H_wfs*H_ma*H_del*H_con;   
   
   ETF = norm(realT(1)/( realT(1) + g*olX));
   NTF = norm(-(H_dm*H_del*g*H_con))*ETF;
   
   //if(std::isnan(NTF)) NTF = 1;
   
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
   
   for(int i=0; i< PSDerr.size(); ++i)
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
   std::vector<realT> re, im;
   
   
   if(ll == 0) ll == _maxFindMin;
   
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
      
      //if( ( re[gi] < 0) && (im[gi+1] < 0 && im[gi] >= 0) ) break;
      
//       if( -1.0/re[gi] > ul && ul != 0 ) 
//       {
//          return ul;
//       }
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

template<typename realT>
realT clGainOpt<realT>::optGainOpenLoop( std::vector<realT> & PSDerr, 
                                         std::vector<realT> & PSDnoise)
{
   realT gmax = 0;
   return optGainOpenLoop(PSDerr, PSDnoise, gmax);
}


template<typename realT>
realT clGainOpt<realT>::optGainOpenLoop( std::vector<realT> & PSDerr, 
                                         std::vector<realT> & PSDnoise, 
                                         realT & gmax )
{
   clGainOptOptGain_OL<realT> olgo;
   olgo.go = this;
   //olgo.freq = &_f;
   olgo.PSDerr = &PSDerr;
   olgo.PSDnoise = &PSDnoise;
   
   if(gmax <= 0) gmax = maxStableGain();
   
   realT gopt;
   
   try
   {
      std::pair<realT,realT> brack;
      
      brack = boost::math::tools::brent_find_minima<clGainOptOptGain_OL<realT>, realT>(olgo, _minFindMin, _minFindMaxFact*gmax, _minFindBits, _minFindMaxIter);

      gopt = brack.first;
   }
   catch(...)
   {
      std::cerr << "optGainOpenLoop: No root found\n";
      gopt = _minFindMaxFact*gmax;
   }
   
   return gopt;
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
   
   for(int f=0; f< _f.size(); ++f)
   {
      etf = g*olXfer(f);// clETF(f, g);
      
      re[f] = real(etf);
      im[f] = imag(etf);
   }
}

//------------ Workers ---------------------

///Bisection worker struct for finding optimum closed loop gain from open loop PSDs
template<typename realT>
struct clGainOptOptGain_OL
{
   clGainOpt<realT> * go;
   //std::vector<realT> * freq;
   std::vector<realT> * PSDerr;
   std::vector<realT> * PSDnoise;
   
   realT operator()(const realT & g)
   {
      return go->clVariance(*PSDerr, *PSDnoise, g);
   }
};

} //namespace analysis
} //namespace AO
} //namespace mx

#endif //clGainOpt_hpp
