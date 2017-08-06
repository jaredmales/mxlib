/** \file clGainOpt.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to manage closed loop gain optimization.
  * \ingroup mxAO_files
  * 
  */

#ifndef __clGainOpt_hpp__
#define __clGainOpt_hpp__




#include <boost/math/tools/roots.hpp>
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
   
//forward declarations
template<typename realT>
struct clMaxStableGainOpt;

template<typename realT>
struct clGainOptOptGain_OL;

template<typename realT>
struct clGainOptOptGain_CL;

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
   int _maxFindBits; ///< The bits of precision to use for finding maximum stable gain. Defaults to 8.
   boost::uintmax_t _maxFindMaxIter; ///< The maximum iterations allowed for finding maximum stable gain.
   realT _dg; ///< Gain stepsize for finite difference derivative calculation. Default = 0.01;
   int _stabULMaxIt; ///< Maximum number of iterations for establishing upper limit.  Default = 10;
   
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
     * \param Ti is the desired loop sampling interval.
     * \param tau is the desired loop delay.
     */
   clGainOpt(realT Ti, realT tau);

   /// Initialize this instance.
   void init();

   /// Get the number of integrations in the (optional) moving average
   /**
     * \returns the current value of _N.
     */ 
   int N(); 
   
   /// Set the number of integrations in the moving average
   /**
     * \param [in] newN is the value of _N.
     */ 
   void N( int newN );
   
   /// Get the loop sampling interval
   /**
     * \returns the current value of _Ti.
     */ 
   realT Ti();
   
   /// Set the loop sampling interval
   /**
     * \param [in] newTi is the new value of _Ti.
     */ 
   void Ti( realT newTi);
   
   /// Get the loop delay
   /**
     * \returns the current value of _tau.
     */ 
   realT tau();
   
   /// Set the loop delay
   /**
     * \param [in] newTau is the new value of _tau.
     */
   void tau(realT newTau);

   /// Set the vector of FIR coefficients
   /**
     * \param [in] newB is a vector of coefficients, which is copied to _b.
     */ 
   void b( const std::vector<realT> & newB); 
   
   /// Set the vector of FIR coefficients
   /**
     * \param [in] newB is a column-vector Eigen::Array of coefficients, which is copied to _b.
     */ 
   void b( const Eigen::Array<realT, -1, -1> & newB);
   
   /// Set the vector of IIR coefficients
   /**
     * \param [in] newA is a vector of coefficients, which is copied to _a.
     */ 
   void a (const std::vector<realT> & newA);
   
   /// Set the vector of IIR coefficients
   /**
     * \param [in] newA is a column-vector Eigen::Array of coefficients, which is copied to _a.
     */    
   void a( const Eigen::Array<realT, -1, -1> & newA);
   
   /// Set the FIR and IIR coefficients so that the control law is a leaky integrator.
   /** Set remember to 1.0 for a pure integrator control law.
     * 
     * \param [in] remember is a number usually close to 1 setting the amount "remembered" from previous iterations.
     */
   void setLeakyIntegrator(realT remember);
   
   /// Set the vector of frequencies
   /**
     * \param [in] newF is a pointer to an array containing the new frequencies
     * \param [in] nf is the number of elements of size(realT) in newF.
     */ 
   void f (realT * newF, size_t nF); 
   
   /// Set the vector of frequencies
   /**
     * \param [in] newF is a vector containing the new frequencies
     */
   void f (std::vector<realT> & newF); 
   
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
     * \param [in] fi is the index of the frequency
     * 
     * \returns the complex value of the open-loop transfer function at f[fi].
     */   
   complexT olXfer(int fi);
   
   /// Calculate the open-loop transfer function
   /** 
     * \param [in] fi is the index of the frequency
     * \param [out] H_dm [optional] the transfer function of the DM
     * \param [out] H_del [optional] the delay transfer function
     * \param [out] H_con [optional] the controller transfer functionf.
     * 
     * \return the complex value of the open-loop transfer function at f.
     */
   complexT olXfer( int fi, 
                    complexT & H_dm, 
                    complexT & H_del, 
                    complexT & H_con );
   
   /// Return the closed loop error transfer function (ETF) at frequency f for gain g.
   /**
     * \param [in] fi the index of the frequency at which to calculate the ETF
     * \param [in] g the loop gain.
     * 
     * \returns the closed loop ETF at f and g.
     */
   complexT clETF( int fi, 
                   realT g );
   
   /// Return the norm of the closed loop error transfer function (ETF) at frequency f for gain g.
   /**
     * \param [in] fi the index of the frequency at which to calculate the ETF
     * \param [in] g the loop gain.
     * 
     * \returns the norm of the closed loop ETF at f and g.
     */
   realT clETF2( int fi, 
                 realT g );
   
   /// Return the norm of the closed loop noise transfer function (NTF) at frequency f for gain g.
   /**
     * \param [in] fi the index of the frequency at which to calculate the NTF
     * \param [in] g the loop gain.
     * 
     * \returns the value of the closed loop NTF at f and g.
     */
   realT clNTF2( int fi, 
                 realT g );
   
   /// Return the norm of the closed loop transfer functions at frequency f for gain g.
   /** Calculates both the error transfer function (ETF) and the noise transfer function (NTF).
     * This minimizes the various complex number operations, compared to calling both clETF2 and clNTF2.
     * 
     * \param [out] ETF is set to the ETF at f and g
     * \param [out] NTF is set to the NTF at f and g
     * \param [in] fi is the index of the frequency at which to calculate the ETF
     * \param [in] g is the loop gain.
     *
     */
   void clTF2( realT & ETF, 
               realT & NTF, 
               int fi, 
               realT g );
   
         
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
     * \param [in] freq [in] is the frequency grid.
     * \param [in] PSDerr [in] is the open-loop process error PSD. 
     * \param [in] PSDnoise [in] is the open-loop measurement noise PSD.
     * \param [in] g [in] is the gain
     * \param [out] varErr is the variance in the residual process error
     * \param [out] varNoise is the variance in the residual measurement noise
     * 
     * \returns the total variance (error + noise) in closed loop
     */ 
   realT clVariance( std::vector<realT> & PSDerr,
                     std::vector<realT> & PSDnoise,
                     realT g,
                     realT & varErr,
                     realT & varNoise );      
   
   /// Calculate the closed loop variance given open-loop PSDs and gain
   /** Overload of clVariance without the varErr and varNoise output parameters.
     * 
     * \param [in] PSDerr  is the open-loop process error PSD. 
     * \param [in] PSDnoise [in] is the open-loop measurement noise PSD.
     * \param [in] g  is the gain
     * 
     * \returns the total variance (error + noise) in closed loop
     */ 
   realT clVariance( std::vector<realT> & PSDerr,
                     std::vector<realT> & PSDnoise,
                     realT g );
   
   /// Calculate the summed derivative of the ETF w.r.t. to the gain.
   /** The derivative of the ETF is calculated by the finite difference approximation, and summing over 
     * the frequency.  The size of the finite difference is 
     * specified by \ref _dg.
     *  
     * \Note To find the true derivative multiply by 1 / (2 * \ref _dg)
     * 
     * \param [in] g is the gain
     * 
     * \returns the summed derivative of the ETF
     */ 
   realT etfDerivative(realT g);

   realT findPeakBelow(realT);
   
   /// Find the maximum stable gain for the loop parameters
   /** Uses bisection to find the root of the derivative of the summed ETF, using the boost bisect routine.
     *
     * \returns the maximum stable gain for the loop parameters
     */
   realT maxStableGain();
   
   /// Find the maximum stable gain for the loop parameters
   /** Uses bisection to find the root of the derivative of the summed ETF, using the boost bisect routine.
     *
     * \param ul [out] [optional] is set to the upper limit used for the bisection.
     *
     * \returns the maximum stable gain for the loop parameters
     */ 
   realT maxStableGain(realT & ll, realT & ul);
      
   realT maxStableGainDo(realT & ll, realT & ul);
      
   ///Return the optimum closed loop gain given an open loop PSD
   /**
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
   int pseudoOpenLoop( std::vector<realT> & PSD, ///< [in/out] input closed loop PSD, on output ocntains the pseudo open loop error PSD ,
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
   _maxFindBits = 8;//std::numeric_limits<realT>::digits;
   _maxFindMaxIter = 1000;
   _dg = 0.01;
   _stabULMaxIt = 10;
   
   
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
   if( newA.size()+1 > _cs.cols() && newA.size()+1 > _cs.cols())
   {
      _fChanged = true;
   }
   
   _a = newA;
   _changed = true;
}

template<typename realT>
void clGainOpt<realT>::a( const Eigen::Array<realT, -1, -1>  & newA)
{
   if( newA.cols() +1 > _cs.cols() && newA.cols() + 1 > _cs.cols())
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
realT clGainOpt<realT>::clVariance( std::vector<realT> & PSDerr,
                                    std::vector<realT> & PSDnoise,
                                    realT g,
                                    realT & varErr,
                                    realT & varNoise )
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
                                    realT g )
{
   realT varErr;
   realT varNoise;
   
   return clVariance(PSDerr, PSDnoise, g, varErr, varNoise);
}

template<typename realT>
realT clGainOpt<realT>::etfDerivative(realT g)
{
   realT ds = 0;

   for(int f=0; f< _f.size(); ++f)
   {
      ds += clETF2(f, g+_dg) - clETF2(f, g-_dg);
   }
   
   return ds;
}

template<typename realT>
realT clGainOpt<realT>::maxStableGain()
{
   realT ul = 0;
   realT ll = _maxFindMin;
   
   return maxStableGain(ll, ul);
}

template<typename realT>
realT clGainOpt<realT>::findPeakBelow(realT st)
{
   realT cur = etfDerivative(st);
   realT last = cur;
   
   st -= 0.001;
   
   cur = etfDerivative(st);
   
   while( cur > last)
   {
      st-= 0.001;
      last = cur;
      cur = etfDerivative(st);
   }
   
   st += 0.001;
   
   return st;
}

template<typename realT>
realT clGainOpt<realT>::maxStableGain(realT & ll, realT & ul)
{
   std::vector<realT> re, im;
   
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
      
      if( ( re[gi] < 0) && (im[gi+1] < 0 && im[gi] >= 0) ) break;
      
      if( -1.0/re[gi] > ul && ul != 0 ) 
      {
         return ul;
      }
   }
   
   return -1.0/ re[gi_c];
   

#if 0   
   realT msg;
   
   realT minG = ll;
   realT maxG = ul;
   
   msg = maxStableGainDo(ll, ul);
   
   //std::cerr << "ETF Deriv: " << etfDerivative(.99*msg) << "\n";
   
   bool again = true;
   int tries = 0;
   
   while(again && tries < _maxFindMaxIter)
   {
      ++tries;
      
      if( findPeakBelow(msg) > 1e8)
      {
         maxG = 0.98*msg;
      
         while( etfDerivative(maxG) > 0 && maxG > ll) maxG -= 0.01;
      
         if(maxG > ll) 
         {
            //std::cerr << "1) Searching: " << ll << " " << maxG << "\n";
            realT msg2 = maxStableGainDo(ll, maxG);
      
            //std::cerr << "Search result: " << msg2 << " " << etfDerivative(0.99*msg2) << "\n";
            if(msg2 > -1 && findPeakBelow(msg2 ) > 1e8)
            {
               //std::cerr << "2nd pole found!\n";
               msg = msg2;
            }
            else again=false;
         }
         else
         {
            again = false;
         }
      }
      else
      {
         //std::cerr << "Probably not real pole found\n";
            
         //again = false;
         minG = 1.02*ll;
         
         while( etfDerivative(minG) < 0 && minG < ul) minG += 0.01;
         
         if(minG < ul)
         {
            ul = 0.0;
            //std::cerr << "2) Searching: " << minG << " " << ul << "\n";
            realT msg2 = maxStableGainDo(minG, ul);
      
            //std::cerr << "Search result: " << msg2 << " " << etfDerivative(0.99*msg2) << "\n";
            if(msg2 > -1 && findPeakBelow(msg2) > 1e8)
            {
               //std::cerr << "2nd pole found!\n";
               msg = msg2;
               ll = minG;
            }
            else
            {
               //std::cerr << "No Luck\n";
               ll = minG; //now we try again.
               msg = msg2;
            }
         }
         else
         {
            again = false;
         } 
      }
   }
   return msg;
#endif
}
   
   
   

template<typename realT>
realT clGainOpt<realT>::maxStableGainDo(realT & ll, realT & ul)
{
   
   
   if(ll ==0) ll = _maxFindMin;
   
   realT s0 = etfDerivative(ll);
   
   if(ul == 0)
   {
      //We pick the upper limit (ul) by increasing by factors of 2 until there is a sign change between it an _maxFindMin.
      ul = 1.0;   
      realT s1 = etfDerivative(ul);
   
      int n = 1;
      while(std::signbit(s1) == std::signbit(s0) && n < _stabULMaxIt)
      {
         ul *= 2;
         s1 = etfDerivative(ul);
         ++n;
      }
   
      //error check
      if(n >= _stabULMaxIt || std::signbit(s1) == std::signbit(s0))
      {
         std::cerr << "maxStableGain failed to find upper limit." << std::endl;
         return -1;
      }
   }
   
   clMaxStableGainOpt<realT> msg;
   msg.go = this;

   try
   {
      //Bracket the root using bisection
      std::pair<realT,realT> brack;
      boost::math::tools::eps_tolerance<realT> tol(_maxFindBits);
      brack = boost::math::tools::bisect< clMaxStableGainOpt<realT>, 
                                          realT, 
                                          boost::math::tools::eps_tolerance<realT> >( msg, 
                                                                                      ll, 
                                                                                      ul, 
                                                                                      tol, 
                                                                                      _maxFindMaxIter);

      //Return the average of the bracket                                 
      return  brack.first; //0.5*(brack.first + brack.second);
   }
   catch(...)
   {
      std::cerr << "maxStableGain failed to find root." << std::endl;
      return -1;
   }
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
   
   try
   {
      std::pair<realT,realT> brack;
   
      brack = boost::math::tools::brent_find_minima<clGainOptOptGain_OL<realT>, realT>(olgo, _minFindMin, _minFindMaxFact*gmax, _minFindBits, _minFindMaxIter);

      return  brack.first;
   }
   catch(...)
   {
      std::cerr << "optGainOpenLoop: No root found\n";
      return _minFindMaxFact*gmax;
   }
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

///Bisection worker struct for maximum stable gain
template<typename realT>
struct clMaxStableGainOpt
{
   clGainOpt<realT> * go;
         
   realT operator()(const realT & g)
   {
      return go->etfDerivative(g);
   }
};      

///Bisection worker struct struct for finding optimum closed loop gain from open loop PSDs
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

#endif //__clGainOpt_hpp__
