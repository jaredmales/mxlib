/** \file clAOLinearPredictor.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to manage closed loop gain linear predictor determination.
  * \ingroup mxAO_files
  * 
  */

#ifndef clAOLinearPredictor_hpp
#define clAOLinearPredictor_hpp

#include <vector>

#include "../../math/geo.hpp"

#include "../../sigproc/psdUtils.hpp"

#include "../../sigproc/autocorrelation.hpp"
#include "../../sigproc/linearPredictor.hpp"

#include "clGainOpt.hpp"

namespace mx
{
namespace AO
{
namespace analysis 
{
   
///Class to manage the calculation of linear predictor coefficients for a closed-loop AO system.
/**
  * \tparam _realT the real floating point type in which to do all arithmetic.
  * 
  * \ingroup mxAOAnalytic
  */
template<typename _realT>
struct clAOLinearPredictor
{
   typedef _realT realT;
   
   std::vector<realT> PSDtn;
   
   std::vector<realT> psd2s;
   
   std::vector<realT> ac;
   
   sigproc::autocorrelationFromPSD<realT> acpsd;

   sigproc::linearPredictor<realT> lp;
   
   realT m_min_var0 {0};
   realT m_min_sc0 {10};
   realT m_precision0 {10};
   realT m_max_sc0 {100};
   realT m_dPrecision {3};
   
   realT m_gmax_lp {5}; ///< The maximum allowable gain for LP.
   
   //Stopping conditions:
   realT m_minPrecision {0.001};
   int m_maxIts {100};

   int extrap;
   
   clAOLinearPredictor()
   {
      extrap = 0;
   }
   
   ///Calculate the LP coefficients for a turbulence PSD and a noise PSD.
   /** This combines the two PSDs, augments to two-sided, and calls the linearPredictor.calcCoefficients method.
     * 
     * A regularization constant can be added to the PSD as well.
     * 
     */
   int calcCoefficients( std::vector<realT> & PSDt, ///< [in] the turbulence PSD
                         std::vector<realT> & PSDn, ///< [in] the WFS noise PSD
                         realT PSDreg,              ///< [in] the regularizing constant.  Set to 0 to not use.
                         int Nc,                     ///< [in] the number of LP coefficients.
                         realT condition = 0
                       )
   {
      PSDtn.resize(PSDt.size());

      for(size_t i=0; i< PSDt.size(); ++i)
      {
         PSDtn[i] = PSDt[i] +  PSDn[i] + PSDreg;
      }
      
      sigproc::augment1SidedPSD( psd2s, PSDtn,1);

      ac.resize(psd2s.size());
   
      //#pragma omp critical
      acpsd(ac, psd2s);

      return lp.calcCoefficients(ac, Nc, condition, extrap);
   }
   

   ///Worker function for regularizing the PSD for coefficient calculation.
   /**
     * \tparam printout if true then the results are printed to stdout as they are calculated.
     * 
     * On first call (min_var = 0):
     *     loop over scale factors from min_sc to max_sc (<=) in steps of precision.
     * 
     * On subsequent calls, when min_var and min_sc are passed back in
     *     loop over scale factors from min_sc-precision to max_sc in steps of 
     */ 
   template< bool printout=false>
   int _regularizeCoefficients( realT & min_var,           ///< [in/out] the minimum variance found.  Set to 0 on initial call
                                realT & min_sc,            ///< [in/out] the scale factor at the minimum variance.
                                realT precision,           ///< [in] the step-size for the scale factor
                                realT max_sc,              ///< [in] the maximum scale factor to test
                                clGainOpt<realT> & go_lp,  ///< [in] the gain optimization object
                                std::vector<realT> & PSDt, ///< [in] the turbulence PSD
                                std::vector<realT> & PSDn, ///< [in] the WFS noise PSD
                                int Nc                     ///< [in] the number of coefficients
                              )
   {
      realT gmax_lp;
      realT gopt_lp;
      realT var_lp;
      
      realT sc0;
      
      realT last_var;
      
      
      
      //min_var == 0 indicates first call
      if( min_var == 0)
      {
         sc0 = min_sc;  
         min_var = std::numeric_limits<realT>::max();
         last_var = std::numeric_limits<realT>::max();
      }
      else
      {
         sc0 = min_sc - precision*m_dPrecision;
         //sc0 = min_sc;  
         last_var = min_var;
      }

      //Test from sc0 to max_sc in steps of precision
      for(realT sc = sc0; sc <= max_sc; sc += precision)
      {
         if( calcCoefficients(PSDt, PSDn, PSDt[0]*pow(10, -sc/10), Nc) < 0) return -1;
   
         go_lp.a(lp._c);
         //go_lp.a(std::vector<realT>({1}));
         go_lp.b(lp._c);
         //go_lp.b(std::vector<realT>({lp._c(0,0)}));
         
         realT ll = 0, ul = 0;
         gmax_lp = go_lp.maxStableGain(ll,ul);
         if(gmax_lp> m_gmax_lp) gmax_lp = m_gmax_lp;
         gopt_lp = go_lp.optGainOpenLoop(var_lp, PSDt, PSDn, gmax_lp);
         //var_lp = go_lp.clVariance(PSDt, PSDn, gopt_lp);
      
         if(printout)
         {
            std::cout << -(sc)/10 << " " << gmax_lp << " " << gopt_lp << " " << var_lp << "\n";
         }
         
         if( var_lp < min_var )
         {
            min_var = var_lp;
            min_sc = sc;
         }
   
         //A jump by a factor of 10 indicates the wall
         if( var_lp/last_var > 10 ) return 0;
         
         last_var = var_lp;
      }
      
      return -1;
      
   }
   
   /// Regularize the PSD and calculate the associated LP coefficients.
   /** The PSD is regularized by adding a constant to it.  This constant is found by minimizing the variance of the residual PSD.
     * 
     * \tparam printout if true then the results are printed to stdout as they are calculated.
     */
   template< bool printout=false>
   int regularizeCoefficients( realT & gmax_lp,           ///< [out] the maximum gain calculated for the regularized PSD
                               realT & gopt_lp,           ///< [out] the optimum gain calculated for the regularized PSD
                               realT & var_lp,            ///< [out] the variance at the optimum gain.
                               clGainOpt<realT> & go_lp,  ///< [in] the gain optimization object
                               std::vector<realT> & PSDt, ///< [in] the turbulence PSD
                               std::vector<realT> & PSDn, ///< [in] the WFS noise PSD
                               int Nc                     ///< [in] the number of coefficients
                             )
   {

      realT min_var = m_min_var0;
      realT min_sc = m_min_sc0;
      realT precision = m_precision0;
      realT max_sc = m_max_sc0;
   
      
      int its = 0;
      while(precision > m_minPrecision && its < m_maxIts)
      {
         _regularizeCoefficients<printout>( min_var, min_sc, precision, max_sc, go_lp, PSDt, PSDn, Nc);
         
         if( min_sc == max_sc )
         {
            if( its == 0 )
            {
               min_sc -= precision;
               max_sc = 200;
            }
            else
            {
               //std::cerr << "Error in regularizeCoefficients.\n";
               //return -1;
               break;
            }
         }
         else
         {
            max_sc = min_sc + precision;
            precision /= m_dPrecision;
         }
         
         ++its;
      }
   
      //Now record final values
      if(calcCoefficients(PSDt, PSDn, PSDt[0]*pow(10, -min_sc/10), Nc) < 0) return -1;
   
      go_lp.a(lp._c);
      //go_lp.a(std::vector<realT>({1}));
      go_lp.b(lp._c);
      //go_lp.b(std::vector<realT>({lp._c(0,0)}));
      
      realT ll = 0, ul = 0;
      gmax_lp = go_lp.maxStableGain(ll,ul);
      gopt_lp = go_lp.optGainOpenLoop(var_lp, PSDt, PSDn, gmax_lp);
      
      return 0;
   }
   
   /// Regularize the PSD and calculate the associated LP coefficients.
   /** The PSD is regularized by adding a constant to it.  This constant is found by minimizing the variance of the residual PSD.
     * 
     * \tparam printout if true then the results are printed to stdout as they are calculated.
     */
   template< bool printout=false>
   int optimizeNc( realT & gmax_lp,           ///< [out] the maximum gain calculated for the regularized PSD
                               realT & gopt_lp,           ///< [out] the optimum gain calculated for the regularized PSD
                               int & Nc,
                               realT & var_lp,            ///< [out] the variance at the optimum gain.
                               clGainOpt<realT> & go_lp,  ///< [in] the gain optimization object
                               std::vector<realT> & PSDt, ///< [in] the turbulence PSD
                               std::vector<realT> & PSDn, ///< [in] the WFS noise PSD
                               int minNc,                     ///< [in] the number of coefficients
                               int maxNc
                             )
   {
      realT minVar = std::numeric_limits<realT>::max();
      
      for(int n = minNc; n <= maxNc; ++n)
      {
         realT _gmax_lp;
         realT _gopt_lp;
         realT _var_lp;
         regularizeCoefficients<printout>(_gmax_lp, _gopt_lp, _var_lp, go_lp, PSDt, PSDn, n);
         
         if(_var_lp < minVar)
         {
            gmax_lp = _gmax_lp;
            gopt_lp = _gopt_lp;
            var_lp = _var_lp;
            Nc = n;
            
            minVar = var_lp;
         }
      }
      
      return 0;
   }
   
};

}//namespace analysis
}//namespace ao
}//namespace mx

#endif // clAOLinearPredictor_hpp
