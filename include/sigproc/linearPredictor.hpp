/** \file linearPredictor.hpp
  * \brief Working with linear prediction.
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef linearPredictor_hpp
#define linearPredictor_hpp

#include <vector>
#include <complex>

#include "../math/constants.hpp"
#include "../math/eigenLapack.hpp"

#include "levinsonRecursion.hpp"

namespace mx 
{
namespace sigproc 
{
   
/// A class to support linear prediction.   
/** \ingroup signal_processing 
  *
  * \todo document linearPredictor 
  */ 
template<typename _realT>
struct linearPredictor
{
   typedef _realT realT;
   
   Eigen::Array<realT, -1, -1> _c;
   
   realT _setCondition {0};
   realT _actCondition {0};
   int _nRejected {0};
   
   /// Calculate the LP coefficients given an autocorrelation.
   /** If condition==0 then the levinson recursion is used.
     * Otherwise, SVD pseudo-inversion is used with the given condition number.
     */ 
   int calcCoefficients( std::vector<realT> & ac,
                         size_t Nc,
                         realT condition = 0,
                         size_t extrap = 0
                       )
   {
      
      if(condition == 0)
      {
         return calcCoefficientsLevinson( ac, Nc, extrap );
      }
      
      
      Eigen::Array<realT,-1,-1> Rmat, Rvec, PInv, LPcoeff;   
      
      Rmat.resize(Nc, Nc);
      Rvec.resize(1, Nc);
   
      for(size_t i=0; i<Nc; ++i)
      {
         for(size_t j=0; j<Nc; ++j)
         {
            Rmat(i,j) = ac[ fabs(i-j)];
         }
         
         Rvec(0,i) = ac[i+1];
      }  

      realT tmpCond = condition;
      
      _setCondition = condition;
      math::eigenPseudoInverse(PInv, tmpCond, _nRejected,  Rmat, condition);
      
      
      _actCondition = tmpCond;
      
      _c = Rvec.matrix()*PInv.matrix();

      
      
      return 0;
   }

   int calcCoefficientsLevinson( std::vector<realT> & ac,
                                 size_t Nc,
                                 size_t extrap = 0
                               )
   {
      std::vector<realT> r, x, y;
      
      r.resize(2.*Nc-1);
      x.resize(Nc);
      y.resize(Nc);
      
      for(size_t i=0; i< Nc; ++i) r[i] = ac[Nc-i - 1];
      for(size_t i=Nc; i< 2*Nc-1; ++i) r[i] = ac[i-Nc+1];
      
      for(size_t i=0;i<Nc; ++i) y[i] = ac[i+1];
      
      levinsonRecursion(r.data(), x.data(), y.data(), Nc);
      
      _c.resize(1, Nc);
      for(size_t i=0; i< Nc; ++i) _c(0,i) = x[i];
    
      
      if(extrap == 1)
      {
         Eigen::Array<realT, -1, -1> ex(_c.rows(), _c.cols());
         for(size_t j=0; j < extrap; ++j)
         {
            for(size_t i=0; i< Nc-1; ++i)
            {
               ex(0,i) = _c(0,0)*_c(0,i) + _c(0,i+1);
            }
            ex(0,Nc-1) = _c(0,0)*_c(0,Nc-1);
         
            _c = ex;
         }
      }
      
      
      
      return 0;
   }
   
   realT c(size_t i)
   {
      return _c(0,i);
   }
   
   realT predict( std::vector<realT> & hist,
                  int idx )
   {
      realT x = 0;
      
      if(idx < _c.cols()) return 0;
      
      for(int i=0; i< _c.cols(); ++i)
      {
         x += _c(0,i)*hist[idx-i];
      }
      
      return x;
   }
   
   realT spectralResponse( realT f, realT fs)
   {
      int n = _c.cols();
      
      std::complex<realT> He = 0;
      for(int j=0; j < n; ++j)
      {
         realT s = (j+1.0)* math::two_pi<realT>();
         He += _c(0,j) * exp(s*std::complex<realT>(0,-1.0)*f/fs);
      }

      realT one = 1.0;
      return std::norm( one/ (one - He));
   }
   
};

} //namespace sigproc 
} //namespace mx
#endif //linearPredictor_hpp

