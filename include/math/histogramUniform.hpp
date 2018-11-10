/** \file histogramUniform.hpp
  * \author Jared R. Males
  * \brief  Header for the std::vector utilities
  * \ingroup gen_math_files
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

#ifndef math_histogramUniform_hpp
#define math_histogramUniform_hpp

namespace mx
{
namespace math 
{
   
///A histogram with uniform bin spacing
/** Calculates the frequency in bins with uniform spacing.
  * 
  * \tparam realT the real data type
  * 
  * \ingroup gen_math
  */
template<typename realT>
class histogramUniform
{
public:

   realT _min; ///<The mininum bin location
   realT _max; ///<The maximum bin location
   realT _width; ///<The bin width

   std::vector<realT> _freqs; ///<The frequencies, one for each bin.

   ///Setup the histogram, performing allocations.
   void setup( realT mn, ///< [in] the new minimum bin location
               realT mx, ///< [in] the new maximum bin location
               realT w   ///< [in] the bin width
             )
   {
      _min = mn;
      _max = mx;
      _width = w;
      
      reset();
   }
   
   ///Resize and 0 the frequency vector.  Assumes _min, _max, and _width are set.
   void reset()
   {
      _freqs.resize( (_max - _min)/_width + 1, 0);
   }
   
   ///Accumulate a value in the appropriate bin.
   void accum( const realT & val /**< [in] The value to accumulate */)
   {
      int i = (val - _min) / _width;
      if( i < 0) i = 0;
      if( i >= _freqs.size()) i = _freqs.size()-1;
      
      
      ++_freqs[i];
   }
   
   ///Accumulate a vector of values.
   void accum( const std::vector<realT> & vals /**< [in] The vector of values to accumulate */)
   {
      for(int i=0; i< vals.size(); ++i) accum(vals[i]);
   }
   
   ///Get the frequency in the i-th bin.
   realT freq(int i /**< [in] the bin number */)
   {
      return _freqs[i]; ///\returns the current value of _freqs[i].
   }
   
   ///Get the number of bins
   int bins()
   {
      return _freqs.size(); ///\returns the size of the frequency vector.
   }
   
   ///Get the value of the left-edge of the i-th bin.
   realT binLeft(int i /**< [in] the bin number */)
   {
      return _min + i*_width; ///\returns the value of the left-edge of the i-th bin.
   }
   
   ///Get the value of the middle of the i-th bin.
   realT binMid(int i /**< [in] the bin number */)
   {
      return _min + i*_width + 0.5*_width; ///\returns the value of the middle of the i-th bin.
   }
   
   ///Get the value of the right edge of the i-th bin.
   realT binRight(int i /**< [in] the bin number */)
   {
      return _min + i*_width + 1.0*_width; ///\returns the value of the right edge of the i-th bin. 
   }
   
   ///Normalize the current frequencies so that the integral over all bins is 1.
   /** This normalizes the histogram so that it is a probability distribution, such that the sum
    * \f$ \sum_i P_i \Delta x = 1 \f$ where \f$ \Delta x \f$ is the bin width.
    */
   void normalize( int excludeTop = 0 /**< [in] [optional] specifies a number of bins at the top of the range to exclude from the sum */)
   {
      realT sum = 0;
      
      for(int i=0; i< _freqs.size()-excludeTop; ++i)
      {
         sum += _freqs[i];
      }
      
      for(int i=0; i< _freqs.size(); ++i)
      {
         _freqs[i] /= (sum*_width);
      }
   }
   
};

}//namespace math
}//namespace mx 

#endif //math_histogramUniform_hpp 
