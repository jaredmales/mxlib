/** \file sourceFinder.hpp
  * \brief Declares and defines a class for finding sources in images
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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

#ifndef improc_sourceFinder_hpp
#define improc_sourceFinder_hpp


#include "eigenImage.hpp"

#include "../math/vectorUtils.hpp"

namespace mx
{
namespace improc 
{
   
/// Find stars in an image by SNR thresholding
/**
  * \ingroup image_utils
  */ 
template<typename _realT>
class sourceFinder
{

public:
   typedef _realT realT; ///< The real valued type in which to do calculations

   typedef std::pair<int,int> pixelT; ///< Pixel coordinates

protected:
   
   realT m_threshold {5}; ///< The SNR threshold to use

   realT m_minsep {2}; ///< The minimum separation between stars.  Closer SNR detections are treated as the same star.

   bool m_useMedian {true}; ///< If true, then the median is used for standard deviation calculation. The mean is used otherwise.

   //Working memory, retained to help with repeated calls.
   std::vector<float> m_vals;
   std::vector<pixelT> m_pixs;
   std::vector<size_t> m_snr;

public:
   
   /// Default c'tor
   sourceFinder();

   /// Constructor to set up the algorithm
   sourceFinder( realT thresh, ///< [in]
               realT minsep  ///< [in]
             );
   
   /// Clear the working memory
   void clear();

   /// Set the SNR threshold
   void treshold( const realT & thresh /**< [in] the new threshold*/ );

   /// Get the SNR threshold
   /** \returns the current value of m_threshold;
     */
   realT threshold();

   /// Set the minimum separation
   void minsep( const realT & msep /**< [in] the new minimum separation*/ );

   /// Get the minimum separation between stars
   /** \returns the current value of m_minsep;
     */
   realT minsep();

   /// Set the useMedian flag
   void useMedian( const bool & useMed /**< [in] the new useMedian flag value*/ );

   /// Get the useMedian flag
   /** \returns the current value of m_useMedian;
     */
   bool useMedian();

   /// Find stars in an image
   /** Returns the coordinates as (row, col) of the highest SNR pixel in each unique star.
     *
     * \returns a vector of pixel coordinates
     */ 
   std::vector<pixelT> operator()( const eigenImage<realT> & im /**< [in] the image to search for stars*/ );

};

template<typename realT>
sourceFinder<realT>::sourceFinder()
{
}

template<typename realT>
sourceFinder<realT>::sourceFinder( realT thresh,
                        realT minsep 
                      ) : m_threshold(thresh), m_minsep(minsep)
{
}

template<typename realT>
void sourceFinder<realT>::treshold(const realT & thresh)
{
   m_threshold = thresh;
}

template<typename realT>
void sourceFinder<realT>::clear()
{
   m_vals.clear();
   m_pixs.clear();
   m_snr.clear();
}

template<typename realT>
realT sourceFinder<realT>::threshold()
{
   return m_threshold;
}

template<typename realT>
void sourceFinder<realT>::minsep(const realT & msep)
{
   m_minsep = msep;
}

template<typename realT>
realT sourceFinder<realT>::minsep()
{
   return m_minsep;
}

template<typename realT>
void sourceFinder<realT>::useMedian(const bool & useMed)
{
   m_useMedian = useMed;
}

template<typename realT>
bool sourceFinder<realT>::useMedian()
{
   return m_useMedian;
}

template<typename realT>
std::vector<typename sourceFinder<realT>::pixelT> sourceFinder<realT>::operator()( const eigenImage<realT> & im )
{
   m_vals.resize(im.rows()*im.cols());
   m_pixs.resize(im.rows()*im.cols());

   int p =0;
   for(int cc=0;cc<im.cols();++cc)
   {
      for(int rr=0;rr<im.rows();++rr)
      {
         m_vals[p] = im(rr,cc);
         m_pixs[p] = {rr,cc};
         ++p;
      }
   }
   
   realT mm;

   if(m_useMedian) mm = math::vectorMedian(m_vals);
   else mm = math::vectorMean(m_vals);
   realT std = sqrt(math::vectorVariance(m_vals, mm));

   //Get the SNR for any pixels above the threshold
   m_snr.clear();
   for(size_t n=0; n < m_vals.size(); ++n)
   {
      if( (m_vals[n]-mm) / std >= m_threshold )
      {
         m_snr.push_back(n);
      }
   }

   //Delete the pixels above threshold which are with minsep of a higher pixel
   realT mms2 = m_minsep*m_minsep; //save an op
   for(long n =0; n < m_snr.size(); ++n)
   {
      for(long m=0; m < m_snr.size(); ++m)
      {
         if(m == n) continue;

         //Squared distance:
         realT d2 = pow( m_pixs[m_snr[n]].first - m_pixs[m_snr[m]].first,2) + pow( m_pixs[m_snr[n]].second - m_pixs[m_snr[m]].second,2);

         if(d2 < mms2)
         {
            if( m_vals[m_snr[n]] < m_vals[m_snr[m]] ) //If this is the higher pixel, delete the lower
            {
               m_snr.erase(m_snr.begin() + n);
               --n;
               break; //So we're done comparing to the one we just deleted.
            }
            else //This is lower, so delete it.
            {
               m_snr.erase(m_snr.begin() + m);
               if(m < n) --n; //changing size behind
               --m;
               continue;
            }
         }
      }
   }

   std::vector<pixelT> found(m_snr.size());

   for(size_t n = 0; n < m_snr.size(); ++n)
   {
      found[n] = m_pixs[m_snr[n]];
   }

   return found;
}

} //namespace improc 
} //namespace mx

#endif // improc_sourceFinder_hpp
