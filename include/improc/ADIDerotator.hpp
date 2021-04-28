/** \file ADIDerotator.hpp
  * \author Jared R. Males
  * \brief Defines a generic ADI derotator class.
  * \ingroup hc_imaging_files
  *
  */

#include "../ioutils/fits/fitsHeader.hpp"

#ifndef ADIDerotator_hpp
#define ADIDerotator_hpp

#include "../math/geo.hpp"

namespace mx
{
   
namespace improc
{

///A generic ADI derotator class.  
/** This class is used to calculate the derotation angle for angular differential imaging.
  * 
  * \ingroup hc_imaging
  * 
  */
template<typename _realT>
struct ADIDerotator
{
   typedef  _realT realT;

   ///Vector of keywords to extract from the fits headers
   std::vector<std::string> m_keywords;
   
   ///Vector(s) to hold the keyword values
   std::vector<realT> m_angles;
   

   std::string m_angleKeyword; ///<The keyword for the angle attribute.  Do not set this directly.
   
   ///Set the angle keyword
   /** Populates the kewords vector appropriately.
     */
   void angleKeyword( const std::string & akw /**< [in] The angle keyword */)
   {
      m_angleKeyword = akw;
      m_keywords = {akw};
   }
   
   realT m_angleScale {0}; ///< The scale to multiply the angle by
   realT m_angleConstant {0}; ///< The constant to add to the scaled-angle.
   
   ADIDerotator()
   {
   }
   
   ///To allow ADIobservation to check for errors.
   bool isSetup()
   {
      if( (m_angleKeyword == "" || m_keywords.size() == 0) || (m_angleScale == 0 && m_angleConstant == 0)) return false;
      return true;
   }
         
   ///Method called by ADIobservation to get keyword-values
   void extractKeywords(std::vector<fits::fitsHeader> & heads /**< [in] The headers from the images being reduced.*/)
   {
      m_angles = fits::headersToValues<realT>(heads, m_angleKeyword);
   }
   
   ///Calculate the derotation angle for a given image number
   /**
     * \returns the angle in radians by which to de-rotate the image c.c.w.
     */ 
   realT derotAngle(size_t imno /**< [in] the image number */) const
   {
      realT derot = m_angleScale*m_angles[imno] + m_angleConstant;
      derot = math::angleMod(derot);
      return math::dtor(derot );
   }
};


///@}

extern template struct ADIDerotator<float>;
extern template struct ADIDerotator<double>;

} //namespace improc
} //namespace mx

#endif // ADIDerotator_hpp

