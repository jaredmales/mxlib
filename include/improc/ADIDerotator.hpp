/** \file ADIDerotator.hpp
  * \author Jared R. Males
  * \brief Defines a generic ADI derotator class.
  * \ingroup hc_imaging_files
  *
  */

#include "fitsHeader.hpp"

#ifndef __ADIDerotator_hpp__
#define __ADIDerotator_hpp__

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
   std::vector<std::string> keywords;
   
   ///Vector(s) to hold the keyword values
   std::vector<realT> angles;
   

   std::string _angleKeyword; ///<The keyword for the angle attribute.  Do not set this directly.
   
   ///Set the angle keyword
   /** Populates the kewords vector appropriately.
     */
   void angleKeyword( const std::string & akw /**< The angle keyword */)
   {
      _angleKeyword = akw;
      keywords = {akw};
   }
   
   realT angleScale; ///< The scale to multiply the angle by
   realT angleConstant; ///< The constant to add to the scaled-angle.
   
   ADIDerotator()
   {
      angleScale = 0;
      angleConstant = 0;
   }
   
   ///To allow ADIobservation to check for errors.
   bool isSetup()
   {
      if( (_angleKeyword == "" || keywords.size() == 0) || (angleScale == 0 && angleConstant == 0)) return false;
      return true;
   }
         
   ///Method called by ADIobservation to get keyword-values
   void extractKeywords(std::vector<fitsHeader> & heads)
   {
      angles = headersToValues<realT>(heads, _angleKeyword);
   }
   
   ///Calculate the derotation angle for a given image number
   realT derotAngle(size_t imno) const
   {
      realT derot = angleScale*angles[imno]+ angleConstant;
      //while(derot < 0) derot += 360.0;
      //derot = fmod(derot, 360.0);
      derot = math::angleMod(derot);
      return math::dtor(derot );
   }
};


///@}

} //namespace improc
} //namespace mx

#endif // __ADIDerotator_hpp__

