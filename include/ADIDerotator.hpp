/** \file ADIDerotator.hpp
  * \author Jared R. Males
  * \brief Defines the ADI derotator base class.
  * \ingroup hc_imaging
  *
  */

#include "fitsHeader.hpp"

#ifndef __ADIDerotator_hpp__
#define __ADIDerotator_hpp__

namespace mx
{
   
/** \addtogroup hc_imaging
  * @{
  */

template<typename _floatT>
struct ADIDerotator
{
   typedef  _floatT floatT;

   ///Vector of keywords to extract from the fits headers
   std::vector<std::string> keywords;
   
   ///Vector(s) to hold the keyword values
   //std::vector<floatT> <<<<derived class has something here>>>>;
   
   ///Constructor should populate keywords
   ADIderotator()
   {
      //keywords.push_back(<<<<derived class has something here>>>>);
   }
   
   ///Method called by HCIobservation to get keyword-values
   virtual void extractKeywords(std::vector<fitsHeader> & heads)
   {
      //<<<<derived class converts fits headers to vector here>>>>
   }
   
   ///Calculate the derotation angle for a given image number
   virtual floatT derotAngle(size_t imno) const
   {
      //<<<<derived class converts to derotation angle here>>>>
   }
};

template<typename _floatT>
struct VisAODerotator : public ADIDerotator<_floatT>
{
   
   ///Vector(s) to hold the keyword values
   std::vector<floatT> rotoff;
   
   ///Constructor should populate keywords
   VisAOderotator()
   {
      keywords.push_back("ROTOFF");
   }
   
   ///Method called by HCIobservation to get keyword-values
   virtual void extractKeywords(std::vector<fitsHeader> & heads)
   {
      rotoff = headersToValues<floatT>(heads, "ROTOFF");
   }
   
   ///Calculate the derotation angle for a given image number
   virtual floatT derotAngle(size_t imno) const
   {
      return DTOR(rotoff[imno]+90-0.6);
   }
};


struct ClioDerotator : public ADIDerotator<_floatT>
{
   
   ///Vector(s) to hold the keyword values
   std::vector<floatT> rotoff;
   
   ///Constructor should populate keywords
   VisAOderotator()
   {
      keywords.push_back("ROTOFF");
   }
   
   ///Method called by HCIobservation to get keyword-values
   virtual void extractKeywords(std::vector<fitsHeader> & heads)
   {
      rotoff = headersToValues<floatT>(heads, "ROTOFF");
   }
   
   ///Calculate the derotation angle for a given image number
   virtual floatT derotAngle(size_t imno) const
   {
      return DTOR(rotoff[imno]-180-1.8);
   }
};

///@}

} //namespace mx

#endif // __ADIDerotator_hpp__

