/** \file ADIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the ADI high contrast imaging data type.
  * \ingroup hc_imaging
  *
  */

#include "HCIobservation.hpp"

#ifndef __ADIobservation_hpp__
#define __ADIobservation_hpp__

namespace mx
{
   
template<typename _floatT>
struct derotVisAO
{
   typedef  _floatT floatT;

   ///Vector of keywords to extract from the fits headers
   std::vector<std::string> keywords;
   
   ///Vector(s) to hold the keyword values
   std::vector<floatT> rotoff;
   
   
   ///Constructor should populate keywords
   derotVisAO()
   {
      keywords.push_back("ROTOFF");
   }
   
   ///Method called by DIobservation to get keyword-values
   void extractKeywords(vector<fitsHeader> & heads)
   {
      rotoff = headersToValues<float>(heads, "ROTOFF");
   }
   
   ///Calculate the derotation angle for a given image number
   floatT derotAngle(size_t imno)
   {
      return DTOR(rotoff[imno]+90-0.6);
   }
};

template<typename _floatT>
struct derotClio
{
   typedef  _floatT floatT;

   ///Vector of keywords to extract from the fits headers
   std::vector<std::string> keywords;
   
   ///Vector(s) to hold the keyword values
   std::vector<floatT> rotoff;
   
   
   ///Constructor should populate keywords
   derotClio()
   {
      keywords.push_back("ROTOFF");
   }
   
   ///Method called by DIobservation to get keyword-values
   void extractKeywords(vector<fitsHeader> & heads)
   {
      rotoff = headersToValues<float>(heads, "ROTOFF");
   }
   
   ///Calculate the derotation angle for a given image number
   floatT derotAngle(size_t imno)
   {
      return DTOR(rotoff[imno]-180-1.8);
   }
};

template<typename _floatT>
struct derotODI
{
   typedef  _floatT floatT;

   ///Vector of keywords to extract from the fits headers
   std::vector<std::string> keywords;
   
   ///Vector(s) to hold the keyword values
   std::vector<std::string> rotoff;
   
   
   ///Constructor should populate keywords
   derotODI()
   {
      //no keywords
   }
   
   ///Method called by DIobservation to get keyword-values
   void extractKeywords(vector<fitsHeader> & heads)
   {
//       rotoff.resize(heads.size());
//       for(int i=0; i<heads.size(); ++i) rotoff[i] = i;
   }
   
   ///Calculate the derotation angle for a given image number
   floatT derotAngle(size_t imno)
   {
      return 0;
   }
};

/** \addtogroup hc_imaging
 * @{
 */

///Process an angular differential imaging (ADI) observation
/** Angular differential imaging (ADI) uses sky rotation to differentiate real objects from
  * speckles.
  * 
  * derotFunctObj is a function object with the following minimum interface
  * \code
  * template<typename _arithT>
  * struct derotF
  * {
  *    typedef  _arithT arithT;
  *    std::string angKeyWord;  //This is used to extract the rotation angle from the header
  * 
  *    //Constructor must initialize angKeyWord
  *    derotF()
  *    {
  *       angKeyWord = "ROTOFF";
  *    }
  * 
  *    //operator() calculates the counter-clockwise angle to de-rotate by.
  *    arithT operator()(arithT rotoff)
  *    {
  *       return rotoff+90-0.6;
  *    } 
  * };
  * \endcode
  * 
  * \tparam floatT is the floating point type in which to do calculations
  * \tparam _derotFunctObj is the derotation functor with the interface described above
  */
template<typename _floatT, class _derotFunctObj>
struct ADIobservation : public HCIobservation<_floatT>
{
   typedef _floatT floatT;
   typedef _derotFunctObj derotFunctObj;
   
   derotFunctObj derotF;
   //vector<floatT> derot;
      
   ADIobservation()
   {
   }
   
   ADIobservation(std::string odir, std::string oprefix, std::string oext) : HCIobservation<floatT>(odir,oprefix,oext)
   {
   }
   
   void readFiles()
   {      
      this->keywords.clear();
      for(int i=0;i<derotF.keywords.size();++i)
      {
         this->keywords.push_back(derotF.keywords[i]);
      }
      
      HCIobservation<floatT>::readFiles();
      
      derotF.extractKeywords(this->heads);
      
   }
   
   
   void derotate()
   {
      eigenImagef rotim;
      floatT derot;

      for(int i=0; i<this->imc.planes();i++)
      {
         derot = derotF.derotAngle(i);
         if(derot == 0) 
         {
            this->imc.image(i) = this->psfsub.image(i);
         }
         else
         {
            imageRotate(rotim, this->psfsub.image(i), derot, cubicConvolTransform<floatT>());
            this->imc.image(i) = rotim;
         }
      }
        
   }
   
};

///@}

} //namespace mx

#endif //__ADIobservation_hpp__


