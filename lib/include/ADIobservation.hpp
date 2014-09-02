/** \file ADIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the ADI high contrast imaging data type.
  * \ingroup hc_imaging
  *
  */

#include "HCIobservation.hpp"

#ifndef __ADIobservation_hpp__
#define __ADIobservation_hpp__

template<typename _arithT>
struct derotVisAO
{
   typedef  _arithT arithT;

   std::string angKeyWord;
   
   derotVisAO()
   {
      angKeyWord = "ROTOFF";
   }
   
   arithT operator()(arithT rotoff)
   {
      return rotoff+90-0.6;
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
   vector<floatT> derot;
      
   ADIobservation()
   {
      this->keywords.push_back(derotF.angKeyWord);
   }
   
   ADIobservation(std::string odir, std::string oprefix, std::string oext) : HCIobservation<floatT>(odir,oprefix,oext)
   {
      this->keywords.push_back(derotF.angKeyWord);
   }
   
   void readFiles()
   {      
      HCIobservation<floatT>::readFiles();
      
      derot = headersToValues<float>(this->heads, derotF.angKeyWord);
   
      for(int i=0;i<derot.size(); ++i) derot[i] = derotF(derot[i]);
   }
   
   
   void derotate()
   {
      eigenImagef rotim;

      for(int i=0; i<this->imc.planes();i++)
      {
         imageRotate(rotim, this->psfsub.image(i), DTOR(derot[i]), cubicConvolTransform<float>());
         this->imc.image(i) = rotim;
      }
        
   }
   
};

///@}

#endif //__ADIobservation_hpp__


