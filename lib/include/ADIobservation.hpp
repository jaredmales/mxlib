
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
#endif //__ADIobservation_hpp__


