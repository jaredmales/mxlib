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
      rotoff = headersToValues<floatT>(heads, "ROTOFF");
   }
   
   ///Calculate the derotation angle for a given image number
   floatT derotAngle(size_t imno) const
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
      rotoff = headersToValues<floatT>(heads, "ROTOFF");
   }
   
   ///Calculate the derotation angle for a given image number
   floatT derotAngle(size_t imno) const
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
  * \tparam floatT is the floating point type in which to do calculations
  * 
  * \tparam _derotFunctObj 
  * \parblock 
  * is the derotation functor with the following minimum interface: 
  * \code
  * template<typename _floatT>
  * struct derotF
  * {
  *    typedef  _floatT floatT;
  * 
  *    //Vector of keywords to extract from the fits headers
  *    std::vector<std::string> keywords;
  *    
  *    //Vector(s) to hold the keyword values
  *    std::vector<floatT> keyValue1;
  *    
  *    
  *    //Constructor should populate keywords
  *    derotVisAO()
  *    {
  *       keywords.push_back("KEYWORD1");
  *    }
  *    
  *    //Method called by HCIobservation to get keyword-values
  *    void extractKeywords(vector<fitsHeader> & heads)
  *    {
  *       keyValue1 = headersToValues<float>(heads, "KEYWORD1");
  *    }
  *    
  *    //Calculate the derotation angle for a given image number
  *    floatT derotAngle(size_t imno) const
  *    {
  *       return DTOR(keyValue1[imno]+90-0.6);
  *    }
  * };
  * \endcode
  * \endparblock
  */
template<typename _floatT, class _derotFunctObj>
struct ADIobservation : public HCIobservation<_floatT>
{
   typedef _floatT floatT;
   typedef _derotFunctObj derotFunctObj;
   typedef Array<floatT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   derotFunctObj derotF;
   //vector<floatT> derot;
      
   ADIobservation()
   {
      doFake = 0;
   }
   
   ADIobservation( const std::string & dir, 
                   const std::string & prefix, 
                   const std::string & ext) : HCIobservation<floatT>(dir,prefix,ext)
   {
      doFake = 0;
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
      
      if(doFake) injectFake();
      
   }
   
   int doFake;
   std::string fakeFileName;
   floatT fakeSep;
   floatT fakePA;
   floatT fakeContrast;
   
   
   
   void injectFake()
   {
      typedef Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> imT;
      imT fakePSF;
      fitsFile<floatT> ff;
      
      ff.read(fakeFileName, fakePSF);

      //Check for correct sizing
      if( (fakePSF.rows() < this->imc.rows() && fakePSF.cols() >= this->imc.cols()) || 
                           (fakePSF.rows() >= this->imc.rows() && fakePSF.cols() < this->imc.cols()))
      {
         throw mxException("mxlib:high contrast imaging", -1, "image wrong size",  __FILE__, __LINE__, "fake PSF has different dimensions and can't be sized properly");
      }
      
      //Check if fake needs to be padded out
      if(fakePSF.rows() < this->imc.rows() && fakePSF.cols() < this->imc.cols())
      {
         imT pfake(this->imc.rows(), this->imc.cols());
         padImage(pfake, fakePSF, this->imc.rows(), this->imc.cols());
         fakePSF = pfake;
      }
      
      //Check if fake needs to be cut down
      if(fakePSF.rows() > this->imc.rows() && fakePSF.cols() > this->imc.cols())
      {
         imT cfake(this->imc.rows(), this->imc.cols());
         cutImage(cfake, fakePSF, this->imc.rows(), this->imc.cols());
         fakePSF = cfake;
      }
      
      
      //allocate shifted fake psf
      imT shiftFake(fakePSF.rows(), fakePSF.cols());
      
      floatT ang, dx, dy;
      
      for(int i=0; i<this->imc.planes(); ++i)
      {

         ang = DTOR(-fakePA) + derotF.derotAngle(i);
         
         dx = fakeSep * sin(ang);
         dy = fakeSep * cos(ang);
                  
         imageShift(shiftFake, fakePSF, dx, dy, cubicConvolTransform<floatT>());
      
         this->imc.image(i) = this->imc.image(i) + shiftFake*fakeContrast;
      }
      
      pout("fake injected");
   }
   
   
   
   void derotate()
   {
      
      
      #pragma omp parallel for num_threads(4) schedule(static, 1)
      eigenImageT rotim;
      floatT derot;
      for(int n=0; n<this->psfsub.size(); ++n)
      {
         for(int i=0; i<this->psfsub[n].planes();++i)
         {
            derot = derotF.derotAngle(i);
            if(derot != 0) 
            {
               imageRotate(rotim, this->psfsub[n].image(i), derot, cubicConvolTransform<floatT>());
               this->psfsub[n].image(i) = rotim;
            }
         }
      }
        
   }
   
};

///@}

} //namespace mx

#endif //__ADIobservation_hpp__


