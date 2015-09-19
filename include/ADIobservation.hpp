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
   std::vector<floatT> dateobs;
   
   ///The period of the orbit
   _floatT period;
   
   ///Constructor should populate keywords
   derotODI()
   {
      period = 365.25;
      keywords.push_back("DATEOBS");
   }
   
   ///Method called by DIobservation to get keyword-values
   void extractKeywords(vector<fitsHeader> & heads)
   {
      dateobs = headersToValues<floatT>(heads, "DATEOBS");
   }
   
   ///Calculate the derotation angle for a given image number
   floatT derotAngle(size_t imno) const
   {
      return D2PI-(fmod(dateobs[imno], period)/period)*D2PI;
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
      
   bool doDerotate;
   
   ADIobservation();
   
   ADIobservation( const std::string & dir, 
                   const std::string & prefix, 
                   const std::string & ext) ;

   ///Read in the files
   /** First sets up the keywords, then calls HCIobservation readFiles
     */
   void readFiles();
   
   ///Post read actions, including fake injection
   virtual void postReadFiles();
   
   /** \name Fake Planets
     * @{ 
     */
   int doFake; ///<Flag controlling whether or not fake planets are injected
   std::string fakeFileName; ///<FITS file containing the fake planet PSF to inject
   bool doFakeScale; ///<Flag controlling whether or not a separate scale is used at each point in time
   std::string fakeScaleFileName; ///< One-column text file containing a scale factor for each point in time.
   
   std::vector<floatT> fakeSep; ///< Separation(s) of the fake planet(s)
   std::vector<floatT> fakePA; ///< Position angles(s) of the fake planet(s)
   std::vector<floatT> fakeContrast; ///< Contrast(s) of the fake planet(s)
   
   ///Inect the fake plants
   void injectFake();
   
   /// @}
   
   ///De-rotate the PSF subtracted images
   void derotate();

   
};

template<typename _floatT, class _derotFunctObj>
ADIobservation<_floatT, _derotFunctObj>::ADIobservation()
{
   doDerotate = true;
   doFake = 0;
   doFakeScale = 0;
}

template<typename _floatT, class _derotFunctObj>
ADIobservation<_floatT, _derotFunctObj>::ADIobservation( const std::string & dir, 
                                                         const std::string & prefix, 
                                                         const std::string & ext) : HCIobservation<floatT>(dir,prefix,ext)
{
   doDerotate = true;
   doFake = 0;
   doFakeScale = 0;
}

template<typename _floatT, class _derotFunctObj>
void ADIobservation<_floatT, _derotFunctObj>::readFiles()
{      
   this->keywords.clear();
   for(int i=0;i<derotF.keywords.size();++i)
   {
      this->keywords.push_back(derotF.keywords[i]);
   }
   
   HCIobservation<floatT>::readFiles();
   
}

template<typename _floatT, class _derotFunctObj>
void ADIobservation<_floatT, _derotFunctObj>::postReadFiles()
{
   derotF.extractKeywords(this->heads);
   
   if(doFake) injectFake();
}

template<typename _floatT, class _derotFunctObj>
void ADIobservation<_floatT, _derotFunctObj>::injectFake()
{
   std::cout << "injecting fake planets\n";
   
   typedef Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> imT;
   imT fakePSF;
   fitsFile<floatT> ff;
   std::ifstream scaleFin; //for reading the scale file.
      
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

   //Fake Scale -- default to 1, read from file otherwise
   std::vector<floatT> fakeScale(this->imc.planes(), 1.0);
   if(doFakeScale)
   {
      scaleFin.open(fakeScaleFileName.c_str());
      for(int i=0; i<this->imc.planes(); ++i)
      {
         scaleFin >> fakeScale[i];
      }
      scaleFin.close();
   }
      
      
   for(int i=0; i<this->imc.planes(); ++i)
   {
      for(int j=0;j<fakeSep.size(); ++j)
      {
         ang = DTOR(-fakePA[j]) + derotF.derotAngle(i);
      
         dx = fakeSep[j] * sin(ang);
         dy = fakeSep[j] * cos(ang);
               
         imageShift(shiftFake, fakePSF, dx, dy, cubicConvolTransform<floatT>());
   
         this->imc.image(i) = this->imc.image(i) + shiftFake*fakeScale[i]*fakeContrast[j];
      }
   }
   
   
   pout("fake injected");
}


template<typename _floatT, class _derotFunctObj>
void ADIobservation<_floatT, _derotFunctObj>::derotate()
{
   //On magaoarx it doesn't seem worth it to use more than 4 threads
   #pragma omp parallel num_threads(4)
   {
      eigenImageT rotim;
      floatT derot;
      #pragma omp for schedule(static, 1)
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
     
}
///@}

} //namespace mx

#endif //__ADIobservation_hpp__


