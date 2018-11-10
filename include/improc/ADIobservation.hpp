/** \file ADIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the ADI high contrast imaging data type.
  * \ingroup hc_imaging_files
  * \ingroup image_processing_files
  *
  */

#include "HCIobservation.hpp"
#include "fitsHeader.hpp"

#include "imagePads.hpp"

#ifndef __ADIobservation_hpp__
#define __ADIobservation_hpp__

namespace mx
{
namespace improc
{
   
#if 0
   
 VisAO:     return dtor(rotoff[imno]+90-0.6);
   
 Clio:      return DTOR(rotoff[imno]-180-1.8);
   

template<typename _realT>
struct derotODI
{
   typedef  _realT realT;

   ///Vector of keywords to extract from the fits headers
   std::vector<std::string> keywords;
   
   ///Vector(s) to hold the keyword values
   std::vector<realT> dateobs;
   
   ///The period of the orbit
   _realT period;
   
   ///Constructor should populate keywords
   derotODI()
   {
      period = 365.25;
      keywords.push_back("DATEOBS");
   }
   

   ///Method called by DIobservation to get keyword-values
   void extractKeywords(std::vector<fitsHeader> & heads)
   {
      dateobs = headersToValues<realT>(heads, "DATEOBS");
   }
   
   ///Calculate the derotation angle for a given image number
   realT derotAngle(size_t imno) const
   {
      return D2PI-(fmod(dateobs[imno], period)/period)*D2PI;
   }
};
#endif

///Process an angular differential imaging (ADI) observation
/** Angular differential imaging (ADI) uses sky rotation to differentiate real objects from
  * speckles.
  * 
  * \tparam realT is the floating point type in which to do calculations
  * 
  * \tparam _derotFunctObj 
  * \parblock 
  * is the derotation object with the following minimum interface: 
  * \code
  * template<typename _realT>
  * struct derotF
  * {
  *    typedef  _realT realT;
  * 
  *    //Vector of keywords to extract from the fits headers
  *    std::vector<std::string> keywords;
  *    
  *    //Vector(s) to hold the keyword values
  *    std::vector<realT> keyValue1;
  *    
  *    ///To allow ADIobservation to check for errors.
  *    bool isSetup()
  *    {
  *      if( <any condition indicatint not set up>) return false;
  *      return true;
  *    }
  * 
  *    //Method called by HCIobservation to get keyword-values
  *    void extractKeywords(vector<fitsHeader> & heads)
  *    {
  *       keyValue1 = headersToValues<float>(heads, "KEYWORD1");
  *    }
  *    
  *    //Calculate the derotation angle for a given image number
  *    realT derotAngle(size_t imno) const
  *    {
  *       return some_function_of(keyValue1[imno]); //This function uses keyValue1[imno] to produce the derotation angle in radians.
  *    }
  * };
  * \endcode
  * \endparblock
  * 
  * \ingroup hc_imaging
  */
template<typename _realT, class _derotFunctObj>
struct ADIobservation : public HCIobservation<_realT>
{
   typedef _realT realT;
   typedef _derotFunctObj derotFunctObj;
   typedef Array<realT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   derotFunctObj derotF;

   bool doDerotate;
   
   ADIobservation();
   
   explicit ADIobservation( const std::string & fileListFile) ;
                   
   ADIobservation( const std::string & dir, 
                   const std::string & prefix, 
                   const std::string & ext = ".fits") ;

   void initialize();
   
   
   /** \name Rotation Setup
     * Configuration of the rotation system.
     * @{ 
     */
   
   std::string angleKeyword;
   
   realT angleScale;
   realT angleConstant;
   
   ///@}
   
   ///Read in the files
   /** First sets up the keywords, then calls HCIobservation readFiles
     */
   int readFiles();
   
   ///Post read actions, including fake injection
   virtual int postReadFiles();
   
   /// Read in already PSF-subtracted files
   /** Used to take up final processing after applying some non-klipReduce processing steps to
     * PSF-subtracted images.
     */ 
   int readPSFSub( const std::string & dir,
                   const std::string & prefix,
                   const std::string & ext,
                   size_t nReductions 
                 );
   
   /** \name Fake Planets
     * @{ 
     */
//   int doFake; ///<Flag controlling whether or not fake planets are injected
   int fakeMethod; // 0 == 1 file total, 1 == 1 file per image 
   std::string fakeFileName; ///<FITS file containing the fake planet PSF to inject or a list of fake images
   
//   bool doFakeScale; ///<Flag controlling whether or not a separate scale is used at each point in time
   std::string fakeScaleFileName; ///< One-column text file containing a scale factor for each point in time.
   
   std::vector<realT> fakeSep; ///< Separation(s) of the fake planet(s)
   std::vector<realT> fakePA; ///< Position angles(s) of the fake planet(s)
   std::vector<realT> fakeContrast; ///< Contrast(s) of the fake planet(s)
   

   
   ///Inect the fake plants
   int injectFake();
   
   int injectFake( eigenImageT & fakePSF,
                    int image_i,
                    realT PA,
                    realT sep,
                    realT contrast,
                    realT scale = 1.0 );

   /// @}
   
   void fitsHeader(fitsHeader * head);
   
   virtual void makeMaskCube();
   
   ///De-rotate the PSF subtracted images
   void derotate();

   double t_fake_begin;
   double t_fake_end;
   
   double t_derotate_begin;
   double t_derotate_end;
   
};

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation()
{
   initialize();
}

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation( const std::string & fileListFile) : HCIobservation<realT>(fileListFile)
{
   initialize();
}

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation( const std::string & dir, 
                                                         const std::string & prefix, 
                                                         const std::string & ext) : HCIobservation<realT>(dir,prefix,ext)
{
   initialize();
}

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::initialize()
{
   doDerotate = true;
   
   fakeMethod = 0;
   
   t_fake_begin = 0;
   t_fake_end = 0;
   
   t_derotate_begin = 0;
   t_derotate_end = 0;
   
//    derotF.angleKeyword("ROTOFF");
//    derotF.angleScale = 1.0;
//    derotF.angleConstant = 89.4;
   
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::readFiles()
{      
   this->keywords.clear();
   
   if(!derotF.isSetup())
   {
      mxError("ADIobservation::readFiles", MXE_PARAMNOTSET, "Derotator is not configured.");
      return -1;
   }
   
   for(size_t i=0;i<derotF.keywords.size();++i)
   {
      this->keywords.push_back(derotF.keywords[i]);
   }
   
   /*----- Append the ADI keywords to propagate them if needed -----*/
      
   if( HCIobservation<realT>::readFiles() < 0) return -1;
   
   derotF.extractKeywords(this->heads);
   
   return 0;
   
   /*---- Check for ADI keywords -----*/
   //Maybe fill in a structure of values, in case things are overwritten by new settings
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::postReadFiles()
{
   derotF.extractKeywords(this->heads);
   
   if(fakeFileName != ""  && !this->skipPreProcess) 
   {
      if( injectFake() < 0) return -1;
   }
   
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::readPSFSub( const std::string & dir,
                                                        const std::string & prefix,
                                                        const std::string & ext,
                                                        size_t nReductions 
                                                      )
{
   this->keywords.clear();
   
   if(!derotF.isSetup())
   {
      mxError("ADIobservation::readFiles", MXE_PARAMNOTSET, "Derotator is not configured.");
      return -1;
   }
   
   for(size_t i=0;i<derotF.keywords.size();++i)
   {
      this->keywords.push_back(derotF.keywords[i]);
   }
   
   /*----- Append the ADI keywords to propagate them if needed -----*/
      
   if( HCIobservation<realT>::readPSFSub(dir,prefix,ext, nReductions) < 0) return -1;
   
   derotF.extractKeywords(this->heads);
   
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::injectFake()
{
   std::cerr << "injecting fake planets\n";

   t_fake_begin = get_curr_time();
   
   //typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imT;
   eigenImageT fakePSF;
   std::vector<std::string> fakeFiles; //used if fakeMethod == 1
   
   fitsFile<realT> ff;
   std::ifstream scaleFin; //for reading the scale file.
      

   //Fake Scale -- default to 1, read from file otherwise
   std::vector<realT> fakeScale(this->imc.planes(), 1.0);
   if(fakeScaleFileName != "")
   {      
      std::vector<std::string> sfileNames;
      std::vector<realT> imS;
      
      //Read the quality file and load it into a map
      if( ioutils::readColumns(fakeScaleFileName, sfileNames, imS) < 0) return -1;
      
      std::map<std::string, realT> scales;     
      for(size_t i=0;i<sfileNames.size();++i) scales[basename(sfileNames[i])] = imS[i];
      
      for(size_t i=0; i<this->fileList.size(); ++i)
      {
         if(scales.count(basename(this->fileList[i].c_str())) > 0)
         {
            fakeScale[i] = scales[basename(this->fileList[i].c_str())];
         }
         else
         {
            std::cerr << "File name not found in fakeScaleFile:\n";
            std::cerr << basename(this->fileList[i].c_str()) << "\n";
            exit(-1);
         }
      }
   } //if(fakeScaleFileName != "")
      
   if(fakeMethod == 0)
   {
      if( ff.read( fakePSF, fakeFileName ) < 0) return -1;
   }

   if(fakeMethod == 1)
   {
      if( ioutils::readColumns(fakeFileName, fakeFiles) < 0) return -1;
   }
   
   for(int i=0; i<this->imc.planes(); ++i)
   {
      if(fakeMethod == 1)
      {
         ff.read(fakePSF, fakeFiles[i]);
      }

      for(size_t j=0;j<fakeSep.size(); ++j)
      {
         if( injectFake(fakePSF, i, fakePA[j], fakeSep[j], fakeContrast[j], fakeScale[j]) < 0) return -1;
      }
   }
   
   
   std::cerr << "fake injected\n";
   
   t_fake_end = get_curr_time();
   
   return 0;
}


template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::injectFake( eigenImageT & fakePSF,
                                                          int image_i,
                                                          _realT PA,
                                                          _realT sep,
                                                          _realT contrast,
                                                          _realT scale)
{
   //typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imT;

   //Check for correct sizing
   if( (fakePSF.rows() < this->imc.rows() && fakePSF.cols() >= this->imc.cols()) || 
                        (fakePSF.rows() >= this->imc.rows() && fakePSF.cols() < this->imc.cols()))
   {
      throw mxException("mxlib:high contrast imaging", -1, "image wrong size",  __FILE__, __LINE__, "fake PSF has different dimensions and can't be sized properly");
   }
   
   //Check if fake needs to be padded out
   if(fakePSF.rows() < this->imc.rows() && fakePSF.cols() < this->imc.cols())
   {
      eigenImageT pfake(this->imc.rows(), this->imc.cols());
      padImage(pfake, fakePSF, 0.5*(this->imc.rows()-fakePSF.rows()), 0);
      fakePSF = pfake;
   }
   
   //Check if fake needs to be cut down
   if(fakePSF.rows() > this->imc.rows() && fakePSF.cols() > this->imc.cols())
   {
      eigenImageT cfake(this->imc.rows(), this->imc.cols());
      cutPaddedImage(cfake, fakePSF, 0.5*(fakePSF.rows() - this->imc.rows()));
      fakePSF = cfake;
   }

   /*** Now shift to the separation and PA, scale, apply contrast, and inject ***/
   //allocate shifted fake psf
   eigenImageT shiftFake(fakePSF.rows(), fakePSF.cols());
   
   realT ang, dx, dy;

   ang = math::dtor(-1*PA) + derotF.derotAngle(image_i);
      
   dx = sep * sin(ang);
   dy = sep * cos(ang);
               
   imageShift(shiftFake, fakePSF, dx, dy, cubicConvolTransform<realT>());
   
   this->imc.image(image_i) = this->imc.image(image_i) + shiftFake*scale*contrast;

   return 0;
   
}

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::makeMaskCube()
{
   this->maskCube.resize( this->Nrows, this->Ncols, this->Nims);
   eigenImageT rm;
   
   for(int i=0; i< this->Nims; ++i)
   {
      rotateMask( rm, this->mask, derotF.derotAngle(i));
      
      this->maskCube.image(i) = rm;
   }
   
   fitsFile<realT> ff; 
   ff.write("maskCube.fits", this->maskCube);
}

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::derotate()
{
   t_derotate_begin = get_curr_time();
   
   //On magaoarx it doesn't seem worth it to use more than 4 threads
   #pragma omp parallel num_threads(4)
   {
      eigenImageT rotim;
      realT derot;
      
      #pragma omp for schedule(static, 1)
      for(size_t n=0; n<this->psfsub.size(); ++n)
      {
         for(int i=0; i<this->psfsub[n].planes();++i)
         {
            derot = derotF.derotAngle(i);
            if(derot != 0) 
            {
               imageRotate(rotim, this->psfsub[n].image(i), derot, cubicConvolTransform<realT>());
               this->psfsub[n].image(i) = rotim;
            }
         }
      }
   }
   
   t_derotate_end = get_curr_time();
}


//If fakeFileName == "" or skipPreProcess == true then use the structure of propagated values

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::fitsHeader( mx::improc::fitsHeader * head)
{
   if(head == 0) return;
   
   head->append("", fitsCommentType(), "----------------------------------------");
   head->append("", fitsCommentType(), "mx::ADIobservation parameters:");
   head->append("", fitsCommentType(), "----------------------------------------");

   if(fakeFileName != "")
   head->append("FAKEFILE", fakeFileName, "name of fake planet PSF file");
   
   if(fakeScaleFileName != "")
   head->append("FAKESCFL", fakeScaleFileName, "name of fake planet scale file name");

   std::stringstream str;
   
   if(fakeSep.size() > 0)
   {
      for(size_t nm=0;nm < fakeSep.size()-1; ++nm) str << fakeSep[nm] << ",";
      str << fakeSep[fakeSep.size()-1];      
      head->append<char *>("FAKESEP", (char *)str.str().c_str(), "separation of fake planets");
   }
   
   if(fakePA.size() > 0 )
   {
      str.str("");
      for(size_t nm=0;nm < fakePA.size()-1; ++nm) str << fakePA[nm] << ",";
      str << fakePA[fakePA.size()-1];      
      head->append<char *>("FAKEPA", (char *)str.str().c_str(), "PA of fake planets");
   }
   
   if( fakeContrast.size() > 0)
   {
      str.str("");
      for(size_t nm=0;nm < fakeContrast.size()-1; ++nm) str << fakeContrast[nm] << ",";
      str << fakeContrast[fakeContrast.size()-1];      
      head->append<char *>("FAKECONT", (char *)str.str().c_str(), "Contrast of fake planets");
   }
}



} //namespace improc
} //namespace mx

#endif //__ADIobservation_hpp__


