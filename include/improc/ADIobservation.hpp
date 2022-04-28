/** \file ADIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the ADI high contrast imaging data type.
  * \ingroup hc_imaging_files
  * \ingroup image_processing_files
  *
  */

#ifndef __ADIobservation_hpp__
#define __ADIobservation_hpp__
#include "../ioutils/fileUtils.hpp"

#include "HCIobservation.hpp"
#include "../ioutils/fits/fitsHeader.hpp"

#include "imagePads.hpp"



namespace mx
{
namespace improc
{
namespace HCI
{
   /// Fake injection PSF file specification methods
   /** \ingroup hc_imaging_enums
     */
   enum fakeMethods{ single, ///< A single PSF is used
                     list ///< A list of PSF files, one per input image, is used.
                   };
                
   /// Get the string name of a fake injection method 
   /**
     * \returns the string name corresponding to the fake injection method
     */ 
   std::string fakeMethodsStr( int method /**< [in] the fake injection method */);
   
   /// Get the fake injection method from its string name
   /**
     * \returns the corresponding member of the fakeMethods enum
     */ 
   int fakeMethodFmStr( const std::string & method  /**< [in] the fake injection method name*/);

}

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
  *      if( <any condition indicating not set up>) return false;
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
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   derotFunctObj m_derotF;

   derotFunctObj m_RDIderotF;
   
   bool m_doDerotate {true};
   
   bool m_postMedSub {false};
   
   ADIobservation();
   
   ADIobservation( const std::string &dir,     ///< [in] the directory to search.
                   const std::string &prefix,  ///< [in] the initial part of the file name.  Can be empty "".
                   const std::string &ext      ///< [in] the extension to append to the file name, must include the '.'.
                 );
   
   explicit ADIobservation( const std::string & fileListFile /**< [in] a file name path to read.*/);

   ADIobservation( const std::string &dir,       ///< [in] the directory to search.
                   const std::string &prefix,    ///< [in] the initial part of the file name.  Can be empty "".
                   const std::string &ext,       ///< [in] the extension to append to the file name, must include the '.'.
                   const std::string &RDIdir,    ///< [in] the directory to search for the reference files.
                   const std::string &RDIprefix, ///< [in] the initial part of the file name for the reference files.  Can be empty "".
                   const std::string &RDIext=""  ///< [in] [optional] the extension to append to the RDI file name, must include the '.'.  If empty "" then same extension as target files is used.
                 );
   
   ADIobservation( const std::string & fileListFile,   ///< [in] a file name path to read for the target file names.
                   const std::string & RDIfileListFile ///< [in] a file name path to read for the reference file names.
                 );
   
   
   ///Read in the target files
   /** First sets up the keywords, then calls HCIobservation readFiles
     */
   int readFiles();
   
   ///Post target read actions, including fake injection
   virtual int postReadFiles();
   
   ///Post target coadd actions.
   /** Here updates derotation for new average values.
     */
   virtual int postCoadd();
   
   ///Read in the RDI files
   /** First sets up the keywords, then calls HCIobservation readRDIFiles
     */
   int readRDIFiles();
   
   ///Post reference read actions, including fake injection
   virtual int postRDIReadFiles();
   
   ///Post reference coadd actions.
   /** Here updates derotation for new average values.
     */
   virtual int postRDICoadd();
   
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
   int m_fakeMethod {HCI::single}; ///< Method for reading fake files, either HCI::single or HCI::list.
   
   std::string m_fakeFileName; ///<FITS file containing the fake planet PSF to inject or a list of fake images
   
   std::string m_fakeScaleFileName; ///< One-column text file containing a scale factor for each point in time.
   
   std::vector<realT> m_fakeSep; ///< Separation(s) of the fake planet(s)
   std::vector<realT> m_fakePA; ///< Position angles(s) of the fake planet(s)
   std::vector<realT> m_fakeContrast; ///< Contrast(s) of the fake planet(s)
   
   realT m_RDIFluxScale {1}; ///< Flux scaling to apply to fake planets injected in RDI.  Would depend on the assumed spectrum in SDI.
   realT m_RDISepScale {1}; ///< Scaling to apply to fake planet separation in RDI.  Would be ratio of wavelengths for SDI.

   
   ///Inect the fake plants
   int injectFake( eigenCube<realT> & ims,              ///< [in/out] the image cube in which to inject the fakes.
                   std::vector<std::string> & fileList, ///< [in] a list of file paths used for per-image fake PSFs.  If empty, then m_fakeFileName is used.
                   derotFunctObj & derotF,
                   realT RDIfluxScale,                  ///< [in] the flux scaling for RDI.  In SDI, this is from the planet spectrum.
                   realT RDISepScale                    ///< [in] the separation scale for RDI.  In SDI, this is the ratio of wavlengths after lambda/D scaling.
                 );
   
   int injectFake( eigenImageT & fakePSF,
                   eigenCube<realT> & ims,
                   int image_i,
                   realT derotAngle,
                   realT PA,
                   realT sep,
                   realT contrast,
                   realT scale,
                   realT RDIfluxScale,
                   realT RDISepScale
                 );

   /// @}
   
   void stdFitsHeader(fits::fitsHeader * head);
   
   virtual void makeMaskCube();
   
   ///De-rotate the PSF subtracted images
   void derotate();

   double t_fake_begin {0};
   double t_fake_end {0};
   
   double t_derotate_begin {0};
   double t_derotate_end {0};
   
};

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation()
{
}

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation( const std::string & dir, 
                                                        const std::string & prefix, 
                                                        const std::string & ext
                                                      ) : HCIobservation<realT>(dir,prefix,ext)
{
}

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation( const std::string & fileListFile) : HCIobservation<realT>(fileListFile)
{
}

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation( const std::string & dir, 
                                                        const std::string & prefix, 
                                                        const std::string & ext,
                                                        const std::string & RDIdir,
                                                        const std::string & RDIprefix,
                                                        const std::string & RDIext
                                                      ) : HCIobservation<realT>(dir,prefix,ext, RDIdir, RDIprefix, RDIext)
{
}

template<typename _realT, class _derotFunctObj>
ADIobservation<_realT, _derotFunctObj>::ADIobservation( const std::string & fileListFile,
                                                        const std::string & RDIfileListFile
                                                      ) : HCIobservation<realT>(fileListFile, RDIfileListFile)
{
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::readFiles()
{      
   this->m_keywords.clear();
   
   if(!m_derotF.isSetup())
   {
      mxError("ADIobservation::readFiles", MXE_PARAMNOTSET, "Derotator is not configured.");
      return -1;
   }

   /*----- Append the ADI keywords to propagate them if needed -----*/

   for(size_t i=0;i<m_derotF.m_keywords.size();++i)
   {
      this->m_keywords.push_back(m_derotF.m_keywords[i]);
   }
         
   if( HCIobservation<realT>::readFiles() < 0) return -1;
   
   return 0;   
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::postReadFiles()
{
   m_derotF.extractKeywords(this->m_heads);
   
   if(m_fakeFileName != ""  && !this->m_skipPreProcess) 
   {
      std::cerr << "Injecting fakes in target images...\n";
      if( injectFake(this->m_tgtIms, this->m_fileList, m_derotF, 1, 1) < 0) return -1;
   }
   
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::postCoadd()
{
   m_derotF.extractKeywords(this->m_heads);
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::readRDIFiles()
{      
   this->m_RDIkeywords.clear();
   
   if(!m_RDIderotF.isSetup())
   {
      mxError("ADIobservation::readRDIFiles", MXE_PARAMNOTSET, "Derotator is not configured.");
      return -1;
   }

   /*----- Append the ADI keywords to propagate them if needed -----*/

   for(size_t i=0;i<m_RDIderotF.m_keywords.size();++i)
   {
      this->m_RDIkeywords.push_back(m_RDIderotF.m_keywords[i]);
   }
         
   if( HCIobservation<realT>::readRDIFiles() < 0) return -1;
   
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::postRDIReadFiles()
{
   m_RDIderotF.extractKeywords(this->m_RDIheads);
   
   /*if(m_fakeFileName != ""  && !this->m_skipPreProcess) 
   {
      if( injectFake(this->m_refIms, this->m_RDIfileList, m_RDIderotF, m_RDIFluxScale, m_RDISepScale) < 0) return -1;
   }*/
   
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::postRDICoadd()
{
   m_RDIderotF.extractKeywords(this->m_RDIheads);
   return 0;
}




template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::readPSFSub( const std::string & dir,
                                                        const std::string & prefix,
                                                        const std::string & ext,
                                                        size_t nReductions 
                                                      )
{
   //Load first file to condigure based on its header.
   std::vector<std::string> flist = ioutils::getFileNames(dir, prefix, "000", ext);
   
   fits::fitsHeader fh;
   eigenImage<realT> im;
   fits::fitsFile<realT> ff;
   
   ff.read(im, fh, flist[0]);
   
   if(!m_derotF.isSetup())
   {
      mxError("ADIobservation::readFiles", MXE_PARAMNOTSET, "Derotator is not configured.");
      return -1;
   }
 
   if(fh.count("POSTMEDS") != 0)
   {
      m_postMedSub = fh["POSTMEDS"].Int();
      std::cerr << "postMedSub: " << m_postMedSub << "\n";
   }
   
   if(fh.count("FAKEFILE") != 0)
   {
      m_fakeFileName = fh["FAKEFILE"].String();
      std::cerr << "fakeFileName: " << m_fakeFileName << "\n";
   }
   
   if(fh.count("FAKESCFL") != 0)
   {
      m_fakeScaleFileName = fh["FAKESCFL"].String();
      std::cerr << "fakeScaleFileName: " << m_fakeScaleFileName << "\n";
   }
   
   if(fh.count("FAKESEP") != 0)
   {
      ioutils::parseStringVector(m_fakeSep, fh["FAKESEP"].String(), ",");
      
      if(m_fakeSep.size() == 0)
      {
         mxError("KLIPReduction", MXE_PARSEERR, "FAKESEP vector did not parse correctly.");
         return -1;
      }
      std::cerr << "fakeSep: " << fh["FAKESEP"].String() << "\n";
   }
   
   if(fh.count("FAKEPA") != 0)
   {
      ioutils::parseStringVector(m_fakePA, fh["FAKEPA"].String(), ",");
      
      if(m_fakePA.size() == 0)
      {
         mxError("KLIPReduction", MXE_PARSEERR, "FAKEPA vector did not parse correctly.");
         return -1;
      }
      std::cerr << "fakePA: " << fh["FAKEPA"].String() << "\n";
   }
   
   if(fh.count("FAKECONT") != 0)
   {
      ioutils::parseStringVector(m_fakeContrast, fh["FAKECONT"].String(), ",");
      
      if(m_fakeContrast.size() == 0)
      {
         mxError("KLIPReduction", MXE_PARSEERR, "FAKECONT vector did not parse correctly.");
         return -1;
      }
      std::cerr << "fakeContrast: " << fh["FAKECONT"].String() << "\n";
   }
   
   
   /*----- Append the ADI keywords to propagate them if needed -----*/

   for(size_t i=0;i<m_derotF.m_keywords.size();++i)
   {
      this->m_keywords.push_back(m_derotF.m_keywords[i]);
   }
   
      
   if( HCIobservation<realT>::readPSFSub(dir,prefix,ext, nReductions) < 0) return -1;
   
   m_derotF.extractKeywords(this->m_heads);
   
   return 0;
}

template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::injectFake( eigenCube<realT> & ims,
                                                        std::vector<std::string> & fileList,
                                                        derotFunctObj & derotF,
                                                        realT RDIFluxScale,
                                                        realT RDISepScale
                                                      )
{
   t_fake_begin = sys::get_curr_time();
   
   //typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imT;
   eigenImageT fakePSF;
   std::vector<std::string> fakeFiles; //used if m_fakeMethod == HCI::list
   
   fits::fitsFile<realT> ff;
   std::ifstream scaleFin; //for reading the scale file.
      

   //Fake Scale -- default to 1, read from file otherwise
   std::vector<realT> fakeScale(ims.planes(), 1.0);
   if(m_fakeScaleFileName != "")
   {      
      std::vector<std::string> sfileNames;
      std::vector<realT> imS;
      
      //Read the scale file and load it into a map
      if( ioutils::readColumns(m_fakeScaleFileName, sfileNames, imS) < 0) return -1;
      
      std::map<std::string, realT> scales;     
      for(size_t i=0;i<sfileNames.size();++i) scales[ioutils::pathFilename(sfileNames[i].c_str())] = imS[i];
      
      for(size_t i=0; i<fileList.size(); ++i)
      {
         if(scales.count(ioutils::pathFilename(fileList[i].c_str())) > 0)
         {
            fakeScale[i] = scales[ioutils::pathFilename(fileList[i].c_str())];
         }
         else
         {
            std::cerr << "File name not found in fakeScaleFile:\n";
            std::cerr << ioutils::pathFilename(fileList[i].c_str()) << "\n";
            exit(-1);
         }
      }
   } //if(fakeScaleFileName != "")
      
   if(m_fakeMethod == HCI::single)
   {
      if( ff.read( fakePSF, m_fakeFileName ) < 0) return -1;
   }

   if(m_fakeMethod == HCI::list)
   {
      if( ioutils::readColumns(m_fakeFileName, fakeFiles) < 0) return -1;
   }
   
   for(int i=0; i<ims.planes(); ++i)
   {
      if(m_fakeMethod == HCI::list)
      {
         ff.read(fakePSF, fakeFiles[i]);
      }

      for(size_t j=0;j<m_fakeSep.size(); ++j)
      {
         if( injectFake(fakePSF, ims, i, derotF.derotAngle(i), m_fakePA[j], m_fakeSep[j], m_fakeContrast[j], fakeScale[j], RDIFluxScale, RDISepScale) < 0) return -1;
      }
   }
   
   
   t_fake_end = sys::get_curr_time();
   
   return 0;
}


template<typename _realT, class _derotFunctObj>
int ADIobservation<_realT, _derotFunctObj>::injectFake( eigenImageT & fakePSF,
                                                        eigenCube<realT> & ims,
                                                        int image_i,
                                                        realT derotAngle,
                                                        realT PA,
                                                        realT sep,
                                                        realT contrast,
                                                        realT scale,
                                                        realT RDIFluxScale,
                                                        realT RDISepScale
                                                      )
{

   //Check for correct sizing
   if( (fakePSF.rows() < ims.rows() && fakePSF.cols() >= ims.cols()) || 
                        (fakePSF.rows() >= ims.rows() && fakePSF.cols() < ims.cols()))
   {
      mxThrowException(err::sizeerr, "ADIobservation::injectFake", "fake PSF has different dimensions and can't be sized properly");
   }
   
   //Check if fake needs to be padded out
   if(fakePSF.rows() < ims.rows() && fakePSF.cols() < ims.cols())
   {
      eigenImageT pfake(ims.rows(), ims.cols());
      padImage(pfake, fakePSF, 0.5*(ims.rows()-fakePSF.rows()), 0);
      fakePSF = pfake;
   }
   
   //Check if fake needs to be cut down
   if(fakePSF.rows() > ims.rows() && fakePSF.cols() > ims.cols())
   {
      eigenImageT cfake(ims.rows(), ims.cols());
      cutPaddedImage(cfake, fakePSF, 0.5*(fakePSF.rows() - ims.rows()));
      fakePSF = cfake;
   }

   /*** Now shift to the separation and PA, scale, apply contrast, and inject ***/
   //allocate shifted fake psf
   eigenImageT shiftFake(fakePSF.rows(), fakePSF.cols());
   
   realT ang, dx, dy;

   ang = math::dtor(-1*PA) + derotAngle;
      
   dx = sep * RDISepScale * sin(ang);
   dy = sep * RDISepScale * cos(ang);
               
   imageShift(shiftFake, fakePSF, dx, dy, cubicConvolTransform<realT>());
   
   ims.image(image_i) = ims.image(image_i) + shiftFake*scale*RDIFluxScale*contrast;

   return 0;
   
}

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::makeMaskCube()
{
   if( this->m_mask.rows() != this->m_Nrows || this->m_mask.cols() != this->m_Ncols)
   {
      std::cerr << "\nMask is not the same size as images.\n\n";
      exit(-1);
   }
   
   this->m_maskCube.resize( this->m_Nrows, this->m_Ncols, this->m_Nims);
   
   #pragma omp parallel
   {
      eigenImageT rm;
      
      #pragma omp for
      for(int i=0; i< this->m_Nims; ++i)
      {
         rotateMask( rm, this->m_mask, m_derotF.derotAngle(i));
         this->m_maskCube.image(i) = rm;
      }
   }
   
   fits::fitsFile<realT> ff; 
   ff.write("maskCube.fits", this->m_maskCube);
}

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::derotate()
{
   t_derotate_begin = sys::get_curr_time();
   
   for(size_t n=0; n<this->m_psfsub.size(); ++n)
   {
      #pragma omp parallel
      {
         eigenImageT rotim;
         realT derot;

         #pragma omp for
         for(int i=0; i<this->m_psfsub[n].planes();++i)
         {
            derot = m_derotF.derotAngle(i);
            if(derot != 0) 
            {
               imageRotate(rotim, this->m_psfsub[n].image(i), derot, cubicConvolTransform<realT>());
               this->m_psfsub[n].image(i) = rotim;
            }
         }
      }
   }
   
   t_derotate_end = sys::get_curr_time();
}


//If fakeFileName == "" or skipPreProcess == true then use the structure of propagated values

template<typename _realT, class _derotFunctObj>
void ADIobservation<_realT, _derotFunctObj>::stdFitsHeader( fits::fitsHeader * head)
{
   if(head == 0) return;
   
   head->append("", fits::fitsCommentType(), "----------------------------------------");
   head->append("", fits::fitsCommentType(), "mx::ADIobservation parameters:");
   head->append("", fits::fitsCommentType(), "----------------------------------------");

   head->append("POSTMEDS", m_postMedSub, "median subtraction after processing");
   
   if(m_fakeFileName != "")
   head->append("FAKEFILE", m_fakeFileName, "name of fake planet PSF file");
   
   if(m_fakeScaleFileName != "")
   head->append("FAKESCFL", m_fakeScaleFileName, "name of fake planet scale file name");

   std::stringstream str;
   
   if(m_fakeSep.size() > 0)
   {
      for(size_t nm=0;nm < m_fakeSep.size()-1; ++nm) str << m_fakeSep[nm] << ",";
      str << m_fakeSep[m_fakeSep.size()-1];      
      head->append<char *>("FAKESEP", (char *)str.str().c_str(), "separation of fake planets");
   }
   
   if(m_fakePA.size() > 0 )
   {
      str.str("");
      for(size_t nm=0;nm < m_fakePA.size()-1; ++nm) str << m_fakePA[nm] << ",";
      str << m_fakePA[m_fakePA.size()-1];      
      head->append<char *>("FAKEPA", (char *)str.str().c_str(), "PA of fake planets");
   }
   
   if( m_fakeContrast.size() > 0)
   {
      str.str("");
      for(size_t nm=0;nm < m_fakeContrast.size()-1; ++nm) str << m_fakeContrast[nm] << ",";
      str << m_fakeContrast[m_fakeContrast.size()-1];      
      head->append<char *>("FAKECONT", (char *)str.str().c_str(), "Contrast of fake planets");
   }
}

template<typename realT> class ADIDerotator;

extern template class ADIobservation<float, ADIDerotator<float>>;
extern template class ADIobservation<double, ADIDerotator<double>>;

} //namespace improc
} //namespace mx

#endif //__ADIobservation_hpp__


