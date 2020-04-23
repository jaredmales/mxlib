/** \file HCIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the basic high contrast imaging data type.
  * \ingroup hc_imaging_files
  * \ingroup image_processing_files
  *
  */

#ifndef __HCIobservation_hpp__
#define __HCIobservation_hpp__

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include "../mxlib.hpp"

#include "../mxException.hpp"

#include "../math/templateBLAS.hpp"
#include "../ioutils/fileUtils.hpp"
#include "../timeUtils.hpp"
#include "../ioutils/readColumns.hpp"

#include "eigenImage.hpp"
#include "eigenCube.hpp"
#include "imageFilters.hpp"
#include "imageMasks.hpp"
#include "imageTransforms.hpp"
#include "fitsFile.hpp"

namespace mx
{

namespace improc
{

/*double t_begin;
double t_end;

double t_load_begin;
double t_load_end;

double t_coadd_begin;
double t_coadd_end;

double t_preproc_begin;
double t_preproc_end;

double t_combo_begin;
double t_combo_end;*/


///Namespace for high contrast imaging enums.
/** \ingroup hc_imaging_enums
  */
namespace HCI
{
   ///Possible combination methods
   /** \ingroup hc_imaging_enums
     */
   enum combineMethods{ noCombine, ///< Do not combine the images.
                        medianCombine, ///< Combine with the median.
                        meanCombine, ///< Combine with the mean.
                        sigmaMeanCombine, ///< Combine with the sigma clipped mean.
                        debug};
}

/// The basic high contrast imaging data type
/** This class manages file reading, resizing, co-adding, pre-processing (masking and filtering),
  * and final image combination.
  *
  * \tparam _realT is the floating point type in which to do all arithmetic.
  *
  * \ingroup hc_imaging
  */
template<typename _realT>
struct HCIobservation
{

   ///The arithmetic type used for calculations.  Does not have to match the type in images on disk.
   typedef _realT realT;

   ///The Eigen image array type basted on realT
   typedef Array<realT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;

   /** \name File Reading
     * Options to control which files are read, how they are read, what meta data is extracted
     * from FITS headers, sizing and masking, etc.
     * @{
     */

   ///The list of files to read in.
   /** This can be set on construction or by calling \ref loadFileList
     */
   std::vector<std::string> fileList;

   ///Specify how many files from fileList to delete from the front of the list
   int deleteFront {0};

   ///Specify how many files from fileList to delete from the back of the list
   int deleteBack {0};


   ///File containing 2 space-delimited columns of fileVame qualityValue pairs.
   /** If this is not empty and \ref qualityThreshold is > 0, then only images where
     * qualityValue >= qualityThreshold.
     *
     * The only restriction on qualityThreshold is that it is > 0.  It is intendend to be
     * something like Strehl ratio.
     */
   std::string qualityFile;

   ///Threshold to apply to qualityValues read from \ref qualityFile.
   /** If <= 0, then thresholding is not performed.
    */
   realT qualityThreshold {0};

   ///Just prints the names and qualities of the files which pass threshold, and stop.
   bool thresholdOnly {false};

   ///Name of the keyword to use for the image date.
   /** Specifies the keyword corresponding to the date.  This is
     * the "DATE" keyword for file write time, and usually "DATE-OBS" for actual observation time.
     *
     * Default is "DATE-OBS".
     *
     * If empty "", then image date is not read.
     */
   std::string MJDKeyword {"DATE-OBS"};

   ///Whether or not the date is in ISO 8601 format
   bool MJDisISO8601 {true};

   ///If the date is not ISO 8601, this specifies the conversion to Julian Days (i.e. seconds to days)
   realT MJDUnits {1.0};

   ///Vector of FITS header keywords to read from the files in fileList.
   std::vector<std::string> keywords;

   ///Set the image size.  Images are cut down to this size after reading.
   /** Set to \<= 0 to use images uncut.
     *
     * Image sizes are not increased if this is larger than their size on disk.
     */
   int imSize {0};


   ///@}

   /** \name Coadding
     * These parameters control whether and how the images are coadded after being read.  Coadding can
     * be done up to a given number of images, and/or a given elapsed time.
     *
     * Averages the values of given Keywords as well.
     * @{
     */

   ///Determine how to coadd the raw images.
   /** Possibilities are
     * - HCI::noCombine -- [default] do not combine.  This turns off coadding.
     * - HCI::medianCombine --  coadded image is the median
     * - HCI::meanCombine -- coadded image is the simple mean
     *
     * No other types of combination are currently supported for coadding.
     */
   int coaddCombineMethod {HCI::noCombine};

   ///Maximum number of images to coadd
   int coaddMaxImno {0};

   ///Maximum elapsed time over which to coadd the images.
   int coaddMaxTime {0};

   ///The values of these keywords will be averaged and replaced.
   std::vector<std::string> coaddKeywords;

   ///@} -- coadding


   /** \name Masking
     * A 1/0 mask can be supplied, which is used in pre-processing and in image combination.
     * @{
     */

   ///Specify a mask file to apply
   /**No mask is applied if this is empty.
     */
   std::string maskFile;

   ///@}



public:

   /** \name Pre-Processing
     * These options control the pre-processing masking and filtering.
     * They are performed in the following order:
     * -# mask applied (enabled by preProcess_mask
     * -# radial profile subtraction (enabled by preProcess_subradprof)
     * -# mask applied (enabled by preProcess_mask
     * -# symmetric Gaussian unsharp mask (preProcess_gaussUSM_fwhm)
     * -# mask applied (enabled by preProcess_mask
     * -# azimuthal unsharp mask (preProcess_azUSM_azW, and preProcess_azUSM_radW)
     * -# mask applied (enabled by preProcess_mask)
     * @{
     */

   bool skipPreProcess {false}; ///<Don't do any of the pre-processing steps (including coadding).

   bool preProcess_beforeCoadd {false}; ///<controls whether pre-processing takes place before or after coadding

   bool preProcess_mask {true}; ///<If true, the mask is applied during each pre-processing step.

   bool preProcess_subradprof {false}; ///<If true, a radial profile is subtracted from each image.

   ///Azimuthal boxcar width for azimuthal unsharp mask
   /** If this is 0 then azimuthal-USM is not performed.
     */
   realT preProcess_azUSM_azW {0};

   ///Radial boxcar width for azimuthal unsharp mask
   /** If this is 0 then azimuthal-USM is not performed.
     */
   realT preProcess_azUSM_radW {0};


   ///Kernel FWHM for symmetric unsharp mask (USM)
   /** USM is not performed if this is 0.
    */
   realT preProcess_gaussUSM_fwhm {0};

   ///Set path and file prefix to output the pre-processed images.
   /** If empty, then pre-processed images are not output.
     */
   std::string preProcess_outputPrefix;

   /// If true, then we stop after pre-processing.
   bool preProcess_only {false};


   ///@}

   /** \name Image Combination
     * These options control how the final image combination is performed.
     * @{
     */


   ///Determine how to combine the PSF subtracted images
   /** Possibilities are
     * - HCI::noCombine -- do not combine
     * - HCI::medianCombine -- [default] final image is the median
     * - HCI::meanCombine -- final image is the simple mean
     * - HCI::weightedMeanCombine -- final image is the weighted mean.  weightFile must be provided.
     * - HCI::sigmaMeanCombine -- final image is sigma clipped mean.  If sigmaThreshold \<= 0, then it reverts to meanCombine.
     */
   int combineMethod {HCI::meanCombine};

   ///Specifies a file containing the image weights, for combining with weighted mean.
   /** This 2-column space-delimited ASCII file containing  filenames and weights. It must be specified before readFiles()
     * is executed.  In readFiles this is loaded after the fileList is cutdown and matched to the remaining files.
     */
   std::string weightFile;



   /// The standard deviation threshold used if combineMethod == HCI::sigmaMeanCombine.
   realT sigmaThreshold {0};

   /// The minimum fraction of good (un-masked) pixels to include in the final combination (0.0 to 1.0). If not met, then the pixel will be NaN-ed.
   realT m_minGoodFract {0.0};

   ///@}

   /** \name Output
     * These options control the ouput fo the final combined images and the individual PSF subtracted images.
     * @{
     */

   ///Set whether the final combined image is written to disk
   int doWriteFinim {1};

   ///The directory where to write output files.
   std::string outputDir;

   ///The base file name of the output final image
   /** The complete name is formed by combining with a sequential number and the ".fits" extension.
     * that is: finimName0000.fits.  This behavior can be modified with exactFinimName.
     */
   std::string finimName {"finim_"};

   ///Use finimName exactly as specified, without appending a number or an extension.
   /** Output is still FITS format, regardless of extension.  This will overwrite
     * an existing file without asking.
     */
   bool exactFinimName {false};

   ///Controls whether or not the individual PSF subtracted images are written to disk.
   /** - true -- write to disk
     * - false -- [default] don't write to disk
     */
   bool doOutputPSFSub {false};

   ///Prefix of the FITS file names used to write individual PSF subtracted images to disk if doOutputPSFSub is true.
   std::string PSFSubPrefix;

   ///@}

   ///\name The Raw Data
   /** @{
    */

   ///Vector image times, in MJD.
   std::vector<double> imageMJD;

   ///Vector of FITS headers,one per file, populated with the values for the keywords.
   std::vector<fitsHeader> heads;

   ///Whether or not the fileList has been read.
   bool filesRead {false};

   ///Whether or not the specified files have been deleted from fileList
   bool filesDeleted {false};


   eigenImageT mask; ///< The mask
   eigenCube<realT> maskCube; ///< A cube of masks, one for each input image, which may be modified versions (e.g. rotated) of mask.

   ///Vector to hold the imaged weights read from the weightFile.
   /** After readWeights is executed by readFiles, this will contain the normalized weights.
     */
   std::vector<realT> comboWeights;


   ///The image cube
   eigenCube<realT> imc;

   int Nims {0};  ///<Number of images in imc
   int Nrows {0}; ///<Number of rows of the images in imc
   int Ncols {0}; ///<Number of columns of the images in imc
   int Npix {0}; ///<Pixels per image, that is Nrows*Ncols

   ///@}

   ///\name The Reduced Data
   /** @{
     */
   ///The PSF subtracted images
   /** This is a vector of cubes so that it can contain results from different reductions,
     * e.g. different modes when using KLIP.
     */
   std::vector<eigenCube<realT> > psfsub;

   ///The final combined images, one for each cube in psfsub.
   eigenCube<realT> finim;

   ///@}

   ///\name Construction and Initialization
   /** @{
     */
   ///Sets defaults for all parameters.  Called by each constructor.
   void initialize();

   ///Default c'tor
   HCIobservation();

   ///Construct and load the file list.
   /** Populates the \ref fileList vector by searching on disk for files which match
     * "dir/prefix*.ext".  See \ref loadFileList
     *
     * \param [in] dir is the directory to search.
     * \param [in] prefix is the initial part of the file name.  Can be empty "".
     * \param [in] ext is the extension to append to the file name. If ommitted defaults to ".fits"
     */
   HCIobservation(const std::string &dir, const std::string &prefix, const std::string & ext=".fits");


   ///Construct using a file containing the file lise
   /** Populates the \ref fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     *
     * \param [in] fileListFile is a file name to read.
     */
   explicit HCIobservation(const std::string &fileListFile);

   ///Load the file list from a file
   /** Populates the \ref fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     *
     * \param [in] fileListFile is a file name to read.
     */
   void loadFileList(const std::string &fileListFile);

   ///Load the file list
   /** Populates the \ref fileList vector by searching on disk for files which match the given parameters.
     * Uses \ref mx::getFileNames to search for all files which match "dir/prefix*.ext".
     *
     * \param [in] dir is the directory to search.
     * \param [in] prefix is the initial part of the file name.  Can be empty "".
     * \param [in] ext is the extension to append to the file name. Can be empty "".
     */
   void loadFileList(const std::string &dir, const std::string &prefix, const std::string &ext);

   ///@}

   ///Read the list of files, cut to size, and apply the mask.
   /**
     * \returns 0 on success, -1 on  error.
     */
   int readFiles();

   ///Read the image weights from weightFile
   /** This is called by readFiles().
     *
     * \returns 0 on success, -1 on  error.
     */
   int readWeights();

   ///Read the image qualities from qualityFile and apply the threshold
   /** This is called by readFiles().
     *
     * \returns 0 on success, -1 on  error.
     */
   int threshold();

   ///Perform post-read actions, for use by derived classes
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
   
   ///Populate the mask cube.  Derived classes can do this as appropriate, e.g. by rotating the mask.
   virtual void makeMaskCube();

   ///Coadd the images
   void coaddImages();

   ///Perform post-coadd actions, for use by derived classes.
   /** A key example is to update keywords after the averaging occurs in coaddImages().
     *
     * \returns 0 on success
     * \returns \<0 on error.
     */
   virtual int postCoadd();
   
   ///Do the pre-processing
   void preProcess();

   ///Output the pre-processed images
   void outputPreProcessed();


   ///Combine the images into a single final image.
   /** Images are combined by the method specified in \ref combineMethod
     */
   void combineFinim();

   ///Write the final combined image to disk
   /**
     */
   void writeFinim(fitsHeader * addHead = 0);


   ///Write the PSF subtracted images to disk
   /**
    */
   void outputPSFSub(fitsHeader * addHead = 0);


   double t_begin {0};
   double t_end  {0};

   double t_load_begin  {0};
   double t_load_end  {0};

   double t_coadd_begin  {0};
   double t_coadd_end  {0};

   double t_preproc_begin  {0};
   double t_preproc_end  {0};

   double t_azusm_begin  {0};
   double t_azusm_end  {0};

   double t_gaussusm_begin  {0};
   double t_gaussusm_end  {0};


   double t_combo_begin  {0};
   double t_combo_end  {0};

};

template<typename _realT>
void HCIobservation<_realT>::initialize()
{
   deleteFront = 0;
   deleteBack = 0;
   filesDeleted = false;

   thresholdOnly = false;

   MJDKeyword = "DATE-OBS";
   MJDisISO8601 = true;
   MJDUnits = 1.0;


   coaddCombineMethod = HCI::noCombine;
   coaddMaxImno = 0;
   coaddMaxTime = 0;

   skipPreProcess = false;

   preProcess_mask = false;

   preProcess_beforeCoadd = false;
   preProcess_subradprof = false;
   preProcess_azUSM_azW = 0;
   preProcess_azUSM_radW = 0;
   preProcess_gaussUSM_fwhm = 0;

   preProcess_only = false;



   filesRead = false;

   combineMethod = HCI::meanCombine;
   sigmaThreshold = 0;

   doWriteFinim = 1;
   finimName = "finim_";
   exactFinimName = false;

   doOutputPSFSub = false;


   t_begin = 0;
   t_end = 0;

   t_load_begin = 0;
   t_load_end = 0;

   t_coadd_begin = 0;
   t_coadd_end = 0;

   t_preproc_begin = 0;
   t_preproc_end = 0;

   t_azusm_begin =0;
   t_azusm_end = 0;

   t_gaussusm_begin = 0;
   t_gaussusm_end = 0;

   t_combo_begin = 0;
   t_combo_end = 0;
}

template<typename _realT>
HCIobservation<_realT>::HCIobservation()
{
   //initialize();
}

template<typename _realT>
HCIobservation<_realT>::HCIobservation(const std::string & dir, const std::string & prefix, const std::string & ext)
{
   //initialize();
   loadFileList(dir, prefix, ext);
}

template<typename _realT>
HCIobservation<_realT>::HCIobservation(const std::string & fileListFile)
{
   //initialize();
   loadFileList(fileListFile);
}

template<typename _realT>
inline void HCIobservation<_realT>::loadFileList(const std::string & fileListFile)
{
   ioutils::readColumns(fileListFile, fileList);
   filesDeleted = false;
}

template<typename _realT>
inline void HCIobservation<_realT>::loadFileList(const std::string & dir, const std::string & prefix, const std::string & ext)
{
   fileList = ioutils::getFileNames(dir, prefix, "", ext);
   filesDeleted = false;
}

template<typename _realT>
int HCIobservation<_realT>::readFiles()
{
   if(fileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The fileList has 0 length, there are no files to be read.");
      return -1;
   }

   //First make the list deletions
   if(!filesDeleted)
   {
      if(deleteFront > 0)
      {
         fileList.erase(fileList.begin(), fileList.begin()+deleteFront);
      }

      if(deleteBack > 0)
      {
         fileList.erase(fileList.end()-deleteBack, fileList.end());
      }
      filesDeleted = true;
   }

   if( fileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The fileList has 0 length, there are no files to be read after deletions.");
      return -1;
   }

   if(qualityFile != "")
   {
      if( threshold() < 0) return -1;

      if(thresholdOnly)
      {
         std::cout << "#Files which passed thresholding:\n";
         for(size_t i=0; i<fileList.size(); ++i)
         {
            std::cout << fileList[i] << "\n";
         }

         exit(0);
      }
   }

   if( fileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The fileList has 0 length, there are no files to be read after thresholding.");
      return -1;

   }

   if(weightFile != "")
   {
      if(readWeights() < 0) return -1;
   }

   Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

   fitsFile<realT> f(fileList[0]);

   f.read(im);

   fitsHeader head;

   if(MJDKeyword != "") head.append(MJDKeyword);

   for(size_t i=0;i<keywords.size();++i)
   {
      head.append(keywords[i]);
   }

   /*----- Append the HCI keywords to propagate them if needed -----*/

   heads.clear(); //This is necessary to make sure heads.resize() copies head on a 2nd call
   heads.resize(fileList.size(), head);

   if(imSize > 0)
   {
      //Now make sure we don't read too much.
      if(imSize > im.rows()) imSize = im.rows();
      if(imSize > im.cols()) imSize = im.cols();

      //the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
      f.setReadSize( floor(0.5*(im.rows()-1) - 0.5*(imSize-1) +0.1), floor(0.5*(im.cols()-1.0) - 0.5*(imSize-1.0)+0.1), imSize, imSize);
      im.resize(imSize, imSize);
   }

   imc.resize(im.rows(), im.cols(), fileList.size());

   t_load_begin = get_curr_time();

  /** \bug It is this step that is really slow during file reads *sometimes*
   */
   f.read(imc.data(), heads, fileList);

   f.setReadSize();

   if(MJDKeyword != "")
   {
      imageMJD.resize(heads.size());

      if(MJDisISO8601)
      {
         for(size_t i=0;i<imageMJD.size();++i)
         {
            imageMJD[i] =  mx::ISO8601date2mjd(heads[i][MJDKeyword].String());
         }
      }
      else
      {
         for(size_t i=0;i<imageMJD.size();++i)
         {
            imageMJD[i] =  heads[i][MJDKeyword].Value<realT>()*MJDUnits;
         }
      }
   }

   t_load_end = get_curr_time();

   Nims =  imc.planes();
   Nrows = imc.rows();
   Ncols = imc.cols();
   Npix =  imc.rows()*imc.cols();

   std::cerr << "loading complete\n";
   
   ///\todo zeroNaNs should be a member and be configurable.  Probably should be true by default.
   bool zeroNaNs = true;

   if( zeroNaNs )
   {
      std::cerr << "zero-ing NaNs\n";
      for(int k=0; k<Nims; ++k)
      {
         for(int i=0; i< Nrows; ++i)
         {
            for(int j=0; j<Ncols; ++j)
            {
               if( !std::isnormal( imc.image(k)(i,j)) )
               {
                  imc.image(k)(i,j) = 0;
               }
            }
         }
      }
   }

   /*** Now do the post-read actions ***/
   if( postReadFiles() < 0) return -1;

   /*** Load the mask ***/
   if( maskFile != "")
   {
      std::cerr << "creating mask cube\n";
      fitsFile<realT> ff;
      ff.read(mask, maskFile);
      makeMaskCube();
   }

   if(!skipPreProcess)
   {
      /*** Now do any pre-processing ***/
      if(preProcess_beforeCoadd) preProcess();

      if(coaddCombineMethod != HCI::noCombine)
      {
         coaddImages();
         if( postCoadd() < 0)
         {
            std::cerr << "Post coadd error " << __FILE__ << " " << __LINE__ << "\n";
            return -1;
         }
         
         //Re-make the mask cube if we coadded...
         if( maskFile != "")
         {
            fitsFile<realT> ff;
            ff.read(mask, maskFile);

            makeMaskCube();
            
         }
      }

      /*** Now do any pre-processing if not done already***/
      if(!preProcess_beforeCoadd) preProcess();

      outputPreProcessed();
   }
   
   
   filesRead = true;
   
   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::readPSFSub( const std::string & dir,
                                        const std::string & prefix,
                                        const std::string & ext,
                                        size_t nReductions 
                                      )
{

   psfsub.resize(nReductions);
   
   for(size_t n =0; n<nReductions; ++n)
   {
      char nstr[5];
      snprintf(nstr, sizeof(nstr), "%03zu", n);
      
      std::string nprefix = prefix + "_" + nstr + "_";
      loadFileList(dir, nprefix, ext);
      
      if(fileList.size() == 0)
      {
         mxError("HCIobservation", MXE_FILENOTFOUND, "The fileList has 0 length, there are no files to be read.");
         return -1;
      }

      Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

      fitsFile<realT> f(fileList[0]);

      f.read(im);

      fitsHeader head;

      if(MJDKeyword != "") head.append(MJDKeyword);

      for(size_t i=0;i<keywords.size();++i)
      {
         head.append(keywords[i]);
      }

      /*----- Append the HCI keywords to propagate them if needed -----*/

      heads.clear(); //This is necessary to make sure heads.resize() copies head on a 2nd call
      heads.resize(fileList.size(), head);

      if(imSize > 0)
      {
         //Now make sure we don't read too much.
         if(imSize > im.rows()) imSize = im.rows();
         if(imSize > im.cols()) imSize = im.cols();

         //the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
         f.setReadSize( floor(0.5*(im.rows()-1) - 0.5*(imSize-1) +0.1), floor(0.5*(im.cols()-1.0) - 0.5*(imSize-1.0)+0.1), imSize, imSize);
         im.resize(imSize, imSize);
      }

      
      if(n > 0)
      {
         if(fileList.size() != (size_t) Nims)
         {
            mxError("HCIobservation", MXE_INVALIDARG, "Different number of images in reductions.");
            return -1;
         }
         if(Nrows != im.rows())
         {
            mxError("HCIobservation", MXE_INVALIDARG, "Different number of rows in reductions.");
            return -1;
         }
         if(Ncols != im.cols())
         {
            mxError("HCIobservation", MXE_INVALIDARG, "Different number of cols in reductions.");
            return -1;
         }
      }
      else
      {
         std::cerr << "found " << nReductions << " sets of " << fileList.size() << " " << im.rows() << " x " << im.cols() << " files\n";
      }    
      Nims =  fileList.size();
      Nrows = im.rows();
      Ncols = im.cols();
      Npix =  im.rows()*im.cols();
   
      psfsub[n].resize(Nrows, Ncols, Nims);
      
      t_load_begin = get_curr_time();

      f.read(psfsub[n].data(), heads, fileList);

      f.setReadSize(); 

      if(MJDKeyword != "")
      {
         imageMJD.resize(heads.size());

         if(MJDisISO8601)
         {
            for(size_t i=0;i<imageMJD.size();++i)
            {
               imageMJD[i] =  mx::ISO8601date2mjd(heads[i][MJDKeyword].String());
            }
         }
         else
         {
            for(size_t i=0;i<imageMJD.size();++i)
            {
               imageMJD[i] =  heads[i][MJDKeyword].Value<realT>()*MJDUnits;
            }
         }
      }

      t_load_end = get_curr_time();

      ///\todo zeroNaNs should be a member and be configurable.  Probably should be true by default.
      bool zeroNaNs = true;

      if( zeroNaNs )
      {
         for(int k=0; k<Nims; ++k)
         {
            for(int i=0; i< Nrows; ++i)
            {
               for(int j=0; j<Ncols; ++j)
               {
                  if( !std::isnormal( psfsub[n].image(k)(i,j)) )
                  {
                     psfsub[n].image(k)(i,j) = 0;
                  }
               }
            }
         }
      }

   }
   if(weightFile != "")
   {
      std::vector<std::string> fn;
      ioutils::readColumns(weightFile, fn, comboWeights);

      std::cerr << "read: " << weightFile << " (" << comboWeights.size() << ")\n";
   }
         
   /*** Now do the post-read actions ***/
   if( postReadFiles() < 0) return -1;
  
   /*** Load the mask ***/
   if( maskFile != "")
   {
      fitsFile<realT> ff;
      ff.read(mask, maskFile);

      makeMaskCube();
   }

   filesRead = true;
   
   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::threshold()
{
   std::cerr  << "Thresholding . . .\n";
   if(qualityFile == "")
   {
      mxError("HCIobservation::threshold", MXE_PARAMNOTSET, "qualityFile not set");
      return -1;
   }

   int origsize = fileList.size();

   std::vector<std::string> qfileNames;
   std::vector<realT> imQ;

   //Read the quality file and load it into a map
   ioutils::readColumns(qualityFile, qfileNames, imQ);

   std::map<std::string, realT> quality;
   for(size_t i=0;i<qfileNames.size();++i) quality[basename(qfileNames[i].c_str())] = imQ[i];

   realT q;

   for(size_t i=0; i<fileList.size(); ++i)
   {
      try
      {
         q = quality.at(basename(fileList[i].c_str()));
      }
      catch(...)
      {
         q = qualityThreshold - 1; //Cause it to be erased
      }

      if (q < qualityThreshold)
      {
         fileList.erase(fileList.begin()+i);
         --i;
      }
   }
   std::cerr << "Done.  Selected " << fileList.size() << " out of " << origsize << "\n";

   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::readWeights()
{
   std::ifstream fin;
   std::string str;

   if(weightFile == "")
   {
      mxError("HCIobservation::readWeights", MXE_PARAMNOTSET, "weightFile not set");
      return -1;
   }

   //Read the weight file and load it into a map
   std::vector<std::string> wfileNames;
   std::vector<realT> imW;

   if( ioutils::readColumns(weightFile, wfileNames, imW) < 0) return -1;

   if(imW.size() < fileList.size())
   {
      mxError("HCIobservation::readWeights", MXE_SIZEERR, "not enough weights specified");
      return -1;
   }

   std::map<std::string, realT> weights;
   for(size_t i=0;i<wfileNames.size();++i) weights[basename(wfileNames[i])] = imW[i];

   comboWeights.resize(fileList.size());

   realT wi;
   realT weightSum = 0;
   for(size_t i=0; i<fileList.size(); ++i)
   {
      try
      {
         wi = weights.at(basename(fileList[i]));
      }
      catch(...)
      {
         mxError("HCIobservation::readWeights", MXE_NOTFOUND, "Weight for a file in fileList not found.");
         return -1;
      }
      comboWeights[i] = wi;
      weightSum += wi;
   }

   //Finally normalize the weights
   for(size_t i=0; i< comboWeights.size(); ++i)
   {
      comboWeights[i] /= weightSum;
   }

   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::postReadFiles()
{
   return 0;
}

template<typename _realT>
void HCIobservation<_realT>::makeMaskCube()
{
   if( this->mask.rows() != this->Nrows || this->mask.cols() != this->Ncols)
   {
      std::cerr << "\nMask is not the same size as images.\n\n";
      exit(-1);
   }
   maskCube.resize( Nrows, Ncols, Nims);
   
   for(int i=0; i< Nims; ++i)
   {
      maskCube.image(i) = mask;
   }

}

template<typename _realT>
void HCIobservation<_realT>::coaddImages()
{
   //Validate setup
   if(coaddMaxImno <=0 && coaddMaxTime <= 0) return;
 
   //Validate combine method
   if(coaddCombineMethod == HCI::noCombine) return;

   t_coadd_begin = get_curr_time();

   std::cerr << "coadding raw images\n";

   std::vector<eigenImageT> coadds;

   //We do all math here in double, to avoid losing precision
   std::vector<double> avgMJD;
   std::vector<std::vector<double> > avgVals;

   int combineMethod =  HCI::medianCombine;
   if(coaddCombineMethod == HCI::meanCombine) combineMethod = HCI::meanCombine;

   //Index range of images for next coadd
   int im0, imF;
   im0 = 0;
   imF = 1;

   //Accumulate images to coadd into a cube
   eigenCube<realT> imsToCoadd;

   //Temporary for combination.
   eigenImageT coadd;

   //Accumulate values
   double initMJD;
   std::vector<double> initVals;
   initVals.resize(coaddKeywords.size());

   //For averaging diagnostics
/*   std::ofstream fout;
   fout.open("coadd_diag.txt");
   fout.precision(20);*/
   
   while(im0 < Nims)
   {
      //fout << "*************************************************\nStarting image accumulation:\n*************************************************\n";
      
      //Initialize the accumulators
      initMJD = imageMJD[im0];

      //fout << im0 << " " << imageMJD[im0] << " ";
      for(size_t i=0;i<coaddKeywords.size(); ++i)
      {
         initVals[i] = heads[im0][coaddKeywords[i]].Value<double>();
         //fout << initVals[i] << " ";
      }

      //fout << "\n";
      
      //Now increment imF, then test whether each variable is now outside the range
      bool increment = true;
      while(increment == true)
      {
         ++imF;

         if(imF >= Nims)
         {
            imF = Nims;
            increment = false;
            break;
         }

         if(imF-im0 > coaddMaxImno && coaddMaxImno > 0)
         {
            --imF;
            increment = false;
            break;
         }

         ///\todo this should really include end of exposure too.
         if( (imageMJD[imF] - imageMJD[im0])*86400.0 > coaddMaxTime && coaddMaxTime > 0)
         {
            --imF;
            increment = false;
            break;
         }

      }//while(increment == true)
      //At this point, imF is the first image NOT to include in the next coadd.

      //Now actually do the accumulation
      ///\todo this needs to handle averaging of angles
      for(int imno = im0+1; imno < imF; ++imno)
      {
         initMJD += imageMJD[imno];
         
         //fout << imno << " " << imageMJD[imno] << " ";
         
         for(size_t i=0;i<coaddKeywords.size(); ++i)
         {
            initVals[i] += heads[imno][coaddKeywords[i]].Value<double>();
            //fout << heads[imno][coaddKeywords[i]].Value<double>() << " ";
         }
         //fout << "\n";
      }

      
      
      //And then turn them into an average
      initMJD /= (imF - im0);
      //fout << "--------------------------------------------------------------\n";
      //fout << "avg: " << initMJD;
      for(size_t i=0;i<coaddKeywords.size(); ++i)
      {
         initVals[i] /= (imF-im0);
         //fout << " " << initVals[i];
      }
      //fout << "\n--------------------------------------------------------------\n";

      //Extract the images into the temporary
      imsToCoadd.resize(Nrows, Ncols, imF-im0);
      for(int i =0; i < (imF-im0); ++i)
      {
         imsToCoadd.image(i) = imc.image(im0 + i);
         //fout << "adding image " << im0 + i << "\n";
      }

      
      
      //Here do the combine and insert into the vector
      if(combineMethod == HCI::medianCombine)
      {
         imsToCoadd.median(coadd);
      }
      if(combineMethod == HCI::meanCombine)
      {
         imsToCoadd.mean(coadd);
      }
      coadds.push_back(coadd);

      //Insert the new averages
      avgMJD.push_back(initMJD);
      avgVals.push_back(initVals);

      im0 = imF;
      imF = im0 + 1;
   }//while(im0 < Nims)

   //Now resize imc and copy the coadds to the new cube
   imc.resize(Nrows, Ncols, coadds.size());
   Nims = coadds.size();

   for(int i=0;i<Nims;++i) imc.image(i) = coadds[i];

   Nims =  imc.planes();
   Nrows = imc.rows();
   Ncols = imc.cols();
   Npix =  imc.rows()*imc.cols();

   //Now deal with imageMJD and headers
   imageMJD.erase(imageMJD.begin()+Nims, imageMJD.end());
   heads.erase(heads.begin()+Nims, heads.end());
   
   //fout << "\nWriting averages to header: \n";
   for(int i=0;i<Nims;++i)
   {
      imageMJD[i] = avgMJD[i];
      //fout << i << " " << imageMJD[i];
      for(size_t j=0;j<coaddKeywords.size(); ++j)
      {
         heads[i][coaddKeywords[j]].setValue(avgVals[i][j]);
         //fout << " " << coaddKeywords[j] << ":" << heads[i][coaddKeywords[j]].Double() << "\n";
      }
   }
   //fout << "\n";
   
   t_coadd_end = get_curr_time();

}//void HCIobservation<_realT>::coaddImages()

template<typename _realT>
int HCIobservation<_realT>::postCoadd()
{
   return 0;
}

template<typename _realT>
void HCIobservation<_realT>::preProcess()
{
   t_preproc_begin = get_curr_time();

   //The mask is applied first, and then after each subsequent P.P. step.
   if( maskFile != "" && preProcess_mask)
   {
      std::cerr << "Masking . . .\n";
      #pragma omp parallel for
      for(int i=0;i<imc.planes(); ++i)
      {
         imc.image(i) *= mask;
      }
      std::cerr << "Done\n";
   }

   if( preProcess_subradprof )
   {
      std::cerr << "subtracting radial profile . . .\n";
      eigenImageT rp;

      for(int i=0;i<imc.planes(); ++i)
      {
         Eigen::Map<Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> > imRef(imc.image(i));
         radprofim(rp, imRef, true);
      }

      std::cerr << "done\n";

      if( maskFile != "" && preProcess_mask)
      {
         std::cerr << "Masking . . .\n";
         #pragma omp parallel for
         for(int i=0;i<imc.planes(); ++i)
         {
            imc.image(i) *= mask;
         }
         std::cerr << "Done\n";
      }
   }

   if( preProcess_gaussUSM_fwhm > 0 && preProcess_mask)
   {
      std::cerr << "Applying Gauss USM . . .\n";
      t_gaussusm_begin = get_curr_time();

      #pragma omp parallel for
      for(int i=0;i<imc.planes(); ++i)
      {
         eigenImageT fim, im;
         im = imc.image(i);
         filterImage(fim, im, gaussKernel<eigenImage<float>,2>(preProcess_gaussUSM_fwhm), 0.5*(imc.cols()-1) - preProcess_gaussUSM_fwhm*4);
         im = (im-fim);
         imc.image(i) = im;
      }

      if( maskFile != "")
      {
         std::cerr << "Masking . . .\n";
         #pragma omp parallel for
         for(int i=0;i<imc.planes(); ++i)
         {
            imc.image(i) *= mask;
         }

      }
      t_gaussusm_end = get_curr_time();
      std::cerr << "Done\n";
   }

   if( preProcess_azUSM_azW && preProcess_azUSM_radW )
   {
      std::cerr << "Applying azimuthal USM . . .\n";
      t_azusm_begin = get_curr_time();
      #pragma omp parallel for
      for(int i=0;i<imc.planes(); ++i)
      {
         eigenImageT fim, im;
         im = imc.image(i);
         filterImage(fim, im, azBoxKernel<eigenImage<float>>(preProcess_azUSM_radW, preProcess_azUSM_azW), 0.5*(imc.cols()-1) - preProcess_azUSM_radW);
         im = (im-fim);
         imc.image(i) = im;
      }

      if( maskFile != "" && preProcess_mask)
      {
         std::cerr << "Masking . . .\n";
         #pragma omp parallel for
         for(int i=0;i<imc.planes(); ++i)
         {
            imc.image(i) *= mask;
         }

      }

      t_azusm_end = get_curr_time();
      std::cerr  << "Done\n";
   }





   t_preproc_end = get_curr_time();
}

template<typename _realT>
void HCIobservation<_realT>::outputPreProcessed()
{
   if(preProcess_outputPrefix == "") return;

   std::string bname, fname;

   fitsFile<_realT> ff;

   /** \todo Should add a HISTORY card here */
   for(int i=0; i< Nims; ++i)
   {
      bname = fileList[i];
      fname = preProcess_outputPrefix + basename(bname.c_str()) + ".fits";
      ff.write(fname, imc.image(i).data(), Ncols, Nrows, 1, heads[i]);
   }


}

template<typename _realT>
void HCIobservation<_realT>::combineFinim()
{
   if(combineMethod == HCI::noCombine) return;

   t_combo_begin = get_curr_time();

   //Validate the combineMethod setting
   int method = HCI::medianCombine;

   if(combineMethod == HCI::medianCombine)
   {
      method = HCI::medianCombine;
   }
   else if(combineMethod == HCI::meanCombine)
   {
      method = HCI::meanCombine;
   }
   else if(combineMethod == HCI::sigmaMeanCombine)
   {
      if(sigmaThreshold > 0)
      {
         method = HCI::sigmaMeanCombine;
      }
      else
      {
         method = HCI::meanCombine;
      }
   }
   else if(combineMethod == HCI::debug)
   {
      method = HCI::debug;
   }


   //Create and size temporary image for averaging
   eigenImageT tfinim;

   finim.resize(psfsub[0].rows(), psfsub[0].cols(), psfsub.size());

   //Now cycle through each set of psf subtractions
   for(size_t n= 0; n < psfsub.size(); ++n)
   {
      if(method == HCI::medianCombine)
      {
         psfsub[n].median(tfinim);
         finim.image(n) = tfinim;
      }
      else if(method == HCI::meanCombine)
      {
         if(comboWeights.size() == (size_t) Nims)
         {
            if( maskFile != "" )
            {
               psfsub[n].mean(tfinim, comboWeights, maskCube, m_minGoodFract);
            }
            else
            {
               psfsub[n].mean(tfinim, comboWeights);
            }
         }
         else
         {
            if( maskFile != "" )
            {
               psfsub[n].mean(tfinim, maskCube, m_minGoodFract);
            }
            else
            {
               psfsub[n].mean(tfinim);
            }
         }
         finim.image(n) = tfinim;
      }
      else if(method == HCI::sigmaMeanCombine)
      {
         if(comboWeights.size() == (size_t) Nims)
         {
            if( maskFile != "" )
            {
               psfsub[n].sigmaMean(tfinim, comboWeights, maskCube, sigmaThreshold, m_minGoodFract);
            }
            else
            {
               psfsub[n].sigmaMean(tfinim, comboWeights, sigmaThreshold);
            }
         }
         else
         {
            if( maskFile != "" )
            {
               psfsub[n].sigmaMean(tfinim, maskCube, sigmaThreshold, m_minGoodFract);
            }
            else
            {
               psfsub[n].sigmaMean(tfinim, sigmaThreshold);
            }
         }
         finim.image(n) = tfinim;
      }
      else if(method == HCI::debug)
      {
         finim.image(n) = psfsub[n].image(0);
      }
   }

   t_combo_end = get_curr_time();
}


template<typename _realT>
inline void HCIobservation<_realT>::writeFinim(fitsHeader * addHead)
{
   std::string fname = finimName;

   if(outputDir != "")
   {
      fname = outputDir + "/" + fname;
   }

   if(!exactFinimName)
   {
      fname = ioutils::getSequentialFilename(fname, ".fits");
   }

   fitsHeader head;

   head.append("", fitsCommentType(), "----------------------------------------");
   head.append("", fitsCommentType(), "mx::HCIobservation parameters:");
   head.append("", fitsCommentType(), "----------------------------------------");

   head.append<int>("FDELFRNT", deleteFront, "images deleted from front of file list");
   head.append<int>("FDELBACK", deleteBack, "images deleted from back of file list");


   head.append("QFILE", basename(qualityFile.c_str()), "quality file for thresholding");
   head.append<realT>("QTHRESH", qualityThreshold, "quality threshold");
   head.append<int>("NUMIMS", Nims, "number of images processed");

   head.append<int>("IMSIZE", imSize, "image size after reading");

   head.append<int>("COADMTHD", coaddCombineMethod, "coadd combination method");
   if(coaddCombineMethod != HCI::noCombine)
   {
      head.append<int>("COADIMNO", coaddMaxImno, "max number of images in each coadd");
      head.append<int>("COADTIME", coaddMaxTime, "max time  in each coadd");
   }

   head.append("MASKFILE", maskFile, "mask file");

   head.append<int>("PPBEFORE", preProcess_beforeCoadd, "pre-process before coadd flag");
   head.append<int>("PPMASK", preProcess_mask, "pre-process mask flag");
   head.append<int>("PPSUBRAD", preProcess_subradprof, "pre-process subtract radial profile flag");
   head.append<realT>("PPAUSMAW", preProcess_azUSM_azW, "pre-process azimuthal USM azimuthal width");
   head.append<realT>("PPAUSMRW", preProcess_azUSM_radW, "pre-process azimuthal USM radial width");

   head.append<realT>("PPGUSMFW", preProcess_gaussUSM_fwhm, "pre-process Gaussian USM fwhm");

   head.append<int>("COMBMTHD", combineMethod, "combination method");

   if(weightFile != "")
      head.append("WEIGHTF", weightFile, "file containing weights for combination");


   if(combineMethod == HCI::sigmaMeanCombine)
      head.append<realT>("SIGMAT", sigmaThreshold, "threshold for sigma clipping");

   if(addHead)
   {
      head.append(*addHead);
   }

   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);

   fitsFile<realT> f;

   f.write(fname, finim, head);

   std::cerr << "Final image written to: " <<  fname << "\n";
}


template<typename _realT>
inline void HCIobservation<_realT>::outputPSFSub(fitsHeader * addHead)
{
   std::string fname;

   fitsHeader head;

   head.append("", fitsCommentType(), "----------------------------------------");
   head.append("", fitsCommentType(), "mx::HCIobservation parameters:");
   head.append("", fitsCommentType(), "----------------------------------------");

   head.append<int>("FDELFRNT", deleteFront, "image deleted from front of file list");
   head.append<int>("FDELBACK", deleteBack, "image deleted from back of file list");

   head.append<int>("IMSIZE", imSize, "image size after reading");

   head.append<int>("COADMTHD", coaddCombineMethod, "coadd combination method");
   if(coaddCombineMethod != HCI::noCombine)
   {
      head.append<int>("COADIMNO", coaddMaxImno, "max number of images in each coadd");
      head.append<int>("COADTIME", coaddMaxTime, "max time  in each coadd");
   }

   if(addHead)
   {
      head.append(*addHead);
   }

   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);

   fitsFile<realT> f;

   std::ofstream wout;
   
   if(comboWeights.size() > 0)
   {
      fname = PSFSubPrefix + "weights.dat";
      if(outputDir != "")
      {
         fname = outputDir + "/" + fname;
      }
      wout.open(fname);
      std::cerr << "writing comboWeights: " << fname << "\n";
   }
   
   char num[256];
   for(size_t n=0; n<psfsub.size(); ++n)
   {
      for(int p=0; p< psfsub[n].planes(); ++p)
      {
         snprintf(num, 256, "_%03zu_%05d.fits",n,p);
         fname = PSFSubPrefix + num;

         if(outputDir != "")
         {
            fname = outputDir + "/" + fname;
         }
   
         fitsHeader h = head;

         h.append(heads[p]);
         f.write(fname, psfsub[n].image(p).data(), psfsub[n].rows(), psfsub[n].cols(), 1, h);
         
         if(comboWeights.size() > 0 && n == 0)
         {
            wout << fname << " " << comboWeights[p] << "\n";
         }
      }
   }

   if(comboWeights.size() > 0)
   {
      wout.close();
   }
}

///@}

} //namespace improc
} //namespace mx

#endif //__HCIobservation_hpp__
