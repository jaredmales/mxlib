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

#include <sys/stat.h>

#include "../mxlib.hpp"

#include "../mxException.hpp"

#include "../math/templateBLAS.hpp"
#include "../sys/timeUtils.hpp"
#include "../ioutils/fileUtils.hpp"
#include "../ioutils/readColumns.hpp"
#include "../ioutils/fits/fitsFile.hpp"
#include "../ipc/ompLoopWatcher.hpp"

#include "eigenImage.hpp"
#include "eigenCube.hpp"
#include "imageFilters.hpp"
#include "imageMasks.hpp"
#include "imageTransforms.hpp"
#include "imageUtils.hpp"

namespace mx
{

namespace improc
{



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
                        
   /// Get the string name of the combineMethod integer
   /**
     * \returns a string with the name of the combineMethod 
     */ 
   std::string combineMethodStr( int method /**< [in] one of the combineMethods enum members */);
 
   /// Get the combineMethod from the corresponding string name
   /**
     * \returns the combineMethods enum member corresponding to the string name.
     */ 
   int combineMethodFmStr( const std::string & method /**< [in] the string name of the combineMethod */);
   
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
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;

   ///\name Construction and Initialization
   /** @{
     */
   ///Default c'tor
   HCIobservation();

   ///Construct and load the target file list.
   /** Populates the \ref m_fileList vector by searching on disk for files which match
     * "dir/prefix*.ext".  See \ref load_fileList
     *
     */
   HCIobservation( const std::string &dir,     ///< [in] the directory to search.
                   const std::string &prefix,  ///< [in] the initial part of the file name.  Can be empty "".
                   const std::string &ext      ///< [in] the extension to append to the file name, must include the '.'.
                 );

   /// Construct using a file containing the target file list
   /** Populates the \ref m_fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     */
   explicit HCIobservation( const std::string & fileListFile /**< [in] a file name path to read.*/ );

   ///Construct and load the target file list and the RDI file list
   /** Populates the \ref m_fileList vector by searching on disk for files which match
     * "dir/prefix*.ext".  See \ref load_fileList
     *
     *  Populates the \ref m_RDIfileList vector by searching on disk for files which match
     * "RDIdir/RDIprefix*.RDIext".  See \ref load_RDIfileList
     */
   HCIobservation( const std::string &dir,       ///< [in] the directory to search.
                   const std::string &prefix,    ///< [in] the initial part of the file name.  Can be empty "".
                   const std::string &ext,       ///< [in] the extension to append to the file name, must include the '.'.
                   const std::string &RDIdir,    ///< [in] the directory to search for the reference files.
                   const std::string &RDIprefix, ///< [in] the initial part of the file name for the reference files.  Can be empty "".
                   const std::string &RDIext=""  ///< [in] [optional] the extension to append to the RDI file name, must include the '.'.  If empty "" then same extension as target files is used.
                 );

   /// Construct using a file containing the target file list and a file containing the RDI target file list
   /** Populates the \ref m_fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     * 
     * Populates the \ref m_RDIfileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     */
   HCIobservation( const std::string & fileListFile,   ///< [in] a file name path to read for the target file names.
                   const std::string & RDIfileListFile ///< [in] a file name path to read for the reference file names.
                 );
   
   /// Load the file list
   /** Populates the \ref m_fileList vector by searching on disk for files which match the given parameters.
     * Uses \ref mx::getFileNames to search for all files which match "dir/prefix*.ext".
     *
     */
   void load_fileList( const std::string &dir,    ///< [in] the directory to search.
                       const std::string &prefix, ///< [in] the initial part of the file name.  Can be empty "".
                       const std::string &ext     ///< [in] the extension to append to the file name, which must include the '.'. Can be empty "".
                     );

   /// Load the file list from a file
   /** Populates the \ref m_fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     *
     */
   void load_fileList( const std::string & fileListFile  /**< [in] a file name path to read.*/ );

      /// Load the RDI basis file list
   /** Populates the \ref m_RDIfileList vector by searching on disk for files which match the given parameters.
     * Uses \ref mx::getFileNames to search for all files which match "dir/prefix*.ext".
     *
     */
   void load_RDIfileList( const std::string &dir,    ///< [in] the directory to search.
                          const std::string &prefix, ///< [in] the initial part of the file name.  Can be empty "".
                          const std::string &ext     ///< [in] the extension to append to the file name, which must include the '.'. Can be empty "".
                        );
   
   /// Load the RDI basis file list from a file
   /** Populates the \ref m_fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     *
     */
   void load_RDIfileList( const std::string & fileListFile  /**< [in] a file name path to read.*/ );
   
   
   ///@}
   
   
   ///\name The Input Target Images
   /** @{
    */

   ///The target image cube
   eigenCube<realT> m_tgtIms;

   int m_Nims {0};  ///<Number of images in m_tgtIms
   int m_Nrows {0}; ///<Number of rows of the images in m_tgtIms
   int m_Ncols {0}; ///<Number of columns of the images in m_tgtIms
   int m_Npix {0};  ///<Pixels per image, that is Nrows*Ncols

   ///Vector of target image times, in MJD.
   std::vector<double> m_imageMJD;

   ///Vector of FITS headers, one per file, populated with the values for the keywords.
   std::vector<fits::fitsHeader> m_heads;

   ///Whether or not the m_fileList has been read.
   bool m_filesRead {false};

   ///Whether or not the specified files have been deleted from m_fileList
   bool m_filesDeleted {false};

   ///Vector to hold the image weights read from the m_weightFile.
   /** After readWeights is executed by readFiles, this will contain the normalized weights.
     * \todo check how comboWeights are handled in coadding 
     */
   std::vector<realT> m_comboWeights;

   ///@}
   
   ///\name The Input Reference Images
   /** @{
    */

   ///The optional reference image cube
   eigenCube<realT> m_refIms;

   ///Vector of reference image times, in MJD.
   std::vector<double> m_RDIimageMJD;

   ///Vector of FITS headers, one per reference file, populated with the values for the keywords.
   std::vector<fits::fitsHeader> m_RDIheads;

   ///Whether or not the reference files have been read.
   bool m_RDIfilesRead {false};

   ///Whether or not the specified files have been deleted from m_RDIfileList
   bool m_RDIfilesDeleted {false};

   ///@}

   ///\name The Reduced Data
   /** @{
     */
   ///The PSF subtracted images
   /** This is a vector of cubes so that it can contain results from different reductions,
     * e.g. different modes when using KLIP.
     */
   std::vector<eigenCube<realT> > m_psfsub;

   ///The final combined images, one for each cube in psfsub.
   eigenCube<realT> m_finim;

   ///@}
   
   
   /** \name Target File Reading
     * Options to control which files are read, how they are read, what meta data is extracted
     * from FITS headers, sizing and masking, etc.
     * @{
     */

   ///The list of files to read in.
   /** This can be set on construction or by calling \ref load_fileList
     */
   std::vector<std::string> m_fileList;

   ///Specify how many files from m_fileList to delete from the front of the list
   int m_deleteFront {0};

   ///Specify how many files from m_fileList to delete from the back of the list
   int m_deleteBack {0};

   ///File containing 2 space-delimited columns of fileVame qualityValue pairs.
   /** If this is not empty and \ref qualityThreshold is > 0, then only images where
     * qualityValue >= qualityThreshold are read and processed.
     *
     * The only restriction on qualityThreshold is that it is > 0.  It is intendend to be
     * something like Strehl ratio.
     */
   std::string m_qualityFile;

   ///Threshold to apply to qualityValues read from \ref qualityFile.
   /** If <= 0, then thresholding is not performed.
    */
   realT m_qualityThreshold {0};

   ///Just prints the names and qualities of the files which pass threshold, and stop.
   /** Useful mainly for debugging.
     */
   bool m_thresholdOnly {false};

   ///Name of the keyword to use for the image date.
   /** Specifies the keyword corresponding to the date.  This is
     * the "DATE" keyword for file write time, and usually "DATE-OBS" for actual observation time.
     *
     * Default is "DATE-OBS".
     *
     * If empty "", then image date is not read.
     */
   std::string m_MJDKeyword {"DATE-OBS"};

   ///Whether or not the date is in ISO 8601 format
   bool m_MJDisISO8601 {true};

   ///If the date is not ISO 8601, this specifies the conversion to Julian Days (i.e. seconds to days)
   realT m_MJDUnits {1.0};

   ///Vector of FITS header keywords to read from the files in m_fileList.
   std::vector<std::string> m_keywords;

   ///Set the image size.  Images are cut down to this size after reading.
   /** Set to \<= 0 to use images uncut.
     *
     * Image sizes are not increased if this is larger than their size on disk.
     */
   int m_imSize {0};

   ///Read the list of files, cut to size, and preprocess.
   /**
     * \returns 0 on success, -1 on  error.
     */
   int readFiles();

   /// Perform post-read actions for the target images, for use by derived classes
   virtual int postReadFiles();

   /// Perform post-coadd actions for the target images, for use by derived classes.
   /** 
     * \returns 0 on success
     * \returns \<0 on error.
     */
   virtual int postCoadd();
   
   ///@}

   /** \name Reference File Reading
     * For RDI, Options to control which files are read, how they are read, what meta data is extracted
     * from FITS headers, sizing and masking, etc.
     * @{
     */

   ///The list of files to read in for the RDI basis.
   /** This is set by calling \ref load_RDIfileList
     */
   std::vector<std::string> m_RDIfileList;

   ///Specify how many files from m_RDIfileList to delete from the front of the list
   int m_RDIdeleteFront {0};

   ///Specify how many files from m_RDIfileList to delete from the back of the list
   int m_RDIdeleteBack {0};

   ///File containing 2 space-delimited columns of fileMame qualityValue pairs for the reference images.
   /** If this is not empty and \ref m_RDIqualityThreshold is > 0, then only images where
     * qualityValue >= qualityThreshold are read and processed.
     *
     * The only restriction on m_RDIqualityThreshold is that it is > 0.  It is intendend to be
     * something like Strehl ratio.
     */
   std::string m_RDIqualityFile;

   ///Threshold to apply to qualityValues read from \ref qualityFile.
   /** If <= 0, then thresholding is not performed.
    */
   realT m_RDIqualityThreshold {0};

   ///Vector of FITS header keywords to read from the files in m_fileList.
   std::vector<std::string> m_RDIkeywords;

   /// Read the list of reference files, cut to size, and preprocess.
   /** The target files must be read with \ref readFiles() before calling this method.
     * 
     * \returns 0 on success, -1 on  error.
     */
   int readRDIFiles();

   ///Perform post-read actions for the RDI images, for use by derived classes
   virtual int postRDIReadFiles();

   ///Perform post-coadd actions, for use by derived classes.
   /** A key example is to update keywords after the averaging occurs in coaddImages().
     *
     * \returns 0 on success
     * \returns \<0 on error.
     */
   virtual int postRDICoadd();
   
   
   ///@}
   //--RDI
   
   /** \name Thresholding
     * Thresholds are applied to a list of files before it is read, based on the image qualities supplied.
     * @{
     */
   
   /// Read the image qualities from a qualityFile and apply the threshold to a fileList
   /** This is called by readFiles().
     *
     * \returns 0 on success, -1 on  error.
     */
   int threshold( std::vector<std::string> & fileList, ///< [in/out] the fileList to threshold
                  const std::string & qualityFile,     ///< [in] the path to the file containing qualities, one per file.
                  realT qualityThreshold               ///< [in] the quality threshold to apply
                );
   
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
   int m_coaddCombineMethod {HCI::noCombine};

   ///Maximum number of images to coadd
   int m_coaddMaxImno {0};

   ///Maximum elapsed time over which to coadd the images.
   realT m_coaddMaxTime {0};

   ///The values of these keywords will be averaged and replaced.
   std::vector<std::string> m_coaddKeywords;

   ///Coadd the images
   void coaddImages( int coaddCombineMethod,
                     int coaddMaxImno,
                     int coaddMaxTime,
                     std::vector<std::string> & coaddKeywords,
                     std::vector<double> & imageMJD,
                     std::vector<fits::fitsHeader> & heads,
                     eigenCube<realT> & ims
                   );

   ///@} -- coadding


   /** \name Masking
     * A 1/0 mask can be supplied, which is used in pre-processing and in image combination.
     * @{
     */

   ///Specify a mask file to apply
   /**No mask is applied if this is empty.
     */
   std::string m_maskFile;
   
   eigenImageT m_mask; ///< The mask
   
   eigenCube<realT> m_maskCube; ///< A cube of masks, one for each input image, which may be modified versions (e.g. rotated) of mask.
   
   ///Read the mask file, resizing to imSize if needed.
   void readMask();
   
   ///Populate the mask cube which is used for post-processing.  Derived classes can do this as appropriate, e.g. by rotating the mask.
   virtual void makeMaskCube();


   ///@}



public:

   /** \name Pre-Processing
     * These options control the pre-processing masking and filtering.
     * They are performed in the following order:
     * -# mask applied (enabled by m_preProcess_mask
     * -# radial profile subtraction (enabled by m_preProcess_subradprof)
     * -# mask applied (enabled by m_preProcess_mask
     * -# symmetric Gaussian unsharp mask (m_preProcess_gaussUSM_fwhm)
     * -# mask applied (enabled by m_preProcess_mask
     * -# azimuthal unsharp mask (m_preProcess_azUSM_azW, and m_preProcess_azUSM_radW)
     * -# mask applied (enabled by m_preProcess_mask)
     * @{
     */

   bool m_skipPreProcess {false}; ///<Don't do any of the pre-processing steps (including coadding).

   bool m_preProcess_beforeCoadd {false}; ///<controls whether pre-processing takes place before or after coadding

   bool m_preProcess_mask {true}; ///<If true, the mask is applied during each pre-processing step.

   bool m_preProcess_subradprof {false}; ///<If true, a radial profile is subtracted from each image.

   /// Azimuthal boxcar width for azimuthal unsharp mask [pixels]
   /** If this is 0 then azimuthal-USM is not performed.
     */
   realT m_preProcess_azUSM_azW {0};

   /// Mazimum azimuthal boxcar width for azimuthal unsharp mask [degrees]
   /** Limits width close to center, preventing wrap-around.  Default is 45 degrees.  Set to 0 for no maximum.
     */
   realT m_preProcess_azUSM_maxAz {45};

   /// Radial boxcar width for azimuthal unsharp mask [pixels]
   /** If this is 0 then azimuthal-USM is not performed.
     */
   realT m_preProcess_azUSM_radW {0};

   ///Kernel FWHM for symmetric unsharp mask (USM)
   /** USM is not performed if this is 0.
    */
   realT m_preProcess_gaussUSM_fwhm {0};

   ///Set path and file prefix to output the pre-processed images.
   /** If empty, then pre-processed images are not output.
     */
   std::string m_preProcess_outputPrefix;

   /// If true, then we stop after pre-processing.
   bool m_preProcess_only {false};

   ///Do the pre-processing
   void preProcess( eigenCube<realT> & ims /**< [in] the image cube, should be either m_tgtIms or m_refIms */);

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
     * - HCI::weightedMeanCombine -- final image is the weighted mean.  m_weightFile must be provided.
     * - HCI::sigmaMeanCombine -- final image is sigma clipped mean.  If m_sigmaThreshold \<= 0, then it reverts to meanCombine.
     */
   int m_combineMethod {HCI::meanCombine};

   ///Specifies a file containing the image weights, for combining with weighted mean.
   /** This 2-column space-delimited ASCII file containing  filenames and weights. It must be specified before readFiles()
     * is executed.  In readFiles this is loaded after the m_fileList is cutdown and matched to the remaining files.
     */
   std::string m_weightFile;



   /// The standard deviation threshold used if combineMethod == HCI::sigmaMeanCombine.
   realT m_sigmaThreshold {0};

   /// The minimum fraction of good (un-masked) pixels to include in the final combination (0.0 to 1.0). If not met, then the pixel will be NaN-ed.
   realT m_minGoodFract {0.0};

   ///Read the image weights from m_weightFile
   /** This is called by readFiles().
     *
     * \returns 0 on success, -1 on  error.
     */
   int readWeights();

   ///Combine the images into a single final image.
   /** Images are combined by the method specified in \ref combineMethod
     */
   void combineFinim();
   
   ///@}

   /** \name Output
     * These options control the ouput of the final combined images and the individual PSF subtracted images.
     * @{
     */

   ///Set whether the final combined image is written to disk
   int m_doWriteFinim {1};

   ///The directory where to write output files.
   std::string m_outputDir;

   ///The base file name of the output final image
   /** The complete name is formed by combining with a sequential number and the ".fits" extension.
     * that is: m_finimName0000.fits.  This behavior can be modified with m_exactFinimName.
     */
   std::string m_finimName {"finim_"};

   ///Use m_finimName exactly as specified, without appending a number or an extension.
   /** Output is still FITS format, regardless of extension.  This will overwrite
     * an existing file without asking.
     */
   bool m_exactFinimName {false};

   ///Controls whether or not the individual PSF subtracted images are written to disk.
   /** - true -- write to disk
     * - false -- [default] don't write to disk
     */
   bool m_doOutputPSFSub {false};

   ///Prefix of the FITS file names used to write individual PSF subtracted images to disk if m_doOutputPSFSub is true.
   std::string m_PSFSubPrefix;
   
   /// Output the pre-processed target images
   void outputPreProcessed();

   ///Output the pre-processed reference images
   void outputRDIPreProcessed();
   
   /// Fill in the HCIobservation standard FITS header 
   /**
     */
   void stdFitsHeader(fits::fitsHeader & head /**< [in/out] the fistHeader structure which will have cards appended to it. */);
   
   ///Write the final combined image to disk
   /**
     */
   void writeFinim(fits::fitsHeader * addHead = 0);
   
   ///Write the PSF subtracted images to disk
   /**
    */
   void outputPSFSub(fits::fitsHeader * addHead = 0);


   ///@}



   /// Read in already PSF-subtracted files
   /** Used to take up final processing after applying some non-klipReduce processing steps to
     * PSF-subtracted images.
     */ 
   int readPSFSub( const std::string & dir,
                   const std::string & prefix,
                   const std::string & ext,
                   size_t nReductions 
                 );
   


   









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

// -- construction and initialization

template<typename _realT>
HCIobservation<_realT>::HCIobservation()
{
}

template<typename _realT>
HCIobservation<_realT>::HCIobservation( const std::string & dir, 
                                        const std::string & prefix, 
                                        const std::string & ext
                                      )
{
   load_fileList(dir, prefix, ext);
}

template<typename _realT>
HCIobservation<_realT>::HCIobservation(const std::string & fileListFile)
{
   load_fileList(fileListFile);
}


template<typename _realT>
HCIobservation<_realT>::HCIobservation( const std::string & dir, 
                                        const std::string & prefix, 
                                        const std::string & ext,
                                        const std::string & RDIdir, 
                                        const std::string & RDIprefix, 
                                        const std::string & RDIext
                                      )
{
   std::cerr << "HCI 6\n";
   
   load_fileList(dir, prefix, ext);
   
   std::string re;
   if(RDIext == "") re = ext;
   else re = RDIext;
   
   load_RDIfileList(RDIdir, RDIprefix, re);
}



template<typename _realT>
HCIobservation<_realT>::HCIobservation( const std::string & fileListFile,
                                        const std::string & RDIfileListFile
                                      )
{
   load_fileList( fileListFile);
   load_RDIfileList( RDIfileListFile);
}

template<typename _realT>
void HCIobservation<_realT>::load_fileList( const std::string & dir, 
                                            const std::string & prefix, 
                                            const std::string & ext
                                          )
{
   m_fileList = ioutils::getFileNames(dir, prefix, "", ext);
   m_filesDeleted = false;
}

template<typename _realT>
void HCIobservation<_realT>::load_fileList(const std::string & fileListFile)
{
   ioutils::readColumns(fileListFile, m_fileList);
   m_filesDeleted = false;
}

template<typename _realT>
void HCIobservation<_realT>::load_RDIfileList( const std::string & dir, 
                                                      const std::string & prefix, 
                                                      const std::string & ext
                                                    )
{
   m_RDIfileList = ioutils::getFileNames(dir, prefix, "", ext);
   m_RDIfilesDeleted = false;
}

template<typename _realT>
void HCIobservation<_realT>::load_RDIfileList( const std::string & fileListFile )
{
   ioutils::readColumns(fileListFile, m_RDIfileList);
   m_RDIfilesDeleted = false;
}

// --< construction and initialization


template<typename _realT>
int HCIobservation<_realT>::readFiles()
{
   if(m_fileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The target fileList has 0 length, there are no files to be read.");
      return -1;
   }

   //First make the list deletions
   if(!m_filesDeleted)
   {
      if(m_deleteFront > 0)
      {
         m_fileList.erase(m_fileList.begin(), m_fileList.begin()+m_deleteFront);
      }

      if(m_deleteBack > 0)
      {
         m_fileList.erase(m_fileList.end()-m_deleteBack, m_fileList.end());
      }
      m_filesDeleted = true;
   }

   if( m_fileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The target fileList has 0 length, there are no files to be read after deletions.");
      return -1;
   }

   if(m_qualityFile != "")
   {
      std::cerr << "Thresholding target images...";
      size_t origsize = m_fileList.size();
      
      if( threshold(m_fileList, m_qualityFile, m_qualityThreshold) < 0) return -1;

      if( m_fileList.size() == 0)
      {  
         mxError("HCIobservation", MXE_FILENOTFOUND, "The fileList has 0 length, there are no files to be read after thresholding.");
         return -1;
      }
      
      std::cerr << "Done.  Selected " << m_fileList.size() << " out of " << origsize << "\n";

      
      if(m_thresholdOnly)
      {
         std::cout << "#Files which passed thresholding:\n";
         for(size_t i=0; i<m_fileList.size(); ++i)
         {
            std::cout << m_fileList[i] << "\n";
         }

         exit(0);
      }
   }


   if(m_weightFile != "")
   {
      if(readWeights() < 0) return -1;
   }

   /*----- Append the HCI keywords to propagate them if needed -----*/
   
   fits::fitsHeader head;

   if(m_MJDKeyword != "") head.append(m_MJDKeyword);

   for(size_t i=0;i<m_keywords.size();++i)
   {
      if(head.count(m_keywords[i]) == 0) head.append(m_keywords[i]);
   }

   m_heads.clear(); //This is necessary to make sure heads.resize() copies head on a 2nd call
   m_heads.resize(m_fileList.size(), head);

   /*----- Read in first image to test size -----*/

   Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

   fits::fitsFile<realT> f(m_fileList[0]);

   f.read(im);

   //We set imSize to match the first image, but we make it a square.
   if(m_imSize == 0)
   {
      m_imSize = im.rows();
      if(m_imSize > im.cols()) m_imSize = im.cols();
   }
   else
   {
      //Now make sure we don't read too much.
      if(m_imSize > im.rows()) m_imSize = im.rows();
      if(m_imSize > im.cols()) m_imSize = im.cols();
   }
      
   //And now set the read size so we only read what we want.   
   //the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
   f.setReadSize( floor(0.5*(im.rows()-1) - 0.5*(m_imSize-1) +0.1), floor(0.5*(im.cols()-1.0) - 0.5*(m_imSize-1.0)+0.1), m_imSize, m_imSize);
   im.resize(m_imSize, m_imSize);
   
   
   m_tgtIms.resize(im.rows(), im.cols(), m_fileList.size());

   t_load_begin = sys::get_curr_time();

   f.read(m_tgtIms.data(), m_heads, m_fileList);

   f.setReadSize();

   /* read in the image timestamps, depending on how MJD is stored in the header */
   if(m_MJDKeyword != "")
   {
      m_imageMJD.resize(m_heads.size());

      if(m_MJDisISO8601)
      {
         for(size_t i=0;i<m_imageMJD.size();++i)
         {
            m_imageMJD[i] =  sys::ISO8601date2mjd(m_heads[i][m_MJDKeyword].String());
         }
      }
      else
      {
         for(size_t i=0;i<m_imageMJD.size();++i)
         {
            m_imageMJD[i] = m_heads[i][m_MJDKeyword].template value<realT>()*m_MJDUnits;
         }
      }
   }

   t_load_end = sys::get_curr_time();

   m_Nims =  m_tgtIms.planes();
   m_Nrows = m_tgtIms.rows();
   m_Ncols = m_tgtIms.cols();
   m_Npix =  m_tgtIms.rows()*m_tgtIms.cols();

   std::cerr << "loading complete\n";
   
   std::cerr << "zero-ing NaNs\n";
   zeroNaNCube(m_tgtIms);
   std::cerr << "done\n";
   
   /*** Now do the post-read actions ***/
   if( postReadFiles() < 0) return -1;

   /*** Read in the mask if present ***/
   readMask();
   

   /*** Now begin processing ***/
   if(!m_skipPreProcess)
   {
      /*** Now do any pre-processing ***/
      if(m_preProcess_beforeCoadd) preProcess(m_tgtIms);

      if(m_coaddCombineMethod != HCI::noCombine)
      {
         std::cerr << "Coadding target images...\n";
         coaddImages(m_coaddCombineMethod, m_coaddMaxImno, m_coaddMaxTime, m_coaddKeywords, m_imageMJD, m_heads, m_tgtIms);
         
         m_Nims =  m_tgtIms.planes();
         m_Nrows = m_tgtIms.rows();
         m_Ncols = m_tgtIms.cols();
         m_Npix =  m_tgtIms.rows()*m_tgtIms.cols();
   
         if( postCoadd() < 0)
         {
            std::cerr << "Post coadd error " << __FILE__ << " " << __LINE__ << "\n";
            return -1;
         }
         std::cerr << "Done.\n";
         
         //Re-make the mask cube if we coadded...
         if( m_maskFile != "")
         {
            makeMaskCube();            
         }
      }

      /*** Now do any pre-processing if not done already***/
      if(!m_preProcess_beforeCoadd) preProcess(m_tgtIms);

      outputPreProcessed();
   }
   
   
   m_filesRead = true;
   
   return 0;
} //readFiles()

template<typename _realT>
int HCIobservation<_realT>::postReadFiles()
{
   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::postCoadd()
{
   return 0;
}

//------------------- readRDIFiles
template<typename _realT>
int HCIobservation<_realT>::readRDIFiles()
{
   
   /* First check if the target files have been read */
   if(m_Nrows == 0 || m_Ncols == 0)
   {
      mxError("HCIobservation", MXE_PARAMNOTSET, "The target image size must be set before reading RDI files.");
      return -1;
   }
   
   if(m_RDIfileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The RDI fileList has 0 length, there are no files to be read.");
      return -1;
   }

   //First make the list deletions
   if(!m_RDIfilesDeleted)
   {
      if(m_RDIdeleteFront > 0)
      {
         m_RDIfileList.erase(m_RDIfileList.begin(), m_RDIfileList.begin()+m_RDIdeleteFront);
      }

      if(m_RDIdeleteBack > 0)
      {
         m_RDIfileList.erase(m_RDIfileList.end()-m_RDIdeleteBack, m_RDIfileList.end());
      }
      m_RDIfilesDeleted = true;
   }

   if( m_RDIfileList.size() == 0)
   {
      mxError("HCIobservation", MXE_FILENOTFOUND, "The RDI fileList has 0 length, there are no files to be read after deletions.");
      return -1;
   }

   if(m_RDIqualityFile != "")
   {
      std::cerr << "Thresholding RDI images...";
      size_t origsize = m_RDIfileList.size();
      
      if( threshold(m_RDIfileList, m_RDIqualityFile, m_RDIqualityThreshold) < 0) return -1;
      
      if( m_RDIfileList.size() == 0)
      {
         mxError("HCIobservation", MXE_FILENOTFOUND, "The fileList has 0 length, there are no files to be read after thresholding.");
         return -1;
      }
   
      std::cerr << "Done.  Selected " << m_RDIfileList.size() << " out of " << origsize << "\n";

   }
   

   /*----- Append the HCI keywords to propagate them if needed -----*/

   fits::fitsHeader head;

   if(m_MJDKeyword != "") head.append(m_MJDKeyword); //Currently assuming the MJD keyword will be the same

   for(size_t i=0;i<m_RDIkeywords.size();++i)
   {
      head.append(m_RDIkeywords[i]);
   }

   m_RDIheads.clear(); //This is necessary to make sure heads.resize() copies head on a 2nd call
   m_RDIheads.resize(m_RDIfileList.size(), head);


   /*----- Read in first image to adjust size ----*/
   Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

   fits::fitsFile<realT> f(m_RDIfileList[0]);

   f.read(im);

   if(im.rows() < m_imSize || im.cols() < m_imSize)
   {
      mxError("HCIobservation", MXE_SIZEERR, "The reference images are too small, do not match the target images.");
      return -1;
   }

   //And now set the read size so we only read what we want.   
   //the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
   f.setReadSize( floor(0.5*(im.rows()-1) - 0.5*(m_imSize-1) +0.1), floor(0.5*(im.cols()-1.0) - 0.5*(m_imSize-1.0)+0.1), m_imSize, m_imSize);
   
   m_refIms.resize(m_imSize, m_imSize, m_RDIfileList.size());

   t_load_begin = sys::get_curr_time();

   f.read(m_refIms.data(), m_RDIheads, m_RDIfileList);

   f.setReadSize();

   if(m_MJDKeyword != "")
   {
      m_RDIimageMJD.resize(m_RDIheads.size());

      if(m_MJDisISO8601)
      {
         for(size_t i=0;i<m_RDIimageMJD.size();++i)
         {
            m_RDIimageMJD[i] =  sys::ISO8601date2mjd(m_RDIheads[i][m_MJDKeyword].String());
         }
      }
      else
      {
         for(size_t i=0;i<m_RDIimageMJD.size();++i)
         {
            m_RDIimageMJD[i] = m_RDIheads[i][m_MJDKeyword].template value<realT>()*m_MJDUnits;
         }
      }
   }

   t_load_end = sys::get_curr_time();

   std::cerr << "loading complete\n";
   
   std::cerr << "zero-ing NaNs\n";
   zeroNaNCube(m_refIms);
   std::cerr << "Done.\n";
   
   /*** Now do the post-read actions ***/
   if( postRDIReadFiles() < 0) return -1;

   /*** Now begin processing ***/
   if(!m_skipPreProcess)
   {
      /*** Now do any pre-processing ***/
      if(m_preProcess_beforeCoadd) preProcess(m_refIms);

      if(m_coaddCombineMethod != HCI::noCombine)
      {
         std::cerr << "Coadding reference images...\n";
         coaddImages(m_coaddCombineMethod, m_coaddMaxImno, m_coaddMaxTime, m_coaddKeywords, m_RDIimageMJD, m_RDIheads, m_refIms);
         
         if( postRDICoadd() < 0)
         {
            std::cerr << "Post coadd error " << __FILE__ << " " << __LINE__ << "\n";
            return -1;
         }
         std::cerr << "Done.\n";         
      }

      /*** Now do any pre-processing if not done already***/
      if(!m_preProcess_beforeCoadd) preProcess(m_refIms);

      //outputRDIPreProcessed();
   }
   
   
   m_RDIfilesRead = true;
   
   return 0;
} //readRDIFiles()

template<typename _realT>
int HCIobservation<_realT>::postRDIReadFiles()
{
   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::postRDICoadd()
{
   return 0;
}

template<typename _realT>
int HCIobservation<_realT>::threshold( std::vector<std::string> & fileList,
                                       const std::string & qualityFile,
                                       realT qualityThreshold
                                     )
{
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
   for(size_t i=0;i<qfileNames.size();++i) quality[ioutils::pathFilename(qfileNames[i].c_str())] = imQ[i];

   realT q;

   for(size_t i=0; i<fileList.size(); ++i)
   {
      try
      {
         q = quality.at(ioutils::pathFilename(fileList[i].c_str()));
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

   return 0;
}

template<typename _realT>
void HCIobservation<_realT>::coaddImages( int coaddCombineMethod,
                                          int coaddMaxImno,
                                          int coaddMaxTime,
                                          std::vector<std::string> & coaddKeywords,
                                          std::vector<double> & imageMJD,
                                          std::vector<fits::fitsHeader> & heads,
                                          eigenCube<realT> & ims
                                        )
{
   std::cerr << "***************************************************************\n";
   std::cerr << "                       *** WARNING ***                         \n";
   std::cerr << "       coadding is poorly tested.  proceed with caution.       \n";
   std::cerr << "***************************************************************\n";
   
   
   //Validate setup
   if(coaddMaxImno <=0 && coaddMaxTime <= 0) return;
 
   //Validate combine method
   if(coaddCombineMethod == HCI::noCombine) return;

   int Nims = ims.planes(); 
   int Nrows = ims.rows();
   int Ncols = ims.cols();
   
   t_coadd_begin = sys::get_curr_time();

   std::vector<eigenImageT> coadds;

   //We do all math here in double, to avoid losing precision
   std::vector<double> avgMJD;
   std::vector<std::vector<double> > avgVals;

   int combineMethod =  HCI::medianCombine;
   if( coaddCombineMethod == HCI::meanCombine) combineMethod = HCI::meanCombine;

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

   while(im0 < Nims)
   {
      //Initialize the accumulators
      initMJD = imageMJD[im0];

      for(size_t i=0;i< coaddKeywords.size(); ++i)
      {
         initVals[i] = heads[im0][coaddKeywords[i]].value<double>();
      }

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
         
         for(size_t i=0;i<coaddKeywords.size(); ++i)
         {
            initVals[i] += heads[imno][coaddKeywords[i]].value<double>();
         }
      }

      
      
      //And then turn them into an average
      initMJD /= (imF - im0);
      for(size_t i=0;i<coaddKeywords.size(); ++i)
      {
         initVals[i] /= (imF-im0);
      }

      //Extract the images into the temporary
      imsToCoadd.resize(Nrows, Ncols, imF-im0);
      for(int i =0; i < (imF-im0); ++i)
      {
         imsToCoadd.image(i) = ims.image(im0 + i);
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

   //Now resize ims and copy the coadds to the new cube
   ims.resize(Nrows, Ncols, coadds.size());
   Nims = coadds.size();

   for(int i=0;i<Nims;++i) ims.image(i) = coadds[i];

   //Now deal with imageMJD and headers
   imageMJD.erase(imageMJD.begin()+Nims, imageMJD.end());
   heads.erase(heads.begin()+Nims, heads.end());
   
   for(int i=0;i<Nims;++i)
   {
      imageMJD[i] = avgMJD[i];
      for(size_t j=0;j<coaddKeywords.size(); ++j)
      {
         heads[i][coaddKeywords[j]].value(avgVals[i][j]);
      }
   }
   
   t_coadd_end = sys::get_curr_time();

   
}//void HCIobservation<_realT>::coaddImages()

template<typename _realT>
void HCIobservation<_realT>::readMask()
{
   /*** Load the mask ***/
   if( m_maskFile != "")
   {
      std::cerr << "creating mask cube\n";
      fits::fitsFile<realT> ff;
      ff.read(m_mask, m_maskFile);
      
      ///\todo here re-size mask if needed to match imSize
      if(m_mask.rows() > m_imSize || m_mask.cols() > m_imSize)
      {
         eigenImageT tmask = m_mask.block( (int)(0.5*(m_mask.rows()-1) - 0.5*(m_imSize-1)), (int)(0.5*(m_mask.rows()-1) - 0.5*(m_imSize-1)), m_imSize, m_imSize);
         m_mask = tmask;
      }
      
      makeMaskCube();
   }

}

template<typename _realT>
void HCIobservation<_realT>::makeMaskCube()
{
   if( m_mask.rows() != m_Nrows || m_mask.cols() != m_Ncols)
   {
      std::cerr << "\nMask is not the same size as images.\n\n";
      exit(-1);
   }
   m_maskCube.resize( m_Nrows, m_Ncols, m_Nims);
   
   for(int i=0; i< m_Nims; ++i)
   {
      m_maskCube.image(i) = m_mask;
   }

}

template<typename _realT>
void HCIobservation<_realT>::preProcess( eigenCube<realT> & ims )
{
   t_preproc_begin = sys::get_curr_time();

   //The mask is applied first, and then after each subsequent P.P. step.
   if( m_maskFile != "" && m_preProcess_mask)
   {
      std::cerr << "Masking . . .\n";
      #pragma omp parallel for
      for(int i=0;i<ims.planes(); ++i)
      {
         ims.image(i) *= m_mask;
      }
      std::cerr << "done\n";
   }

   if( m_preProcess_subradprof )
   {
      std::cerr << "subtracting radial profile . . .\n";
      eigenImageT rp;

      for(int i=0;i<ims.planes(); ++i)
      {
         Eigen::Map<Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> > imRef(ims.image(i));
         radprofim(rp, imRef, true);
      }

      std::cerr << "done\n";

      if( m_maskFile != "" && m_preProcess_mask)
      {
         std::cerr << "Masking . . .\n";
         #pragma omp parallel for
         for(int i=0;i<ims.planes(); ++i)
         {
            ims.image(i) *= m_mask;
         }
         std::cerr << "done\n";
      }
   }

   if( m_preProcess_gaussUSM_fwhm > 0 && m_preProcess_mask)
   {
      std::cerr << "Applying Gauss USM . . .\n";
      t_gaussusm_begin = sys::get_curr_time();

      #pragma omp parallel for
      for(int i=0;i<ims.planes(); ++i)
      {
         eigenImageT fim, im;
         im = ims.image(i);
         filterImage(fim, im, gaussKernel<eigenImage<_realT>,2>(m_preProcess_gaussUSM_fwhm), 0.5*(ims.cols()-1) - m_preProcess_gaussUSM_fwhm*4);
         im = (im-fim);
         ims.image(i) = im;
      }

      if( m_maskFile != "")
      {
         std::cerr << "Masking . . .\n";
         #pragma omp parallel for
         for(int i=0;i<ims.planes(); ++i)
         {
            ims.image(i) *= m_mask;
         }

      }
      t_gaussusm_end = sys::get_curr_time();
      std::cerr << "done\n";
   }

   if( m_preProcess_azUSM_azW && m_preProcess_azUSM_radW )
   {
      ipc::ompLoopWatcher<> status( ims.planes(), std::cerr);

      std::cerr << "Applying azimuthal USM . . .\n";
      t_azusm_begin = sys::get_curr_time();
      #pragma omp parallel for
      for(int i=0;i<ims.planes(); ++i)
      {
         eigenImageT fim, im;
         im = ims.image(i);
         filterImage(fim, im, azBoxKernel<eigenImage<realT>>(m_preProcess_azUSM_radW, m_preProcess_azUSM_azW, m_preProcess_azUSM_maxAz));
         im = (im-fim);
         ims.image(i) = im;
         status.incrementAndOutputStatus();
      }

      if( m_maskFile != "" && m_preProcess_mask)
      {
         std::cerr << "Masking . . .\n";
         #pragma omp parallel for
         for(int i=0;i<ims.planes(); ++i)
         {
            ims.image(i) *= m_mask;
         }
      }

      t_azusm_end = sys::get_curr_time();
      std::cerr  << "done (" << t_azusm_end - t_azusm_begin << " sec)                                \n";
   }

   t_preproc_end = sys::get_curr_time();
   
} //void HCIobservation<_realT>::preProcess()

template<typename _realT>
int HCIobservation<_realT>::readWeights()
{
   std::ifstream fin;
   std::string str;

   if(m_weightFile == "")
   {
      mxError("HCIobservation::readWeights", MXE_PARAMNOTSET, "m_weightFile not set");
      return -1;
   }

   //Read the weight file and load it into a map
   std::vector<std::string> wfileNames;
   std::vector<realT> imW;

   if( ioutils::readColumns(m_weightFile, wfileNames, imW) < 0) return -1;

   if(imW.size() < m_fileList.size())
   {
      mxError("HCIobservation::readWeights", MXE_SIZEERR, "not enough weights specified");
      return -1;
   }

   std::map<std::string, realT> weights;
   for(size_t i=0;i<wfileNames.size();++i) weights[ioutils::pathFilename(wfileNames[i].c_str())] = imW[i];

   m_comboWeights.resize(m_fileList.size());

   realT wi;
   realT weightSum = 0;
   for(size_t i=0; i<m_fileList.size(); ++i)
   {
      try
      {
         wi = weights.at(ioutils::pathFilename(m_fileList[i].c_str()));
      }
      catch(...)
      {
         mxError("HCIobservation::readWeights", MXE_NOTFOUND, "Weight for a file in m_fileList not found.");
         return -1;
      }
      m_comboWeights[i] = wi;
      weightSum += wi;
   }

   //Finally normalize the weights
   for(size_t i=0; i< m_comboWeights.size(); ++i)
   {
      m_comboWeights[i] /= weightSum;
   }

   return 0;
} //int HCIobservation<_realT>::readWeights()

template<typename _realT>
void HCIobservation<_realT>::combineFinim()
{
   if(m_combineMethod == HCI::noCombine) return;

   t_combo_begin = sys::get_curr_time();

   //Validate the combineMethod setting
   int method = HCI::medianCombine;

   if(m_combineMethod == HCI::medianCombine)
   {
      method = HCI::medianCombine;
   }
   else if(m_combineMethod == HCI::meanCombine)
   {
      method = HCI::meanCombine;
   }
   else if(m_combineMethod == HCI::sigmaMeanCombine)
   {
      if(m_sigmaThreshold > 0)
      {
         method = HCI::sigmaMeanCombine;
      }
      else
      {
         method = HCI::meanCombine;
      }
   }
   else if(m_combineMethod == HCI::debug)
   {
      method = HCI::debug;
   }

   //Create and size temporary image for averaging
   eigenImageT tfinim;

   m_finim.resize(m_psfsub[0].rows(), m_psfsub[0].cols(), m_psfsub.size());

   //Now cycle through each set of psf subtractions
   for(size_t n= 0; n < m_psfsub.size(); ++n)
   {
      if(method == HCI::medianCombine)
      {
         m_psfsub[n].median(tfinim);
         m_finim.image(n) = tfinim;
      }
      else if(method == HCI::meanCombine)
      {
         if(m_comboWeights.size() == (size_t) m_Nims)
         {
            if( m_maskFile != "" )
            {
               m_psfsub[n].mean(tfinim, m_comboWeights, m_maskCube, m_minGoodFract);
            }
            else
            {
               m_psfsub[n].mean(tfinim, m_comboWeights);
            }
         }
         else
         {
            if( m_maskFile != "" )
            {
               m_psfsub[n].mean(tfinim, m_maskCube, m_minGoodFract);
            }
            else
            {
               m_psfsub[n].mean(tfinim);
            }
         }
         m_finim.image(n) = tfinim;
      }
      else if(method == HCI::sigmaMeanCombine)
      {
         if(m_comboWeights.size() == (size_t) m_Nims)
         {
            if( m_maskFile != "" )
            {
               m_psfsub[n].sigmaMean(tfinim, m_comboWeights, m_maskCube, m_sigmaThreshold, m_minGoodFract);
            }
            else
            {
               m_psfsub[n].sigmaMean(tfinim, m_comboWeights, m_sigmaThreshold);
            }
         }
         else
         {
            if( m_maskFile != "" )
            {
               m_psfsub[n].sigmaMean(tfinim, m_maskCube, m_sigmaThreshold, m_minGoodFract);
            }
            else
            {
               m_psfsub[n].sigmaMean(tfinim, m_sigmaThreshold);
            }
         }
         m_finim.image(n) = tfinim;
      }
      else if(method == HCI::debug)
      {
         m_finim.image(n) = m_psfsub[n].image(0);
      }
   }

   t_combo_end = sys::get_curr_time();
} //void HCIobservation<_realT>::combineFinim()

template<typename _realT>
void HCIobservation<_realT>::outputPreProcessed()
{
   if(m_preProcess_outputPrefix == "") return;

   std::string bname, fname;

   fits::fitsFile<_realT> ff;

   /** \todo Should add a HISTORY card here */
   for(int i=0; i< m_Nims; ++i)
   {
      bname = m_fileList[i];
      fname = m_preProcess_outputPrefix + ioutils::pathFilename(bname.c_str()) + ".fits";
      ff.write(fname, m_tgtIms.image(i).data(), m_Ncols, m_Nrows, 1, m_heads[i]);
   }
} //void HCIobservation<_realT>::outputPreProcessed()

template<typename _realT>
void HCIobservation<_realT>::stdFitsHeader(fits::fitsHeader & head)
{
   head.append("", fits::fitsCommentType(), "----------------------------------------");
   head.append("", fits::fitsCommentType(), "mx::HCIobservation parameters:");
   head.append("", fits::fitsCommentType(), "----------------------------------------");

   head.append<int>("FDELFRNT", m_deleteFront, "images deleted from front of file list");
   head.append<int>("FDELBACK", m_deleteBack, "images deleted from back of file list");


   head.append("QFILE", ioutils::pathFilename(m_qualityFile.c_str()), "quality file for thresholding");
   head.append<realT>("QTHRESH", m_qualityThreshold, "quality threshold");
   head.append<int>("NUMIMS", m_Nims, "number of images processed");

   head.append<int>("IMSIZE", m_imSize, "image size after reading");

   head.append<std::string>("COADMTHD", HCI::combineMethodStr(m_coaddCombineMethod), "coadd combination method");
   if(m_coaddCombineMethod != HCI::noCombine)
   {
      head.append<int>("COADIMNO", m_coaddMaxImno, "max number of images in each coadd");
      head.append<realT>("COADTIME", m_coaddMaxTime, "max time in each coadd");
   }
   else
   {
      head.append<int>("COADIMNO", 0, "max number of images in each coadd");
      head.append<realT>("COADTIME", 0, "max time in each coadd");
   }
   
   head.append("MASKFILE", m_maskFile, "mask file");

   head.append<int>("PREPROC BEFORE", m_preProcess_beforeCoadd, "pre-process before coadd flag");
   head.append<int>("PREPROC MASK", m_preProcess_mask, "pre-process mask flag");
   head.append<int>("PREPROC SUBRADPROF", m_preProcess_subradprof, "pre-process subtract radial profile flag");
   head.append<realT>("PREPROC AZUSM AZWIDTH", m_preProcess_azUSM_azW, "pre-process azimuthal USM azimuthal width [pixels]");
   head.append<realT>("PREPROC AZUSM MAXAZ", m_preProcess_azUSM_maxAz, "pre-process azimuthal USM maximum azimuthal width [degrees]");
   head.append<realT>("PREPROC AZUSM RADWIDTH", m_preProcess_azUSM_radW, "pre-process azimuthal USM radial width [pixels]");

   head.append<realT>("PREPROC GAUSSUSM FWHM", m_preProcess_gaussUSM_fwhm, "pre-process Gaussian USM fwhm [pixels]");

}

template<typename _realT>
void HCIobservation<_realT>::writeFinim(fits::fitsHeader * addHead)
{
   std::string fname = m_finimName;

   if(m_outputDir != "")
   {
      mkdir(m_outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   
      fname = m_outputDir + "/" + fname;
   }

   if(!m_exactFinimName)
   {
      fname = ioutils::getSequentialFilename(fname, ".fits");
   }

   fits::fitsHeader head;

   //Add HCIobservation standard header:
   stdFitsHeader(head);
   
   //Now add the final combination details:
   head.append<std::string>("COMBMTHD", HCI::combineMethodStr(m_combineMethod), "combination method");

   if(m_weightFile != "")
      head.append("WEIGHTF", m_weightFile, "file containing weights for combination");


   if(m_combineMethod == HCI::sigmaMeanCombine)
      head.append<realT>("SIGMAT", m_sigmaThreshold, "threshold for sigma clipping");

   if(addHead)
   {
      head.append(*addHead);
   }

   fits::fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);

   fits::fitsFile<realT> f;

   f.write(fname, m_finim, head);

   std::cerr << "Final image written to: " <<  fname << "\n";
} //void HCIobservation<_realT>::writeFinim(fits::fitsHeader * addHead)

template<typename _realT>
void HCIobservation<_realT>::outputPSFSub(fits::fitsHeader * addHead)
{
   
   
   std::string fname;

   fits::fitsHeader head;

   //Add the HCIobservation standard fits header
   stdFitsHeader(head);
   

   if(addHead)
   {
      head.append(*addHead);
   }

   fits::fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);

   fits::fitsFile<realT> f;

   std::ofstream wout;
   
   if(m_comboWeights.size() > 0)
   {
      fname = m_PSFSubPrefix + "weights.dat";
      if(m_outputDir != "")
      {
         mkdir(m_outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
         fname = m_outputDir + "/" + fname;
      }
      wout.open(fname);
      std::cerr << "writing comboWeights: " << fname << "\n";
   }
   
   char num[256];
   for(size_t n=0; n<m_psfsub.size(); ++n)
   {
      for(int p=0; p< m_psfsub[n].planes(); ++p)
      {
         snprintf(num, 256, "_%03zu_%05d.fits",n,p);
         fname = m_PSFSubPrefix + num;

         if(m_outputDir != "")
         {
            mkdir(m_outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            fname = m_outputDir + "/" + fname;
         }
   
         fits::fitsHeader h = head;

         h.append(m_heads[p]);
         
         
         f.write(fname, m_psfsub[n].image(p).data(), m_psfsub[n].rows(), m_psfsub[n].cols(), 1, h);
         
         if(m_comboWeights.size() > 0 && n == 0)
         {
            wout << fname << " " << m_comboWeights[p] << "\n";
         }
      }
   }

   if(m_comboWeights.size() > 0)
   {
      wout.close();
   }
} //void HCIobservation<_realT>::outputPSFSub(fits::fitsHeader * addHead)





template<typename _realT>
int HCIobservation<_realT>::readPSFSub( const std::string & dir,
                                        const std::string & prefix,
                                        const std::string & ext,
                                        size_t nReductions 
                                      )
{

   
   m_psfsub.resize(nReductions);
   
   //Load first file to condigure based on its header.
   std::vector<std::string> flist = ioutils::getFileNames(dir, prefix, "000", ext);
   fits::fitsHeader fh;
   eigenImage<realT> im;
   fits::fitsFile<realT> ff;
   
   ff.read(im, fh, flist[0]);
   
   if(fh.count("FDELFRNT") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "FDELFRNT not found in FITS header.");
      return -1;
   }
   m_deleteFront = fh["FDELFRNT"].value<int>();
   std::cerr << "deleteFront: " << m_deleteFront << "\n";
   
   if(fh.count("FDELBACK") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "FDELBACK not found in FITS header.");
      return -1;
   }
   m_deleteBack = fh["FDELBACK"].value<int>();
   std::cerr << "deleteBack: " << m_deleteBack << "\n";
   
   if(fh.count("QFILE") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "QFILE not found in FITS header.");
      return -1;
   }
   m_qualityFile = fh["QFILE"].String();
   std::cerr << "qualityFile: " << m_qualityFile << "\n";
   
   if(fh.count("QTHRESH") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "QTHRESH not found in FITS header.");
      return -1;
   }
   m_qualityThreshold = fh["QTHRESH"].value<realT>();
   std::cerr << "qualityThreshold: " << m_qualityThreshold << "\n";
   
   if(fh.count("COADMTHD") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "COADMTHD not found in FITS header.");
      return -1;
   }
   m_coaddCombineMethod = HCI::combineMethodFmStr(fh["COADMTHD"].String());
   std::cerr << "coaddCombineMethod: " << m_coaddCombineMethod << "\n";
   
   if(fh.count("COADIMNO") != 0)
   {
      m_coaddMaxImno = fh["COADIMNO"].value<int>();
      std::cerr << "coaddMaxImno: " << m_coaddMaxImno << "\n";
   }
   
   if(fh.count("COADTIME") != 0)
   {
      m_coaddMaxImno = fh["COADTIME"].value<realT>();
      std::cerr << "coaddMaxtime: " << m_coaddMaxTime << "\n";
   }
   
   if(m_maskFile == "")
   {
      if(fh.count("MASKFILE") == 0)
      {
         mxError("KLIPReduction", MXE_PARAMNOTSET, "MASKFILE not found in FITS header and not set in configuration.");
         return -1;
      }
      m_maskFile = fh["MASKFILE"].String();
   }   
   std::cerr << "maskFile: " << m_maskFile << "\n";

   if(fh.count("PPBEFORE") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "PPBEFORE not found in FITS header.");
      return -1;
   }
   m_preProcess_beforeCoadd = fh["PPBEFORE"].value<int>();
   std::cerr << "preProcess_beforeCoadd: " << m_preProcess_beforeCoadd << "\n";
   
   if(fh.count("PPMASK") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "PPMASK not found in FITS header.");
      return -1;
   }
   m_preProcess_mask = fh["PPMASK"].value<int>();
   std::cerr << "preProcess_mask: " << m_preProcess_mask << "\n";
   
   if(fh.count("PPSUBRAD") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "PPSUBRAD not found in FITS header.");
      return -1;
   }
   m_preProcess_subradprof = fh["PPSUBRAD"].value<int>();
   std::cerr << "preProcess_subradprof: " << m_preProcess_subradprof << "\n";
   
   if(fh.count("PPAUSMAW") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "PPAUSMAW not found in FITS header.");
      return -1;
   }
   m_preProcess_azUSM_azW = fh["PPAUSMAW"].value<realT>();
   std::cerr << "preProcess_azUSM_azW: " << m_preProcess_azUSM_azW << "\n";
   
   if(fh.count("PPAUSMRW") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "PPAUSMRW not found in FITS header.");
      return -1;
   }
   m_preProcess_azUSM_radW = fh["PPAUSMRW"].value<realT>();
   std::cerr << "preProcess_azUSM_radW: " << m_preProcess_azUSM_radW << "\n";
   
   if(fh.count("PPGUSMFW") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "PPGUSMFW not found in FITS header.");
      return -1;
   }
   m_preProcess_gaussUSM_fwhm = fh["PPGUSMFW"].value<realT>();
   std::cerr << "preProcess_gaussUSM_fwhm: " << m_preProcess_gaussUSM_fwhm << "\n";
   

   fits::fitsHeader head;

   if(m_MJDKeyword != "") head.append(m_MJDKeyword);

   for(size_t i=0;i<m_keywords.size();++i)
   {
      head.append(m_keywords[i]);
   }
   
   for(size_t n =0; n<nReductions; ++n)
   {
      char nstr[5];
      int nwr = snprintf(nstr, sizeof(nstr), "%03zu", n);
      if(nwr < 0 || n >= sizeof(nstr))
      {
         std::cerr << "possibly bad formatting in filename\n";
      }
         
      
      std::string nprefix = prefix + "_" + nstr + "_";
      load_fileList(dir, nprefix, ext);
      
      if(m_fileList.size() == 0)
      {
         mxError("HCIobservation", MXE_FILENOTFOUND, "The m_fileList has 0 length, there are no files to be read.");
         return -1;
      }

      Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

      fits::fitsFile<realT> f(m_fileList[0]);

      fits::fitsHeader fh = head;
      f.read(im, fh);

      if(n == 0)
      {
         std::cerr << fh["FAKEPA"].String() << "\n";
      }
      
      //We set imSize to match the first image, but we make it a square.
      if(m_imSize == 0)
      {
         m_imSize = im.rows();
         if(m_imSize > im.cols()) m_imSize = im.cols();
      }
      else
      {
         //Now make sure we don't read too much.
         if(m_imSize > im.rows()) m_imSize = im.rows();
         if(m_imSize > im.cols()) m_imSize = im.cols();
      }
      
      //the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
      f.setReadSize( floor(0.5*(im.rows()-1) - 0.5*(m_imSize-1) +0.1), floor(0.5*(im.cols()-1.0) - 0.5*(m_imSize-1.0)+0.1), m_imSize, m_imSize);
      im.resize(m_imSize, m_imSize);

      if(n > 0)
      {
         if(m_fileList.size() != (size_t) m_Nims)
         {
            mxError("HCIobservation", MXE_INVALIDARG, "Different number of images in reductions.");
            return -1;
         }
         if(m_Nrows != im.rows())
         {
            mxError("HCIobservation", MXE_INVALIDARG, "Different number of rows in reductions.");
            return -1;
         }
         if(m_Ncols != im.cols())
         {
            mxError("HCIobservation", MXE_INVALIDARG, "Different number of cols in reductions.");
            return -1;
         }
      }
      else
      {
         std::cerr << "found " << nReductions << " sets of " << m_fileList.size() << " " << im.rows() << " x " << im.cols() << " files\n";
      }    
      m_Nims =  m_fileList.size();
      m_Nrows = im.rows();
      m_Ncols = im.cols();
      m_Npix =  im.rows()*im.cols();
   
      m_psfsub[n].resize(m_Nrows, m_Ncols, m_Nims);
      
      m_heads.clear(); //This is necessary to make sure heads.resize() copies head on a 2nd call
      m_heads.resize(m_fileList.size(), head);

      t_load_begin = sys::get_curr_time();
      
      f.read(m_psfsub[n].data(), m_heads, m_fileList);

      f.setReadSize(); 

      if(m_MJDKeyword != "")
      {
         m_imageMJD.resize(m_heads.size());

         if(m_MJDisISO8601)
         {
            for(size_t i=0;i<m_imageMJD.size();++i)
            {
               m_imageMJD[i] =  sys::ISO8601date2mjd(m_heads[i][m_MJDKeyword].String());
            }
         }
         else
         {
            for(size_t i=0;i<m_imageMJD.size();++i)
            {
               m_imageMJD[i] = m_heads[i][m_MJDKeyword].template value<realT>()*m_MJDUnits;
            }
         }
      }

      t_load_end = sys::get_curr_time();
      
      for(size_t n=0; n<m_psfsub.size(); ++n)
      {
         zeroNaNCube(m_psfsub[n]);
      }

   }
   
   if(m_weightFile != "")
   {
      std::vector<std::string> fn;
      ioutils::readColumns(m_weightFile, fn, m_comboWeights);

      std::cerr << "read: " << m_weightFile << " (" << m_comboWeights.size() << ")\n";
   }
         
   /*** Now do the post-read actions ***/
   if( postReadFiles() < 0) return -1;
  
   readMask();

   m_filesRead = true;
   
   return 0;
}


///@}

extern template class  HCIobservation<float>;
extern template class  HCIobservation<double>;

} //namespace improc
} //namespace mx

#endif //__HCIobservation_hpp__
