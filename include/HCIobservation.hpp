/** \file HCIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the basic high contrast imaging data type.
  * \ingroup hc_imaging
  *
  */

#ifndef __HCIobservation_hpp__
#define __HCIobservation_hpp__ 

#include <vector>
#include <string>
#include <fstream>

#include "mxlib.h"

#include "mxException.hpp"

#include "templateBLAS.hpp"
#include "fileUtils.hpp"
#include "eigenImage.hpp"
#include "eigenCube.hpp"
#include "imageTransforms.hpp"
#include "fitsFile.hpp"
#include "timeUtils.hpp"
#include "pout.hpp"
#include "imageMasks.hpp"
#include "mxException.hpp"

namespace mx
{
   
/** \addtogroup hc_imaging
  * @{
  */

namespace HCI 
{
   ///Possible combination methods
   enum combineMethods{ noCombine, medianCombine, meanCombine, weightedMeanCombine, debug};
}

/// The basic high contrast imaging data type
/** 
  */
template<typename _floatT>
struct HCIobservation
{

   ///The arithmetic type used for calculations.  Does not have to match the type in images on disk.
   typedef _floatT floatT;
   
   ///The Eigen image array type basted on floatT
   typedef Array<floatT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
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
   int deleteFront;
   
   ///Specify how many files from fileList to delete from the back of the list
   int deleteBack;
   
   ///Name of the keyword to use for the image date. 
   /** Specifies the keyword corresponding to .  This is
     * the "DATE" keyword for file write time, and usually "DATE-OBS" for actual observation time.
     *  
     * Default is "DATE-OBS".  
     * 
     * If empty "", then image date is not read.
     */
   std::string MJDKeyword;
   
   ///Whether or not the date is in ISO 8601 format
   bool MJDisISO8601;
   
   ///If the date is not ISO 8601, this specifies the conversion to Julian Days (i.e. seconds to days)
   floatT MJDUnits;
   
   ///Vector of FITS header keywords to read from the files in fileList.
   std::vector<std::string> keywords;

   ///Set the image size.  Images are cut down to this size after reading.
   /** Set to <= 0 to use images uncut.
     */
   int imSize;
   
   ///Controls whether the mask is applied.
   /** If true, then the mask described by \ref maskIdx and \ref maskVal is applied after the file read.
     */
   bool applyMask;
   
   ///Indices of the mask to apply to each image if \ref applyMask is true
   std::vector<size_t> maskIdx;
   
   ///Value to insert as the mask according to \ref maskIdx.
   floatT maskVal;
   
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
   int coaddCombineMethod;
  
   ///Maximum number of images to coadd
   int coaddMaxImno;
   
   ///Maximum elapsed time over which to coadd the images.
   int coaddMaxTime;
  
   ///The values of these keywords will be averaged and replaced.
   std::vector<std::string> coaddKeywords;
   
   ///@} -- coadding
   
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
     */
   int combineMethod;
   
   ///Specifies a file containing the image weights, for combining with weighted mean.
   /** This simple ASCII file must contain one weight per image, and must be specified before readFiles()
     * is executed.  Weights should be white-space separated.
     */
   std::string weightFile;
   
   ///Vector to hold the weights read from the weightFile.
   vector<floatT> comboWeights;
   
   
   ///@}
   
   /** \name Output
     * These options control the ouput fo the final combined images and the individual PSF subtracted images.
     * @{
     */
   
   ///Set whether the final combined image is written to disk
   int doWriteFinim;
   
   ///The base file name of the output final image
   /** The complete name is formed by combining with a sequential number and the fits extension.
     * that is: finimName0000.fits 
     */
   std::string finimName;
   
   
   ///Controls whether or not the individual PSF subtracted images are written to disk.  
   /** - true -- write to disk
     * - false -- [default] don't write to disk 
     */
   bool doOutputPSFSub;
   
   ///Prefix of the FITS file names used to write individual PSF subtracted images to disk if doOutputPSFSub is true.
   std::string PSFSubPrefix;
         
   ///@}
   
   ///\name The Raw Data
   /** @{
    */
   
   ///Vector image times, in MJD.
   std::vector<double> imageMJD;
   
   ///Vector of FITS headers,one per file, populated with the values for the keywords.
   vector<fitsHeader> heads;
   
   ///Whether or not the fileList has been read.
   bool filesRead;
   
   ///Whether or not the specified files have been deleted from fileList
   bool filesDeleted;
   
   ///The image cube
   eigenCube<floatT> imc;

   int Nims;  ///<Number of images in imc
   int Nrows; ///<Number of rows of the images in imc
   int Ncols; ///<Number of columns of the images in imc
   int Npix; ///<Pixels per image, that is Nrows*Ncols

   ///@}
   
   ///\name The Reduced Data
   /** @{
     */
   ///The PSF subtracted images
   /** This is a vector of cubes so that it can contain results from different reductions,
     * e.g. different modes when using KLIP.
     */ 
   std::vector<eigenCube<floatT> > psfsub;
   
   ///The final combined images, one for each cube in psfsub.
   eigenCube<floatT> finim;
   
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
     * "./prefix*".  See \ref loadFileList
     * 
     * \param [in] prefix is the initial part of the file name.  Can be empty "".
     */
   HCIobservation(const std::string &prefix);
   
   ///Construct and load the file list.
   /** Populates the \ref fileList vector by searching on disk for files which match
     * "dir/prefix*".  See \ref loadFileList
     * 
     * \param [in] dir is the directory to search.
     * \param [in] prefix is the initial part of the file name.  Can be empty "".
     */
   HCIobservation(const std::string &dir, const std::string &prefix);
   
   ///Construct and load the file list.
   /** Populates the \ref fileList vector by searching on disk for files which match
     * "dir/prefix*.ext".  See \ref loadFileList
     * 
     * \param [in] dir is the directory to search.
     * \param [in] prefix is the initial part of the file name.  Can be empty "".
     * \param [in] ext is the extension to append to the file name. Can be empty "".
     */
   HCIobservation(const std::string &dir, const std::string &prefix, const std::string &ext);

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
   void readFiles();
      
   ///Perform post-read actions, for use by derived classes
   virtual void postReadFiles();

   ///Read the image weights from a single column text file.
   void readWeights();

   ///Coadd the images
   void coaddImages();
   
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
   
};

template<typename _floatT>
void HCIobservation<_floatT>::initialize()
{   
   deleteFront = 0;
   deleteBack = 0;
   filesDeleted = false;
   
   MJDKeyword = "DATE-OBS";
   MJDisISO8601 = true;   
   MJDUnits = 1.0;
   
   imSize = 0;
   applyMask = false;
  
   coaddCombineMethod = HCI::noCombine;
   coaddMaxImno = 0;
   coaddMaxTime = 0;
 
   filesRead = false;
   
   combineMethod = HCI::medianCombine;
   
   doWriteFinim = 1;
   finimName = "finim_";
      
   doOutputPSFSub = false;
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation()
{
   initialize();
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation(const std::string & prefix)
{
   initialize();
   loadFileList("./", prefix, ".fits"); 
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation(const std::string & dir, const std::string & prefix)
{
   initialize();   
   loadFileList(dir, prefix, ".fits");
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation(const std::string & dir, const std::string & prefix, const std::string & ext)
{
   initialize();
   loadFileList(dir, prefix, ext);
}
   

template<typename _floatT>
inline void HCIobservation<_floatT>::loadFileList(const std::string & dir, const std::string & prefix, const std::string & ext)
{
   fileList = getFileNames(dir, prefix, ext);
   filesDeleted = false;
}

template<typename _floatT>
inline void HCIobservation<_floatT>::readFiles()
{
   
   if(fileList.size() == 0)
   {
      throw mxException("mx::HCIobservation",1, "no files to read",__FILE__, __LINE__-2, 
                      "The fileList has 0 length, there are no files to be read.");
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

   Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> im;
      
   fitsFile<floatT> f(fileList[0]);

   f.read(im);

        
   fitsHeader head;

   if(MJDKeyword != "") head.append(MJDKeyword);
   
   for(int i=0;i<keywords.size();++i)
   {
      head.append(keywords[i]);
   }
      
   heads.clear(); //This is necessary to make sure heads.resize copies head on a 2nd call
   heads.resize(fileList.size(), head);
      
   imc.resize(im.rows(), im.cols(), fileList.size());
   
   f.read(fileList, imc.data(), heads);

   if(MJDKeyword != "")
   {
      imageMJD.resize(heads.size());

      if(MJDisISO8601)
      {
         for(int i=0;i<imageMJD.size();++i)
         {
            imageMJD[i] =  mx::ISO8601date2mjd(heads[i][MJDKeyword].String());
         }
      }
      else
      {
         for(int i=0;i<imageMJD.size();++i)
         {
            imageMJD[i] =  heads[i][MJDKeyword].Value<floatT>()*MJDUnits;
         }
      }
   }

   //Re-size the image
   if(imSize > 0)
   {
      eigenCube<floatT> timc;
   
      timc.shallowCopy(imc, true);
     
      double xc = 0.5*(timc.rows()-1);
      double yc = 0.5*(timc.cols()-1);

      int min_x = floor(xc - (0.5*imSize-0.5) + 0.01);
      if(min_x < 0) min_x = 0;
   
      int max_x = floor(xc + (0.5*imSize-0.5) + 0.51);
      if(max_x >= timc.cols()) max_x = timc.rows()-1;
   
      int min_y = floor(yc - (0.5*imSize-0.5) + 0.01);
      if(min_y < 0) min_y = 0;
   
      int max_y = floor(yc + (0.5*imSize-0.5) + 0.51);
      if(max_y >= timc.rows()) max_y = timc.cols() - 1;
   
      imc.resize( max_y-min_y + 1, max_x-min_x+1, timc.planes());
   
      for(int n=0;n<timc.planes();++n)
      {
         imc.image(n) = timc.image(n).block(min_x, min_y, max_x-min_x + 1, max_y-min_y+1);
      }

   }
   
   Nims =  imc.planes();
   Nrows = imc.rows();
   Ncols = imc.cols();
   Npix =  imc.rows()*imc.cols();
   
   
   /*** Now do the post-read actions ***/
   postReadFiles();
   
   if(weightFile != "")
   {
      readWeights();
   }
   
   if(coaddCombineMethod != HCI::noCombine)
   {
      coaddImages();
   }
   
   if(applyMask)
   {
      for(int n=0;n<Nims;++n)
      {
         typename eigenCube<floatT>::imageRef im = imc.image(n);
         mx::applyMask(im, maskIdx, maskVal);
      }
   }  
   
   eigenImageT kernel;
   gaussKernel(kernel, 15);
   
   for(int n=0;n<Nims;++n)
   {
      eigenImageT imbg;
      typename eigenCube<floatT>::imageRef im = imc.image(n);
      
      smoothImage(imbg, im, kernel);
      im = im - imbg;
      mx::applyMask(im, maskIdx, maskVal);
   }
   
   
   
   filesRead = true;
}
 
template<typename _floatT>
void HCIobservation<_floatT>::postReadFiles()
{
   return;
}

template<typename _floatT>
void HCIobservation<_floatT>::readWeights()
{
   std::ifstream fin;
   std::string str;
   
   fin.open(weightFile.c_str());
   
   comboWeights.resize(Nims);
   
   for(int i=0; i<Nims; ++i)
   {
      fin >> str;
      comboWeights[i] = convertFromString<floatT>(str);
   }
   
   fin.close();
}

template<typename _floatT>
void HCIobservation<_floatT>::coaddImages()
{
   //Validate setup
   if(coaddMaxImno <=0 && coaddMaxTime <= 0) return;
   //Validate combine method
   if(coaddCombineMethod == HCI::noCombine) return;

   pout("coadding raw images\n");
   
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
   eigenCube<floatT> imsToCoadd;  
   
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

      for(int i=0;i<coaddKeywords.size(); ++i)
      {
         initVals[i] = heads[im0][coaddKeywords[i]].Value<double>();
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
         for(int i=0;i<coaddKeywords.size(); ++i)
         {
            initVals[i] += heads[imno][coaddKeywords[i]].Value<double>();
         }
      }
      
      //And then turn them into an average
      initMJD /= (imF - im0);
      for(int i=0;i<coaddKeywords.size(); ++i)
      {
         initVals[i] /= (imF-im0);
      }

      //Extract the images into the temporary
      imsToCoadd.resize(Nrows, Ncols, imF-im0);
      for(int i =0; i < (imF-im0); ++i)
      {
         imsToCoadd.image(i) = imc.image(im0 + i);
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
   for(int i=0;i<Nims;++i)
   {
      imageMJD[i] = avgMJD[i];
      for(int j=0;j<coaddKeywords.size(); ++j)
      {
         heads[i][coaddKeywords[j]].setValue(avgVals[i][j]);
      }
   }
   
   
}//void HCIobservation<_floatT>::coaddImages()

template<typename _floatT>
void HCIobservation<_floatT>::combineFinim()
{
   if(combineMethod == HCI::noCombine) return;
   
   
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
   else if(combineMethod == HCI::weightedMeanCombine && comboWeights.size() == Nims)
   {
      method = HCI::weightedMeanCombine;
   }
   else if(combineMethod == HCI::weightedMeanCombine && comboWeights.size() != Nims)
   {
      method = HCI::meanCombine;
   }
   else if(combineMethod == HCI::debug)
   {
      method = HCI::debug;
   }
   
   
   //Create and size temporary image for averaging
   eigenImageT tfinim;
   
   finim.resize(psfsub[0].rows(), psfsub[0].cols(), psfsub.size());
   
   //Now cycle through each set of psf subtractions
   for(int n= 0; n < psfsub.size(); ++n)
   {
      if(method == HCI::weightedMeanCombine)
      {
         for(int i=0;i<psfsub[n].planes();++i)
         {
            psfsub[n].image(i) = comboWeights[i]*psfsub[n].image(i);
         }
         psfsub[n].mean(tfinim);
         finim.image(n) = tfinim; // /wsum*psfsub[n].planes();
      }
      else if(method == HCI::medianCombine)
      {
         psfsub[n].median(tfinim);
         finim.image(n) = tfinim;
      }
      else if(method == HCI::meanCombine)
      {
         psfsub[n].mean(tfinim);
         finim.image(n) = tfinim;
      }
      else if(method == HCI::debug)
      {
         finim.image(n) = psfsub[n].image(0);
      }
   }
}
   
template<typename _floatT>
inline void HCIobservation<_floatT>::writeFinim(fitsHeader * addHead)
{
   std::string fname;
   
   
   fname = getSequentialFilename(finimName, ".fits");
   
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

   
   head.append<int>("COMBMTHD", combineMethod, "combination method");
                    
   
   if(addHead)
   {
      head.append(*addHead);
   }
   
   fitsHeaderGitStatus(head, "mxlib_comp",  mxlib_compiled_git_sha1(), mxlib_compiled_git_repo_modified());
   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   
   fitsFile<floatT> f;
      
   f.write(fname, finim.data(), finim.rows(), finim.cols(), finim.planes(), &head);
   
   pout("Final image written to: ", fname);
}


template<typename _floatT>
inline void HCIobservation<_floatT>::outputPSFSub(fitsHeader * addHead)
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
   
   fitsHeaderGitStatus(head, "mxlib_comp",  mxlib_compiled_git_sha1(), mxlib_compiled_git_repo_modified());
   fitsHeaderGitStatus(head, "mxlib_uncomp",  MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED);
   
   fitsFile<floatT> f;

   char num[16];
   for(int n=0; n<psfsub.size(); ++n)
   {
      snprintf(num, 16, "%05d.fits", n);
      fname = PSFSubPrefix + num;
   
      f.write(fname, psfsub[n].data(), psfsub[n].rows(), psfsub[n].cols(), psfsub[n].planes(), &head);
   }
   
}

///@} 

} //namespace mx

#endif //__HCIobservation_hpp__
