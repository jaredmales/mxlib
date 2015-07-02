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

#include "templateBLAS.hpp"
#include "fileUtils.hpp"
#include "eigenImage.hpp"
#include "eigenCube.hpp"
#include "imageTransforms.hpp"
#include "fitsFile.hpp"
#include "timeUtils.hpp"
#include "pout.hpp"
#include "imageMasks.hpp"

namespace mx
{
   
/** \addtogroup hc_imaging
  * @{
  */

namespace HCI 
{
   ///Possible combination methods
   enum combineMethods{ noCombine, medianCombine, meanCombine, weightedMeanCombine};
   
   
}

/// The basic high contrast imaging data type
/** 
  */
template<typename _floatT>
struct HCIobservation
{

   typedef _floatT floatT;
   
   std::string dir; ///<The directory where to search for files
   std::string prefix; ///<The prefix of the files
   std::string ext; ///<The extension of the files, default is ".fits"
   
   ///Set the image size.  Images are cut down to this size after reading.
   /** Set to <= 0 to use images uncut.
     */
   int imSize;
   
   ///Controls whether the mask is applied applied.
   /** If true, then the mask described by \ref maskIdx and \ref maskVal is applied after the file read.
     */
   bool applyMask;
   
   ///Indices of the mask to apply to each image if \ref applyMask is true
   std::vector<size_t> maskIdx;
   
   ///Value to insert as the mask according to \ref maskIdx.
   floatT maskVal;
   
   ///Determine how to combine the images
   int combineMethod;
   
   
   ///Set whether the final combined image is written to disk
   int doWriteFinim;
   
   ///The base file name of the output final image
   /** The complete name is formed by combining with a sequential number
     */
   std::string finimName;
   
   
   
   bool doWeightedCombo;
   std::string weightFile;
   vector<floatT> comboWeights;
   
   bool doOutputPsfsub;
   std::string psfsubPrefix;
      
   std::vector<std::string> keywords;
   vector<fitsHeader> heads;
   
   bool filesRead;
   
   eigenCube<floatT> imc;
   std::vector<eigenCube<floatT> > psfsub;
   
   eigenCube<floatT> finim;
   
   int Nims;
   int Nrows;
   int Ncols;
   int Npix;
   
   void initialize();

   HCIobservation();
   
   HCIobservation(std::string oprefix);
   HCIobservation(std::string odir, std::string oprefix);
   HCIobservation(std::string odir, std::string oprefix, std::string oext);

   
   void readFiles(const std::vector<std::string> &flist);
   
   ///Read the list of files, cut to size, and apply the mask.
   void readFiles();
   
   ///Read the image weights from a single column text file.
   void readWeights();

   ///Combine the images into a single final image.
   /** Images can be combined by either average, weighted-average, or median.
     */
   void combineFinim();
   
   ///Write the final combined image to disk
   /** 
    */
   void writeFinim(fitsHeader * addHead = 0);
};

template<typename _floatT>
void HCIobservation<_floatT>::initialize()
{
   dir = "./";
   prefix = "";
   ext = ".fits";
   
   imSize = 0;
   applyMask = false;
   
   filesRead = false;
   
   combineMethod = HCI::medianCombine;
   
   doWriteFinim = 1;
   finimName = "finim_";
   
   doWeightedCombo = false;
   
   doOutputPsfsub = false;
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation()
{
   initialize();
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation(std::string oprefix)
{
   initialize();
   prefix = oprefix;
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation(std::string odir, std::string oprefix)
{
   initialize();
   
   dir = odir;
   prefix = oprefix;
}

template<typename _floatT>
HCIobservation<_floatT>::HCIobservation(std::string odir, std::string oprefix, std::string oext)
{
   initialize();
   
   dir = odir;
   prefix = oprefix;
   ext = oext;      
}
   

template<typename _floatT>
inline void HCIobservation<_floatT>::readFiles()
{
   vector<string> flist = getFileNames(dir, prefix, ext);
      
   sort(flist.begin(), flist.end());

   readFiles(flist);
}

template<typename _floatT>
inline void HCIobservation<_floatT>::readFiles(const std::vector<std::string> & flist)
{      
   Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> im;
      
   fitsFile<floatT> f(flist[0]);

   f.read(im);

   fitsHeader head;

   for(int i=0;i<keywords.size();++i)
   {
      head.append(keywords[i], 0);
   }
   
   heads.resize(flist.size(), head);
   head.clear();
   
   imc.resize(im.rows(), im.cols(), flist.size());
   
   f.read(flist, imc.data(), heads);

   
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
   
      std::cout << min_x << " " << max_x << " " << min_y << " " << max_y << "\n";
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
   
   
   if(applyMask)
   {
      for(int n=0;n<Nims;++n)
      {
         typename eigenCube<floatT>::imageRef im = imc.image(n);
         mx::applyMask(im, maskIdx, maskVal);
      }
   }  
   
   if(weightFile != "")
   {
      readWeights();
   }
   
   filesRead = true;
}
 
template<typename _floatT>
inline void HCIobservation<_floatT>::readWeights()
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

   doWeightedCombo = true;
}



template<typename _floatT>
inline void HCIobservation<_floatT>::combineFinim()
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
   
   
   //Create and size temporary image for averaging
   eigenImagef tfinim;
   
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
   head.append<int>("IMSIZE", imSize, "image size after reading");
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

///@} 

} //namespace mx

#endif //__HCIobservation_hpp__
