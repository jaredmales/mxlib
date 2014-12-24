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

/// The basic high contrast imaging data type
/** 
  */
template<typename _floatT>
struct HCIobservation
{

   typedef _floatT floatT;
   
   std::string dir;
   std::string prefix;
   std::string ext;
   
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
   
   bool doFinimCombine;
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
   
   void initialize()
   {
      imSize = 0;
      applyMask = false;
      
      filesRead = false;
      
      doFinimCombine = true;
      finimName = "finim.fits";
      
      doWeightedCombo = false;
      
      doOutputPsfsub = false;
   }
   
   HCIobservation()
   {
      initialize();
   }
   
   HCIobservation(std::string odir, std::string oprefix, std::string oext)
   {
      initialize();
      
      dir = odir;
      prefix = oprefix;
      ext = oext;      
   }
   
   ///Read the list of files, cut to size, and apply the mask.
   void readFiles();
   
   ///Read the image weights from a single column text file.
   void readWeights();

   ///Combine the images into a single final image.
   /** Images can be combined by either average, weighted-average, or median.
     */
   void combineFinim();
   
};



template<typename _floatT>
inline void HCIobservation<_floatT>::readFiles()
{      
   eigenImagef im;

   vector<string> flist = getFileNames(dir, prefix, ext);
   sort(flist.begin(), flist.end());

   fitsFilef f(flist[0]);

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
   eigenImagef tfinim;
   
   finim.resize(psfsub[0].rows(), psfsub[0].cols(), psfsub.size());
   
   for(int n= 0; n < psfsub.size(); ++n)
   {
      if(doWeightedCombo)
      {
         floatT wsum = 0.0;
         for(int i=0;i<psfsub[n].planes();++i)
         {
            psfsub[n].image(i) = comboWeights[i]*psfsub[n].image(i);
            wsum += comboWeights[i];
         }
      
         psfsub[n].mean(tfinim);
         finim.image(n) = tfinim;// /wsum*psfsub[n].planes();
      }
      else
      {
         psfsub[n].median(tfinim);
         finim.image(n) = tfinim;
      }
   }
}
   
///@} 

} //namespace mx

#endif //__HCIobservation_hpp__
