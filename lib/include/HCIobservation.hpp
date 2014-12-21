/** \file HCIobservation.hpp
  * \author Jared R. Males
  * \brief Defines the basic high contrast imaging data type.
  * \ingroup hc_imaging
  *
  */

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

#ifndef __HCIobservation_hpp__
#define __HCIobservation_hpp__

namespace mx
{
   
/** \addtogroup hc_imaging
  * @{
  */

///The basic high contrast imaging data type
/** Why you no doc?
 */
template<typename _floatT>
struct HCIobservation
{
   
   typedef _floatT floatT;
   
   std::string dir;
   std::string prefix;
   std::string ext;
   
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
   
   void readFiles()
   {      
      eigenImagef im;
   
      vector<string> flist = getFileNames(dir, prefix, ext);
      sort(flist.begin(), flist.end());

      //flist.erase(flist.begin(), flist.end()-100);
   
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
   
      Nims =  imc.planes();
      Nrows = imc.rows();
      Ncols = imc.cols();
      Npix =  imc.rows()*imc.cols();
      
      
      if(weightFile != "")
      {
         readWeights();
      }
      
      filesRead = true;
   }
    
   void readWeights()
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
   
   void combineFinim()
   {
      eigenImagef tfinim;
      
      finim.resize(psfsub[0].rows(), psfsub[0].cols(), psfsub.size());
      
      for(int n= 0; n < psfsub.size(); ++n)
      {
         if(doWeightedCombo)
         {
            for(int i=0;i<psfsub[n].planes();++i)
            {
               psfsub[n].image(i) = comboWeights[i]*psfsub[n].image(i);
            }
         
            psfsub[n].mean(tfinim);
            finim.image(n) = tfinim;
         }
         else
         {
            psfsub[n].median(tfinim);
            finim.image(n) = tfinim;
         }
      }
   }
   
//    void outputPsfsub()
//    {
//       fitsFile<floatT> ff;
//       std::string fname;
//       
//       for(int i=0;i<psfsub.planes(); ++i)
//       {
//          fname = psfsubPrefix + convertToString(i) + ".fits";
//          ff.write(fname, psfsub.image(i).data(), psfsub.rows(), psfsub.cols());
//       }
//    }
   
};


/// @}

} //namespace mx

#endif //__HCIobservation_hpp__
