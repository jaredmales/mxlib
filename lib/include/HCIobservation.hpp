
#include <vector>
#include <string>

#include "templateBLAS.hpp"
#include "fileUtils.hpp"
#include "eigenImage.hpp"
#include "eigenCube.hpp"
#include "imageTransforms.hpp"
#include "fitsFile.hpp"
#include "timeUtils.hpp"
#include "pout.hpp"


using namespace mx;

#ifndef __HCIobservation_hpp__
#define __HCIobservation_hpp__



template<typename _floatT>
struct HCIobservation
{
   typedef _floatT floatT;
   
   std::string dir;
   std::string prefix;
   std::string ext;
   
   
   std::vector<std::string> keywords;
   vector<fitsHeader> heads;
   
   bool filesRead;
   
   eigenCube<floatT> imc;
   eigenCube<floatT> psfsub;
   
   Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> finim;
   
   int Nims;
   int Nrows;
   int Ncols;
   int Npix;
   
   HCIobservation()
   {
      filesRead = false;
   }
   
   HCIobservation(std::string odir, std::string oprefix, std::string oext)
   {
      dir = odir;
      prefix = oprefix;
      ext = oext;
      
      filesRead = false;
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
         head.append(keywords[i], "");
      }
      
      heads.resize(flist.size(), head);

      
      imc.resize(im.rows(), im.cols(), flist.size());
      
      f.read(flist, imc.data(), heads);
   
      Nims =  imc.planes();
      Nrows = imc.rows();
      Ncols = imc.cols();
      Npix =  imc.rows()*imc.cols();
      
      
      filesRead = true;
   }
    
   void combine()
   {
      imc.median(finim);
   }
   
};



#endif //__HCIobservation_hpp__
