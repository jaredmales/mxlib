/** \file KLIPreduction.hpp
  * \author Jared R. Males
  * \brief Declarations and definitions for an implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm. 
  * \ingroup hc_imaging
  *
  */



#ifndef __KLIPreduction_hpp__
#define __KLIPreduction_hpp__

#include "ADIobservation.hpp"

namespace mx
{
   
double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, tf;
double dcut, dcv, dklims, dgemm, dsyevr, dcfs, drot, dcombo, dread;

/** \addtogroup hc_imaging
  * @{
  */

namespace HCI 
{
   enum excludeMethods{excludeNone, excludePixel, excludeAngle, excludeImno};   
}

/// An implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm.
/**
  * 
  */ 
template<typename _floatT, class _derotFunctObj>
struct KLIPreduction : public ADIobservation<_floatT, _derotFunctObj>
{
   typedef _floatT floatT;
   
   typedef Array<floatT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   int padSize;
   
   /// Specifies the number of modes to include in the PSF.
   /** The output image is a cube with a plane for each entry in Nmodes.
     * Only the number of eigenvalues required for the maximum value of Nmodes
     * are calculated, so this can have an important impact on speed.
     * 
     * Can be initialized as:
     * \code 
     * red.Nmodes={5,10,15,20};
     * \endcode
     * 
     */ 
   std::vector<int> Nmodes;
   
   int maxNmodes;
   
   floatT mindpx;
   
   
   /// Controls how reference images are excluded, if at all, from the covariance matrix for each target image.
   /** Can have the following values:
     *  - <b>HCI::excludeNone</b> = no exclusion, all images included
     *  - <b>HCI::excludePixel</b> = exclude based on pixels of rotation at the inner edge of the region
     *  - <b>HCI::excludeAngle</b> = exclude based on degrees of rotation at the inner edge of the region
     *  - <b>HCI::excludeImno</b> = exclude based on number of images
     */
   int excludeMethod;
   
   ///Number of reference images to include in the covariance matrix
   /** If > 0, then at most this many images, determined by highest cross-correlation, are included.
     * This is determined after rotational/image-number exclusion. 
     * If == 0, then all reference images are included. 
     */
   int includeRefNum; 
   
   
   
   KLIPreduction()
   {
      initialize();
   }
   
   KLIPreduction( const std::string & odir, 
                  const std::string & oprefix, 
                  const std::string & oext) : ADIobservation<_floatT, _derotFunctObj>(odir, oprefix, oext)
   {
      initialize();
   }
   
   void initialize()
   {
      excludeMethod = HCI::excludeNone;
      includeRefNum = 0;
      padSize = 4;
   }
   
   void meanSubtract(eigenCube<floatT> & ims, std::vector<floatT> & sds);
   void medianSubtract(eigenCube<floatT> & ims, std::vector<floatT> & sds);
   void getStdDevs(std::vector<floatT> sd, eigenCube<floatT> & ims);
   
   void regions( vector<floatT> minr, 
                 vector<floatT> maxr, 
                 vector<floatT> minq, 
                 vector<floatT> maxq);
   
   void regions( floatT minr, 
                 floatT maxr, 
                 floatT minq, 
                 floatT maxq)
   {
      std::vector<floatT> vminr(1, minr);
      std::vector<floatT> vmaxr(1, maxr);
      std::vector<floatT> vminq(1, minq);
      std::vector<floatT> vmaxq(1, maxq);
      
      regions(vminr, vmaxr, vminq, vmaxq);
   }
   
   void worker(eigenCube<floatT> & rims, vector<size_t> & idx, floatT dang);
   
   template<typename eigenT, typename eigenT1>
   void calcKLIms( eigenT & klims, 
                   eigenT & cv, 
                   const eigenT1 & Rims, 
                   int n_modes = 0 );   
   
   template<typename eigenT, typename eigenTv>
   void collapseCovar( eigenT & cutCV, 
                       const eigenT & CV,
                       const std::vector<floatT> & sds,
                       eigenT & rimsCut,
                       const eigenTv & rims,
                       int imno,
                       double dang );
                       
   
};

template<typename _floatT, class _derotFunctObj>
inline
void KLIPreduction<_floatT, _derotFunctObj>::meanSubtract(eigenCube<floatT> & ims, std::vector<_floatT> & norms)
{

   norms.resize(ims.planes());

   if(this->applyMask)
   {
      //S#pragma omp parallel for schedule(static, 1)
      for(int n=0;n<ims.planes(); ++n)
      {
         _floatT mn = 0, Navg = 0;
      
         for(int j=0;j<ims.rows()*ims.cols();++j)
         {
            if( ims.image(n)(j) != this->maskVal)
            {
               mn += ims.image(n)(j);
               ++Navg;
            }
         }
         mn /= Navg;

         for(int j=0;j<ims.rows()*ims.cols();++j)
         {
            if( ims.image(n)(j) != this->maskVal)
            {      
               ims.image(n)(j) -= mn;//ims.image(i).mean();
            }
            else 
            {
               ims.image(n)(j) = 0.;
            }
         }
         norms[n] = ims.image(n).matrix().norm();
      }
   }
   else
   {
      for(int n=0;n<ims.planes(); ++n)
      {
         ims.image(n) -= ims.image(n).mean();
         norms[n] = ims.image(n).matrix().norm();
      }
   }
}
 
template<typename _floatT, class _derotFunctObj>
inline
void KLIPreduction<_floatT, _derotFunctObj>::medianSubtract(eigenCube<floatT> & ims, std::vector<_floatT> & sds)
{
         
   sds.resize(ims.planes());
   //#pragma omp parallel for schedule(static, 1)
   for(int i=0;i<ims.planes(); ++i)
   {
      _floatT med = eigenMedian(ims.image(i));
      ims.image(i) -= med;
      sds[i] = ims.image(i).matrix().norm();//This is the standard deviation relative to the median.
   }
} 

template<typename floatT, class derotFunctObj>
inline
void KLIPreduction<floatT, derotFunctObj>::regions( vector<floatT> minr, 
                                                    vector<floatT> maxr, 
                                                    vector<floatT> minq, 
                                                    vector<floatT> maxq)
{   
   t0 = get_curr_time();
      
   dklims = 0;
   dgemm = 0;
   dsyevr = 0;
   dcfs=0;
   dcut = 0;
   dcv = 0;
   dread = 0;
   
   maxNmodes = Nmodes[0];
   for(int i = 1; i<Nmodes.size(); ++i)
   {
      if(Nmodes[i] > maxNmodes) maxNmodes = Nmodes[i];
   }
   
   pout("Beginning\n");
      

   t1 = get_curr_time();
   
   if(!this->filesRead) 
   {
      this->imSize = 2*(*std::max_element(maxr.begin(),maxr.end()) + padSize);
            
      this->readFiles();
      
   }
   
   
   
   pout("Files read\n");
   dread = get_curr_time()-t1;
   
   this->psfsub.resize(Nmodes.size());
   for(int n=0;n<Nmodes.size(); ++n)
   {
      this->psfsub[n].resize(this->Nrows, this->Ncols, this->Nims);
      this->psfsub[n].cube().setZero();
   }
   
   //Make radius and angle images
   eigenImagef rIm(this->Nrows,this->Ncols);
   eigenImagef qIm(this->Nrows,this->Ncols);
   
   radAngImage(rIm, qIm, .5*(this->Nrows-1), .5*(this->Ncols-1));

   pout("starting regions", minr.size());
   //******** For each region do this:
   for(int regno = 0; regno < minr.size(); ++regno)
   {
      vector<size_t> idx = imageRegionIndices(rIm, qIm, .5*(this->Nrows-1), .5*(this->Ncols-1), 
                                                    minr[regno], maxr[regno], minq[regno], maxq[regno]);
   
      //Create storage for the R-ims and psf-subbed Ims
      eigenCube<float> rims(idx.size(), 1, this->Nims);
   
      t3 = get_curr_time();
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i< this->Nims; ++i)
      {
         auto rim = rims.image(i);
         cutImageRegion(rim, this->imc.image(i), idx, false);
      }
      dcut += get_curr_time() - t3;

      floatT dang = 0;
      
      if(excludeMethod == HCI::excludePixel)
      {
         dang = fabs(atan(mindpx/minr[regno]));
      }
      else if(excludeMethod == HCI::excludeAngle)
      {
         dang = DTOR(mindpx);
      }
      else if(excludeMethod == HCI::excludeImno)
      {
         dang = mindpx;
      }
      
      //*** Dispatch the work
      worker(rims, idx, dang);
      pout("worker done");
      
   }
   pout("deroting");
   drot = get_curr_time();
   this->derotate();
   drot= get_curr_time()-drot;
   
   if(this->combineMethod > 0)
   {
      pout("combining");
      dcombo = get_curr_time();
      this->combineFinim();
      dcombo = get_curr_time() - dcombo;
      
      pout("writing");
      
      fitsHeader head;
      
      head.append("", fitsCommentType(), "----------------------------------------");
      head.append("", fitsCommentType(), "mx::KLIPreduction parameters:");
      head.append("", fitsCommentType(), "----------------------------------------");
   
      std::stringstream str;
      for(int nm=0;nm < Nmodes.size()-1; ++nm) str << Nmodes[nm] << ",";
      str << Nmodes[Nmodes.size()-1];      
      head.append<char *>("NMODES", (char *)str.str().c_str(), "number of modes");
      
      str.str("");
      for(int nm=0;nm < minr.size()-1; ++nm) str << minr[nm] << ",";
      str << minr[minr.size()-1];      
      head.append<char *>("REGMINR", (char *)str.str().c_str(), "region inner edge(s)");
      
      str.str("");
      for(int nm=0;nm < maxr.size()-1; ++nm) str << maxr[nm] << ",";
      str << maxr[maxr.size()-1];      
      head.append<char *>("REGMAXR", (char *)str.str().c_str(), "region outer edge(s)");
      
      str.str("");
      for(int nm=0;nm < minq.size()-1; ++nm) str << minq[nm] << ",";
      str << minq[minq.size()-1];      
      head.append<char *>("REGMINQ", (char *)str.str().c_str(), "region minimum angle(s)");
      
      str.str("");
      for(int nm=0;nm < maxq.size()-1; ++nm) str << maxq[nm] << ",";
      str << maxq[maxq.size()-1];      
      head.append<char *>("REGMAXQ", (char *)str.str().c_str(), "region maximum angle(s)");
      
      
      head.append<int>("EXCLMTHD", excludeMethod, "value of excludeMethod");
      head.append<floatT>("MINDPX", mindpx, "minimum pixel delta");
      head.append<int>("INCLREFN", includeRefNum, "value of includeRefNum");

      this->writeFinim(&head);
   }
   
//    if(this->doOutputPsfsub)
//    {
//       pout("outputting");
//       this->outputPsfsub();
//    }
   
   tf = get_curr_time();
   pout("total time: ", tf-t0);
   pout("reading: " , dread);
   pout("Non-reading: ", tf-t0 - (dread));
   pout("rim cut: ", dcut);
   
   pout("covar: ", dcv);
   pout("klims: ", dklims);
   pout("  evecs: ", dsyevr);
   pout("  gemm: ", dgemm);
   pout("cfs: ", dcfs);
   pout("derotation: ", drot);
   pout("combination: ", dcombo);
   
}

template<typename floatT, class derotFunctObj>
inline
void KLIPreduction<floatT, derotFunctObj>::worker(eigenCube<floatT> & rims, vector<size_t> & idx, floatT dang)
{
   pout("beginning worker");
   
//   eigenImagef evecs;//, evals;
   eigenImagef cfs; //The coefficients
   eigenImagef psf;
   eigenImagef rims_cut;
   eigenImagef klims;
   eigenImagef cv;
   eigenImagef cv_cut;
   std::vector<floatT> sds;
      
   //*** First mean subtract ***//
   pout("Median subtracting\n");
   
   //meanSubtract(rims);
   meanSubtract(rims, sds);

   //*** Form lower-triangle covariance matrix
   pout("calculating covariance matrix");
   t5 = get_curr_time();   
   eigenSYRK(cv, rims.cube());
   dcv += get_curr_time() - t5;

   if( excludeMethod == HCI::excludeNone )
   {
      pout("calculating K-L images");
      /**** Now calculate the K-L Images ****/
      t7 = get_curr_time();
      calcKLIms(klims, cv, rims.cube(), maxNmodes);
      t9 = get_curr_time();
      dklims += t9 - t7;
   }

   for(int imno = 0; imno < this->Nims; ++imno)
   {
      pout("image:", imno, "/", this->Nims);

       if( excludeMethod != HCI::excludeNone )
       {
         collapseCovar( cv_cut,  cv, sds, rims_cut, rims.asVectors(), imno, dang);
         /**** Now calculate the K-L Images ****/
         t7 = get_curr_time();
         calcKLIms(klims, cv_cut, rims_cut, maxNmodes);
         t9 = get_curr_time();
         dklims += t9 - t7;
   
      } //if(mindpx != 0)
      cfs.resize(1, klims.rows());
   
      t10 = get_curr_time();
   
      #pragma omp parallel for schedule(static, 1)
      for(int j=0; j<cfs.size(); ++j)
      {
         cfs(j) = klims.row(j).matrix().dot(rims.cube().col(imno).matrix());
      }
      dcfs += get_curr_time()-t10;
  
      for(int mode_i =0; mode_i < Nmodes.size(); ++mode_i)
      {
         psf = cfs(maxNmodes-1)*klims.row(maxNmodes-1);

         //Count down, since eigenvalues are returned in increasing order
         for(int j=maxNmodes-2; j>=maxNmodes-Nmodes[mode_i]; --j)
         {
            psf += cfs(j)*klims.row(j);
         }    
         insertImageRegion(this->psfsub[mode_i].cube().col(imno), rims.cube().col(imno) - psf.transpose(), idx);
      }
   }

}


template<typename floatT, class derotFunctObj>
template<typename eigenT, typename eigenT1>
inline
void KLIPreduction<floatT, derotFunctObj>::calcKLIms( eigenT & klims, 
                                                      eigenT & cv, 
                                                      const eigenT1 & Rims, 
                                                      int n_modes )
{

   eigenT evecs, evals;
   
   if(cv.rows() != cv.cols())
   {
      std::cerr << "Non-square covariance matrix input to klip_klims\n";
      return;
   }

   if(cv.rows() != Rims.cols())
   {
      std::cerr << "Covariance matrix - reference image size mismatch in klip_klims\n";
      return;
   }


   int tNims = cv.rows();
   int tNpix = Rims.rows();

   if(n_modes <= 0 || n_modes > tNims) n_modes = tNims;

   t8 = get_curr_time();

   //Calculate eigenvectors and eigenvalues
   /* SYEVR sorts eigenvalues in ascending order, so we specifiy the top n_modes
    */   
   eigenSYEVR(cv, evecs, evals, tNims - n_modes, tNims);

   
   dsyevr += get_curr_time() - t8;

   //Normalize the eigenvectors
   //evals = (1./evals.sqrt());

    for(int i=0;i< n_modes; ++i)
    {
       evecs.col(i) = evecs.col(i)/sqrt(evals(i));
    }

   
   klims.resize(n_modes, tNpix);

   t11 = get_curr_time();
   //Now calculate KL images
   /*
    *  KL = E^T * R  ==> C = A^T * B
    */
   gemm<typename eigenT::Scalar>(CblasColMajor, CblasTrans, CblasTrans, n_modes, tNpix,
                              tNims, 1., evecs.data(), cv.rows(), Rims.data(), Rims.rows(),
                                 0., klims.data(), klims.rows());

   dgemm += get_curr_time() - t11;
} //calcKLIms


struct cvEntry
{
   int index;
   double cvVal;
};

//This sorts with greater-than
bool cvEntryComp( const cvEntry & cvE1, const cvEntry & cvE2)
{
   return ( fabs(cvE1.cvVal) > fabs(cvE2.cvVal) );
}

//This sorts with less-than
bool cvEntryCompIndex( const cvEntry & cvE1, const cvEntry & cvE2)
{
   return ( cvE1.index < cvE2.index );
}

template<typename eigenT, typename eigenTin>
void extractRowsAndCols(eigenT & out, const eigenTin & in, const std::vector<cvEntry> & idx)
{
   
   out.resize(idx.size(), idx.size());
   
   for(int i=0; i< idx.size(); ++i)
   {
      for(int j=0; j < idx.size(); ++j)
      {
         out(i,j) = in(idx[i].index, idx[j].index);
      }
   }
   
}

template<typename eigenT, typename eigenTin>
void extractCols(eigenT & out, const eigenTin & in, const std::vector<cvEntry> & idx)
{
   
   out.resize(in.rows(), idx.size()); 
   
   for(int i=0; i< idx.size(); ++i)
   {
      out.col(i) = in.col(idx[i].index); //it1->index);
   }
   
}

template<typename floatT, class derotFunctObj>
template<typename eigenT, typename eigenTv>
void KLIPreduction<floatT, derotFunctObj>::collapseCovar( eigenT & cutCV, 
                                                          const eigenT & CV,
                                                          const std::vector<floatT> & sds,
                                                          eigenT & rimsCut,
                                                          const eigenTv & rims,
                                                          int imno,
                                                          double dang )
{
   std::vector<cvEntry> allidx(this->Nims);
   
   
   //Initialize the vector cvEntries
   for(int i=0; i<this->Nims; ++i)
   {
      allidx[i].index = i;
      
      //CV is lower-triangular
      if(i <= imno)
      {
         allidx[i].cvVal = CV(imno,i)/(sds[i]*sds[imno]);
      }
      else
      {
         allidx[i].cvVal = CV(i,imno)/(sds[i]*sds[imno]);
      }
   }
   
   int rotoff0 = 0;
   int rotoff1 = 0;
   
   if(excludeMethod == HCI::excludePixel || excludeMethod == HCI::excludeAngle )
   {
      rotoff1 = this->Nims;
      
      //Find first rotoff within dang
      int j;
      for(j=0; j< this->Nims; ++j)
      {
         if( fabs(angleDiff<1>(this->derotF.derotAngle(j), this->derotF.derotAngle(imno))) <= dang )
         {
            rotoff0 = j;
            ++j;
            break;
         }
      }
      //Find first rotoff outside dang --> this is the first image that will be included again
      for(; j< this->Nims; ++j)
      {
         if( fabs(angleDiff<1>(this->derotF.derotAngle(j), this->derotF.derotAngle(imno))) > dang )
         {
            rotoff1 = j;
            break;
         }
      }
   }
   else if(excludeMethod == HCI::excludeImno)
   {
      rotoff1 = this->Nims;
      //Find first imno within dang
      int j;
      for(j=0; j< this->Nims; ++j)
      {
         if( fabs( j - imno ) <= dang )
         {
            rotoff0 = j;
            ++j;
            break;
         }
      }

      //Find last imno outside dang
      for(; j< this->Nims; ++j)
      {
         if( fabs( j -  imno ) > dang)
         {
            rotoff1 = j;
            break;
         }
      }
   }
   
   pout("rejecting", rotoff1-rotoff0, "images");
   if(rotoff1-rotoff0 > 0)
   {
      allidx.erase(allidx.begin()+rotoff0, allidx.begin()+rotoff1);
   }
   
   
   if( includeRefNum > 0 && includeRefNum < allidx.size())
   {
      //First partially sort the correlation values
      std::nth_element(allidx.begin(), allidx.begin()+ includeRefNum, allidx.end(), cvEntryComp);
      std::cout << "Minimum correlation: " << allidx[includeRefNum].cvVal << "\n";
      
      //Now delete the lower correlation values
      allidx.erase(allidx.begin()+includeRefNum, allidx.end());
      
      //Now we have to re-sort the remaining indices so that the CV matrix will still be U.T.
      std::sort(allidx.begin(), allidx.end(), cvEntryCompIndex);
   }
    
   std::cout << "Keeping " << allidx.size() << " reference images\n";
   

   extractRowsAndCols(cutCV, CV, allidx);
   extractCols(rimsCut, rims, allidx);
}

///@}

} //namespace mx

#endif //__KLIPreduction_hpp__
