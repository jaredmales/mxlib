

#include "ADIobservation.hpp"

#ifndef __KLIPreduction_hpp__
#define __KLIPreduction_hpp__

namespace mx
{
   
double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, tf;
double dcut, dcv, dklims, dgemm, dsyevr, dcfs, drot, dread;

template<typename _floatT, class _derotFunctObj>
struct KLIPreduction : public ADIobservation<_floatT, _derotFunctObj>
{
   typedef _floatT floatT;
   
   typedef Array<floatT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   int Nmodes;
   floatT dang;
   
   int dang_method;
   
   KLIPreduction()
   {
      Nmodes = 0;
      dang_method = 0;
   }
   
   KLIPreduction( const std::string & odir, 
                  const std::string & oprefix, 
                  const std::string & oext) : ADIobservation<_floatT, _derotFunctObj>(odir, oprefix, oext)
   {
      Nmodes = 0;
   }
   
   void meanSubtract();
   
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
   
   void worker(eigenCube<floatT> & rims, vector<size_t> & idx);
   
   template<typename eigenT, typename eigenT1>
   void calcKLIms( eigenT & klims, 
                   eigenT & cv, 
                   const eigenT1 & Rims, 
                   eigenT & evecs, 
                   int Nmodes = 0 );   
   
};

template<typename _floatT, class _derotFunctObj>
inline
void KLIPreduction<_floatT, _derotFunctObj>::meanSubtract()
{
   
   eigenImageT meanIm;
   
   this->imc.mean(meanIm);
      
   #pragma omp parallel for schedule(static, 1)
   for(int i=0;i<this->imc.planes(); ++i)
   {
      this->imc.image(i) -= meanIm;
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
   
   pout("Beginning\n");
      

   t1 = get_curr_time();
   
   if(!this->filesRead) this->readFiles();
   
   pout("Files read\n");
   dread = get_curr_time()-t1;
   
   pout("Sizes:", this->Nrows, this->Ncols, this->Nims);
   this->psfsub.resize(this->Nrows, this->Ncols, this->Nims);   
   this->psfsub.cube().setZero();

   
   pout("mean subtracting\n");
   //Subtract mean image
   meanSubtract();

   //Make radius and angle images
   eigenImagef rIm(this->Nrows,this->Ncols);
   eigenImagef qIm(this->Nrows,this->Ncols);
   radAngImage(rIm, qIm, .5*(this->Nrows-1), .5*(this->Ncols-1));

      pout("starting regions", minr.size());
   //******** For each region do this:
   for(int regno = 0; regno < minr.size(); ++regno)
   {
      pout("getting indices");
      vector<size_t> idx = imageRegionIndices(rIm, qIm, .5*(this->Nrows-1), .5*(this->Ncols-1), 
                                                    minr[regno], maxr[regno], minq[regno], maxq[regno]);
   
      //Create storage for the R-ims and psf-subbed Ims
      eigenCube<float> rims(idx.size(), 1, this->Nims);
   
      pout("cutting regions \n");
      t3 = get_curr_time();
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i< this->Nims; ++i)
      {
         auto rim = rims.image(i);
         cutImageRegion(rim, this->imc.image(i), idx, false);
      }
      dcut += get_curr_time() - t3;
   
    
      //*** Dispatch the work
      worker(rims, idx);
      pout("worker done");
      
   }
   pout("deroting");
   
   this->derotate();
   
   if(this->doFinimCombine)
   {
      this->combineFinim();
   
      fitsFile<floatT> f;
      f.write(this->finimName, this->finim.data(), this->finim.rows(), this->finim.cols());
   }
   
   if(this->doOutputPsfsub)
   {
      pout("outputting");
      this->outputPsfsub();
   }
   
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
   
}

template<typename floatT, class derotFunctObj>
inline
void KLIPreduction<floatT, derotFunctObj>::worker(eigenCube<floatT> & rims, vector<size_t> & idx)
{
   pout("beginning worker");
   
   int rotoff0, rotoff1;
   eigenImagef evecs;//, evals;
   eigenImagef cfs; //The coefficients
   eigenImagef psf;
   eigenImagef rims_cut;
   eigenImagef klims;
   eigenImagef cv;
   eigenImagef cv_cut;
   

   t5 = get_curr_time();

   //*** Form lower-triangle covariance matrix
   eigenSYRK(cv, rims.cube());
   
   dcv += get_curr_time() - t5;

   if(dang == 0)
   {
      pout("calculating K-L images");
      /**** Now calculate the K-L Images ****/
      t7 = get_curr_time();
      calcKLIms(klims, cv, rims.cube(), evecs, Nmodes);
      t13 = get_curr_time();
      dklims += t13 - t7;
   }
   
   for(int imno = 0; imno < this->Nims; ++imno)
   {
      pout("image:", imno, "/", this->Nims);
   
      if(dang != 0)
      {
         rotoff0 = 0;
         rotoff1 = this->Nims;
   
         if(dang_method == 0)
         {
            //Find first rotoff within dang
            int j;
            for(j=0; j< this->Nims; ++j)
            {
               if( fabs(angleDiff(this->derotF.derotAngle(j), this->derotF.derotAngle(imno))) <= DTOR(dang) )
               {
                  rotoff0 = j;
                  break;
               }
            }

            //Find last rotoff outside dang
            for(; j< this->Nims; ++j)
            {
               if( fabs(angleDiff(this->derotF.derotAngle(j), this->derotF.derotAngle(imno))) >= DTOR(dang) )
               {
                  rotoff1 = j;
                  break;
               }
            }
         }
         
         if(dang_method == 1)
         {
            //Find first rotoff within dang
            int j;
            for(j=0; j< this->Nims; ++j)
            {
               if( fabs( j - imno ) <= dang )
               {
                  rotoff0 = j;
                  break;
               }
            }

            //Find last rotoff outside dang
            for(; j< this->Nims; ++j)
            {
               if( fabs( j -  imno ) >= dang)
               {
                  rotoff1 = j;
                  break;
               }
            }
         }
   
         pout("rejecting", rotoff1-rotoff0, "images");

   
         /**** Now remove the rejected images ****/
 
         removeRowsAndCols(cv_cut, cv, rotoff0, rotoff1-rotoff0);
         removeCols(rims_cut, rims.asVectors(), rotoff0, rotoff1-rotoff0);
   
         /**** Now calculate the K-L Images ****/
         t7 = get_curr_time();
         calcKLIms(klims, cv_cut, rims_cut, evecs, Nmodes);
         t13 = get_curr_time();
         dklims += t13 - t7;

         pout("klims:", t13-t7, "secs");
   
      } //if(dang != 0)
      
      cfs.resize(1, klims.rows());
   
      t14 = get_curr_time();
   
      #pragma omp parallel for schedule(static, 1)
      for(int j=0; j<cfs.size(); ++j)
      {
         cfs(j) = klims.row(j).matrix().dot(rims.cube().col(imno).matrix());
      }
      dcfs += get_curr_time()-t14;
  
      psf = cfs(Nmodes-1)*klims.row(Nmodes-1);

      
      //Count down, since eigenvalues are returned in increasing order
      for(int j=Nmodes-2; j>=0; --j)
      {
         psf += cfs(j)*klims.row(j);
      }    
            
      insertImageRegion(this->psfsub.cube().col(imno), rims.cube().col(imno) - psf.transpose(), idx);
   }

}





template<typename floatT, class derotFunctObj>
template<typename eigenT, typename eigenT1>
inline
void KLIPreduction<floatT, derotFunctObj>::calcKLIms( eigenT & klims, 
                                                      eigenT & cv, 
                                                      const eigenT1 & Rims, 
                                                      eigenT & evecs,  
                                                      int Nmodes )
{

   eigenT evals;
   
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

   if(Nmodes <= 0 || Nmodes > tNims) Nmodes = tNims;

   t8 = get_curr_time();

   //Calculate eigenvectors and eigenvalues
   /* SYEVR sorts eigenvalues in ascending order, so we specifiy the top Nmodes
    */   
   eigenSYEVR(cv, evecs, evals, tNims - Nmodes, tNims);

   
   dsyevr += get_curr_time() - t8;

   //Normalize the eigenvectors
   //evals = (1./evals.sqrt());

   for(int i=0;i< Nmodes; ++i)
   {
      evecs.col(i) = evecs.col(i)/sqrt(evals(i));
   }

   
   klims.resize(Nmodes, tNpix);


   t11 = get_curr_time();
   //Now calculate KL images
   /*
    *  KL = E^T * R  ==> C = A^T * B
    */

   pout(evecs.rows(), evecs.cols());
   pout(cv.rows(), cv.cols());
   pout(Rims.rows(), Rims.cols());
   pout(klims.rows(), klims.cols());

   gemm<typename eigenT::Scalar>(CblasColMajor, CblasTrans, CblasTrans, Nmodes, tNpix,
                              tNims, 1., evecs.data(), cv.rows(), Rims.data(), Rims.rows(),
                                 0., klims.data(), klims.rows());

   dgemm += get_curr_time() - t11;
} //calcKLIms

} //namespace mx

#endif //__KLIPreduction_hpp__
