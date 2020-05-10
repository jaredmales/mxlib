/** \file KLIPreduction.hpp
  * \author Jared R. Males
  * \brief Declarations and definitions for an implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm. 
  * \ingroup hc_imaging_files
  * \ingroup image_processing_files
  *
  */



#ifndef __KLIPreduction_hpp__
#define __KLIPreduction_hpp__


#include <vector>
#include <map>

#include <omp.h>

#include "../ompLoopWatcher.hpp"
#include "../math/geo.hpp"
#include "../math/eigenLapack.hpp"

#include "ADIobservation.hpp"

namespace mx
{
namespace improc
{
   
//double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, tf;
//double dcut, dscv, dklims, dgemm, dsyevr, dcfs, drot, dcombo, dread;

   



namespace HCI 
{
   ///Image exclusion methods
   /** \ingroup hc_imaging_enums
     */
   enum excludeMethods{ excludeNone, ///< Exclude no images
                        excludePixel, ///< Exclude by pixels of rotation
                        excludeAngle, ///< Exclude by angle of roration
                        excludeImno ///< Exclude by number of images
                      };   

   ///Image inclusion methods
   /** \ingroup hc_imaging_enums
     */
   enum includeMethods{ includeAll,   ///< include all images
                        includeCorr,  ///< include images which are most correlated with the target
                        includeTime,  ///< include images which are closest in time to the target
                        includeAngle, ///< include images which are closest in angle to the target
                        includeImno   ///< include images which are closest in imno to the target
                      };

}

/// An implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm.
/** KLIP\cite soummer_2012 is a principle components analysis (PCA) based technique for PSF estimation.
  * 
  * 
  * \tparam _realT  is the floating point type in which to do calculations
  * \tparam _derotFunctObj the ADIobservation derotator class.
  * \tparam _evCalcT the real type in which to do eigen-decomposition.  Should generally be double for stable results.  
  * \ingroup hc_imaging 
  */ 
template<typename _realT, class _derotFunctObj, typename _evCalcT = double>
struct KLIPreduction : public ADIobservation<_realT, _derotFunctObj>
{
   typedef _realT realT;
   
   typedef Array<realT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   typedef _evCalcT evCalcT;
   
   int padSize;
   
   /// Specifies the number of modes to include in the PSF.
   /** The output image is a cube with a plane for each entry in m_Nmodes.
     * Only the number of eigenvalues required for the maximum value of m_Nmodes
     * are calculated, so this can have an important impact on speed.
     * 
     * Can be initialized as:
     * \code 
     * red.m_Nmodes={5,10,15,20};
     * \endcode
     * 
     */ 
   std::vector<int> m_Nmodes;
   
   int maxNmodes {0};
   
   /// Specify the minimum pixel difference at the inner edge of the search region 
   realT m_minDPx {0};
   
   /// Specify the maximum pixel difference at the inner edge of the search region 
   realT m_maxDPx {0};
   
   /// Controls how reference images are excluded, if at all, from the covariance matrix for each target image based on a minimum criterion.
   /** Can have the following values:
     *  - <b>HCI::excludeNone</b> = no exclusion, all images included [default]
     *  - <b>HCI::excludePixel</b> = exclude based on pixels of rotation at the inner edge of the region
     *  - <b>HCI::excludeAngle</b> = exclude based on degrees of rotation at the inner edge of the region
     *  - <b>HCI::excludeImno</b> = exclude based on number of images
     */
   int m_excludeMethod {HCI::excludeNone};
   
   /// Controls how reference images are excluded, if at all, from the covariance matrix for each target image based on a maximum criterion.
   /** Can have the following values:
     *  - <b>HCI::excludeNone</b> = no exclusion, all images included [default]
     *  - <b>HCI::excludePixel</b> = exclude based on pixels of rotation at the inner edge of the region
     *  - <b>HCI::excludeAngle</b> = exclude based on degrees of rotation at the inner edge of the region
     *  - <b>HCI::excludeImno</b> = exclude based on number of images
     */
   int m_excludeMethodMax  {HCI::excludeNone};
   
   ///Number of reference images to include in the covariance matrix
   /** If > 0, then at most this many images, determined by highest cross-correlation, are included.
     * This is determined after rotational/image-number exclusion. 
     * If == 0, then all reference images are included. 
     */
   int m_includeRefNum {0}; 
   
   /// Controls how number of included images is calculated.
   /** The number of included images is calculated after exclusion is complete.
     * Can have the following values:
     * - <b>HCI::includeAll</b> = all remaining images are included [default]
     * - <b>HCI::includeCorr</b> = the m_includeRefNum of the remaining images which are most correlated with the target are included
     * - <b>HCI::includeTime</b> = the m_includeRefNum of the remaining images which are closest in time to the target are included
     * - <b>HCI::includeAngle</b> = the m_includeRefNum of the remaining images which are closest in angle to the target are included
     * - <b>HCI::includeImno</b> = the m_includeRefNum of the remaining images which are closest in image number to the target are included
     */
   int m_includeMethod {HCI::includeAll};
   
   eigenImage<int> m_imsIncluded;
   
   KLIPreduction()
   {
      initialize();
   }
   
   KLIPreduction( const std::string & dir, 
                  const std::string & prefix, 
                  const std::string & ext = ".fits") : ADIobservation<_realT, _derotFunctObj>(dir, prefix, ext)
   {
      initialize();
   }
   
   explicit KLIPreduction( const std::string & fileListFile ) : ADIobservation<_realT, _derotFunctObj>(fileListFile)
   {
      initialize();
   }
   
   virtual ~KLIPreduction() {}
   
   void initialize()
   {
      padSize = 4;
      
      t_worker_begin =0;
      t_worker_end = 0;
      
      t_eigenv = 0;
      
      t_klim = 0;
      
      t_psf = 0;
      
   }
   
   void meanSubtract(eigenCube<realT> & ims, std::vector<realT> & sds);
   void medianSubtract(eigenCube<realT> & ims, std::vector<realT> & sds);
   void getStdDevs(std::vector<realT> sd, eigenCube<realT> & ims);
   
   /// Run KLIP in a set of geometric search regions.
   /** The arguments are 4 vectors, where each entry defines one component of the  search region.
     *
     * \returns 0 on success
     * \returns -1 on error
     * 
     */
   int regions( std::vector<realT> minr, ///< [in]
                std::vector<realT> maxr, ///< [in]
                std::vector<realT> minq, ///< [in]
                std::vector<realT> maxq  ///< [in]
              );
   
   /// Run KLIP in a geometric search region.
   /** \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     * 
     */
   int regions( realT minr, ///< [in]
                realT maxr, ///< [in]
                realT minq, ///< [in]
                realT maxq  ///< [in]
              )
   {
      std::vector<realT> vminr(1, minr);
      std::vector<realT> vmaxr(1, maxr);
      std::vector<realT> vminq(1, minq);
      std::vector<realT> vmaxq(1, maxq);
      
      return regions(vminr, vmaxr, vminq, vmaxq);
   }
   
   void worker(eigenCube<realT> & rims, std::vector<size_t> & idx, realT dang, realT dangMax);

   double t_worker_begin;
   double t_worker_end;
   
   double t_eigenv;
   double t_klim;
   double t_psf;
   
   void dump_times()
   {
      printf("KLIP reduction times: \n");
      printf("  Total time: %f sec\n", this->t_end-this->t_begin);
      printf("    Loading: %f sec\n", this->t_load_end-this->t_load_begin);
      printf("    Fake Injection: %f sec\n", this->t_fake_end-this->t_fake_begin);
      printf("    Coadding: %f sec\n", this->t_coadd_end-this->t_coadd_begin);
      printf("    Preprocessing: %f sec\n", this->t_preproc_end - this->t_preproc_begin);
      printf("      Az USM: %f sec\n", this->t_azusm_end - this->t_azusm_begin);
      printf("      Gauss USM: %f sec\n", this->t_gaussusm_end - this->t_gaussusm_begin);
      printf("    KLIP algorithm: %f elapsed real sec\n", this->t_worker_end - this->t_worker_begin);
      double klip_cpu = this->t_eigenv + this->t_klim + this->t_psf;
      printf("      EigenDecomposition %f cpu sec (%f%%)\n", this->t_eigenv, this->t_eigenv/klip_cpu*100);
      printf("      KL image calc %f cpu sec (%f%%)\n", this->t_klim, this->t_klim/klip_cpu*100);
      printf("      PSF calc/sub %f cpu sec (%f%%)\n", this->t_psf, this->t_psf/klip_cpu*100);
      printf("    Derotation: %f sec\n", this->t_derotate_end-this->t_derotate_begin);
      printf("    Combination: %f sec\n", this->t_combo_end-this->t_combo_begin);
   }

   int processPSFSub( const std::string & dir,
                      const std::string & prefix,
                      const std::string & ext
                    );

};

template<typename _realT, class _derotFunctObj, typename _evCalcT>
inline
void KLIPreduction<_realT, _derotFunctObj, _evCalcT>::meanSubtract(eigenCube<realT> & ims, std::vector<_realT> & norms)
{

   norms.resize(ims.planes());

   for(int n=0;n<ims.planes(); ++n)
   {
      ims.image(n) -= ims.image(n).mean();
      norms[n] = ims.image(n).matrix().norm();
   }
}
 
template<typename _realT, class _derotFunctObj, typename _evCalcT>
inline
void KLIPreduction<_realT, _derotFunctObj, _evCalcT>::medianSubtract(eigenCube<realT> & ims, std::vector<_realT> & sds)
{
         
   sds.resize(ims.planes());
   //#pragma omp parallel for schedule(static, 1)
   for(int i=0;i<ims.planes(); ++i)
   {
      _realT med = eigenMedian(ims.image(i));
      ims.image(i) -= med;
      sds[i] = ims.image(i).matrix().norm();//This is the standard deviation relative to the median.
   }
} 

template<typename _realT, class _derotFunctObj, typename _evCalcT>
inline
int KLIPreduction<_realT, _derotFunctObj, _evCalcT>::regions( std::vector<_realT> minr, 
                                                              std::vector<_realT> maxr, 
                                                              std::vector<_realT> minq, 
                                                              std::vector<_realT> maxq)
{   
   this->t_begin = get_curr_time();
   
   maxNmodes = m_Nmodes[0];
   for(size_t i = 1; i < m_Nmodes.size(); ++i)
   {
      if( m_Nmodes[i] > maxNmodes) maxNmodes = m_Nmodes[i];
   }
   
   std::cerr << "Beginning\n";
      
   if(this->imSize == 0)
   {
      this->imSize = 2*(*std::max_element(maxr.begin(),maxr.end()) + padSize);
   }
   
   if(!this->filesRead) 
   {         
      if( this->readFiles() < 0) return -1;
   }
   
   if(this->preProcess_only && !this->skipPreProcess)
   {
      std::cerr << "Pre-processing complete, stopping.\n";
      return 0;
   }

   std::cerr << "allocating psf subtracted cubes\n";
   this->psfsub.resize(m_Nmodes.size());
   for(size_t n=0;n<m_Nmodes.size(); ++n)
   {
      this->psfsub[n].resize(this->Nrows, this->Ncols, this->Nims);
      this->psfsub[n].cube().setZero();
   }
   
   //Make radius and angle images
   eigenImageT rIm(this->Nrows,this->Ncols);
   eigenImageT qIm(this->Nrows,this->Ncols);
   
   radAngImage(rIm, qIm, .5*(this->Nrows-1), .5*(this->Ncols-1));

   
   m_imsIncluded.resize(this->Nims,this->Nims);
   m_imsIncluded.setConstant(1);
   
   std::cerr << "starting regions " << minr.size() << "\n";
   //******** For each region do this:
   for(size_t regno = 0; regno < minr.size(); ++regno)
   {
      eigenImageT * maskPtr = 0;
      if( this->mask.rows() == this->Nrows && this->mask.cols() == this->Ncols) maskPtr = &this->mask;
      
      std::vector<size_t> idx = annulusIndices(rIm, qIm, .5*(this->Nrows-1), .5*(this->Ncols-1), 
                                                    minr[regno], maxr[regno], minq[regno], maxq[regno], maskPtr);
   
      //Create storage for the R-ims and psf-subbed Ims
      eigenCube<realT> rims(idx.size(), 1, this->Nims);
   
      //#pragma omp parallel for schedule(static, 1)
      for(int i=0;i< this->Nims; ++i)
      {
         auto rim = rims.image(i);
         cutImageRegion(rim, this->imc.image(i), idx, false);
      }

      realT dang = 0;
      realT dangMax = 0;
      
      if(m_minDPx < 0) m_excludeMethod = HCI::excludeNone;
      if(m_maxDPx < 0) m_excludeMethodMax = HCI::excludeNone;
      
      if(m_excludeMethod == HCI::excludePixel)
      {
         dang = fabs(atan(m_minDPx/minr[regno]));
      }
      else if(m_excludeMethod == HCI::excludeAngle)
      {
         dang = math::dtor(m_minDPx);
      }
      else if(m_excludeMethod == HCI::excludeImno)
      {
         dang = m_minDPx;
      }
      
      
      if(m_excludeMethodMax == HCI::excludePixel)
      {
         dangMax = fabs(atan(m_maxDPx/minr[regno]));
      }
      else if(m_excludeMethodMax == HCI::excludeAngle)
      {
         dangMax = math::dtor(m_maxDPx);
      }
      else if(m_excludeMethodMax == HCI::excludeImno)
      {
         dangMax = m_maxDPx;
      }
      
      //*** Dispatch the work
      worker(rims, idx, dang, dangMax);
      std::cerr << "worker done\n";
      
   }
   
   fitsFile<int> ffii;
   ffii.write("imsIncluded.fits", m_imsIncluded);
   
   if(this->doDerotate)
   {
      std::cerr << "derotating\n";
      this->derotate();
   }
   
   
   if(this->combineMethod > 0)
   {
      std::cerr << "combining\n";
      this->combineFinim();
      
   }
   
   if(this->doWriteFinim == true || this->doOutputPSFSub == true)
   {
      std::cerr << "writing\n";
      
      fitsHeader head;
      
      this->ADIobservation<_realT, _derotFunctObj>::fitsHeader(&head);
      
      head.append("", fitsCommentType(), "----------------------------------------");
      head.append("", fitsCommentType(), "mx::KLIPreduction parameters:");
      head.append("", fitsCommentType(), "----------------------------------------");
   
      std::stringstream str;
      
      if(m_Nmodes.size() > 0)
      {
         for(size_t nm=0;nm < m_Nmodes.size()-1; ++nm) str << m_Nmodes[nm] << ",";
         str << m_Nmodes[m_Nmodes.size()-1];      
         head.append<char *>("NMODES", (char *)str.str().c_str(), "number of modes");
      }
      
      if(minr.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < minr.size()-1; ++nm) str << minr[nm] << ",";
         str << minr[minr.size()-1];      
         head.append<char *>("REGMINR", (char *)str.str().c_str(), "region inner edge(s)");
      }
      
      if(maxr.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < maxr.size()-1; ++nm) str << maxr[nm] << ",";
         str << maxr[maxr.size()-1];      
         head.append<char *>("REGMAXR", (char *)str.str().c_str(), "region outer edge(s)");
      }
      
      if(minq.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < minq.size()-1; ++nm) str << minq[nm] << ",";
         str << minq[minq.size()-1];      
         head.append<char *>("REGMINQ", (char *)str.str().c_str(), "region minimum angle(s)");
      }
      
      if(maxq.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < maxq.size()-1; ++nm) str << maxq[nm] << ",";
         str << maxq[maxq.size()-1];      
         head.append<char *>("REGMAXQ", (char *)str.str().c_str(), "region maximum angle(s)");
      }
      
      head.append<int>("EXCLMTHD", m_excludeMethod, "value of excludeMethod");
      head.append<realT>("MINDPX", m_minDPx, "minimum pixel delta");
      head.append<realT>("MAXDPX", m_minDPx, "maximum pixel delta");
      head.append<int>("INCLREFN", m_includeRefNum, "value of includeRefNum");

      if(this->doWriteFinim == true && this->combineMethod > 0)
      {
         this->writeFinim(&head);
      }
      
      if(this->doOutputPSFSub)
      {
         this->outputPSFSub(&head);
      }
   }
   
   this->t_end = get_curr_time();
 
   dump_times();
   
   return 0;
}

struct cvEntry
{
   int index;
   double cvVal;
   double angle;
   bool included {true};
};


template<typename eigenT, typename eigenTin>
void extractRowsAndCols(eigenT & out, const eigenTin & in, const std::vector<size_t> & idx)
{
   
   out.resize(idx.size(), idx.size());
   
   for(size_t i=0; i< idx.size(); ++i)
   {
      for(size_t j=0; j < idx.size(); ++j)
      {
         out(i,j) = in(idx[i], idx[j]);
      }
   }
   
}

template<typename eigenT, typename eigenTin>
void extractCols(eigenT & out, const eigenTin & in, const std::vector<size_t> & idx)
{
   
   out.resize(in.rows(), idx.size()); 
   
   for(size_t i=0; i< idx.size(); ++i)
   {
      out.col(i) = in.col(idx[i]); //it1->index);
   }
   
}

template<typename realT, typename eigenT, typename eigenTv, class derotFunctObj>
void collapseCovar( eigenT & cutCV, 
                    const eigenT & CV,
                    const std::vector<realT> & sds,
                    eigenT & rimsCut,
                    const eigenTv & rims,
                    int imno,
                    double dang,
                    double dangMax,
                    int Nims,
                    int excludeMethod,
                    int excludeMethodMax,
                    int includeRefNum,
                    const derotFunctObj & derotF,
                    eigenImage<int> & imsIncluded
                  )
{   
   std::vector<cvEntry> allidx(Nims);
   
   std::cerr << "dangs: " << dang << " " << dangMax << "\n";
   
   
   //Initialize the vector of cvEntries
   for(int i=0; i < Nims; ++i)
   {
      allidx[i].index = i;
      allidx[i].angle = derotF.derotAngle(i);
      
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
   
   if(excludeMethod == HCI::excludePixel || excludeMethod == HCI::excludeAngle )
   {
      for(size_t j=0; j < Nims; ++j)
      {
         if( fabs(math::angleDiff<math::radians>( derotF.derotAngle(j), derotF.derotAngle(imno))) <= dang ) allidx[j].included = false;
      }
   }
   else if(excludeMethod == HCI::excludeImno)
   {
      for(size_t j=0; j < Nims; ++j) 
      {
         if( fabs( (long) j - imno)  <= dang ) allidx[j].included = false;
      }      
   }
      
   if(excludeMethodMax == HCI::excludePixel || excludeMethodMax == HCI::excludeAngle )
   {
      for(size_t j=0; j < Nims; ++j)
      {
         if( fabs(math::angleDiff<math::radians>( derotF.derotAngle(j), derotF.derotAngle(imno))) > dangMax ) allidx[j].included = false;
      }
   }
   else if(excludeMethodMax == HCI::excludeImno)
   {
      for(size_t j=0; j < Nims; ++j)
      {
         if( fabs( (long) j - imno)  > dangMax ) allidx[j].included = false;
      }  
   }      
   
     
   if( includeRefNum > 0 && (size_t) includeRefNum < allidx.size())
   {
      long kept = 0;
      for(size_t j=0; j < Nims; ++j)
      {
         if(allidx[j].included == true) ++kept;
      }
      
      //Get a vector for sorting
      std::vector<realT> cvVal;
      cvVal.resize(kept);
      size_t k = 0;
      for(size_t j=0; j < Nims; ++j)
      {
         if(allidx[j].included == true) 
         {  
            cvVal[k] = allidx[j].cvVal;
            ++k;
         }
      }
      
      //Partially sort the correlation values
      std::nth_element(cvVal.begin(), cvVal.begin()+ (kept-includeRefNum), cvVal.end());
      
      realT mincorr = cvVal[kept-includeRefNum];
      std::cerr << "    Minimum correlation: " << mincorr << "\n";
      
      
      for(size_t j=0; j < Nims; ++j)
      {
         if( allidx[j].cvVal < mincorr ) allidx[j].included = false;
      }
         
   }

   std::vector<size_t> keepidx;
   for(size_t j=0;j<Nims;++j)
   {
      imsIncluded(imno,j) = allidx[j].included;
      
      if(allidx[j].included) keepidx.push_back(j);
   }
   
   std::cerr << "  Keeping " << keepidx.size() << " reference images out of " << Nims << " (" << Nims-keepidx.size() << " rejected)\n";
   
   if(keepidx.size() == 0)
   {
      std::cerr << "\n\n" << imno << "\n\n";
   }

   extractRowsAndCols(cutCV, CV, keepidx);
   extractCols(rimsCut, rims, keepidx);
   
}


template<typename _realT, class _derotFunctObj, typename _evCalcT>
void KLIPreduction<_realT, _derotFunctObj, _evCalcT>::worker(eigenCube<_realT> & rims, std::vector<size_t> & idx, realT dang, realT dangMax)
{
   std::cerr << "beginning worker\n";

   t_worker_begin = get_curr_time();
   
   std::vector<realT> sds;

   //*** First mean subtract ***//   
   meanSubtract(rims, sds);  

   //*** Form lower-triangle covariance matrix      
   eigenImageT cv;
 
   math::eigenSYRK(cv, rims.cube());
      
   fitsFile<realT> ff;
   ff.write("cv.fits", cv);
   ompLoopWatcher<> status( this->Nims, std::cerr);
   
   //Pre-calculate KL images once if we are exclude none
   eigenImageT master_klims;
   if( m_excludeMethod == HCI::excludeNone && m_excludeMethodMax == HCI::excludeNone && m_includeRefNum == 0)
   {
      double teigenv;
      double tklim;
      math::calcKLModes<double>(master_klims, cv, rims.cube(), maxNmodes, nullptr, &teigenv, &tklim);
      t_eigenv += teigenv;
      t_klim += tklim;
   }
      
   //int nTh = 0;
   #pragma omp parallel //num_threads(20) 
   {
      //We need local copies for each thread.  Only way this works, for whatever reason.
      eigenImageT cfs; //The coefficients
      eigenImageT psf;
      eigenImageT rims_cut;
      eigenImageT cv_cut;
      eigenImageT klims;
      
      math::syevrMem<evCalcT> mem;

      if( m_excludeMethod == HCI::excludeNone && m_excludeMethodMax == HCI::excludeNone && m_includeRefNum == 0 )
      {
         klims = master_klims;
      }
   
      #pragma omp for 
      for(int imno = 0; imno < this->Nims; ++imno)
      {
         status.incrementAndOutputStatus();
         
         if( m_excludeMethod != HCI::excludeNone || m_excludeMethodMax != HCI::excludeNone || m_includeRefNum != 0 )
         {         

            collapseCovar<realT>( cv_cut,  cv, sds, rims_cut, rims.asVectors(), imno, dang, dangMax, this->Nims, this->m_excludeMethod, this->m_excludeMethodMax, this->m_includeRefNum, this->derotF, m_imsIncluded);
            
            /**** Now calculate the K-L Images ****/
            double teigenv, tklim;
            math::calcKLModes(klims, cv_cut, rims_cut, maxNmodes, &mem, &teigenv, &tklim);
            t_eigenv += teigenv;
            t_klim += tklim;
         }
         cfs.resize(1, klims.rows());
   
  
         double t0 = get_curr_time();
         
         for(int j=0; j<cfs.size(); ++j)
         {
            cfs(j) = klims.row(j).matrix().dot(rims.cube().col(imno).matrix());    
         }

         for(size_t mode_i =0; mode_i < m_Nmodes.size(); ++mode_i)
         {
            psf = cfs(cfs.size()-1)*klims.row(cfs.size()-1);

            //Count down, since eigenvalues are returned in increasing order
            //  handle case where cfs.size() < m_Nmodes[mode_i], i.e. when more modes than images.
            for(int j=cfs.size()-2; j>=cfs.size()-m_Nmodes[mode_i] && j >= 0; --j)
            {
               psf += cfs(j)*klims.row(j);
            }  
            
            //#pragma omp critical
            insertImageRegion( this->psfsub[mode_i].cube().col(imno), rims.cube().col(imno) - psf.transpose(), idx);
         }
         

         t_psf += (get_curr_time() - t0) ;/// omp_get_num_threads();
         
         
      } //for imno
   }//openmp parrallel  
   
   t_worker_end = get_curr_time();
}


template<typename _realT, class _derotFunctObj, typename _evCalcT>
inline
int KLIPreduction<_realT, _derotFunctObj, _evCalcT>::processPSFSub( const std::string & dir,
                                                                    const std::string & prefix,
                                                                    const std::string & ext
                                                                  )

{   
  
   
   std::cerr << "Beginning\n";
      
   this->skipPreProcess = true;

   this->readPSFSub(dir, prefix, ext, m_Nmodes.size());
   
   
   //This is generic to both regions and this from here on . . .
   
   if(this->doDerotate)
   {
      std::cerr << "derotating\n";
      this->derotate();
   }
   
   
   if(this->combineMethod > 0)
   {
      std::cerr << "combining\n";
      this->combineFinim();
      
   }
   
   if(this->doWriteFinim == true || this->doOutputPSFSub == true)
   {
      std::cerr << "writing\n";
      
      fitsHeader head;
      
      this->ADIobservation<_realT, _derotFunctObj>::fitsHeader(&head);
      
      head.append("", fitsCommentType(), "----------------------------------------");
      head.append("", fitsCommentType(), "mx::KLIPreduction parameters:");
      head.append("", fitsCommentType(), "----------------------------------------");
   
      std::stringstream str;
      
      if(m_Nmodes.size() > 0)
      {
         for(size_t nm=0;nm < m_Nmodes.size()-1; ++nm) str << m_Nmodes[nm] << ",";
         str << m_Nmodes[m_Nmodes.size()-1];      
         head.append<char *>("NMODES", (char *)str.str().c_str(), "number of modes");
      }
            
      head.append<int>("EXCLMTHD", m_excludeMethod, "value of excludeMethod");
      head.append<realT>("MINDPX", m_minDPx, "minimum pixel delta");
      head.append<realT>("MAXDPX", m_maxDPx, "maximum pixel delta");
      head.append<int>("INCLREFN", m_includeRefNum, "value of includeRefNum");

      if(this->doWriteFinim == true && this->combineMethod > 0)
      {
         this->writeFinim(&head);
      }
      
      if(this->doOutputPSFSub)
      {
         this->outputPSFSub(&head);
      }
   }
   
   
   return 0;
}

///@}

} //namespace improc
} //namespace mx

#endif //__KLIPreduction_hpp__
