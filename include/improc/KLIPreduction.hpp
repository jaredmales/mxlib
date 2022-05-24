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

#include "../ipc/ompLoopWatcher.hpp"
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
   
   ///Mean subtraction methods
   /** These control how the data in each search region is centered to meet the PCA requirement.
     * \ingroup hc_imaging_enums 
     */
   enum meansubMethods{ imageMean,   ///< The mean of each image (within the search region) is subtracted from itself
                        imageMedian, ///< The median of each image (within the search region) is subtracted from itself
                        imageMode,   ///< The mode of each image (within the search region) is subtracted from itself
                        meanImage,   ///< The mean image of the data is subtracted from each image
                        medianImage  ///< The median image of the data is subtracted from each image
                      };
   
   std::string meansubMethodStr( int method );
   
   int meansubMethodFmStr( const std::string & method );
   
   ///Image exclusion methods
   /** \ingroup hc_imaging_enums
     */
   enum excludeMethods{ excludeNone, ///< Exclude no images
                        excludePixel, ///< Exclude by pixels of rotation
                        excludeAngle, ///< Exclude by angle of roration
                        excludeImno ///< Exclude by number of images
                      };   
                      
   std::string excludeMethodStr(int method);
                      
   int excludeMethodFmStr(const std::string & method);
   
   ///Image inclusion methods
   /** \ingroup hc_imaging_enums
     */
   enum includeMethods{ includeAll,   ///< include all images
                        includeCorr,  ///< include images which are most correlated with the target
                        includeTime,  ///< include images which are closest in time to the target
                        includeAngle, ///< include images which are closest in angle to the target
                        includeImno   ///< include images which are closest in imno to the target
                      };
                      
   std::string includeMethodStr( int method );
                     
   int includeMethodFmStr( const std::string & method );
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
   
   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> eigenImageT;
   
   typedef _evCalcT evCalcT;
   
   /// Default c'tor
   KLIPreduction();
   
   /// Construct and load the target file list.
   /** Populates the \ref m_fileList vector by searching on disk for files which match
     * "dir/prefix*.ext".  See \ref load_fileList
     *
     */
   KLIPreduction( const std::string & dir,    ///< [in] the directory to search for files.
                  const std::string & prefix, ///< [in] the prefix of the files
                  const std::string & ext     ///< [in] the extension of the files
                );
   
   /// Construct using a file containing the target file list
   /** Populates the \ref m_fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     */
   explicit KLIPreduction( const std::string & fileListFile /**< [in] the full path to a file containing a list of files */ );
   
   //Construct and load the target file list and the RDI file list
   /** Populates the \ref m_fileList vector by searching on disk for files which match
     * "dir/prefix*.ext".  See \ref load_fileList
     *
     *  Populates the \ref m_RDIfileList vector by searching on disk for files which match
     * "RDIdir/RDIprefix*.RDIext".  See \ref load_RDIfileList
     */   
   KLIPreduction( const std::string & dir,       ///< [in] the directory to search for files.
                  const std::string & prefix,    ///< [in] the prefix of the files
                  const std::string & ext,       ///< [in] the extension of the files, default is .fits
                  const std::string & RDIdir,    ///< [in] the directory to search for the reference files.
                  const std::string & RDIprefix, ///< [in] the initial part of the file name for the reference files.  Can be empty "".
                  const std::string & RDIext=""  ///< [in] [optional] the extension to append to the RDI file name, must include the '.'.  If empty "" then same extension as target files is used.
                );
   
   /// Construct using a file containing the target file list and a file containing the RDI target file list
   /** Populates the \ref m_fileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     * 
     * Populates the \ref m_RDIfileList vector by reading the file, which should be a single
     * column of new-line delimited file names.
     */
   explicit KLIPreduction( const std::string & fileListFile,   ///< [in] a file name path to read for the target file names.
                           const std::string & RDIfileListFile ///< [in] a file name path to read for the reference file names.
                        );
   
   
   virtual ~KLIPreduction();
   
   
   
   
   
   
   
   /// Specify how the data are centered for PCA within each search region
   /** Can have the following values:
     * - <b>HCI::imageMean</b> = the mean of each image (within the search region) is subtracted from itself
     * - <b>HCI::imageMedian</b> = the median of each image (within the search region) is subtracted from itself
     * - <b>HCI::imageMode</b>  = the mode of each image (within the search region) is subtracted from itself
     * - <b>HCI::meanImage</b> = the mean image of the data is subtracted from each image
     * - <b>HCI::medianImage</b> = the median image of the data is subtracted from each image
     */
   int m_meanSubMethod {HCI::imageMean};
   
   
   int m_padSize {4};
   
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
   
   int m_maxNmodes {0};
   
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

      
   /// Subtract the basis mean from each of the images
   /** The mean is subtracted according to m_meanSubMethod.
    */
   void meanSubtract( eigenCube<realT> & rims, ///< [in/out] The reference images.  These are mean subtracted on output.
                      eigenCube<realT> & tims, ///< [in/out] The target images, which can be the same cube as rims (tested by pointer comparison), in which case they will be ignored.  Mean subtractedon output.
                      std::vector<realT> & sds ///< [out] The standard deviation of the mean subtracted refernce images.
                    );

   std::vector<_realT> m_minr; 
   std::vector<_realT> m_maxr;
   std::vector<_realT> m_minq;
   std::vector<_realT> m_maxq;
                                                              
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
   
   void worker( eigenCube<realT> & rims,
                eigenCube<realT> & tims,
                std::vector<size_t> & idx, 
                realT dang, 
                realT dangMax
              );

   int finalProcess();
   
   int processPSFSub( const std::string & dir,
                      const std::string & prefix,
                      const std::string & ext
                    );
   
   double t_worker_begin {0};
   double t_worker_end {0};
   
   double t_eigenv {0};
   double t_klim {0};
   double t_psf {0};
   
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



};

template<typename _realT, class _derotFunctObj, typename _evCalcT>
KLIPreduction<_realT, _derotFunctObj, _evCalcT>::KLIPreduction()
{
}

template<typename _realT, class _derotFunctObj, typename _evCalcT>
KLIPreduction<_realT, _derotFunctObj, _evCalcT>::KLIPreduction( const std::string & dir, 
                                                                const std::string & prefix, 
                                                                const std::string & ext
                                                              ) : ADIobservation<_realT, _derotFunctObj>(dir, prefix, ext)
{
}
   
template<typename _realT, class _derotFunctObj, typename _evCalcT>
KLIPreduction<_realT, _derotFunctObj, _evCalcT>::KLIPreduction( const std::string & fileListFile
                                                              ) : ADIobservation<_realT, _derotFunctObj>(fileListFile) 
{
}

template<typename _realT, class _derotFunctObj, typename _evCalcT>
KLIPreduction<_realT, _derotFunctObj, _evCalcT>::KLIPreduction( const std::string & dir, 
                                                                const std::string & prefix, 
                                                                const std::string & ext,
                                                                const std::string & RDIdir, 
                                                                const std::string & RDIprefix, 
                                                                const std::string & RDIext
                                                              ) : ADIobservation<_realT, _derotFunctObj>(dir, prefix, ext, RDIdir, RDIprefix, RDIext)
{
   std::cerr << "KLIP 6\n";
}
   
template<typename _realT, class _derotFunctObj, typename _evCalcT>
KLIPreduction<_realT, _derotFunctObj, _evCalcT>::KLIPreduction( const std::string & fileListFile,
                                                                const std::string & RDIfileListFile
                                                              ) : ADIobservation<_realT, _derotFunctObj>(fileListFile, RDIfileListFile) 
{
}

template<typename _realT, class _derotFunctObj, typename _evCalcT>
KLIPreduction<_realT, _derotFunctObj, _evCalcT>::~KLIPreduction() 
{
}
   
template<typename _realT, class _derotFunctObj, typename _evCalcT>
void KLIPreduction<_realT, _derotFunctObj, _evCalcT>::meanSubtract( eigenCube<realT> & ims,
                                                                    eigenCube<realT> & tims,
                                                                    std::vector<_realT> & norms
                                                                  )
{

   norms.resize(ims.planes());

   if(m_meanSubMethod == HCI::meanImage || m_meanSubMethod == HCI::medianImage)
   {
      eigenImageT mean;
      
      if(m_meanSubMethod == HCI::meanImage)
      {
         ims.mean(mean);
      }
      else if(m_meanSubMethod == HCI::medianImage)
      {
         ims.median(mean);
      }  
      
      for(int n=0;n<ims.planes(); ++n)
      {
         ims.image(n) -= mean;
         
         realT immean = ims.image(n).mean();
         norms[n] = (ims.image(n)-immean).matrix().norm();
      }
      
      if(&tims != &ims)
      {
         for(int n=0;n<tims.planes(); ++n)
         {
            tims.image(n) -= mean;
         }
      }
   }
   else
   {
      realT mean;
      std::vector<realT> work; //Working memmory for median calc
      
      for(int n=0;n<ims.planes(); ++n)
      {
         if( m_meanSubMethod == HCI::imageMean )
         {
            mean = ims.image(n).mean();
         }
         else if(m_meanSubMethod == HCI::imageMedian)
         {
            mean = imageMedian(ims.image(n), &work);
         }
         
         ims.image(n) -= mean;
         
         //Because we might not have used the mean, we need to re-mean to make this the standard deviation
         realT immean = ims.image(n).mean();
         norms[n] = (ims.image(n)-immean).matrix().norm();
         
      }
      
      if(&tims != &ims)
      {
         for(int n=0;n<tims.planes(); ++n)
         {
            if( m_meanSubMethod == HCI::imageMean )
            {  
               mean = tims.image(n).mean();
            }
            else if(m_meanSubMethod == HCI::imageMedian)
            {
               mean = imageMedian(tims.image(n), &work);
            }
         
            tims.image(n) -= mean;
         }
      }
   }
}
 


template<typename _realT, class _derotFunctObj, typename _evCalcT>
int KLIPreduction<_realT, _derotFunctObj, _evCalcT>::regions( std::vector<_realT> minr, 
                                                              std::vector<_realT> maxr, 
                                                              std::vector<_realT> minq, 
                                                              std::vector<_realT> maxq
                                                            )
{   
   this->t_begin = sys::get_curr_time();
   
   m_minr = minr;
   m_maxr = maxr;
   m_minq = minq;
   m_maxq = maxq;
   
   
   m_maxNmodes = m_Nmodes[0];
   for(size_t i = 1; i < m_Nmodes.size(); ++i)
   {
      if( m_Nmodes[i] > m_maxNmodes) m_maxNmodes = m_Nmodes[i];
   }
   
   std::cerr << "Beginning\n";
      
   if(this->m_imSize == 0)
   {
      this->m_imSize = 2*(*std::max_element(maxr.begin(),maxr.end()) + m_padSize);
   }
   
   if(!this->m_filesRead) 
   {         
      if( this->readFiles() < 0) return -1;
   }
   
   //CHECK IF RDI HERE
   if(!this->m_RDIfilesRead && this->m_RDIfileList.size() != 0 )
   {
      if( this->readRDIFiles() < 0) return -1;
   }
   
   if(this->m_preProcess_only && !this->m_skipPreProcess)
   {
      std::cerr << "Pre-processing complete, stopping.\n";
      return 0;
   }

   std::cerr << "allocating psf subtracted cubes\n";
   this->m_psfsub.resize(m_Nmodes.size());
   for(size_t n=0;n<m_Nmodes.size(); ++n)
   {
      this->m_psfsub[n].resize(this->m_Nrows, this->m_Ncols, this->m_Nims);
      this->m_psfsub[n].cube().setZero();
   }
   
   //Make radius and angle images
   eigenImageT rIm(this->m_Nrows,this->m_Ncols);
   eigenImageT qIm(this->m_Nrows,this->m_Ncols);
   
   radAngImage(rIm, qIm, .5*(this->m_Nrows-1), .5*(this->m_Ncols-1));

   m_imsIncluded.resize(this->m_Nims,this->m_Nims);
   m_imsIncluded.setConstant(1);
   
   std::cerr << "starting regions " << minr.size() << "\n";
   
   //******** For each region do this:
   for(size_t regno = 0; regno < minr.size(); ++regno)
   {
      eigenImageT * maskPtr = 0;
      if( this->m_mask.rows() == this->m_Nrows && this->m_mask.cols() == this->m_Ncols) maskPtr = &this->m_mask;
      
      std::vector<size_t> idx = annulusIndices(rIm, qIm, .5*(this->m_Nrows-1), .5*(this->m_Ncols-1), 
                                                    minr[regno], maxr[regno], minq[regno], maxq[regno], maskPtr);
   
      //Create storage for the R-ims and psf-subbed Ims
      eigenCube<realT> tims(idx.size(), 1, this->m_Nims);

      //------If doing RDI, create bims
      eigenCube<realT> rims;
      
      if(this->m_refIms.planes() > 0)
      {
         rims.resize(idx.size(), 1, this->m_Nims);
      }
      
      //#pragma omp parallel for schedule(static, 1)
      for(int i=0;i< this->m_Nims; ++i)
      {
         auto tim = tims.image(i);
         cutImageRegion(tim, this->m_tgtIms.image(i), idx, false);
      }
      
      for(int p=0;p<this->m_refIms.planes();++p)
      {
         auto rim = rims.image(p);
         cutImageRegion(rim, this->m_refIms.image(p), idx, false);
      }

      realT dang = 0;
      realT dangMax = 0;

      if(m_minDPx < 0) m_excludeMethod = HCI::excludeNone;
      if(m_maxDPx < 0) m_excludeMethodMax = HCI::excludeNone;

      //------- If doing RDI, excludeMethod and excludeMethodMax must be none!
      if(this->m_refIms.planes() > 0)
      {
         m_excludeMethod = HCI::excludeNone;
         m_excludeMethodMax = HCI::excludeNone;
      }

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
      
      //------- If doing RDI, call this with rims, bims
      //*** Dispatch the work
      if(this->m_refIms.planes() > 0) //RDI
      {
         std::cerr << "\n\n******* RDI MODE **********\n\n";
         worker(rims, tims, idx, dang, dangMax);
      }
      else //ADI
      {
         std::cerr << "\n\n******* ADI MODE **********\n\n";
         worker(tims, tims, idx, dang, dangMax);
      }
      std::cerr << "worker done\n";
      
   }
   
   fits::fitsFile<int> ffii;
   ffii.write("imsIncluded.fits", m_imsIncluded);
   

   if(finalProcess() < 0)
   {
      std::cerr << "Error in final processing\n";
   }
   
   this->t_end = sys::get_curr_time();
 
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
   
   //std::cerr << "dangs: " << dang << " " << dangMax << "\n";
   
   
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
   
   //std::cerr << "  Keeping " << keepidx.size() << " reference images out of " << Nims << " (" << Nims-keepidx.size() << " rejected)\n";
   
   if(keepidx.size() == 0)
   {
      std::cerr << "\n\n" << imno << "\n\n";
   }

   extractRowsAndCols(cutCV, CV, keepidx);
   extractCols(rimsCut, rims, keepidx);
   
}


template<typename _realT, class _derotFunctObj, typename _evCalcT>
void KLIPreduction<_realT, _derotFunctObj, _evCalcT>::worker( eigenCube<_realT> & rims,
                                                              eigenCube<_realT> & tims,
                                                              std::vector<size_t> & idx, 
                                                              realT dang, 
                                                              realT dangMax
                                                            )
{
   std::cerr << "beginning worker\n";

   t_worker_begin = sys::get_curr_time();
   
   std::vector<realT> sds;

   eigenImageT meanim;
   
   //*** First mean subtract ***//   
   meanSubtract(rims, tims, sds);  

   //*** Form lower-triangle covariance matrix      
   eigenImageT cv;
 
   math::eigenSYRK(cv, rims.cube());
    
   fits::fitsFile<realT> ff;
   ff.write("cv.fits", cv);
   ipc::ompLoopWatcher<> status( this->m_Nims, std::cerr);
   
   //Pre-calculate KL images once if we are exclude none OR IF RDI
   eigenImageT master_klims;
   if( m_excludeMethod == HCI::excludeNone && m_excludeMethodMax == HCI::excludeNone && m_includeRefNum == 0)
   {
      double teigenv;
      double tklim;
      
      std::cerr << cv.rows() << " " << cv.cols() << " " << rims.rows() << " " << rims.cols() << " " << rims.planes() << " " << m_maxNmodes << "\n";
      math::calcKLModes<double>(master_klims, cv, rims.cube(), m_maxNmodes, nullptr, &teigenv, &tklim);
      
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

      if( m_excludeMethod == HCI::excludeNone && m_excludeMethodMax == HCI::excludeNone && m_includeRefNum == 0 ) //OR RDI
      {
         klims = master_klims;
      }
   
      #pragma omp for 
      for(int imno = 0; imno < this->m_Nims; ++imno)
      {
         status.incrementAndOutputStatus();
         
         if( m_excludeMethod != HCI::excludeNone || m_excludeMethodMax != HCI::excludeNone || m_includeRefNum != 0 )
         {         
            collapseCovar<realT>( cv_cut,  cv, sds, rims_cut, rims.asVectors(), imno, dang, dangMax, this->m_Nims, this->m_excludeMethod, this->m_excludeMethodMax, this->m_includeRefNum, this->m_derotF, m_imsIncluded);
            
            /**** Now calculate the K-L Images ****/
            double teigenv, tklim;
            math::calcKLModes(klims, cv_cut, rims_cut, m_maxNmodes, &mem, &teigenv, &tklim);
            t_eigenv += teigenv;
            t_klim += tklim;
         }
         cfs.resize(1, klims.rows());
   
  
         double t0 = sys::get_curr_time();
         
         for(int j=0; j<cfs.size(); ++j)
         {
            cfs(j) = klims.row(j).matrix().dot(tims.cube().col(imno).matrix());    
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
            insertImageRegion( this->m_psfsub[mode_i].cube().col(imno), tims.cube().col(imno) - psf.transpose(), idx);
         }
         

         t_psf += (sys::get_curr_time() - t0) ;/// omp_get_num_threads();
         
         
      } //for imno
   }//openmp parrallel  
   
   t_worker_end = sys::get_curr_time();
}

template<typename _realT, class _derotFunctObj, typename _evCalcT>
int KLIPreduction<_realT, _derotFunctObj, _evCalcT>::finalProcess()
{
   if(this->m_postMedSub)
   {
      std::cerr << "Subtracting medians in post\n";
      
      for(size_t n=0; n<this->m_psfsub.size(); ++n)
      {
         #pragma omp parallel
         {
            eigenImage<realT> medim;
         
            this->m_psfsub[n].median(medim);
         
            #pragma omp for
            for(int i=0; i<this->m_psfsub[n].planes();++i)
            {
               this->m_psfsub[n].image(i) -= medim;
            }
         }
      }   
   }
   
   if(this->m_doDerotate)
   {
      std::cerr << "derotating\n";
      this->derotate();
   }
   
   
   if(this->m_combineMethod > 0)
   {
      std::cerr << "combining\n";
      this->combineFinim();
      
   }
   
   if(this->m_doWriteFinim == true || this->m_doOutputPSFSub == true)
   {
      std::cerr << "writing\n";
      
      fits::fitsHeader head;
      
      this->ADIobservation<_realT, _derotFunctObj>::stdFitsHeader(&head);
      
      head.append("", fits::fitsCommentType(), "----------------------------------------");
      head.append("", fits::fitsCommentType(), "mx::KLIPreduction parameters:");
      head.append("", fits::fitsCommentType(), "----------------------------------------");
   
      
      head.append("MEANSUBM", HCI::meansubMethodStr(m_meanSubMethod), "PCA mean subtraction method");
      
      
      std::stringstream str;
      
      if(m_Nmodes.size() > 0)
      {
         for(size_t nm=0;nm < m_Nmodes.size()-1; ++nm) str << m_Nmodes[nm] << ",";
         str << m_Nmodes[m_Nmodes.size()-1];      
         head.append<char *>("NMODES", (char *)str.str().c_str(), "number of modes");
      }
      
      if(m_minr.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < m_minr.size()-1; ++nm) str << m_minr[nm] << ",";
         str << m_minr[m_minr.size()-1];      
         head.append<char *>("REGMINR", (char *)str.str().c_str(), "region inner edge(s)");
      }
      
      if(m_maxr.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < m_maxr.size()-1; ++nm) str << m_maxr[nm] << ",";
         str << m_maxr[m_maxr.size()-1];      
         head.append<char *>("REGMAXR", (char *)str.str().c_str(), "region outer edge(s)");
      }
      
      if(m_minq.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < m_minq.size()-1; ++nm) str << m_minq[nm] << ",";
         str << m_minq[m_minq.size()-1];      
         head.append<char *>("REGMINQ", (char *)str.str().c_str(), "region minimum angle(s)");
      }
      
      if(m_maxq.size() > 0)
      {
         str.str("");
         for(size_t nm=0;nm < m_maxq.size()-1; ++nm) str << m_maxq[nm] << ",";
         str << m_maxq[m_maxq.size()-1];      
         head.append<char *>("REGMAXQ", (char *)str.str().c_str(), "region maximum angle(s)");
      }
      
      head.append<std::string>("EXMTHDMN", HCI::excludeMethodStr(m_excludeMethod), "exclusion method (min)");
      head.append<realT>("MINDPX", m_minDPx, "minimum delta (units based on EXMTHDMN)");
      
      head.append<std::string>("EXMTHDMX", HCI::excludeMethodStr(m_excludeMethodMax), "exclusion method (max)");
      head.append<realT>("MAXDPX", m_maxDPx, "maximum delta (units based on EXMTHDMX)");
      
      head.append<std::string>("INMTHDMX", HCI::includeMethodStr(m_includeMethod), "inclusion method");
      head.append<int>("INCLREFN", m_includeRefNum, "number of images included by INMTHDMX");

      if(this->m_doWriteFinim == true && this->m_combineMethod > 0)
      {
         this->writeFinim(&head);
      }
      
      if(this->m_doOutputPSFSub)
      {
         this->outputPSFSub(&head);
      }
   }
   
   return 0;
}
   
template<typename _realT, class _derotFunctObj, typename _evCalcT>
int KLIPreduction<_realT, _derotFunctObj, _evCalcT>::processPSFSub( const std::string & dir,
                                                                    const std::string & prefix,
                                                                    const std::string & ext
                                                                  )

{   
   std::cerr << "Beginning PSF Subtracted Image Processing\n";

   //Load first file to condigure based on its header.
   std::vector<std::string> flist = ioutils::getFileNames(dir, prefix, "000", ext);
   
   fits::fitsHeader fh;
   eigenImage<realT> im;
   fits::fitsFile<realT> ff;
   
   ff.read(im, fh, flist[0]);
   
   if(fh.count("MEANSUBM") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "MEANSUBM not found in FITS header.");
      return -1;
   }
   m_meanSubMethod = HCI::meansubMethodFmStr(fh["MEANSUBM"].String());
   std::cerr << "meanSubMethod: " << HCI::meansubMethodStr(m_meanSubMethod) << "\n";

   if(fh.count("NMODES") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "NMODES not found in FITS header.");
      return -1;
   }
   ioutils::parseStringVector(m_Nmodes, fh["NMODES"].String(), ",");
   if(m_Nmodes.size() == 0)
   {
      mxError("KLIPReduction", MXE_PARSEERR, "NMODES vector did not parse correctly.");
      return -1;
   }
   std::cerr << "nModes: " << fh["NMODES"].String() << "\n";
   
   /* -- REGMINR -- */
   if(fh.count("REGMINR") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "REGMINR not found in FITS header.");
      return -1;
   }
   ioutils::parseStringVector(m_minr, fh["REGMINR"].String(), ",");
   if(m_minr.size() == 0)
   {
      mxError("KLIPReduction", MXE_PARSEERR, "REGMINR vector did not parse correctly.");
      return -1;
   }
   std::cerr << "minr: " << fh["REGMINR"].String() << "\n";
   
   /* -- REGMAXR -- */
   if(fh.count("REGMAXR") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "REGMAXR not found in FITS header.");
      return -1;
   }
   ioutils::parseStringVector(m_maxr, fh["REGMAXR"].String(), ",");
   if(m_maxr.size() == 0)
   {
      mxError("KLIPReduction", MXE_PARSEERR, "REGMAXR vector did not parse correctly.");
      return -1;
   }
   std::cerr << "minr: " << fh["REGMAXR"].String() << "\n";
   
   
   /* -- REGMINQ -- */
   if(fh.count("REGMINQ") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "REGMINQ not found in FITS header.");
      return -1;
   }
   ioutils::parseStringVector(m_minq, fh["REGMINQ"].String(), ",");
   if(m_minq.size() == 0)
   {
      mxError("KLIPReduction", MXE_PARSEERR, "REGMINQ vector did not parse correctly.");
      return -1;
   }
   std::cerr << "minq: " << fh["REGMINQ"].String() << "\n";
   
   /* -- REGMAXR -- */
   if(fh.count("REGMAXR") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "REGMAXR not found in FITS header.");
      return -1;
   }
   ioutils::parseStringVector(m_maxq, fh["REGMAXR"].String(), ",");
   if(m_maxq.size() == 0)
   {
      mxError("KLIPReduction", MXE_PARSEERR, "REGMAXR vector did not parse correctly.");
      return -1;
   }
   std::cerr << "minr: " << fh["REGMAXR"].String() << "\n";
   
   if(fh.count("EXMTHDMN") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "EXMTHDMN not found in FITS header.");
      return -1;
   }
   m_excludeMethod = HCI::excludeMethodFmStr(fh["EXMTHDMN"].String());
   std::cerr << "excludeMethod: " << HCI::excludeMethodStr(m_excludeMethod) << "\n";
   
   if(fh.count("MINDPX") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "MINDPX not found in FITS header.");
      return -1;
   }
   m_minDPx = fh["MINDPX"].value<realT>();
   std::cerr << "minDPx: " << m_minDPx << "\n";
   
   if(fh.count("EXMTHDMX") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "EXMTHDMX not found in FITS header.");
      return -1;
   }
   m_excludeMethodMax = HCI::excludeMethodFmStr(fh["EXMTHDMX"].String());
   std::cerr << "excludeMethodMax: " << HCI::excludeMethodStr(m_excludeMethodMax) << "\n";
   
   if(fh.count("MAXDPX") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "MAXDPX not found in FITS header.");
      return -1;
   }
   m_maxDPx = fh["MAXDPX"].value<realT>();
   std::cerr << "maxDPx: " << m_maxDPx << "\n";

   if(fh.count("INMTHDMX") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "INMTHDMX not found in FITS header.");
      return -1;
   }
   m_includeMethod = HCI::includeMethodFmStr(fh["INMTHDMX"].String());
   std::cerr << "includeMethod: " << HCI::includeMethodStr(m_includeMethod) << "\n";
   
   if(fh.count("INCLREFN") == 0)
   {
      mxError("KLIPReduction", MXE_PARAMNOTSET, "INCLREFN not found in FITS header.");
      return -1;
   }
   m_includeRefNum = fh["INCLREFN"].value<int>();
   std::cerr << "includedRefNum: " << m_includeRefNum << "\n";


   
   this->m_skipPreProcess = true;

   this->m_keywords.clear();
   
   
   this->readPSFSub(dir, prefix, ext, m_Nmodes.size());
   
   
   finalProcess();
   
   return 0;
}

///@}


template<typename realT> class ADIDerotator;

extern template struct KLIPreduction<float, ADIDerotator<float>, float>;
extern template struct KLIPreduction<float, ADIDerotator<float>, double>;
extern template struct KLIPreduction<double, ADIDerotator<double>, double>;

} //namespace improc
} //namespace mx

#endif //__KLIPreduction_hpp__
