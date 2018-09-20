#ifndef __pywfsSlopeReconstructor_hpp__
#define __pywfsSlopeReconstructor_hpp__

#include <mx/improc/eigenCube.hpp>


namespace mx
{
namespace AO
{
namespace sim
{
 
struct pywfsSlopeReconstructorSpec
{
   std::string dmName; 
   std::string basisName; 
   
   std::string rMatId;
};

///A Pyramid Wavefront Sensor slope reconstructor.
/** Calculates slopes, normalized by total flux in the image.
  */ 
template<typename _floatT> 
class pywfsSlopeReconstructor
{
public:
   
   typedef _floatT floatT;
   
   ///The type of the measurement (i.e. the slope vector)
   //typedef Eigen::Array<floatT,-1,-1> measurementT;
   
   ///The type of the WFS image
   typedef Eigen::Array<floatT, -1, -1> imageT;
   
   ///The type of the response matrix
   typedef Eigen::Array<floatT, -1, -1> rmatT;
 
   typedef pywfsSlopeReconstructorSpec specT;
   
protected:
   Eigen::Array<floatT,-1,-1> _recon; ///< The reconstructor matrix.
   
   int _maskType; ///0 is centrally obscured circle, 1 is supplied by fits files
   
   floatT _maskRadius; ///<The radius of the quadrant mask
   floatT _maskObscuration; ///<The central obscuration of the quadrant mask

   std::string _maskFile; ///<The name of the quadrant mask file
   
   bool _maskMade; ///<Whether or not the mask has been made

   Eigen::Array<floatT, -1,-1> _quadMask; ///<The quadrant mask
   void calcMask(); ///<Calculates the quadrant mask
   int _measurementSize; ///<The number of slopes in the measurement
 
      
   floatT _calAmp; ///<The calibration amplitude used for response matrix acquisition
   
   int _nModes; ///<The number of modes to be reconstructed
   
   int _detRows; ///<The size of the WFS image, in rows
   int _detCols;///<The size of the WFS image, in columns
   
   
   imageT _rMat; ///<The response matrix
   improc::eigenCube<floatT> _rImages;
   
public:   
   
   int _binFact; ///<The binning to apply before reconstructing.
   
   ///Default c'tor
   pywfsSlopeReconstructor();
   
   template<typename AOSysT>
   void initialize(AOSysT & AOSys, specT & spec)
   {
      std::string recMatrix = mx::AO::path::sys::cal::iMat( AOSys._sysName, 
                                                            spec.dmName, 
                                                            AOSys._wfsName, 
                                                            AOSys._pupilName, 
                                                            spec.basisName, 
                                                            spec.rMatId)  ;
      loadRecon(recMatrix);
      
      
   }


   ///Get the quadrant mask radius (_maskRadius)
   floatT maskRadius();
   
   ///Set the quadrant mask radius (_maskRadius)
   /** Calling this will cause the quadrant mask to be recalculated next time it is needed.
     * 
     * \param mr [in] the new mask radius
     */ 
   void maskRadius(floatT mr);
   
   ///Get the quadrant mask central obscuration ratio (_maskObscuration)
   floatT maskObscuration();
   
   ///Set the quadrant mask central obscuration ratio (_maskObscuration)
   /** Calling this will cause the quadrant mask to be recalculated next time it is needed.
     * 
     * \param mo [in] the new central obscuration ratio 
     */
   void maskObscuration(floatT mo);
   
   ///Get the quadrant mask file name (_maskFile)
   std::string maskFile();
   
   ///Set the quadrant mask file name (_maskFile)
   /** Calling this will cause the quadrant mask to be reloaded next time it is needed.
     * 
     * \param mf [in] the new mask file
     */
   void maskFile(const std::string & mf);
   
   ///Get the calibration amplitude used in response matrix acquisition (_calAmp)
   floatT calAmp();
   
   ///Set the calibration amplitude used in response matrix acquisition (_calAmp)
   /**
     * \param ca [in] the new calibration amplitude
     */ 
   void calAmp(floatT ca);

   ///Get the number of modes (_nModes)
   int nModes();
   
   ///Get the number of detector rows (_detRows)
   int detRows();
   
   ///Set the number of detector rows (_detRows)
   void detRows(int dr);

   ///Get the number of detector columns (_detCols)   
   int detCols();
   
   ///Set the number of detector columns (_detCols)
   void detCols(int dc);

   ///Load the reconstrutor from the specified FITS file 
   /** 
     * \param fname is the name of the FITS file, including path
     */   
   void loadRecon(std::string fname);
   
   ///Return the size of the unbinned measurement
   int measurementSize();
   
   ///Calculate the slope measurement
   /**
     * \param slopes [out] a (_measurementSize X 2)  array of slopes
     * \param wfsImage [in] the WFS image from which to measure the slopes
     */   
   template<typename measurementT, typename wfsImageT>
   void calcMeasurement(measurementT & slopes, wfsImageT & wfsImage);
         
   ///Reconstruct the wavefront from the input image, producing the modal amplitude vector
   template<typename measurementT, typename wfsImageT>
   void reconstruct(measurementT & commandVect, wfsImageT & wfsImage);

   ///Initialize the response matrix for acquisition
   /** 
     * \param nmodes the number of modes 
     * \param calamp the calibration amplitude
     * \param detrows the number of detector rows
     * \param detcols the number of detector columns
     */ 
   void initializeRMat(int nmodes, floatT calamp, int detrows,int detcols);
   
   ///Accumalte the next measurement in the response matrix
   /** 
     * \param i the measurement index
     * \param measureVec is the i-th measurement vector
     */ 
   template<typename measurementT>
   void accumulateRMat(int i, measurementT &measureVec);
   
   template<typename measurementT, typename wfsImageT>
   void accumulateRMat(int i, measurementT &measureVec, wfsImageT & wfsImage);
   
   ///Write the accumulated response matrix to disk
   /**
     * \param fname the name, including path, of the response matrix
     */ 
   void saveRMat(std::string fname);
      
   void saveRImages(std::string fname);
      
   
   
};


template<typename floatT> 
pywfsSlopeReconstructor<floatT>::pywfsSlopeReconstructor()
{
   _maskRadius = 0;
   _maskObscuration = 0;
   
   _maskMade = false;
   
   _binFact = 1;
}

template<typename floatT> 
floatT pywfsSlopeReconstructor<floatT>::maskRadius()
{
   return _maskRadius;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::maskRadius(floatT mr)
{
   _maskRadius = mr;
   _maskType = 0;
   _maskMade = false;
   
}
   
template<typename floatT> 
floatT pywfsSlopeReconstructor<floatT>::maskObscuration()
{
   return _maskObscuration;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::maskObscuration(floatT mo)
{
   _maskObscuration = mo;
   _maskType = 0;
   _maskMade = false;
}
 
template<typename floatT> 
std::string pywfsSlopeReconstructor<floatT>::maskFile()
{
   return _maskFile;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::maskFile(const std::string & mf)
{
   _maskFile = mf;
   _maskType = 1;
   _maskMade = false;
}

template<typename floatT> 
floatT pywfsSlopeReconstructor<floatT>::calAmp()
{
   return _calAmp;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::calAmp(floatT ca)
{
   _calAmp = ca;
}

template<typename floatT> 
int pywfsSlopeReconstructor<floatT>::nModes()
{
   return _nModes;
}

template<typename floatT> 
int pywfsSlopeReconstructor<floatT>::detRows()
{
   return _detRows;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::detRows(int dr)
{
   _detRows = dr;
   _maskMade = false;
}

template<typename floatT> 
int pywfsSlopeReconstructor<floatT>::detCols()
{
   return _detCols;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::detCols(int dc)
{
   _detCols = dc;
   _maskMade = false;
}




template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::loadRecon(std::string fname)
{
   improc::fitsFile<floatT> ff;
   improc::fitsHeader head;
   
   ff.read(fname, _recon, head);

   _maskRadius = head["MASKRAD"].Value<floatT>();
   _maskObscuration = head["MASKOBS"].Value<floatT>();
   _maskFile = head["MASKFILE"].String();
   
   if( _maskFile != "") _maskType = 1;
   
   _calAmp = head["CALAMP"].Value<floatT>();

   _nModes = head["NMODES"].Value<int>();

   _detRows = head["DETROWS"].Int();
   _detCols = head["DETCOLS"].Int();
   
   _maskMade = false;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::calcMask()
{
   if(_maskType == 1)
   {
      improc::fitsFile<floatT> ff;
            
      std::cerr << "Loading Mask: " << _maskFile << "\n";
      ff.read(_maskFile, _quadMask);
      

   }
   {
      _quadMask.resize(0.5*_detRows/_binFact,0.5*_detCols/_binFact);
      wfp::circularPupil( _quadMask, _maskObscuration,_maskRadius/_binFact);
   }
      
   _measurementSize = 2* _quadMask.sum();// + 64*64;
   
   _maskMade = true;
   
    improc::fitsFile<floatT> ff;
    ff.write("quadMask.fits", _quadMask);
}

template<typename floatT> 
int pywfsSlopeReconstructor<floatT>::measurementSize()
{
   if(!_maskMade) calcMask();
   return _measurementSize;

}

template<typename floatT> 
template<typename measurementT, typename wfsImageT>
void pywfsSlopeReconstructor<floatT>::calcMeasurement(measurementT & slopes, wfsImageT & wfsImage)
{
   if(!_maskMade) calcMask();
   
   imageT _wfsImage;
   
   if(_binFact > 1)
   {
      std::cout << "rebinning" << "\n";
      _wfsImage.resize( wfsImage.image.rows()/_binFact, wfsImage.image.cols()/_binFact);
      imageDownSample( _wfsImage, wfsImage.image);
   }
   else
   {
      _wfsImage = wfsImage.image;
   }

   int nsz = _wfsImage.rows();
   
   int nPix = _quadMask.sum();
 
   std::vector<int> x(nPix), y(nPix);

   int k = 0;
   for(int i=0;i<_quadMask.rows(); ++i)
   {
      for(int j=0;j<_quadMask.rows(); ++j)
      {
         if(_quadMask(i,j) == 1)
         {
            x[k] = i;
            y[k] = j;
            ++k;
         }
      }
   }

   slopes.measurement.resize(1, 2.*nPix);// + 64*64); //wfsImage.tipImage.rows()*wfsImage.tipImage.cols());
   
   floatT I0, I1, I2, I3;
   
   floatT norm = _wfsImage.sum();
   
   for(int i=0;i<nPix; ++i)
   {
      I0 = _wfsImage(x[i], y[i]);
      I1 = _wfsImage(x[i] + 0.5*nsz, y[i]);
      I2 = _wfsImage(x[i], y[i]+0.5*nsz);
      I3 = _wfsImage(x[i]+0.5*nsz, y[i]+0.5*nsz);
 
      if(norm == 0)
      {
         slopes.measurement(0,i) =0;;
         slopes.measurement(0,i+nPix) = 0;//
      }
      else
      {
         slopes.measurement(0,i)          = (I0 + I1 - I2 - I3)/norm;//(I0+I1+I2+I3);
         slopes.measurement(0,i+nPix) = (I0 + I2 - I1 - I3)/norm; //(I0+I1+I2+I3);
      }
   }
   
   
   /*
   for(int i=0; i< wfsImage.tipImage.rows(); ++i)
   {
      for(int j=0;j<wfsImage.tipImage.cols(); ++j)
      {
         slopes.measurement(0, 2*nPix + i*wfsImage.tipImage.cols() + j) = wfsImage.tipImage(i,j);
      }
   }*/
   
}
     
template<typename floatT>
template<typename measurementT, typename wfsImageT>
void pywfsSlopeReconstructor<floatT>::reconstruct(measurementT & commandVect, wfsImageT & wfsImage)
{
   measurementT slopes;
   
   calcMeasurement(slopes, wfsImage);
   
   commandVect.measurement = slopes.measurement.matrix()*_recon.matrix();
   
   commandVect.iterNo = wfsImage.iterNo;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::initializeRMat(int nModes, floatT calamp, int detRows, int detCols)
{
   
   _nModes = nModes;
   
   calAmp(calamp);
   
   _detRows = detRows;
   _detCols = detCols;
   
   _rMat.resize(measurementSize(), nModes);
   _rMat.setZero();
   
   _rImages.resize(_detRows, _detCols, _nModes);
   
}

template<typename floatT>   
template<typename measurementT>
void pywfsSlopeReconstructor<floatT>::accumulateRMat(int i, measurementT &measureVec)
{
   _rMat.col(i) = measureVec.measurement.row(0);
}

template<typename floatT>   
template<typename measurementT, typename wfsImageT>
void pywfsSlopeReconstructor<floatT>::accumulateRMat(int i, measurementT &measureVec, wfsImageT & wfsImage)
{
   accumulateRMat(i, measureVec);
   
   _rImages.image(i) = wfsImage.image;
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::saveRMat(std::string fname)
{
   improc::fitsFile<floatT> ff;
   improc::fitsHeader head;
   
   if(_maskType == 1)
   {
      head.append("MASKFILE", maskFile(), "Name of mask file");
   }
   else
   {
      head.append("MASKRAD", maskRadius(), "Mask radius, in pixels");
      head.append("MASKOBS", maskObscuration(), "Mask fractional central obscuration");
   }
   
   head.append("DETROWS", _detRows, "WFS detector rows");
   head.append("DETCOLS", _detCols, "WFS detector cols");
   head.append("CALAMP", _calAmp, "DM Calibration amplitude");
   head.append("NMODES", _nModes, "Number of modes included in the response matrix.");
   
   //ff.write(fname, _rMat.data(), _rMat.rows(), _rMat.cols(), 1, &head);
   ff.write(fname, _rMat, head);
}

template<typename floatT> 
void pywfsSlopeReconstructor<floatT>::saveRImages(std::string fname)
{
   improc::fitsFile<floatT> ff;
   improc::fitsHeader head;
   
   if(_maskType == 1)
   {
      head.append("MASKFILE", maskFile(), "Name of mask file");
   }
   else
   {
      head.append("MASKRAD", maskRadius(), "Mask radius, in pixels");
      head.append("MASKOBS", maskObscuration(), "Mask fractional central obscuration");
   }
   
   head.append("DETROWS", _detRows, "WFS detector rows");
   head.append("DETCOLS", _detCols, "WFS detector cols");
   head.append("CALAMP", _calAmp, "DM Calibration amplitude");
   head.append("NMODES", _nModes, "Number of modes included in the response matrix.");
   
   //ff.write(fname, _rImages.data(), _rImages.rows(), _rImages.cols(), _rImages.planes(), &head);
   ff.write(fname, _rImages, head);
}

} //namespace sim 
} //namespace AO
} //namespace mx

#endif //__pywfsSlopeReconstructor_hpp__

