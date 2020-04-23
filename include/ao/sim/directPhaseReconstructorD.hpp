#ifndef directPhaseReconstructor_hpp
#define directPhaseReconstructor_hpp

#pragma GCC system_header
#include <Eigen/Dense>

#include <mx/signalWindows.hpp>

#ifdef DEBUG
#define BREAD_CRUMB std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n"; 
#else
#define BREAD_CRUMB 
#endif

namespace mx
{

   
   
namespace AO
{
   
namespace sim
{

struct slopeReconstructorSpec
{
   std::string dmName; 
   std::string basisName; 
   
   std::string rMatId;
};
   
///A Pyramid Wavefront Sensor slope reconstructor.
/** Calculates slopes, normalized by total flux in the image.
  */ 
template<typename realT> 
class directPhaseReconstructor
{
public:
   
   ///The type of the measurement (i.e. the slope vector)
   //typedef Eigen::Array<realT,-1,-1> measurementT;
   
   ///The type of the WFS image
   typedef Eigen::Array<realT, -1, -1> imageT;
   
   ///The type of the response matrix
   typedef Eigen::Array<realT, -1, -1> rmatT;
 
   typedef slopeReconstructorSpec specT;
      
protected:
   Eigen::Array<realT,-1,-1> _recon; ///< The reconstructor matrix.
   
   imageT _rMat; ///<The response matrix
   eigenCube<realT> _rImages;
   
   int _nModes {0}; ///<The number of modes to be reconstructed
   
   int _detRows {0}; ///<The size of the WFS image, in rows
   int _detCols {0};///<The size of the WFS image, in columns
   
   int _measurementSize {0}; ///<The number of values in the measurement
      
      
   //The mirror modes
   mx::eigenCube<realT> *_modes {nullptr};
   
   imageT * _pupil {nullptr};
   
   
   imageT _mask;
   
   realT norm {0};
   
   ds9_interface ds9i;
      
public:
   mx::eigenImaged _spectrum;
   
   imageT * _gains {nullptr}; 
   
public:   
   ///Default c'tor
   directPhaseReconstructor();

   template<typename AOSysT>
   void initialize(AOSysT & AOSys, specT & spec)
   {
         
      _modes = &AOSys.dm._infF;
   
      _nModes = _modes->planes();
      _detRows = _modes->rows();
      _detCols = _modes->cols();
   
      _measurementSize = _detRows*_detCols;
   
      _pupil = &AOSys._pupil;
      _measurementSize = _pupil->sum();
   
      _mask.resize(_pupil->rows(), _pupil->cols());
   
//          tukey2d(realT *filt, int dim, realT N, realT alpha, realT xc, realT yc)
      //mx::tukey2d(_mask.data(), _mask.rows(), (realT) _mask.rows(), (realT) 0.0, (realT) 0.5*(_mask.rows()-1), (realT) 0.5*(_mask.cols()-1));
      _mask = *_pupil;
      
      mx::fitsFile<realT> ff;
      ff.write("dprMask.fits", _mask);
      
      std::string recMatrix = mx::AO::path::sys::cal::iMat( AOSys._sysName, 
                                                            spec.dmName, 
                                                            AOSys._wfsName, 
                                                            AOSys._pupilName, 
                                                            spec.basisName, 
                                                            spec.rMatId)  ;
      loadRecon(recMatrix);
   
   }
   
   template<typename AOSysT>
   void linkSystem(AOSysT & AOSys);
   
   ///Get the calibration amplitude used in response matrix acquisition (_calAmp)
   realT calAmp();
   
   
   ///Set the calibration amplitude used in response matrix acquisition (_calAmp)
   /**
     * \param ca [in] the new calibration amplitude
     */ 
   void calAmp(realT ca);
   
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
   void initializeRMat(int nmodes, realT calamp, int detrows,int detcols);
   
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


template<typename realT> 
directPhaseReconstructor<realT>::directPhaseReconstructor()
{
   
}

template<typename realT> 
template<typename AOSysT>
void directPhaseReconstructor<realT>::linkSystem(AOSysT & AOSys)
{
   _modes = &AOSys.dm._modes;
   
   _nModes = _modes->planes();
   _detRows = _modes->rows();
   _detCols = _modes->cols();
   
   _measurementSize = _detRows*_detCols;
   
   _pupil = &AOSys._pupil;
   _measurementSize = _pupil->sum();
   
   _mask.resize(_pupil->rows(), _pupil->cols());
   mx::tukey2d(_mask.data(), _mask.rows(), (realT) _mask.rows(), 0.0, 0.5*(_mask.rows()-1), 0.5*(_mask.cols()-1));
   
   mx::fitsFile<realT> ff;
   ff.write("dprMask.fits", _mask);
   
}

template<typename realT> 
realT directPhaseReconstructor<realT>::calAmp()
{
   return 0.5*780.0e-9/two_pi<realT>();
}

template<typename realT> 
void directPhaseReconstructor<realT>::calAmp(realT ca)
{
   return;
}

template<typename realT> 
int directPhaseReconstructor<realT>::nModes()
{
   return _nModes;
}

template<typename realT> 
int directPhaseReconstructor<realT>::detRows()
{
   return _detRows;
}

template<typename realT> 
int directPhaseReconstructor<realT>::detCols()
{
   return _detCols;
}

template<typename realT> 
void directPhaseReconstructor<realT>::loadRecon(std::string fname)
{
   mx::fitsFile<realT> ff;
   mx::fitsHeader head;
   
   ff.read(fname, _recon, head);

   
}

template<typename realT> 
int directPhaseReconstructor<realT>::measurementSize()
{
   return _measurementSize;
}

template<typename realT>
template<typename measurementT, typename wfsImageT>
void directPhaseReconstructor<realT>::calcMeasurement(measurementT & slopes, wfsImageT & wfsImage)
{
   
   slopes.measurement.resize(1, _measurementSize );
   
   int k = 0;
   
   for(int i=0; i< wfsImage.image.rows(); ++i)
   {
      for(int j=0; j < wfsImage.image.cols(); ++j)
      {
         if( (*_pupil)(i,j) == 0) continue;
         
         slopes.measurement(0, k) = wfsImage.image(i,j);
         ++k;
      }
   }
   
}
     
template<typename realT> 
template<typename measurementT, typename wfsImageT>
void directPhaseReconstructor<realT>::reconstruct(measurementT & commandVect, wfsImageT & wfsImage)
{
   
   BREAD_CRUMB;
   
   measurementT slopes;
   
   calcMeasurement(slopes, wfsImage);
   
   BREAD_CRUMB;
   commandVect.measurement = slopes.measurement.matrix()*_recon.matrix();

   BREAD_CRUMB;
   commandVect.iterNo = wfsImage.iterNo;
}

template<typename realT> 
void directPhaseReconstructor<realT>::initializeRMat(int nModes, realT calamp, int detRows, int detCols)
{
   _nModes = nModes;
   
   _detRows = detRows;
   _detCols = detCols;
   
   _rMat.resize(measurementSize(), nModes);
   _rMat.setZero();
   
   _rImages.resize(_detRows, _detCols, _nModes);
}

template<typename realT> 
template<typename measurementT>
void directPhaseReconstructor<realT>::accumulateRMat(int i, measurementT &measureVec)
{
   _rMat.col(i) = measureVec.measurement.row(0);
}
  
template<typename realT> 
template<typename measurementT, typename wfsImageT>
void directPhaseReconstructor<realT>::accumulateRMat(int i, measurementT &measureVec, wfsImageT & wfsImage)
{
   BREAD_CRUMB;
   accumulateRMat(i, measureVec);
   
   BREAD_CRUMB;
   
   _rImages.image(i) = wfsImage.image;
}

template<typename realT> 
void directPhaseReconstructor<realT>::saveRMat(std::string fname)
{
   mx::fitsFile<realT> ff;
   mx::fitsHeader head;
   
   ff.write(fname, _rMat);
   
}

template<typename realT> 
void directPhaseReconstructor<realT>::saveRImages(std::string fname)
{
   mx::fitsFile<realT> ff;
   
   ff.write(fname, _rImages);
}


} //namespace sim
} //namespace AO
} //namespace mx

#endif //directPhaseReconstructor_hpp

