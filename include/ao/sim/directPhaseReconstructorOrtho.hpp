#ifndef directPhaseReconstructor_hpp
#define directPhaseReconstructor_hpp


#include "../../improc/eigenImage.hpp"
#include "../../improc/eigenCube.hpp"
#include "../../improc/fitsFile.hpp"
using namespace mx::improc;

#include "../../sigproc/signalWindows.hpp"

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


struct directPhaseReconstructorSpec
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
 
   typedef directPhaseReconstructorSpec specT;
      
protected:
    Eigen::Array<realT,-1,-1> _recon; ///< The reconstructor matrix.
    
   int _nModes; ///<The number of modes to be reconstructed
   
   int _detRows; ///<The size of the WFS image, in rows
   int _detCols;///<The size of the WFS image, in columns
   
   int _measurementSize; ///<The number of values in the measurement
      
   realT _calAmp; ///<The calibration amplitude used for response matrix acquisition
   
   imageT _rMat; ///<The response matrix
      
   eigenCube<realT> _rImages;
   
   //The mirror modes
   improc::eigenCube<realT> *_modes;
   
   
   
   
   imageT * _pupil;

   
   
   realT norm;
   
public:
   
   //The orthogonalized basis
   improc::eigenCube<realT> ortho;
   
   //The orthogonal spectrum 
   improc::eigenImage<realT> spectrum;
   
   std::vector<realT> b;
   
   imageT _mask;
   
   improc::eigenImage<double> _spectrum;
   
   imageT * _gains; 
   
   int _npix;
   
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
   
      _mask = *_pupil;
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
   void loadRecon(const std::string & fname);
   
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
   _nModes = 0;
   
   norm = 0;
   
   _gains =0;
   
   _npix = 0;
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
   
   _mask = *_pupil;
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
void directPhaseReconstructor<realT>::loadRecon(const std::string & fname)
{
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

}
     
template<typename realT> 
template<typename measurementT, typename wfsImageT>
void directPhaseReconstructor<realT>::reconstruct(measurementT & commandVect, wfsImageT & wfsImage)
{
   
   BREAD_CRUMB;
   
   if(_npix == 0)
   {
      _npix = _pupil->sum();
      
   }
   
   BREAD_CRUMB;   

   wfsImage.image *= _mask;
   
   BREAD_CRUMB;
   

   
   b.resize(_nModes);
   
   #pragma omp parallel for
   for(int j=0; j< _nModes; ++j)
   {
      b[j] = (wfsImage.image*ortho.image(j)).sum()/ _npix;
   }

   //commandVect.measurement.setZero();
      
   //#pragma omp parallel for
   for(int j=0; j< modes.planes(); ++j)
   {
      realT amp = b[0]*spectrum(0,j);
      
      for(int i=1;i<modes.planes();++i) 
      {
         amp += b[i]*spectrum(i,j);
      }
      commandVect.measurement(0,j) = amp;
   }

   //std::cerr << "--- " << b[0] << " " << b[0]*spectrum(0,0) << " " << commandVect.measurement(0,0) << "\n";
   for(int k=0; k< _nModes; ++k)
   {
      std::cerr << commandVect.measurement(0,k) << " ";
   }
   std::cerr << "\n";
   
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
   int l = 0;
   for(int j=0; j<measureVec.measurement.rows(); ++j)
   {
      for(int k=0; k<measureVec.measurement.cols(); ++k)
      {
         _rMat(l, i) = measureVec.measurement(j,k);
         ++l;
      }
   }
   
  //_rMat.col(i) = measureVec.measurement.row(0);
}
  
template<typename realT>   
template<typename measurementT,typename wfsImageT>
void directPhaseReconstructor<realT>::accumulateRMat(int i, measurementT &measureVec, wfsImageT & wfsImage)
{
   accumulateRMat(i, measureVec);
   _rImages.image(i) = wfsImage.image;   
}

template<typename realT> 
void directPhaseReconstructor<realT>::saveRMat(std::string fname)
{
   fitsFile<realT> ff;
   fitsHeader head;
   
   head.append("DETROWS", _detRows, "WFS detector rows");
   head.append("DETCOLS", _detCols, "WFS detector cols");
   head.append("CALAMP", _calAmp, "DM Calibration amplitude");
   head.append("NMODES", _nModes, "Number of modes included in the response matrix.");
   
   ff.write(fname, _rMat, head);
}

template<typename realT> 
void directPhaseReconstructor<realT>::saveRImages(std::string fname)
{
   fitsFile<realT> ff;
   fitsHeader head;
   
   
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

#endif //directPhaseReconstructor_hpp


//--- Code for mean and tip/tilt subtraction in reconstruct:
//    realT mean = (wfsImage.image * *_pupil).sum()/npix;
//    
//    wfsImage.image -= mean;
//    wfsImage.image *= *_pupil;

#if 0 
   
   BREAD_CRUMB;
   if(norm == 0)
   {
      for(int i = 0; i < wfsImage.rows(); ++i)
      {
         for(int j=0; j < wfsImage.cols(); ++j)
         {
            norm += pow((i - 0.5*(wfsImage.rows()-1.0))*(*_pupil)(i,j),2);
         }
      }
    BREAD_CRUMB;  
      norm = sqrt(norm);
      
      //std::cout << "NORM: === " << norm << "\n";
   }
   BREAD_CRUMB;

   //std::cout << "NORM: === " << norm << "\n";
   
   //Remove tip/tilt:
   realT ampx = 0, ampy = 0;
   BREAD_CRUMB;   
   for(int i = 0; i < wfsImage.rows(); ++i)
   {
      for(int j=0; j < wfsImage.cols(); ++j)
      {
         ampx += wfsImage(i,j) * (i - 0.5*(wfsImage.rows()-1.0))/norm;
         ampy += wfsImage(i,j) * (j - 0.5*(wfsImage.cols()-1.0))/norm;
      }
   }
   BREAD_CRUMB;
   for(int i = 0; i < wfsImage.rows(); ++i)
   {
      for(int j=0; j < wfsImage.cols(); ++j)
      {
         wfsImage(i,j) -= ampx * (i - 0.5*(wfsImage.rows()-1.0))/norm * (*_pupil)(i,j);
         wfsImage(i,j) -= ampy * (j - 0.5*(wfsImage.cols()-1.0))/norm * (*_pupil)(i,j);
      }
   }
#endif

