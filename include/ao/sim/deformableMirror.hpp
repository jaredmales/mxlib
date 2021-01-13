/** \file
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief
  * \ingroup mxAO_sim_files
  *
  */

#ifndef deformableMirror_hpp
#define deformableMirror_hpp

#include <cmath>

#include "../../wfp/imagingUtils.hpp"
#include "../../sys/timeUtils.hpp"
#include "../../ioutils/fits/fitsFile.hpp"

#include "../../ioutils/readColumns.hpp"

#include "../../math/constants.hpp"

#include "../../math/cuda/cudaPtr.hpp"
#include "../../math/cuda/templateCublas.hpp"

#include "../aoPaths.hpp"
#include "wavefront.hpp"




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

struct deformableMirrorSpec
{
   std::string name;
   std::string basisName;
};

template<typename _realT>
class deformableMirror
{
public:

   typedef _realT realT;

   typedef std::complex<realT> complexT;

   ///The wavefront data type
   typedef wavefront<realT>  wavefrontT;

   ///The pupil image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

   typedef deformableMirrorSpec specT;

public: //<-give these accesors and make protected

   std::string _name;

   std::string _basisName;

   std::string m_pupilName;

   ///Time to move from old shape to new shape, in loop time steps
   int _settleTime;

   int _settling;

   int _settleTime_counter;

   ///The amplitude used when measuring the response matrix of the current basis set.
   float _calAmp;


   imageT m_pupil;    ///The system pupil, possibly apodized, etc.
   realT m_pupilSum; ///The sum of the pupil mask, nominally the number of pixels.
   std::vector<size_t> m_idx; /// The offset coordinates of non-zero pixels in the pupil
   

public:

   //The modes-2-command matrix for the basis
   Eigen::Array<realT, -1, -1> m_m2c;

   //The mirror influence functions
   improc::eigenCube<realT> m_infF;
   
   #ifdef MXAO_USE_GPU
   cublasHandle_t *m_cublasHandle;
   cuda::cudaPtr<realT> m_one;
   cuda::cudaPtr<realT> m_alpha;
   cuda::cudaPtr<realT> m_zero;
   
   cuda::cudaPtr<realT> m_devM2c;
   cuda::cudaPtr<realT> m_devInfF;
   
   cuda::cudaPtr<realT> m_devModeCommands;
   cuda::cudaPtr<realT> m_devActCommands;
   cuda::cudaPtr<realT> m_devShape;
   #endif
   
   size_t m_nActs;
   size_t m_nRows;
   size_t m_nCols;

protected:
   //The current shape of the mirror
   imageT _shape;

   //The shape of the mirror when movement begins
   imageT _oldShape;

   //The next shape of the mirror after movement ends.
   imageT _nextShape;


public:

   ///Default c'tor.
   deformableMirror();

   ~deformableMirror()
   {
      if(_commandFileOpen)
      {
         _commandFout.close();
      }
   }

   template<typename AOSysT>
   void initialize( AOSysT & AOSys, 
                    specT & spec,
                    const std::string & pupil
                  );

   std::string name()
   {
      return _name;
   }

   std::string basisName()
   {
      return _basisName;
   }


   ///Get the settling time of the DM
   int settleTime();

   ///Set the settling time of the DM
   void settleTime(int st);

   ///Get the calibration amplitude.
   /** The modal commands are relative to this value.
     */
   realT calAmp();

   ///Set the calibration amplitude.
   void calAmp(realT ca);

   ///Get the number of modes in the  M2C.
   int nModes();

   ///Apply a single mode.
   void applyMode(wavefrontT & wf, int modeNo, realT amp, realT lambda);

   realT _settlingIter;
   realT _settledIter;

   template<typename commandT>
   void setShape(commandT & commandV);

   void applyShape( wavefrontT & wf,  
                    realT lambda
                  );

   double t0, t1, t_mm, t_sum;

   bool _writeCommands;
   bool _commandFileOpen;
   std::string _commandFile;
   std::ofstream _commandFout;

   realT _commandLimit;

   Eigen::Array<double,-1,-1> _pos, _map;

   sigproc::psdFilter<realT,2> m_filter;

   bool m_applyFilter {false};

   void setFilter( int width );
};


template<typename _realT>
deformableMirror<_realT>::deformableMirror()
{
   _settleTime = 1;
   _settling = 0;
   _settleTime_counter = 0;
   _calAmp = 1e-6;
   
   t_mm = 0;
   t_sum = 0;

   _writeCommands = false;
   _commandFileOpen = false;

   _commandLimit = 0;

   //readColumns("sigma.dat", sigma);

   //std::cerr << "sigma.size() = " << sigma.size() << "\n";
}

template<typename _realT>
template<typename AOSysT>
void deformableMirror<_realT>::initialize( AOSysT & AOSys, 
                                           specT & spec,
                                           const std::string & pupil )
{
   _name = spec.name;
   _basisName = spec.basisName;
   m_pupilName = pupil;

   fits::fitsFile<_realT> ff;

   std::string pName;
   pName = mx::AO::path::pupil::pupilFile(m_pupilName);
   ff.read(m_pupil, pName);
   m_pupilSum = m_pupil.sum();

   m_idx.clear();
   int kk = 0;
   for(int rr=0;rr<m_pupil.rows();++rr)
   {
      for(int cc=0;cc<m_pupil.cols();++cc)
      {
         if(m_pupil(rr,cc) == 1)
         {
            m_idx.push_back(kk);
         }
         ++kk;
      }
   }
   
   if(_name == "modalDM")
   {
      std::cerr << "modalDM\n";

      std::string ifName;
      ifName = mx::AO::path::basis::modes(_basisName);
      
      ff.read(m_infF, ifName);

      m_nActs = m_infF.planes();
      m_nRows = m_infF.rows();
      m_nCols = m_infF.cols();

      m_m2c.resize( m_nActs, m_nActs);
      m_m2c.setZero();

      for(int i=0;i<m_m2c.rows();++i) m_m2c(i,i) = 1.0;
      
      #ifdef MXAO_USE_GPU
      
      m_cublasHandle = &AOSys.m_cublasHandle;
      
      m_one.resize(1);
      realT one = 1.0;
      m_one.upload(&one, 1);
      
      m_alpha.resize(1);
      realT alpha = -1*_calAmp;
      m_alpha.upload(&alpha, 1);
      
      m_zero.resize(1);
      m_zero.initialize();
      
      mx::improc::eigenImage<realT> modft;
      
      modft.resize( m_idx.size(),m_nActs);
   
      for(int pp=0;pp<m_nActs;++pp)
      {
         for(int nn=0; nn < m_idx.size(); ++nn)
         {
            *(modft.data() + pp*m_idx.size() + nn) = *(m_infF.image(pp).data() + m_idx[nn]);
         }
      }
   
      m_devInfF.upload(modft.data(), modft.rows()*modft.cols());
      
      m_devM2c.upload(m_m2c.data(), m_m2c.rows()*m_m2c.cols());
      
      m_devModeCommands.resize(m_nActs);
      m_devActCommands.resize(m_nActs);
      
      #endif
      
   }
   else
   {
      std::cerr << "Non-modal DM is currently not implemented\n";
      exit(-1);
#if 0
      std::string ifName;
      ifName = mx::AO::path::dm::influenceFunctions(_name);


      improc::eigenCube<realT> infFLoad;
      ff.read(infFLoad, ifName);

      realT c = 0.5*(infFLoad.rows()-1);
      realT w = 0.5*(m_pupil.rows()-1);

      m_infF.resize( m_pupil.rows(), m_pupil.cols(), infFLoad.planes());

      for(int i=0;i<infFLoad.planes(); ++i)
      {
         m_infF.image(i) = infFLoad.image(i).block( c-w, c-w, m_pupil.rows(), m_pupil.rows());
      }

      std::string m2cName;

      m2cName = mx::AO::path::dm::M2c( _name, _basisName );

      ff.read(m_m2c, m2cName);


      std::string posName =  mx::AO::path::dm::actuatorPositions(_name, true);

      ff.read(_pos, posName);
#endif

   }

   _shape.resize(m_nRows, m_nCols);

   _shape.setZero();
   _nextShape = _shape;
   _oldShape = _shape;

   return;
}



template<typename _realT>
int deformableMirror<_realT>::settleTime()
{
   return settleTime;
}

template<typename _realT>
void deformableMirror<_realT>::settleTime(int st)
{
   if(st < 1)
   {
      std::cerr << "DM: settleTime must be > 1.  Correcting.\n";
      st = 1;
   }

   _settleTime = st;
}

template<typename _realT>
_realT deformableMirror<_realT>::calAmp()
{
   return _calAmp;
}

template<typename _realT>
void deformableMirror<_realT>::calAmp(realT ca)
{
   _calAmp = ca;
   
    #ifdef MXAO_USE_GPU
   realT alpha = -1*_calAmp;
   m_alpha.upload(&alpha, 1);
   #endif
      
}


template<typename _realT>
int deformableMirror<_realT>::nModes()
{
   return m_m2c.cols();
}


template<typename _realT>
void deformableMirror<_realT>::applyMode(wavefrontT & wf, int modeNo, realT amp, realT lambda)
{

   Eigen::Array<_realT,-1,-1> commandV(1, nModes());
   commandV.setZero();

   commandV(0, modeNo) = 1;

   Eigen::Array<_realT, -1, -1> c;

   t0 = sys::get_curr_time();
   c = m_m2c.matrix() * commandV.matrix().transpose();
   t1 = sys::get_curr_time();
   t_mm += t1-t0;

   imageT shape( m_nRows, m_nCols);


   shape = c(0,0)*m_infF.image(0);


   t0 = sys::get_curr_time();
   #pragma omp parallel
   {
      Eigen::Array<_realT, -1, -1> tmp;
      //tmp.resize(m_nRows, m_nCols);
      //tmp.setZero();
      tmp.setZero(m_nRows, m_nCols);

      #pragma omp for schedule(static)
      for(int i=1;i < m_nActs; ++i)
      {
         tmp +=  c(i,0) * m_infF.image(i);
      }
      #pragma omp critical
      shape += tmp;
   }

   t1 = sys::get_curr_time();
   t_sum += t1-t0;

   wf.phase += 2*amp*shape*m_pupil*math::two_pi<realT>()/lambda;

}

template<typename realT>
void makeMap( Eigen::Array<realT, -1, -1> & map,  Eigen::Array<realT, -1, -1> & pos, Eigen::Array<realT, -1, -1> & act)
{

   realT minx = pos.row(0).minCoeff();
   realT maxx = pos.row(0).maxCoeff();

   int i=0;
   realT dx = 0;
   while(dx == 0)
   {
      dx = fabs(pos(0,i)- pos(0,0));
      ++i;
   }

   realT miny = pos.row(1).minCoeff();
   realT maxy = pos.row(1).maxCoeff();

   i = 0;
   realT dy = 0;

   while(dy == 0)
   {
      dy = fabs(pos(1,i)- pos(1,0));
      ++i;
   }

   int nx = (maxx-minx)/dx + 1;
   int ny = (maxy-miny)/dy + 1;


   map.resize(nx, ny);
   map.setZero();

   realT x, y;

   for(int j=0;j < pos.cols(); ++j)
   {
      x = (pos(0,j) - minx)/dx;

      y = ( pos(1,j) - miny ) / dy;

      map( (int) x, (int) y) = act(j,0);
   }



}


template<typename _realT>
template<typename commandT>
void deformableMirror<_realT>::setShape(commandT & commandV)
{

   if(_settling)
   {
      std::cerr << "DM: new command received while still settling.\n";
      return;
   }


   static Eigen::Array<_realT, -1, -1> c; //static to avoid re-alloc, todo: make class member

   #ifdef MXAO_USE_GPU
   m_devModeCommands.upload(commandV.measurement.data(), commandV.measurement.size());
   realT alpha;
   cublasStatus_t stat = cuda::cublasTgemv<realT>(*m_cublasHandle, CUBLAS_OP_N,  m_nActs, m_nActs, m_alpha(), m_devM2c(), m_devModeCommands(), m_zero(), m_devActCommands());

   if(stat != CUBLAS_STATUS_SUCCESS)
   {
      std::cerr << "cublas error\n";
   }
   

   //c.resize(m_nActs,1);
   //m_devActCommands.download(c.data());
   
   
   #else
   Eigen::Map<Eigen::Array<realT,-1,-1>> commandV_measurement(commandV.measurement.data(), 1, commandV.measurement.size());
   c = -1*_calAmp*m_m2c.matrix() * commandV_measurement.matrix().transpose();
   
   #endif

#if 0
   
   if(_commandLimit > 0)
   {

      realT r1 = sqrt( pow(_pos(0,1) - _pos(0,0),2) + pow(_pos(1,1) - _pos(1,0),2));

      realT r;

      int nLimited = 0;

      for(int i=0; i< _pos.cols(); ++i)
      {
         for(int j=i+1; j< _pos.cols(); ++j)
         {
            r = sqrt( pow(_pos(0,j) - _pos(0,i),2) + pow(_pos(1,j) - _pos(1,i),2));

            if( fabs(r1 - r) < .01 )
            {
               realT iact = fabs( c(i,0) - c(j,0) );
               if( iact > _commandLimit )
               {
                  std::cerr << "Limited Actuators " << i << " " << j << "\n";
                  ++nLimited;
                  c(i,0) *= _commandLimit/iact;
                  c(j,0) *= _commandLimit/iact;
               }
            }
         }
      }
      if(nLimited > 0) std::cerr << nLimited << " strokes limited\n";

   }

#if 0
   if(_commandLimit > 0 )
   {
      for(int i=0; i < c.rows(); ++i)
      {
         if(c(i,0) > _commandLimit ) c(i,0) = _commandLimit;
         if(c(i,0) < -1*_commandLimit ) c(i,0) = -1*_commandLimit;
      }
   }
#endif



   if( _writeCommands )
   {
      if(!_commandFileOpen)
      {
         _commandFout.open( _commandFile );
         _commandFout << std::scientific;
         _commandFileOpen = true;
      }

      _commandFout << commandV.iterNo << "> ";
      for(int i=0; i<c.rows(); ++i)
      {
         _commandFout << c(i,0) << " ";
      }
      _commandFout << std::endl;
   }


   bool skipFrame = 0;

   ///\todo Should check for command limits here.
   for(int i=0; i< m_nActs; ++i)
   {
      if( std::isnan( c(i,0) ) || !std::isfinite(c(i,0)))
      {
         skipFrame = true;
         break;
      }
   }

   if(skipFrame)
   {
      std::cerr << "SKIP FRAME\n";
      return;
   }

#endif

   #ifdef MXAO_USE_GPU
   
   static std::vector<realT> tmp;
   tmp.resize(m_idx.size());
   
   m_devShape.resize(m_idx.size()); //no-op except first time.
   m_devShape.initialize();
   
   for(int i=0; i < m_nActs; ++i)
   {
      //realT alpha = c(i,0);
      cuda::cublasTaxpy<realT>(*m_cublasHandle, m_idx.size(), m_devActCommands()+i, m_devInfF() + i*m_idx.size(), 1, m_devShape(), 1);
   }
   
   m_devShape.download(tmp.data());
   
   _nextShape.setZero();
   
   realT * imP = _nextShape.data();
   
   for(size_t n=0; n<m_idx.size(); ++n)  imP[m_idx[n]] = tmp[n];
         
   #else
   
   _nextShape = c(0,0)*m_infF.image(0);//*sigma[0];
   
   #pragma omp parallel
   {
      Eigen::Array<_realT, -1, -1> tmp ;
      tmp.resize(m_nRows, m_nCols);
      tmp.setZero();

      realT * tmpP = tmp.data();
      #pragma omp for
      for(int i=1;i < m_nActs; ++i)
      {
         realT * imP = m_infF.image(i).data();
         for(size_t n=0; n<m_idx.size(); ++n) *(tmpP + m_idx[n]) += c(i,0) * *(imP + m_idx[n]);
      }
      #pragma omp critical
      {
         realT * imP = _nextShape.data();
         for(size_t n=0; n<m_idx.size(); ++n)  *(imP + m_idx[n]) += *(tmpP + m_idx[n]);
      }
   }
   #endif

    //================ filter here!!
    

   if(m_applyFilter)
   {
      std::cerr << "DM filtering\n";
      m_filter.filter(_nextShape);
   }
      
   _oldShape = _shape;
#if 1
   _settling = 1;
   _settleTime_counter = 0;
   _settlingIter = commandV.iterNo;
#else
//Ignoring settling time.
   _shape = _nextShape;
#endif

}

template<typename _realT>
void deformableMirror<_realT>::applyShape(wavefrontT & wf,  realT lambda)
{

#if 1
   BREAD_CRUMB;

   if(_settling)
   {
      BREAD_CRUMB;
       _shape = _oldShape + (_nextShape-_oldShape)*( (realT) _settleTime_counter+1.0)/( (realT) _settleTime);

       realT mn = (_shape*m_pupil).sum()/m_pupilSum;
       _shape = (_shape - mn)*m_pupil;
       
       ++_settleTime_counter;

       if(_settleTime_counter >= _settleTime)
       {
          _settling = 0;
          _settledIter = _settlingIter;
       }
   }
#endif


   BREAD_CRUMB;

   wf.phase += 2*_shape*math::two_pi<realT>()/lambda;
   
   BREAD_CRUMB;

   #ifdef MXAO_DIAG_LOOPCOUNTS
   std::cerr << "DM " << m_nActs << " Shape applied: " << wf.iterNo << " " << _settledIter << "\n";
   #endif
}

template<typename _realT>
void deformableMirror<_realT>::setFilter( int width )
{
   int nr = _shape.rows();
   int nc = _shape.cols();

   typename wfsImageT<_realT>::imageT filterMask;

   filterMask.resize( nr, nc );
   filterMask.setZero();

   filterMask.block(0,0, width+1, width+1).setConstant(1.0);
   filterMask.block(0, nc - width, width+1, width).setConstant(1.0);
   filterMask.block(nr - width, 0, width, width+1).setConstant(1.0);
   filterMask.block(nr - width, nc - width, width, width).setConstant(1.0);

         
   m_filter.psdSqrt(filterMask, 1.0/_shape.rows(), 1.0/_shape.cols());



}

} //sim
} //AO
} //namespace mx

#endif //deformableMirror_hpp
