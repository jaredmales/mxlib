/** \file ccdDetector.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to simulate a CCD.
  * \ingroup mxAO_files
  * 
  */

#ifndef __ccdDetector_hpp__
#define __ccdDetector_hpp__


//#include <Eigen/Dense>

#include "mx/imageTransforms.hpp"
#include "mx/randomT.hpp"

#include "wavefront.hpp"

#include <fcntl.h>
#include <iostream>

namespace mx
{
namespace AO
{
namespace sim
{

   
   
///A simulated CCD detector
/** A simulated CCD detector, including an optional EMCCD noise model.
  *
  * \ingroup mxAOSim
  */   
template<typename _floatT>
class ccdDetector 
{
   
public:
   
   typedef _floatT floatT;
   
   typedef wavefront<floatT> wavefrontT;
   
   typedef Eigen::Array<floatT, Eigen::Dynamic, Eigen::Dynamic> imageT;
   
   typedef mx::randomT<floatT, std::mt19937_64, std::normal_distribution<floatT> > norm_distT;
      
   typedef mx::randomT<int, std::mt19937_64, std::poisson_distribution<int> > poisson_distT;
   
   
   typedef mx::randomT<floatT, std::mt19937_64, std::gamma_distribution<floatT> > gamma_distT;
   
   
   
   ccdDetector();

protected:
      
   norm_distT normVar; ///<Gets normal-distributed variates
   poisson_distT poissonVar; ///<Gets Poisson distributed variates
   gamma_distT gammaVar; ///<Gets gamma distributed variates
   
   floatT _qe; ///<The quantum efficiency.
   
   floatT _darkCurrent; ///<The dark current, per pixel per second
   floatT _ron; ///<The readout noise, electrons per pixel per read #include "mx/randomT.hpp"
   
   floatT _cic; ///< EMCCD clock induced charge, electrons per pixel per read.
   floatT _gain; ///<Electron multiplication gain.  If >1, then EMCCD is modeled.
   
   floatT _expTime; ///<The exposure time, in seconds.
   
   int _rows; ///<The detector size, in rows.
   int _cols; ///<The detector size, in columns.
   
   bool _noNoise;
   
   long get_seed();
   
public:
   
   floatT qe();
   
   void qe(const floatT & q);
   
   floatT darkCurrent();
   
   void darkCurrent(const floatT & dc);
   
   floatT ron();
   
   void ron(const floatT & r);
   
   floatT cic();
   
   void cic(const floatT & c);
   
   floatT gain();
   
   void gain(const floatT & g);

   floatT expTime();
   
   void expTime(const floatT &dt);
   
   int rows();
   
   void rows(const int & r);
   
   int cols();
   
   void cols(const int & c);
   
   void setSize(const int &r, const int &c);
   
   bool noNoise();
   
   void noNoise(bool nn);
   
   ///Rebin and add noise to the input image, placing the result in the output image.
   /** The output image must be the same size or smaller than the input image.  
     * The output image is resized only if necessary.
     * The input image is multiplied by expTime, so its flux should be in photons/sec.
     * Noise is modeled as follows:
     * -# The number of input photons for pixel (i,j) is \f$ n_{ie} = (in(i,j) + darkCurrent)*expTime + cic \f$
     * -# The Poisson distribution with mean and variance \f$ n_{ie} \f$ is sampled. 
     * -# If gain > 1, the gamma distribution is sampled
     * -# Read noise is applied
     * -# Result is rounded to nearest whole integer.
     */
   void exposeImage( imageT & out, ///< [out] The output image, after all above steps applied. 
                     imageT & in ///< [in] the input image, in photons/sec flux units.
                   );
   
      
   
};



template<typename floatT>
ccdDetector<floatT>::ccdDetector()
{
//   long seed = get_seed();
   normVar.seed();
   poissonVar.seed();
   gammaVar.seed();
   
   _qe = 1;
   
   _darkCurrent = 0;
   _ron = 0;
   
   _noNoise = false;
   
   _cic = 0;
   _gain = 1;
   
   _expTime = 1;
   
   _rows = 0;
   _cols = 0;
   
   
   
   
}


template<typename floatT>
long ccdDetector<floatT>::get_seed()
{
   long int seedval;
     
   int fd;
      
   fd = open("/dev/urandom", O_RDONLY);
      
   int rv =read(fd, &seedval, sizeof(long int));
      
    close(fd);
      
    
   if(rv < 0)
   {
      std::cerr << "read from /dev/urandom returned error\n";
         
      return 0;
   }
     
   return seedval;
}

template<typename floatT>
floatT ccdDetector<floatT>::qe()
{
   return _qe;
}

template<typename floatT>
void ccdDetector<floatT>::qe(const floatT & q)
{
   _qe = q;
}

template<typename floatT>
floatT ccdDetector<floatT>::darkCurrent()
{
   return _darkCurrent;
}

template<typename floatT>
void ccdDetector<floatT>::darkCurrent(const floatT & dc)
{
   _darkCurrent = dc;
}
   
template<typename floatT>
floatT ccdDetector<floatT>::ron()
{
   return _ron;
}
  
template<typename floatT>
void ccdDetector<floatT>::ron(const floatT & r)
{
   _ron = r;
}

template<typename floatT>
floatT ccdDetector<floatT>::cic()
{
   return _cic;
}
   
template<typename floatT>
void ccdDetector<floatT>::cic(const floatT & c)
{
   _cic = c;
}

template<typename floatT>
floatT ccdDetector<floatT>::gain()
{
   return _gain;
}
   
template<typename floatT>
void ccdDetector<floatT>::gain(const floatT & g)
{
   _gain = g;
}

template<typename floatT>
floatT ccdDetector<floatT>::expTime()
{
   return _expTime;
}
   
template<typename floatT>
void ccdDetector<floatT>::expTime(const floatT & dt)
{
   _expTime = dt;
}


template<typename floatT>
int ccdDetector<floatT>::rows()
{
   return _rows;
}
   
template<typename floatT>
void ccdDetector<floatT>::rows(const int & r)
{
   _rows = r;
}

template<typename floatT>
int ccdDetector<floatT>::cols()
{
   return _cols;
}
   
template<typename floatT>
void ccdDetector<floatT>::cols(const int & c)
{
   _cols = c;
}

template<typename floatT>
void ccdDetector<floatT>::setSize(const int & r, const int & c)
{
   _rows = r;
   _cols = c;
}


template<typename floatT>
bool ccdDetector<floatT>::noNoise()
{
   return _noNoise;
}

template<typename floatT>
void ccdDetector<floatT>::noNoise(bool nn)
{
   _noNoise = nn;
}


template<typename floatT>
void ccdDetector<floatT>::exposeImage(imageT & out, imageT & in)
{
   
   using poisson_param_t = typename std::poisson_distribution<int>::param_type;
   using gamma_param_t = typename std::gamma_distribution<floatT>::param_type;
   
   
   out.resize(_rows, _cols);
   
   imageDownSample(out, in);
         
   if(_noNoise) return;
   
   for(int i=0;i<_rows;++i)
   {
      for(int j=0; j<_cols;++j)
      {
         floatT charge = (out(i,j)*_qe + _darkCurrent) * _expTime + _cic;
         
         if(charge > 1000)
         {
            out(i,j) = charge + normVar*sqrt(charge);
         }
         else
         {
            poissonVar.distribution.param(poisson_param_t{charge});         
            out(i,j) = poissonVar;
         }
         
         if(_gain > 1)
         {
            gammaVar.distribution.param(gamma_param_t{out(i,j), _gain});
            
            
            out(i,j) = gammaVar;
            
         }
   
         out(i,j) += normVar*_ron;
         
         if(_gain > 1) out(i,j) /= _gain;
         
         out(i,j) = round(out(i,j));
         
      }
   }
         
         
         
   
   
}

} //namespace sim
} //namespace AO
} //namespace mx

#endif //__ccdDetector_hpp__
