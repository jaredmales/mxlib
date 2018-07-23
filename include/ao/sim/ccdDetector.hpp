/** \file ccdDetector.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Provides a class to simulate a CCD.
  * \ingroup mxAO_files
  *
  */

#ifndef ccdDetector_hpp
#define ccdDetector_hpp


//#include <Eigen/Dense>

#include "../../improc/imageTransforms.hpp"
#include "../../randomT.hpp"

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

   norm_distT m_normVar; ///<Gets normal-distributed variates
   poisson_distT m_poissonVar; ///<Gets Poisson distributed variates
   gamma_distT m_gammaVar; ///<Gets gamma distributed variates

   floatT m_qe {1}; ///<The quantum efficiency.

   floatT m_darkCurrent {0}; ///<The dark current, per pixel per second
   floatT m_ron {0}; ///<The readout noise, electrons per pixel per read #include "mx/randomT.hpp"

   floatT m_cic {0}; ///< EMCCD clock induced charge, electrons per pixel per read.
   floatT m_gain {1}; ///<Electron multiplication gain.  If >1, then EMCCD is modeled.

   floatT m_expTime {1}; ///<The exposure time, in seconds.

   int m_rows {0}; ///<The detector size, in rows.
   int m_cols {0}; ///<The detector size, in columns.

   bool m_noNoise {false};

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
   m_normVar.seed();
   m_poissonVar.seed();
   m_gammaVar.seed();
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
   return m_qe;
}

template<typename floatT>
void ccdDetector<floatT>::qe(const floatT & q)
{
   m_qe = q;
}

template<typename floatT>
floatT ccdDetector<floatT>::darkCurrent()
{
   return m_darkCurrent;
}

template<typename floatT>
void ccdDetector<floatT>::darkCurrent(const floatT & dc)
{
   m_darkCurrent = dc;
}

template<typename floatT>
floatT ccdDetector<floatT>::ron()
{
   return m_ron;
}

template<typename floatT>
void ccdDetector<floatT>::ron(const floatT & r)
{
   m_ron = r;
}

template<typename floatT>
floatT ccdDetector<floatT>::cic()
{
   return m_cic;
}

template<typename floatT>
void ccdDetector<floatT>::cic(const floatT & c)
{
   m_cic = c;
}

template<typename floatT>
floatT ccdDetector<floatT>::gain()
{
   return m_gain;
}

template<typename floatT>
void ccdDetector<floatT>::gain(const floatT & g)
{
   m_gain = g;
}

template<typename floatT>
floatT ccdDetector<floatT>::expTime()
{
   return m_expTime;
}

template<typename floatT>
void ccdDetector<floatT>::expTime(const floatT & dt)
{
   m_expTime = dt;
}


template<typename floatT>
int ccdDetector<floatT>::rows()
{
   return m_rows;
}

template<typename floatT>
void ccdDetector<floatT>::rows(const int & r)
{
   m_rows = r;
}

template<typename floatT>
int ccdDetector<floatT>::cols()
{
   return m_cols;
}

template<typename floatT>
void ccdDetector<floatT>::cols(const int & c)
{
   m_cols = c;
}

template<typename floatT>
void ccdDetector<floatT>::setSize(const int & r, const int & c)
{
   m_rows = r;
   m_cols = c;
}


template<typename floatT>
bool ccdDetector<floatT>::noNoise()
{
   return m_noNoise;
}

template<typename floatT>
void ccdDetector<floatT>::noNoise(bool nn)
{
   m_noNoise = nn;
}


template<typename floatT>
void ccdDetector<floatT>::exposeImage(imageT & out, imageT & in)
{

   using poisson_param_t = typename std::poisson_distribution<int>::param_type;
   using gamma_param_t = typename std::gamma_distribution<floatT>::param_type;


   out.resize(m_rows, m_cols);

   imageDownSample(out, in);

   if(m_noNoise) return;

   for(int i=0;i<m_rows;++i)
   {
      for(int j=0; j<m_cols;++j)
      {
         floatT charge = (out(i,j)*m_qe + m_darkCurrent) * m_expTime + m_cic;

         if(charge > 1000)
         {
            out(i,j) = charge + m_normVar*sqrt(charge);
         }
         else
         {
            m_poissonVar.distribution.param(poisson_param_t{charge});
            out(i,j) = m_poissonVar;
         }

         if(m_gain > 1)
         {
            m_gammaVar.distribution.param(gamma_param_t{out(i,j), m_gain});


            out(i,j) = m_gammaVar;

         }

         out(i,j) += m_normVar*m_ron;

         if(m_gain > 1) out(i,j) /= m_gain;

         out(i,j) = round(out(i,j));

      }
   }





}

} //namespace sim
} //namespace AO
} //namespace mx

#endif //ccdDetector_hpp
