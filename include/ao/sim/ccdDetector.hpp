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
#include "../../math/randomT.hpp"

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
template<typename _realT>
class ccdDetector
{

public:

   typedef _realT realT;

   typedef wavefront<realT> wavefrontT;

   typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

   typedef mx::math::randomT<realT, std::mt19937_64, std::normal_distribution<realT> > norm_distT;

   typedef mx::math::randomT<int, std::mt19937_64, std::poisson_distribution<int> > poisson_distT;


   typedef mx::math::randomT<realT, std::mt19937_64, std::gamma_distribution<realT> > gamma_distT;



   ccdDetector();

protected:

   norm_distT m_normVar; ///<Gets normal-distributed variates
   poisson_distT m_poissonVar; ///<Gets Poisson distributed variates
   gamma_distT m_gammaVar; ///<Gets gamma distributed variates

   realT m_qe {1}; ///<The quantum efficiency.

   realT m_darkCurrent {0}; ///<The dark current, per pixel per second
   realT m_ron {0}; ///<The readout noise, electrons per pixel per read #include "mx/randomT.hpp"

   realT m_cic {0}; ///< EMCCD clock induced charge, electrons per pixel per read.
   realT m_gain {1}; ///<Electron multiplication gain.  If >1, then EMCCD is modeled.

   realT m_expTime {1}; ///<The exposure time, in seconds.

   int m_rows {0}; ///<The detector size, in rows.
   int m_cols {0}; ///<The detector size, in columns.

   bool m_noNoise {false}; ///< If true no noise is added to the exposed image.

   //This is probably vestigial, commented out on 2018-07-23
   //Delete after testing with an aoSim compile.
   //long get_seed();

public:

   /// Get the current value of qe.
   /**
     * \returns the current value of m_qe
     */
   realT qe();

   void qe(const realT & q);

   realT darkCurrent();

   void darkCurrent(const realT & dc);

   realT ron();

   void ron(const realT & r);

   realT cic();

   void cic(const realT & c);

   realT gain();

   void gain(const realT & g);

   realT expTime();

   void expTime(const realT &dt);

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



template<typename realT>
ccdDetector<realT>::ccdDetector()
{
   m_normVar.seed();
   m_poissonVar.seed();
   m_gammaVar.seed();
}


/*template<typename realT>
long ccdDetector<realT>::get_seed()
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
*/

template<typename realT>
realT ccdDetector<realT>::qe()
{
   return m_qe;
}

template<typename realT>
void ccdDetector<realT>::qe(const realT & q)
{
   m_qe = q;
}

template<typename realT>
realT ccdDetector<realT>::darkCurrent()
{
   return m_darkCurrent;
}

template<typename realT>
void ccdDetector<realT>::darkCurrent(const realT & dc)
{
   m_darkCurrent = dc;
}

template<typename realT>
realT ccdDetector<realT>::ron()
{
   return m_ron;
}

template<typename realT>
void ccdDetector<realT>::ron(const realT & r)
{
   m_ron = r;
}

template<typename realT>
realT ccdDetector<realT>::cic()
{
   return m_cic;
}

template<typename realT>
void ccdDetector<realT>::cic(const realT & c)
{
   m_cic = c;
}

template<typename realT>
realT ccdDetector<realT>::gain()
{
   return m_gain;
}

template<typename realT>
void ccdDetector<realT>::gain(const realT & g)
{
   m_gain = g;
}

template<typename realT>
realT ccdDetector<realT>::expTime()
{
   return m_expTime;
}

template<typename realT>
void ccdDetector<realT>::expTime(const realT & dt)
{
   m_expTime = dt;
}


template<typename realT>
int ccdDetector<realT>::rows()
{
   return m_rows;
}

template<typename realT>
void ccdDetector<realT>::rows(const int & r)
{
   m_rows = r;
}

template<typename realT>
int ccdDetector<realT>::cols()
{
   return m_cols;
}

template<typename realT>
void ccdDetector<realT>::cols(const int & c)
{
   m_cols = c;
}

template<typename realT>
void ccdDetector<realT>::setSize(const int & r, const int & c)
{
   m_rows = r;
   m_cols = c;
}


template<typename realT>
bool ccdDetector<realT>::noNoise()
{
   return m_noNoise;
}

template<typename realT>
void ccdDetector<realT>::noNoise(bool nn)
{
   m_noNoise = nn;
}


template<typename realT>
void ccdDetector<realT>::exposeImage(imageT & out, imageT & in)
{

   using poisson_param_t = typename std::poisson_distribution<int>::param_type;
   using gamma_param_t = typename std::gamma_distribution<realT>::param_type;


   out.resize(m_rows, m_cols);

   improc::imageDownSample(out, in);

   if(m_noNoise) return;

   for(int i=0;i<m_rows;++i)
   {
      for(int j=0; j<m_cols;++j)
      {
         realT charge = (out(i,j)*m_qe + m_darkCurrent) * m_expTime + m_cic;

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
