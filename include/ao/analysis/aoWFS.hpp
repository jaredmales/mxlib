/** \file aoWFS.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions of various wavefront sensors
  * \ingroup mxAOAnalytic_files
  * 
  */

#ifndef __aoWFS_hpp__
#define __aoWFS_hpp__

namespace mx
{

namespace AO
{

///The unmodulated pyramid wavefront sensor sensitivity function.
/** Calculates the \f$\beta_p\f$ parameter of Guyon, 2005 \cite guyon_2005
  *
  * \tparam floatT is the floating point type used for calculations
  */ 
template<typename floatT>
struct pywfsUnmod
{
   const char * _id = "Unmodulated Pyramid";
   
   ///Get the sensitivity at a spatial frequency.
   /**
     * \param m the spatial frequency index for u
     * \param n the spatial frequency index for v
     * \param D the diameter
     * 
     * \returns the sensitivity to photon noise
     */ 
   floatT beta_p(int m, int n, floatT D)
   {
      using namespace boost::math::constants;
      return root_two<floatT>();
   }
   
   ///Dump the details of the WFS to a io stream.
   template<typename iosT>
   iosT & dumpWFS(iosT & ios)
   {
      ios << "# WFS Parameters:\n";
      ios << "#    ID = " << _id << '\n';
      
      return ios;
   }
};



} //namespace AO
} //namespace mx

#endif // __aoWFS_hpp__

