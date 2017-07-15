/** \file aoWFS.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions of various wavefront sensors
  * \ingroup mxAO_files
  * 
  */

#ifndef __aoWFS_hpp__
#define __aoWFS_hpp__

namespace mx
{

namespace AO
{

///The ideal wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
  * for the ideal WFS.
  * 
  * This is the base class for all WFS.
  * 
  * \tparam realT is the floating point type used for calculations
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT>
struct wfs 
{
   std::string _id;
   
   wfs()
   {
      _id = "Ideal WFS";
   }
   
   /** The sensitivity of the ideal WFS is 1  at all k.
     * 
     * \returns the sensitivity to photon noise parameter
     */ 
   virtual realT beta_p( int m, ///< [in] the spatial frequency index for u
                         int n,  ///< [in] the spatial frequency index for v
                         realT D ///< [in] the diameter
                       )
   {
      return static_cast<realT>(1);
   }
   
   ///Dump the details of the WFS to an io stream.
   template<typename iosT>
   iosT & dumpWFS(iosT & ios)
   {
      ios << "# WFS Parameters:\n";
      ios << "#    ID = " << _id << '\n';
      
      return ios;
   }
};
   
///The unmodulated pyramid wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
  * for the unmodulated PyWFS.
  * 
  * \tparam realT is the floating point type used for calculations
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT>
struct pywfsUnmod : public wfs<realT>
{
   pywfsUnmod()
   {
      this->_id = "Unmodulated Pyramid";
   }
   
   ///Get the sensitivity at a spatial frequency.
   /** The sensitivity of the unmodulated PyWFS is \f$ \sqrt{2} \f$ at all k.
     * 
     * \returns the sensitivity to photon noise parameter
     */ 
   realT beta_p( int m, ///< [in] the spatial frequency index for u
                 int n,  ///< [in] the spatial frequency index for v
                 realT D ///< [in] the diameter
               )
   {
      using namespace boost::math::constants;
      return root_two<realT>();
   }
   
};


///The asymptotic modulated pyramid wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
  * for the modulated PyWFS in the asymptotic limit.
  * 
  * \tparam realT is the floating point type used for calculations
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT>
struct pywfsModAsymptotic : public wfs<realT>
{
   pywfsModAsymptotic()
   {
      this->_id = "Asymptotic Modulated Pyramid";
   }
   
   ///Get the sensitivity at a spatial frequency.
   /** The sensitivity of the asymptotic modulated PyWFS is \f$ 2 \sqrt{2} \f$ at all k.
     * 
     * \returns the sensitivity to photon noise parameter
     */ 
   realT beta_p( int m, ///< [in] the spatial frequency index for u
                 int n,  ///< [in] the spatial frequency index for v
                 realT D ///< [in] the diameter
               )
   {
      using namespace boost::math::constants;
      return static_cast<realT>(2) * root_two<realT>();
   }
   
};


} //namespace AO
} //namespace mx

#endif // __aoWFS_hpp__

