/** \file aoWFS.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Definitions of various analytic wavefront sensors
  * \ingroup mxAO_files
  * 
  */

#ifndef aoWFS_hpp
#define aoWFS_hpp

#include "../../math/constants.hpp"

namespace mx
{

namespace AO
{

namespace analysis
{
   
///The ideal wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
  * for the ideal WFS.
  * 
  * This is the base class for all WFS.
  * 
  * \tparam realT is the floating point type used for calculations
  * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT, typename iosT = std::ostream>
struct wfs 
{
   std::string _id;
   
   ///Constructor
   /** Only sets the value of _id.
     */
   wfs()
   {
      _id = "Ideal WFS";
   }
   
   ///Destructor
   /** Declared virtual so more complicated derived types can be created.
     */
   virtual ~wfs()
   {
      return;
   }
   
   ///Get the sensitivity at a spatial frequency.
   /** The sensitivity of the ideal WFS is 1 at all k \cite guyon_2005.
     * 
     * \returns the sensitivity to photon noise parameter
     */ 
   virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                         int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                         realT D, ///< [in] the telescope diameter (not used by this WFS)
                         realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                         realT r0 ///< [in] Fried's parameter (not used by this WFS)
                       )
   {
      //Stuff a sock in the compiler's mouth:
      static_cast<void>(m);
      static_cast<void>(n);
      static_cast<void>(D);
      static_cast<void>(d);
      static_cast<void>(r0);
      
      return static_cast<realT>(1);
   }
   
   ///Dump the details of the WFS to an io stream.
   /** Is virtual so that derived types can add parameters.
     */
   virtual iosT & dumpWFS(iosT & ios)
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
  * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT, typename iosT = std::ostream>
struct pywfsUnmod : public wfs<realT, iosT>
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
   virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                         int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                         realT D, ///< [in] the telescope diameter (not used by this WFS)
                         realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                         realT r0 ///< [in] Fried's parameter (not used by this WFS)  
                       )
   {
      //Stuff a sock in the compiler's mouth:
      static_cast<void>(m);
      static_cast<void>(n);
      static_cast<void>(D);
      static_cast<void>(d);
      static_cast<void>(r0);
      
      return math::root_two<realT>();
   }
   
};


///The asymptotic modulated pyramid wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
  * for the modulated PyWFS in the asymptotic limit.
  * 
  * \tparam realT is the floating point type used for calculations
  * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT, typename iosT = std::ostream>
struct pywfsModAsymptotic : public wfs<realT, iosT>
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
   virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u (not used by this WFS)
                         int n,   ///< [in] the spatial frequency index for v (not used by this WFS)
                         realT D, ///< [in] the telescope diameter (not used by this WFS)
                         realT d, ///< [in] the sub-ap spacing (not used by this WFS)
                         realT r0 ///< [in] Fried's parameter (not used by this WFS)
                       )
   {
      //Stuff a sock in the compiler's mouth:
      static_cast<void>(m);
      static_cast<void>(n);
      static_cast<void>(D);
      static_cast<void>(d);
      static_cast<void>(r0);
      
      return static_cast<realT>(2) * math::root_two<realT>();
   }
   
};

///The shack hartmann wavefront sensor sensitivity function.
/** Provides the \f$ \beta_p \f$ parameter of Guyon, 2005 \cite guyon_2005
  * for the shack hartmann WFS.
  * 
  * \tparam realT is the floating point type used for calculations
  * \tparam iosT is an output stream type with operator \<\< defined (default is std::ostream)
  * 
  * \ingroup mxAOAnalytic
  */ 
template<typename realT, typename iosT = std::ostream>
struct shwfs : public wfs<realT, iosT>
{
   shwfs()
   {
      this->_id = "Shack Hartmann";
   }
   
   ///Get the sensitivity at a spatial frequency.
   /** The sensitivity of the shack hartmann WFS
     * 
     * \returns the sensitivity to photon noise parameter
     */ 
   virtual realT beta_p( int m,   ///< [in] the spatial frequency index for u
                         int n,   ///< [in] the spatial frequency index for v
                         realT D, ///< [in] the telescope diameter
                         realT d, ///< [in] the sub-ap spacing
                         realT r0 ///< [in] Fried's parameter
               )
   {
      realT k = sqrt(m*m + n*n)/D;
      
      return 1.48/(k*d)*sqrt(1 + pow(d/r0,2));
      
   }
   
};

} //namespace analysis
} //namespace AO
} //namespace mx

#endif // aoWFS_hpp

