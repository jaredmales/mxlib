/** \file units.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Unit specifications and conversions.
  * \ingroup astrofiles
  * 
  */

#ifndef __mx_astro_units_hpp__
#define __mx_astro_units_hpp__

#include "constants.hpp"

namespace mx
{

namespace astro
{
   

namespace units
{

/** \defgroup astrounits Unit Conversions 
  * \brief Definitions of unit coversions for physical constants.
  * \ingroup phyconstants
  * 
  * These types provide the ratios to convert from SI units to the desired system.
  * 
  * @{
  */

/// International System of Units (SI) units-type
/** Since all constants are specified in SI units, all conversions here are 1.0.
  */
template<typename _realT>
struct si
{
   typedef _realT realT; ///< The real floating point type in which to specify constants.
   static constexpr realT length = static_cast<realT>(1.0); ///< Conversion from SI (m) to SI (m)
   static constexpr realT time = static_cast<realT>(1.0); ///< Conversion from SI (s) to SI (s)
   static constexpr realT mass = static_cast<realT>(1.0); ///< Conversion from SI (kg) to SI (kg)
   static constexpr realT energy = static_cast<realT>(1.0); ///< Conversion from SI (J) to SI (J)
   static constexpr realT temperature = static_cast<realT>(1.0); ///< Conversion from SI (K) to SI (K)
};

/// Centimeter-Gram-Second (cgs) units-type
/** The units are:
  * - Length: centimeter (cm)
  * - Time: second (s)
  * - Mass: gram (g)
  * - Energy: erg (erg)
  * - Temperature: Kelvin (K)
  */
template<typename _realT>
struct cgs
{
   typedef _realT realT; ///< The real floating point type in which to specify constants.
   static constexpr realT length = static_cast<realT>(100.0);  ///< Conversion from SI (m) to cgs (cm)
   static constexpr realT time = static_cast<realT>(1.0); ///< Conversion from SI (s) to cgs (s)
   static constexpr realT mass = static_cast<realT>(1000.0); ///< Conversion from SI (kg) to cgs (g)
   static constexpr realT energy = static_cast<realT>(1e7); ///< Conversion from SI (J) to cgs (erg)
   static constexpr realT temperature = static_cast<realT>(1.0); ///< Conversion from SI (K) to cgs (K)
};

/// Solar units-type
/** The units are:
  * - Length: au
  * - Time: year
  * - Mass: solar mass
  * - Energy: solar-luminosity-year
  * - Temperature: solar effective temperature
  */
template<typename _realT>
struct solar
{
   typedef _realT realT; ///< The real floating point type in which to specify constants.
   static constexpr realT length = static_cast<realT>(1)/constants::au<si<realT>>();  ///< Conversion from SI (m) to solar (au)
   static constexpr realT time = static_cast<realT>(1.0) / constants::year<si<realT>>(); ///< Conversion from SI (s) to solar (yr)
   static constexpr realT mass = constants::G<si<realT>>()/constants::GMSun<si<realT>>(); ///< Conversion from SI (kg) to solar (M_sun)
   static constexpr realT energy = static_cast<realT>(1.0)/( constants::lumSun<si<realT>>()*static_cast<realT>(365.25)*static_cast<realT>(86400.0)); ///< Conversion from SI (J) to  solar (Solar-luminosities X year)
   static constexpr realT temperature = static_cast<realT>(1.0) / constants::TeffSun<si<realT>>(); ///< Conversion from SI (K) to solar (5772 K)
};

/// Earth units-type
/** The units are:
  * - Length: Earth-radii
  * - Time: s
  * - Mass: Earth mass
  * - Energy: J
  * - Temperature: K
  */
template<typename _realT>
struct earth
{
   typedef _realT realT; ///< The real floating point type in which to specify constants.
   static constexpr realT length = static_cast<realT>(1)/constants::radEarth<si<realT>>();  ///< Conversion from SI (m) 
   static constexpr realT time = static_cast<realT>(1.0); ///< Conversion from SI (s) 
   static constexpr realT mass = static_cast<realT>(1.0)/constants::massEarth<si<realT>>(); ///< Conversion from SI (kg) 
   static constexpr realT energy = static_cast<realT>(1.0); ///< Conversion from SI (J)
   static constexpr realT temperature = static_cast<realT>(1.0); ///< Conversion from SI (K)
};

/// Jupiter units-type
/** The units are:
  * - Length: Jupiter-radii
  * - Time: s
  * - Mass: Jupiter mass
  * - Energy: J
  * - Temperature: K
  */
template<typename _realT>
struct jupiter
{
   typedef _realT realT; ///< The real floating point type in which to specify constants.
   static constexpr realT length = static_cast<realT>(1)/constants::radJupiter<si<realT>>();  ///< Conversion from SI (m) 
   static constexpr realT time = static_cast<realT>(1.0); ///< Conversion from SI (s) 
   static constexpr realT mass = static_cast<realT>(1.0)/constants::massJupiter<si<realT>>(); ///< Conversion from SI (kg) 
   static constexpr realT energy = static_cast<realT>(1.0); ///< Conversion from SI (J)
   static constexpr realT temperature = static_cast<realT>(1.0); ///< Conversion from SI (K)
};

///@}

}//namespace units

} //namespace astro
} //namespace mx


#endif //__mx_astro_units_hpp__

