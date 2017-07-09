/** \file constants.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Constants for astronomy.
  * \ingroup astrofiles
  * 
  */

#ifndef __mx_astro_constants_hpp__
#define __mx_astro_constants_hpp__


namespace mx
{

namespace astro
{
   
namespace constants
{


/** \defgroup dimensionless_constants Dimensionless Constants 
  * \ingroup astroconstants
  * Constants without dimensions are provided with a template argument specifying the precision.
  * @{
  */

///Tangent of 1.0 arcsec 
/** Calculated with 100 digits of precision, recorded here with 50 digits.
  */
template<typename realT>
constexpr realT tan_arcsec()
{
   return static_cast<realT>(4.8481368111333441675396429478852851658848753880815e-06);
}

///@} Dimensionless Constants


/** \defgroup dimensioned_constants Constants with Dimensions 
  * \ingroup astroconstants
  * Constants with dimensions are provided with a template argument specifying the units system to use, see \ref astrounits.
  * @{
  */


///Newton's Gravitational Constant
/** Source: 2014 CODATA recommended values \cite codata_2014, http://physics.nist.gov/cuu/Constants/index.html
  */
template<typename units>
constexpr typename units::realT G()
{
   return static_cast<typename units::realT>(6.67408e-11) * (units::length*units::length*units::length)/(units::mass * units::time*units::time);
}

///The speed of light.
/** Source: 2014 CODATA recommended values \cite codata_2014, http://physics.nist.gov/cuu/Constants/index.html
  */
template<typename units>
constexpr typename units::realT c()
{
   return static_cast<typename units::realT>(299792458)* (units::length/units::time);
}

///Boltzmann Constant
/** \f$ k = 1.38064852 \times 10^{-23} \: J \: K^{-1}\f$
  * 
  * Source: 2014 CODATA recommended values \cite codata_2014, http://physics.nist.gov/cuu/Constants/index.html
  */
template<typename units>
constexpr typename units::realT k()
{
   return static_cast<typename units::realT>(1.38064852e-23) * units::energy / units::temperature;
}

///Stefan-Boltzmann Constant
/** \f$ \sigma =  5.670367 \times 10^{-8} \: W \: m^{-2} \: K^{-4} \f$
  * 
  * Source: 2014 CODATA recommended values \cite codata_2014, http://physics.nist.gov/cuu/Constants/index.html
  */
template<typename units>
constexpr typename units::realT sigma()
{
   return static_cast<typename units::realT>(5.670367e-8) * (units::energy/units::time) / (units::length*units::length) / (units::temperature*units::temperature*units::temperature*units::temperature);
}

///Planck Constant
/** Source: 2014 CODATA recommended values \cite codata_2014, http://physics.nist.gov/cuu/Constants/index.html
  */
template<typename units>
constexpr typename units::realT h()
{
   return static_cast<typename units::realT>(6.626070040e-34) * units::energy * units::time;
}

///Length of day 
/** As defined by IAU
  */
template<typename units>
constexpr typename units::realT day()
{
   return static_cast<typename units::realT>(86400.0)*units::time;
}

///Length of year
/** As defined by IAU
  */
template<typename units>
constexpr typename units::realT year()
{
   return (static_cast<typename units::realT>(365.25)*static_cast<typename units::realT>(86400.0))*units::time;
}

///Astronomical Unit
/** As defined by IAU Resolution 2012 B2
  */
template<typename units>
constexpr typename units::realT au()
{
   return static_cast<typename units::realT>(149597870700.0)*units::length;
}

///The parsec
template<typename units>
constexpr typename units::realT parsec()
{
   return au<units>()/tan_arcsec<typename units::realT>();
}


///Radius of the Sun 
/** \f$ R_{sun} = 6.957\times 10^{8} \: m\f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT radSun()
{
   return static_cast<typename units::realT>(6.957e8)*units::length;
}

///Solar Irradiance at 1 au
/** \f$ S_{sun} = 1361 \: W \: m^{-2}\f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT solarIrrad()
{
   return static_cast<typename units::realT>(1361)*units::energy/units::time/(units::length*units::length);
}

///Luminosity of the Sun 
/** \f$ L_{sun} = 3.828\times 10^{26} \: W\f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT lumSun()
{
   return static_cast<typename units::realT>(3.828e26)*units::energy/units::time;
}

///Effective Temperature of the Sun 
/** \f$ T_{sun} = 5772 \: K \f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT TeffSun()
{
   return static_cast<typename units::realT>(5772.0)*units::temperature;
}

/// Solar Mass Parameter
/** \f$ GM_{sun} =  1.3271244 \times 10^20 \: m^3 \: s^{-2} \f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT GMSun()
{
   return static_cast<typename units::realT>(1.3271244e20)*(units::length*units::length*units::length)/(units::time*units::time);
}

/// Solar Mass 
/** The mass of Sun is the \ref GMSun() "Solar mass parameter" divided by the \ref G() "graviational constant".
  * 
  * \f$ M_{Sun} =  GM_{Sun}/G \f$
  * 
  */
template<typename units>
constexpr typename units::realT massSun()
{
   return GMSun<units>()/G<units>();
}

/// Radius of Earth (nominal equatorial)
/** \f$ R_{e} =  6.3781\times 10^6 \: m \f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT radEarth()
{
   return static_cast<typename units::realT>(6.3781e6)*units::length;
}

/// Earth Mass Parameter
/** \f$ GM_{e} =  3.986004 \times 10^14 \: m^3 \: s^{-2} \f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT GMEarth()
{
   return static_cast<typename units::realT>(3.986004e14)*(units::length*units::length*units::length)/(units::time*units::time);
}

/// Earth Mass
/** The mass of Earth is the \ref GMEarth() "mass parameter"  divided by the \ref G() "gravitational constant".
  * 
  * \f$ M_{e} =  GM_{e}/G \f$
  *
  */
template<typename units>
constexpr typename units::realT massEarth()
{
   return GMEarth<units>()/G<units>();
}

/// Radius of Jupiter (nominal equatorial)
/** \f$ R_{J} =  7.1492\times 10^7 \: m \f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT radJupiter()
{
   return static_cast<typename units::realT>(7.1492e7)*units::length;
}

/// Jupiter Mass Parameter
/** \f$ GM_{J} =  1.2668653 \times 10^17 \: m^3 \: s^{-2} \f$
  *
  * Source: IAU Resolution 2015 B2, <a href="https://www.iau.org/static/resolutions/IAU2015_English.pdf">https://www.iau.org/static/resolutions/IAU2015_English.pdf</a> 
  */
template<typename units>
constexpr typename units::realT GMJupiter()
{
   return static_cast<typename units::realT>(1.2668653e17)*(units::length*units::length*units::length)/(units::time*units::time);
}

/// Jupiter Mass
/** The mass of Jupiter is the \ref GMJupiter() "Jupiter mass parameter"  divided by the \ref G()  "gravitational constant".
  * 
  * \f$ M_{J} =  GM_{J}/G \f$
  *
  */
template<typename units>
constexpr typename units::realT massJupiter()
{
   return GMJupiter<units>()/G<units>();
}



/// Mass ratio of Mercury to the Sun.
/** 
  * \f$ M_{Sun}/M_{Me} = 6.023657330 \times 10^{6} \f$
  * 
  * Source: IAU Resolution 2009 B2, <a href="https://www.iau.org/static/resolutions/IAU2009_English.pdf">https://www.iau.org/static/resolutions/IAU2009_English.pdf</a>.
  * see also <a href="http://maia.usno.navy.mil/NSFA/NSFA_cbe.html">http://maia.usno.navy.mil/NSFA/NSFA_cbe.html</a>
  * 
  */
template<typename units>
constexpr typename units::realT mrMercury()
{
   return static_cast<typename units::realT>(6.023657330e6);
}

/// Mass of Mercury.
/** 
  * The mass of Mercury is the \ref massSun() "mass of the Sun" divided by the \ref mrMercury() "mass ratio of Mercury to the Sun".
  * 
  */
template<typename units>
constexpr typename units::realT massMercury()
{
   return massSun<units>()/mrMercury<units>();
}

/// Radius of Mercury.
/** 
  * \f$ = 2.4397 \times 10^{6}
  * 
  * Source: Seidelmann et al. (2007) \cite seidelmann_2007
  * 
  */
template<typename units>
constexpr typename units::realT radMercury()
{
   return static_cast<typename units::realT>(2.4397e6)*units::length;
}

/// Mass ratio of Venus to the Sun.
/** 
  * \f$ M_{Sun}/M_{Ve} = 4.08523719 \times 10^{5}\f$
  * 
  * Source: IAU Resolution 2009 B2, <a href="https://www.iau.org/static/resolutions/IAU2009_English.pdf">https://www.iau.org/static/resolutions/IAU2009_English.pdf</a>.
  * see also <a href="http://maia.usno.navy.mil/NSFA/NSFA_cbe.html">http://maia.usno.navy.mil/NSFA/NSFA_cbe.html</a>
  * 
  */
template<typename units>
constexpr typename units::realT mrVenus()
{
   return static_cast<typename units::realT>(4.08523719e5);
}

/// Mass of Venus.
/** 
  * The mass of Venus is the \ref massSun() "mass of the Sun" divided by the \ref mrVenus() "mass ratio of Venus to the Sun".
  * 
  */
template<typename units>
constexpr typename units::realT massVenus()
{
   return massSun<units>()/mrVenus<units>();
}

/// Radius of Venus.
/** 
  * \f$ = 6.0518 \times 10^{6}
  * 
  * Source: Seidelmann et al. (2007) \cite seidelmann_2007
  * 
  */
template<typename units>
constexpr typename units::realT radVenus()
{
   return static_cast<typename units::realT>(6.0518e6)*units::length;
}

/// Mass ratio of Mars to the Sun.
/** 
  * \f$ M_{Sun}/M_{Ma} = 3.09870359 \times 10^{6}\f$
  * 
  * Source: IAU Resolution 2009 B2, <a href="https://www.iau.org/static/resolutions/IAU2009_English.pdf">https://www.iau.org/static/resolutions/IAU2009_English.pdf</a>.
  * see also <a href="http://maia.usno.navy.mil/NSFA/NSFA_cbe.html">http://maia.usno.navy.mil/NSFA/NSFA_cbe.html</a>
  * 
  */
template<typename units>
constexpr typename units::realT mrMars()
{
   return static_cast<typename units::realT>( 3.09870359e6);
}

/// Mass of Mars.
/** 
  * The mass of Mars is the mass of the Sun divided by the mass ratio of Mars to the Sun.
  * 
  */
template<typename units>
constexpr typename units::realT massMars()
{
   return massSun<units>()/mrMars<units>();
}

/// Radius of Mars.
/** 
  * \f$ = 3.39619 \times 10^{6}
  * 
  * Source: Seidelmann et al. (2007) \cite seidelmann_2007
  * 
  */
template<typename units>
constexpr typename units::realT radMars()
{
   return static_cast<typename units::realT>(3.39619e6)*units::length;
}

/// Mass ratio of Saturn to the Sun.
/** 
  * \f$ M_{Sun}/M_{Sa} = 3.4979018 \times 10^{3}\f$
  * 
  * Source: IAU Resolution 2009 B2, <a href="https://www.iau.org/static/resolutions/IAU2009_English.pdf">https://www.iau.org/static/resolutions/IAU2009_English.pdf</a>.
  * see also <a href="http://maia.usno.navy.mil/NSFA/NSFA_cbe.html">http://maia.usno.navy.mil/NSFA/NSFA_cbe.html</a>
  * 
  */
template<typename units>
constexpr typename units::realT mrSaturn()
{
   return static_cast<typename units::realT>( 3.4979018e3);
}

/// Mass of Saturn.
/** 
  * The mass of Saturn is the mass of the Sun divided by the mass ratio of Saturn to the Sun.
  * 
  */
template<typename units>
constexpr typename units::realT massSaturn()
{
   return massSun<units>()/mrSaturn<units>();
}

/// Radius of Saturn.
/** 
  * \f$ = 6.0268 \times 10^{7}
  * 
  * Source: Seidelmann et al. (2007) \cite seidelmann_2007
  * 
  */
template<typename units>
constexpr typename units::realT radSaturn()
{
   return static_cast<typename units::realT>(6.0268e7)*units::length;
}


/// Mass ratio of Uranus to the Sun.
/** 
  * \f$ M_{Sun}/M_{U} = 2.2902951 \times 10^{4}\f$
  * 
  * Source: IAU Resolution 2009 B2, <a href="https://www.iau.org/static/resolutions/IAU2009_English.pdf">https://www.iau.org/static/resolutions/IAU2009_English.pdf</a>.
  * see also <a href="http://maia.usno.navy.mil/NSFA/NSFA_cbe.html">http://maia.usno.navy.mil/NSFA/NSFA_cbe.html</a>
  * 
  */
template<typename units>
constexpr typename units::realT mrUranus()
{
   return static_cast<typename units::realT>( 2.2902951e4);
}

/// Mass of Uranus.
/** 
  * The mass of Uranus is the mass of the Sun divided by the mass ratio of Uranus to the Sun.
  * 
  */
template<typename units>
constexpr typename units::realT massUranus()
{
   return massSun<units>()/mrUranus<units>();
}

/// Radius of Uranus.
/** 
  * \f$ = 2.5559 \times 10^{7}
  * 
  * Source: Seidelmann et al. (2007) \cite seidelmann_2007
  * 
  */
template<typename units>
constexpr typename units::realT radUranus()
{
   return static_cast<typename units::realT>(2.5559e7)*units::length;
}

/// Mass ratio of Neptune to the Sun.
/** 
  * \f$ M_{Sun}/M_{U} = 1.941226 \times 10^{4}\f$
  * 
  * Source: IAU Resolution 2009 B2, <a href="https://www.iau.org/static/resolutions/IAU2009_English.pdf">https://www.iau.org/static/resolutions/IAU2009_English.pdf</a>.
  * see also <a href="http://maia.usno.navy.mil/NSFA/NSFA_cbe.html">http://maia.usno.navy.mil/NSFA/NSFA_cbe.html</a>
  * 
  */
template<typename units>
constexpr typename units::realT mrNeptune()
{
   return static_cast<typename units::realT>( 1.941226e4);
}

/// Mass of Neptune.
/** 
  * The mass of Neptune is the mass of the Sun divided by the mass ratio of Neptune to the Sun.
  * 
  */
template<typename units>
constexpr typename units::realT massNeptune()
{
   return massSun<units>()/mrNeptune<units>();
}

/// Radius of Neptune.
/** 
  * \f$ = 2.4764 \times 10^{7}
  * 
  * Source: Seidelmann et al. (2007) \cite seidelmann_2007
  * 
  */
template<typename units>
constexpr typename units::realT radNeptune()
{
   return static_cast<typename units::realT>(2.4764e7)*units::length;
}


///@} Dimensioned constants

} //namespace constants
} //namespace astro
} //namespace mx


#endif //__mx_astro_constants_hpp__

