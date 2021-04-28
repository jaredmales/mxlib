/** \file planets.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Various utilities related to planets.
  * \ingroup astrofiles
  * 
  */


#ifndef __mx_astro_planets_hpp__
#define __mx_astro_planets_hpp__

#include "../math/constants.hpp"

#include "units.hpp"


namespace mx
{
namespace astro
{


/** \ingroup planets
  * @{
  */

///An ad-hoc planetary mass-to-radius relationship (old version)
/** The goal of this function is to provide a radius given an exoplanet mass, for lightly-irradiated exoplanets.  By lightly-irradiated we mean (roughly) planet's at 
  * Mercury's separation or further, scaled for stellar luminosity.
  * Here we make use of the transition from rocky to gas-dominated composition at \f$ 1.6  R_e \f$ identified by Rogers \cite rogers_2015 
  * (see also Marcy et al. (2014) \cite marcy_2014).  Below this radius we assume Earth composition and so
  * \f$ R \propto M^{1/3}\f$.  Above this we scale with a power law matched to the mean radius and mass of Uranus and Neptune, which defines the radius between 
  * \f$ 1.6^3 M_\oplus \f$ and this Uranus-Neptune mean point.  Above this point we use
  * a polynomial fit (in log(M)) to points including the Uranus/Neptune mean, Saturn, Jupiter, and above Jupiter's mass the average points from the 4.5 Gyr 1 AU models from 
  * Fortney et al. (2007) \cite fortney_2007.  Above 3591.1 \f$ M_\oplus \f$ (\f$\sim 11 M_{jup}\f$) we scale as \f$ M^{-1/8} \f$ based on the curve shown in Fortney et al. (2011) \cite fortney_2010.
  * 
  * \f$
    \frac{R}{R_\oplus} = 
    \begin{cases}
      \left(\frac{M}{M_\oplus}\right)^{1/3}, & M < 4.1 M_\oplus \\
      0.62\left(\frac{M}{M_\oplus}\right)^{0.67}, & 4.1 M_\oplus \le M < 15.84 M_\oplus \\
      14.0211 - 44.8414 \log_{10}\left(\frac{M}{M_\oplus}\right) + 53.6554 \log_{10}^2\left(\frac{M}{M_\oplus}\right)
        -25.3289\log_{10}^3\left(\frac{M}{M_\oplus}\right) + 5.4920\log_{10}^4\left(\frac{M}{M_\oplus}\right) - 0.4586 \log_{10}^5\left(\frac{M}{M_\oplus}\right), & 15.84 \le M < 3591.1 M_\oplus \\
      32.03 \left(\frac{M}{M_\oplus}\right)^{-1/8}, & 3591.1 M_\oplus \le M
    \end{cases}
     \f$
  * 
  * \image html planet_mrdiag_2017.06.19.png "The ad-hoc mass-to-radius relationship compared to known planets.  Blue circles are the Solar system.  Black circles indicate expoplanets with Teq below that of a blackbody at Mercury's separation (413 K)."  
  *
  * This function makes use of the units type system (\ref astrounits) so it can be used with Earth masses, Jupiter masses, kg (SI units), etc.
  *
  * \returns the estimated radius of the planet.
  * 
  * \tparam units is the units-type specifying the units of mass.  See \ref astrounits.
  */ 
template<typename units>
typename units::realT planetMass2Radius( typename units::realT mass /**< The mass of the planet. */)
{
   typedef typename units::realT realT;
   
   using namespace mx::astro::constants;

   if( mass <  4.1*massEarth<units>())
   {
      return pow( mass/massEarth<units>(), math::third<realT>())*radEarth<units>();
   }
   else if( mass < static_cast<realT>(15.84)*massEarth<units>() )
   {
      return 0.62*pow( mass/massEarth<units>(), static_cast<realT>(0.67))*radEarth<units>();
   }
   else if( mass < 3591.1*massEarth<units>())
   {
      realT logM = log10(mass/massEarth<units>());
      
      return (static_cast<realT>(14.0211)-static_cast<realT>(44.8414)*logM+static_cast<realT>(53.6554)*pow(logM,2) - 
                                   static_cast<realT>(25.3289)*pow(logM,3) + static_cast<realT>(5.4920)*pow(logM,4) 
                                                                                         - 0.4586*pow(logM,5))*radEarth<units>();
   }
   else
   {
      return static_cast<realT>(32.03)*pow(mass/massEarth<units>(), static_cast<realT>(-0.125))* radEarth<units>();
   }
   

}

///An ad-hoc planetary mass-to-radius relationship.
/** The goal of this function is to provide a radius given an exoplanet mass, for cool, lightly-irradiated exoplanets.  Here cool means \f$T_{eq} < 1000 \f$ K.
  * 
  * We make use of the transition from rocky to gas-dominated composition at \f$ 1.6  R_\oplus \f$ identified by Rogers \cite rogers_2015 
  * (see also Marcy et al. (2014) \cite marcy_2014).  Below this radius we assume Earth composition and so
  * \f$ R \propto M^{1/3}\f$.  
  * 
  * Above \f$ 15 M_\oplus\f$ we use the empirical relationship of Thorngren et al (2019) \cite{Thorngren_2019}.
  * 
  * Between these two cases we scale with a power law matched to the endpoints.
  *
  * Above  \f$ 12 M_{Jup} we scale as \f$ M^{-1/8} \f$ based on the curve shown in Fortney et al. (2011) \cite fortney_2010.
  * 
  * \image html planet_mrdiag_2017.06.19.png "The ad-hoc mass-to-radius relationship compared to known planets.  Blue circles are the Solar system.  Points
  * with error bars are from the NASA exoplanet catalog, selected for \f$T_{eq} < 1000\f$ K
  *
  * This function makes use of the units type system (\ref astrounits) so it can be used with Earth masses, Jupiter masses, kg (SI units), etc.
  *
  * \returns the estimated radius of the planet.
  * 
  * \tparam units is the units-type specifying the units of mass.  See \ref astrounits.
  */ 
template<typename units>
typename units::realT planetMass2RadiusWThorngren( typename units::realT mass /**< The mass of the planet. */)
{
   typedef typename units::realT realT;
   
   using namespace mx::astro::constants;

   if( mass <  4.1*massEarth<units>())
   {
      return pow( mass/massEarth<units>(), math::third<realT>())*radEarth<units>();
   }
   else if( mass < static_cast<realT>(15.00)*massEarth<units>() )
   {
      return 0.6413*pow( mass/massEarth<units>(), static_cast<realT>(0.6477))*radEarth<units>();
   }
   else if( mass < 12*massJupiter<units>())
   {
      realT lM = log10(mass/massJupiter<units>());
     
      return (static_cast<realT>(0.96)+static_cast<realT>(0.21)*lM - static_cast<realT>(0.20)*pow(lM,2))*radJupiter<units>();
   }
   else
   {
      return static_cast<realT>(29.97)*pow(mass/massEarth<units>(), static_cast<realT>(-0.125))* radEarth<units>();
   }
   

}

/// Empirical mass-to-radius relationship for cool EGPs from Thorngren et al. (2019)
/** A polynomial in log(M) fit to cool (T < 1000 K) exoplanets by Thorngren et al. (2019) \cite Thorngren_2019.
  * Only valid for \f$ 15 M_\oplus < M < 12 M_{Jup} \f$
  * 
  * \returns the estimated radius of the planet.
  * 
  * \tparam units is the units-type specifying the units of mass.  See \ref astrounits.
  */ 
template<typename units>
typename units::realT planetMass2RadiusThorngren( typename units::realT mass /**< The mass of the planet. */)
{
   typedef typename units::realT realT;
   
   using namespace mx::astro::constants;
   
   realT lM = log10(mass/massJupiter<units>());
   
   return (0.96+0.21*lM - 0.2*pow(lM,2))*radJupiter<units>();
}

///The planetary mass-to-radius of Fabrycky et al. (2014)..
/** A broken power law for low mass planets from Fabrycky et al. (2014)\cite fabrycky_2014 (see also \cite lissauer_2011}. 
  * 
  * \f$
    \frac{R}{R_e} = 
    \begin{cases}
      \left(\frac{M}{M_e}\right)^{1/3}, & M < 1 M_e \\
      \left(\frac{M}{M_e}\right)^{1/2.06}, & 1 M_e \le M     
      \end{cases}
     \f$
  * 
  * This makes use of the units type system (\ref astrounits) so it can be used with Earth masses, Jupiter masses, kg (SI units), etc.
  *
  * \returns the estimated radius of the planet.
  * 
  * \tparam units is the units-type specifying the units of mass.  See \ref astrounits.
  */ 
template<typename units>
typename units::realT planetMass2RadiusFab2014( typename units::realT mass /**< The mass of the planet. */)
{
   typedef typename units::realT realT;
   
   using namespace mx::astro::constants;
      
   if( mass <  massEarth<units>())
   {
      return pow( mass/massEarth<units>(), math::third<realT>())*radEarth<units>();
   }
   else 
   {
      return pow( mass/massEarth<units>(), static_cast<realT>(1)/static_cast<realT>(2.06))*radEarth<units>();
   }

}


///@}
} //namespace astro
} //namespace mx
#endif
