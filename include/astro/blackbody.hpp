#ifndef __mx_astro_blackbody_hpp__
#define __mx_astro_blackbody_hpp__

#include <vector>

#include "units.hpp"
#include "constants.hpp"
#include "../math/constants.hpp"

namespace mx
{
namespace astro 
{
   
template<typename units>
typename units::realT equilibriumTemp( typename units::realT L, ///< [in] Stellar luminosity
                                       typename units::realT r, ///< [in] Distance from star
                                       typename units::realT Ab, ///< [in] Bond Albedo 
                                       typename units::realT f  ///< [in] Redistribution factor
                                     )
{
   typedef typename units::realT realT;
   
   return pow( (L * (1-Ab)/f) /( static_cast<realT>(16.0)*constants::sigma<units>()*math::pi<realT>()*r*r), static_cast<realT>(0.25));
}

///The blackbody spectral energy distribution in the mx::astro::astroSpectrum form.
/** You specify the temperature, and optionally the radius and distance of the blackbody. 
  * 
  * \tparam units specifies the units of the spectrum 
  * \tparam freq if true, then the blackbody is calcuated as a frequency distribution.  If false (default), it is a wavelength distribution.
  */
template< typename units, bool freq=false>
struct blackbody
{
   typedef typename units::realT realT; ///<The real floating point type used for calculations

   realT _temperature; ///< The temperature of the blackbody.  Default value is the effective temperature of the Sun.
   realT _radius; ///< The optional radius of the blackbody.
   realT _distance; ///< The optional distance to the blackbody.

   std::vector<realT> _spectrum; ///< The calculated spectral energy distribution.

   ///Default c'tor.
   blackbody()
   {
      _temperature = constants::TeffSun<units>();
      _radius = 0;
      _distance = 0;
   }

   ///Constructor used to initialize parameters.
   blackbody( realT T,  ///< [in] The effective temperature of the blackbody. Units as specified by the units template-parameter.
              realT R = 0, ///< [in] [optional] The radius of the blackbody. Units as specified by the units template-parameter.
              realT d = 0 ///< [in]  [optional] The distance to the blackbody. Units as specified by the units template-parameter.
            )
   {
      setParameters(T, R, d);
   }
   
   ///Set the parameters of the blackbody.
   void setParameters( realT T,  ///< [in] The effective temperature of the blackbody. Units as specified by the units template-parameter.
                       realT R,  ///< [in] The radius of the blackbody. Units as specified by the units template-parameter.
                       realT d   ///< [in] The distance to the blackbody. Units as specified by the units template-parameter.
                     )
   {
      _temperature = T;
      _radius = R;
      _distance = d;
   }

   template<typename gridT>
   int setSpectrum( gridT & grid)
   {
      constexpr realT h = constants::h<units>();
      constexpr realT c = constants::c<units>();
      constexpr realT k = constants::k<units>();
      
      realT solidang=1.0;

      if( _radius != 0 and _distance != 0)
      {
         solidang = math::pi<realT>()*pow(_radius/_distance,2);
      }


      _spectrum.resize(grid.size());
      
      for(int i=0; i < grid.size(); ++i)
      {
         if(!freq)
         {
            _spectrum[i] = 2*h*pow(c,2)*solidang / pow(grid[i],5)  / expm1(h*c/(grid[i]*k*_temperature));// (exp(h*c/(grid[i]*k*_temperature))-1.0);
         }
         else
         {
            _spectrum[i] = 2*h/pow(c,2)*solidang * pow(grid[i],3)  / expm1(exp(h*grid[i]/(k*_temperature)));// (exp(h*grid[i]/(k*_temperature))-1.0);
         }
      }
      return 0;
   }
   
   realT operator[](int i)
   {
      return _spectrum[i];
   }
   
};

} //namespace astro
} //namespace mx
#endif //__mx_astro_blackbody_hpp__
