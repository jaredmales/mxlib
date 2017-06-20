#ifndef __mx_astro_blackbody_hpp__
#define __mx_astro_blackbody_hpp__

#include <vector>

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include "units.hpp"
#include "constants.hpp"

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
   
   return pow( (L * (1-Ab)/f) /( static_cast<realT>(16.0)*constants::sigma<units>()*pi<realT>()*r*r), static_cast<realT>(0.25));
}

template< typename units, bool freq=false>
struct blackbody
{
   typedef typename units::realT realT;

   realT _temperature;
   realT _radius;
   realT _distance;

   std::vector<realT> _spectrum;

   blackbody()
   {
      _temperature = constants::TeffSun<units>();
      _radius = 0;
      _distance = 0;
   }

   blackbody( realT T,
              realT R = 0,
              realT d = 0
            )
   {
      setParameters(T, R, d);
   }
   
   void setParameters( realT T,
                       realT R,
                       realT d
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
         solidang = pi<realT>()*pow(_radius/_distance,2);
      }


      _spectrum.resize(grid.size());
      
      for(int i=0; i < grid.size(); ++i)
      {
         if(!freq)
         {
            _spectrum[i] = 2*h*pow(c,2)*solidang        / pow(grid[i],5)  / (exp(h*c/(grid[i]*k*_temperature))-1.0);
         }
         else
         {
            _spectrum[i] = 2*h/pow(c,2)*solidang * pow(grid[i],3)  / (exp(h*grid[i]/(k*_temperature))-1.0);
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
