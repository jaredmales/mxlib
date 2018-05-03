/** \file astroDynamics.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Various utilities for astrodynamics.
  * \ingroup astrofiles
  * 
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//


#ifndef astroDynamics_hpp
#define astroDynamics_hpp

#include "sofa.hpp"

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include "../math/geo.hpp"

namespace mx
{
namespace astro 
{



///Breaks down a decimal day into hours, minutes and decimal point seconds.
/** Assumes 86400 seconds per day, ignoring leap seconds.
  * 
  * \tparam realT the floating point type for calculations
  * 
  * \retval 0 on success
  * \retval -1 on failure.
  */
template<typename realT>
int hrMinSecFromDay( int & Dy,    ///< [out] the integer day
                     int & hr,    ///< [out] the integer hour
                     int & min,   ///< [out] the integer minute
                     realT & sec, ///< [out] the integer second
                     realT day    ///< [in] the decimal day to be broken down
                   )
{
   Dy = (int) day;
   hr = (day - (realT) Dy)*24.0;
   min = (day - (realT) Dy - ((realT)hr)/24.0) * 60.*24.0;
   sec = (day - (realT) Dy - ((realT)hr)/24.0 - ((realT) min)/(60.0*24.0)) * 3600.0*24.0;
   
   return 0;
}


///Returns Greenwich Mean Sidereal Time for a given UTC time.
/** Uses the SOFA 2006 GMST function, and a (for now) hardcoded value of dUT1
  *
  * \retval 0 on success
  */
template<typename realT>
int getGMST( realT & GMST, ///< [out] GMST in hours
             int Yr,       ///< [in] integer year
             int Mo,       ///< [in] integer month (1-12)
             int Dy,       ///< [in] integer day (1-31)
             int hr,       ///< [in] integer hour (0-23)
             int min,      ///< [in] integer minute (0-59) 
             realT sec     ///< [in] decimal second (0-59)
           )
{
   int rv;
   double utc0, utc1, tai0, tai1, tt0, tt1, ut10, ut11;
   
   //Get UTC time in SOFA 2-double format
   rv = sofa::iauDtf2d ( "UTC", Yr, Mo, Dy, hr, min, static_cast<double>(sec), &utc0, &utc1);
   
   //Convert UTC to TAI
   rv = sofa::iauUtctai(utc0, utc1, &tai0, &tai1);
   
   //Convert TAI to TT
   rv = sofa::iauTaitt (tai0, tai1, &tt0, &tt1);
   
   //Convert UTC to UT1
   rv = sofa::iauUtcut1(utc0, utc1, -.4, &ut10, &ut11); //The current DUT1 - how to keep updated?
   
   GMST = sofa::iauGmst06(ut10, ut11, tt0, tt1)/(two_pi<realT>())*static_cast<realT>(24);
   
   return rv;
}

///Returns Greenwich Mean Sidereal Time for a given UTC time.
/** Uses the SOFA 2006 GMST function, and a (for now) hardcoded value of dUT1
 * \param GMST [output] GMST in hours
 * \param Yr [input] integer year
 * \param Mo [input] integer month (1-12)
 * \param day [input] decimal day (1-31)
 * \retval 0 on sucess
 */
template<typename realT>
int getGMST( realT & GMST, ///< [out] GMST in hours
             int Yr,       ///< [in] integer year
             int Mo,       ///< [in] integer month (1-12)
             realT day     ///< [in] decimal day (1-31)
           )
{
   int Dy, hr, min;
   double sec;
   
   hrMinSecFromDay<realT>(Dy, hr, min, sec, day);
   return getGMST(GMST, Yr, Mo, Dy, hr, min, sec);
}

///Returns Local Mean Sidereal Time for a given UTC time and longitude.
/** Uses the SOFA 2006 GMST function, and a (for now) hardcoded value of dUT1
  *
  * \retval 0 on success
  * 
  */
template<typename realT>
int getLMST( realT & LMST, ///< [out] LMST in hours
             int Yr,       ///< [in] integer year
             int Mo,       ///< [in] integer month (1-12) 
             int Dy,       ///< [in] integer day (1-31)
             int hr,       ///< [in] integer hour (0-23) 
             int min,      ///< [in] integer minute (0-59) 
             realT sec,    ///< [in] decimal second (0-59) 
             realT lon     ///< [in] longitude in degrees, E+, W-
           )
{
   int rv;
   realT GMST;
   
   rv =  getGMST(GMST, Yr,Mo, Dy, hr, min, sec);
   
   LMST =  GMST + lon/static_cast<realT>(15);
   
   LMST = fmod(LMST, 24);
   if(LMST < 0.0) LMST += 24.; 

   return rv;
}

///Returns Local Mean Sidereal Time for a given UTC time.
/** Uses the SOFA 2006 GMST function, and a (for now) hardcoded value of dUT1
  *
  * \tparam realT the floating point type for all calculations
  * 
  * \retval 0 on sucess
  */
template<typename realT>
int getLMST( realT & LMST, ///< [out] LMST in hours
             int Yr,       ///< [in] integer year 
             int Mo,       ///< [in] integer month (1-12)
             realT day,    ///< [in] decimal day (1-31)
             realT lon     ///< [in] longitude in degrees, E+, W-
           )
{
   int Dy, hr, min;
   realT sec;
   
   hrMinSecFromDay<realT>(Dy, hr, min, sec, day);
   
   return getLMST(LMST, Yr, Mo, Dy, hr, min, sec, lon);
}


///Calculates the Azimuth and Elevation for a given Hour-Angle, Declination, and Latitude.
/** References: J. Meeus "Astronomical Algorithms", 1991; V. Pisacane "Fundamentals of Space Systems", 2nd ed., 2005.
 * 
 * \tparam realT the real floating point type for calculations
 * 
 * \returns 0 on success
 * \returns -1 on error
 */
template<typename realT>
int calcAzEl( realT & az, ///< [out] the calculated azimuth [degrees]
              realT & el, ///< [out] the calculate elevation [degrees]
              realT ha,   ///< [in] the hour ange [degrees]
              realT dec,  ///< [in] the declination [degrees]
              realT lat   ///< [in] the latitude [degrees]
            )
{
   realT ha_rad, dec_rad, lat_rad, az_rad, el_rad;
   
   ha_rad = ha*static_cast<realT>(15)*pi<realT>()/static_cast<realT>(180);
   dec_rad = dec*pi<realT>()/static_cast<realT>(180);
   lat_rad = lat*pi<realT>()/static_cast<realT>(180);
   
   az_rad = atan2(sin(ha_rad), cos(ha_rad)*sin(lat_rad)-tan(dec_rad)*cos(lat_rad));
   
   el_rad = asin(sin(lat_rad)*sin(dec_rad) + cos(lat_rad)*cos(dec_rad)*cos(ha_rad));
   
   az = az_rad*static_cast<realT>(180)/pi<realT>() + static_cast<realT>(180);
   el = el_rad*static_cast<realT>(180)/pi<realT>();
   
   return 0 ;
}

///Calculate the Parallactic angle, with angles in degrees
/** \param lat is the observer latitude
 * \param dec is the object declination
 * \param ha is the hour angle, in degrees, positive to the east
 */
template<typename realT>
realT parAngDeg( realT lat, 
                      realT dec, 
                      realT ha
                  )
{
   return math::rtod(atan2(cos(math::dtor(lat))*sin(math::dtor(ha)), sin(math::dtor(lat))*cos(math::dtor(dec)) - cos(math::dtor(lat))*sin(math::dtor(dec))*cos(math::dtor(ha))));
}


} //namespace astro 
} //namespace mx

#endif //astroDynamics_hpp

