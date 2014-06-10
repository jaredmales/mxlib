//constants.h

#ifndef __ASTROCONSTANTS_H__
#define __ASTROCONSTANTS_H__


#include <sofam.h>

//Physical Constants


/** \addtogroup genphyconstants General Constants
  * @{
  */
  
///Square-root of 2
/** This was calculated with long double precision
 */
#define SQRT_2 (1.414213562373095048763788073031832936976570636034)



///Newtons Gravitational Constant, units = m^3 / (kg s^2)
/** source = Report of the IAU WGAS Sub-group on Numerical Standards (1994 Values)
  */
#define GRAV_G (6.67259e-11)

/// @}


/** \addtogroup astroconstants 
  */
/// @{



   
///Astronomical Unit, units = m
/** source = NASA/JPL Ephemeris DE-405
  */
#define AU_M 149597870691.0

///Astronomical Unit, units = km
#define AU_KM (AU_M / 1000.0)


///Parsec, units=m
#define PC_M (AU_M / DAS2R) 

///Parsec, units=km
#define PC_KM (AU_KM / DAS2R) 

///Light-year, units=m
#define LY_M (CMPS*DAYSEC*DJY)






///Solar Gravitational Constant, units = m^3/s^2
/** source = Report of the IAU WGAS Sub-group on Numerical Standards (1994 Values)
  */
#define GM_SOL 1.327124e+20

///Solar Gravitational Constant, units = AU^3/day^2
#define GM_SOL_AUD (0.01720209895*0.01720209895)

//Solar System Mass Ratios
//units = unitless
//source = Report of the IAU WGAS Sub-group on Numerical Standards (1994 Values)
#define MR_SOL     1.0
#define MR_MERCURY 6023600.0
#define MR_VENUS   408523.71
#define MR_EMB     328900.56
#define MR_EL      81.30059
#define MR_EARTH   332946.05
#define MR_LUNA    (MR_EARTH*MR_EL)
#define MR_MARS    3098708.0
#define MR_JUPITER 1047.3486
#define MR_SATURN  3497.898
#define MR_URANUS  22902.98
#define MR_NEPTUNE 19412.24
//Asteroid source: JPL Horizon output of GM
#define MR_CERES   (GM_SOL/(63.2*1e+9))
#define MR_PALLAS  (GM_SOL/(14.3*1e+9))
#define MR_VESTA   (GM_SOL/(17.8*1e+9))
#define MR_HYGIEA  (GM_SOL/(7.0*1e+9))

//Solar System Masses
//units = kg
//source = derived
#define MASS_SOL     (GM_SOL/GRAV_G)
#define MASS_MERCURY (MASS_SOL / MR_MERCURY)
#define MASS_VENUS   (MASS_SOL / MR_VENUS)
#define MASS_EARTH   (MASS_SOL / MR_EARTH)
#define MASS_LUNA    (MASS_SOL / MR_LUNA)
#define MASS_MARS    (MASS_SOL / MR_MARS)
#define MASS_JUPITER (MASS_SOL / MR_JUPITER)
#define MASS_SATURN  (MASS_SOL / MR_SATURN)
#define MASS_URANUS  (MASS_SOL / MR_URANUS)
#define MASS_NEPTUNE (MASS_SOL / MR_NEPTUNE)

#define MASS_CERES   (MASS_SOL / MR_CERES)
#define MASS_PALLAS  (MASS_SOL / MR_PALLAS)
#define MASS_VESTA   (MASS_SOL / MR_VESTA)
#define MASS_HYGIEA  (MASS_SOL / MR_HYGIEA)

//Solar System Dimensions
///Radius of the earth, units = m
/** source = WGS84
  */
#define RAD_EARTH (6378137.0)

///Radius of the Sun, units = m
/** source = IAU 2009
 */
#define RAD_SUN (696000000.0)

///Schwarzschild radius of the Sun, units = m
/** Based on value in AU given in sofam.h
 */
#define SRS_M (SRS * AU_M)


///Radius of Jupiter, units = m
/** source = IAU 2009
 */
#define RAD_JUPITER (71492000.0)


/// @}

#endif //__ASTROCONSTANTS_H__

