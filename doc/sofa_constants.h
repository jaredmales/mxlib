/** \addtogroup genphyconstants General Constants
  * @{
  * 
  */ 

/// Pi (defined in the \ref astro_sofa "SOFA" library sofam.h)
#define DPI (3.141592653589793238462643)

/// 2Pi (defined in the \ref astro_sofa "SOFA" library sofam.h)  
#define D2PI (6.283185307179586476925287)

/// Radians to degrees (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DR2D (57.29577951308232087679815)

/// Degrees to radians (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DD2R (1.745329251994329576923691e-2)

/// Radians to arcseconds (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DR2AS (206264.8062470963551564734)

/// Arcseconds to radians (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DAS2R (4.848136811095359935899141e-6)

/// Seconds of time to radians (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DS2R (7.272205216643039903848712e-5)

/// Arcseconds in a full circle (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define TURNAS (1296000.0)

/// Milliarcseconds to radians (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DMAS2R (DAS2R / 1e3)

/// Length of tropical year B1900 (days) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DTY (365.242198781)

/// Seconds per day. (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DAYSEC (86400.0)

/// Days per Julian year (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJY (365.25)

/// Days per Julian century (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJC (36525.0)

/// Days per Julian millennium (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJM (365250.0)

/// Reference epoch (J2000.0), Julian Date (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJ00 (2451545.0)

/// Julian Date of Modified Julian Date zero (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJM0 (2400000.5)

/// Reference epoch (J2000.0), Modified Julian Date (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJM00 (51544.5)

/// 1977 Jan 1.0 as MJD (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DJM77 (43144.0)

/// TT minus TAI (s) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define TTMTAI (32.184)

/// Astronomical unit (m) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DAU (149597870e3)

/// Speed of light (m/s) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define CMPS 299792458.0

/// Light time for 1 au (s) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define AULT 499.004782

/// Speed of light (AU per day) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define DC (DAYSEC / AULT)

/// L_G = 1 - d(TT)/d(TCG) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define ELG (6.969290134e-10)

/// L_B = 1 - d(TDB)/d(TCB) (defined in the \ref astro_sofa "SOFA" library sofam.h)  
#define ELB (1.550519768e-8)

/// TDB (s) at TAI 1977/1/1.0 (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define TDB0 (-6.55e-5)

/// Schwarzschild radius of the Sun (au) 
/**  = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11 (defined in the \ref astro_sofa "SOFA" library sofam.h) 
 */
#define SRS 1.97412574336e-8

/// dint(A) - truncate to nearest whole number towards zero (double) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define dint(A) ((A)<0.0?ceil(A):floor(A))

/// dnint(A) - round to nearest whole number (double) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define dnint(A) ((A)<0.0?ceil((A)-0.5):floor((A)+0.5))

/// dsign(A,B) - magnitude of A with sign of B (double) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define dsign(A,B) ((B)<0.0?-fabs(A):fabs(A))

/// max(A,B) - larger (most +ve) of two numbers (generic) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define gmax(A,B) (((A)>(B))?(A):(B))

/// min(A,B) - smaller (least +ve) of two numbers (generic) (defined in the \ref astro_sofa "SOFA" library sofam.h) 
#define gmin(A,B) (((A)<(B))?(A):(B))

///@}
