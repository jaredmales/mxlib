/** \file astrodyn.cpp
 * \author Jared R. Males
 * \brief Definitions of some astrodynamic utilities
 *
 */

#include "astrodyn.hpp"


static char UTCSTR[4] = "UTC";

double get_JD(double Y, double M, double D)
{
   double A, B;
   
   if(M<3)
   {
      Y-=1;
      M+=12;
   }
   
   A = floor(Y/100);
   B = 2.0 - A + floor(A/4);
   
   return floor(365.25*(Y+4716.0))+floor(30.6001*(M+1.0)) + D + B - 1524.5;
}

int get_hrminsec_fm_day(int &Dy, int &hr, int &min, double &sec, double day)
{
   
   Dy = day;
   hr = (day - (double) Dy)*24.;
   min = (day - (double) Dy - ((double)hr)/24.) * 60.*24.;
   sec = (day - (double) Dy - ((double)hr)/24. - ((double) min)/(60.*24.)) * 3600.*24.;
   
   return 0;
}

int get_GMST(double &GMST, int Yr, int Mo, int Dy, int hr, int min, double sec)
{
   int rv;
   double utc0, utc1, tai0, tai1, tt0, tt1, ut10, ut11;
   
   //Get UTC time in SOFA 2-double format
   rv = iauDtf2d ( UTCSTR, Yr, Mo, Dy, hr, min, sec, &utc0, &utc1);
   
   //Convert UTC to TAI
   rv = iauUtctai(utc0, utc1, &tai0, &tai1);
   
   //Convert TAI to TT
   rv = iauTaitt (tai0, tai1, &tt0, &tt1);
   
   //Convert UTC to UT1
   rv = iauUtcut1(utc0, utc1, -.4, &ut10, &ut11); //The current DUT1 - how to keep updated?
   
   GMST = iauGmst06(ut10, ut11, tt0, tt1)/(2.*DPI)*24.;
   
   return rv;
}

int get_GMST(double &GMST, int Yr, int Mo, double day)
{
   int Dy, hr, min;
   double sec;
   
   get_hrminsec_fm_day(Dy, hr, min, sec, day);
   return get_GMST(GMST, Yr, Mo, Dy, hr, min, sec);
}

int get_GAST(double &GAST, int Yr, int Mo, int Dy, int hr, int min, double sec)
{
   int rv;
   double utc0, utc1, tai0, tai1, tt0, tt1, ut10, ut11;
   
   //Get UTC time in SOFA 2-double format
   rv = iauDtf2d (UTCSTR, Yr, Mo, Dy, hr, min, sec, &utc0, &utc1);
   
   //Convert UTC to TAI
   rv = iauUtctai(utc0, utc1, &tai0, &tai1);
   
   //Convert TAI to TT
   rv = iauTaitt (tai0, tai1, &tt0, &tt1);
   
   //Convert UTC to UT1
   rv = iauUtcut1(utc0, utc1, -.4, &ut10, &ut11); //The current DUT1 - how to keep updated?
   
   GAST = iauGst06a(utc0, utc1, tt0, tt1)/(2.*DPI)*24.;
   
   return rv;
}

int get_GAST(double &GAST, int Yr, int Mo, double day)
{
   int Dy, hr, min;
   double sec;
   
   get_hrminsec_fm_day(Dy, hr, min, sec, day);
   return get_GAST(GAST, Yr, Mo, Dy, hr, min, sec);
}

/*double get_GMST(double Y, double M, double D)
 * {
 *   double JD0h, GMST0h, GMST;
 *   double JD = get_JD(Y, M, D);
 *   JD0h = (double)(floor(JD))-0.5;
 * 
 *   //double d = JD0h - 2451545.0;
 *   double T = (JD - 2451545.0) / 36525.0;
 * 
 *   GMST0h = 24110.54841 + 8640184.812866 * T + 0.093104 * T*T - 0.0000062 * T*T*T;
 * 
 *   GMST = GMST0h + (JD - JD0h)*1.00273790935*86400.0;
 * 
 *   GMST = fmod(GMST, 86400);
 *   //if(GMST >= 86400.0) GMST -= 86400.0;
 * 
 *   return GMST;
 }*/

int get_LMST(double &LMST, int Yr, int Mo, int Dy, int hr, int min, double sec, double lon)
{
   int rv;
   double GMST;
   rv =  get_GMST(GMST, Yr,Mo, Dy, hr, min, sec);
   
   LMST =  GMST + lon/15.0;
   
   LMST = fmod(LMST, 24);//23.01669);
   if(LMST < 0.0) LMST += 24.; //23.01669;

   return rv;
}

int get_LMST(double &LMST, int Yr, int Mo, double day, double lon)
{
   int Dy, hr, min;
   double sec;
   
   get_hrminsec_fm_day(Dy, hr, min, sec, day);
   return get_LMST(LMST, Yr, Mo, Dy, hr, min, sec, lon);
}

int get_LAST(double &LAST, int Yr, int Mo, int Dy, int hr, int min, double sec, double lon)
{
   
   int rv;
   double GAST;
   rv = get_GAST( GAST, Yr, Mo, Dy, hr, min, sec);
   
   LAST =  GAST + lon/15.0;
   
   LAST = fmod(LAST, 24.);//23.01669);
   if(LAST < 0.0) LAST += 24.;//23.01669;

   return rv;
}

int get_LAST(double &LAST, int Yr, int Mo, double day, double lon)
{
   int Dy, hr, min;
   double sec;
   
   get_hrminsec_fm_day(Dy, hr, min, sec, day);
   return get_LAST(LAST, Yr, Mo, Dy, hr, min, sec, lon);
}

void calc_AZ_EL(double *az, double *el, double ha, double dec, double lat)
{
   double ha_rad, dec_rad, lat_rad, az_rad, el_rad;
   
   //printf("%f\n", get_ha()*15.0);
   ha_rad = ha*15.0*DPI/180.0;
   dec_rad = dec*DPI/180.0;
   lat_rad = lat*DPI/180.0;
   
   az_rad = atan2(sin(ha_rad), cos(ha_rad)*sin(lat_rad)-tan(dec_rad)*cos(lat_rad));
   
   el_rad = asin(sin(lat_rad)*sin(dec_rad) + cos(lat_rad)*cos(dec_rad)*cos(ha_rad));
   
   *az = az_rad*180.0/DPI + 180.0;
   *el = el_rad*180.0/DPI;
   
   return;
}

void azel_to_hadec(double &ha, double &dec, double az, double el, double lat)
{
   dec = asin( sin(el)*sin(lat)+cos(el)*cos(az)*cos(lat) );
   
   ha = atan2(-sin(az)*cos(el)/cos(dec),  (sin(el)-sin(dec)*sin(lat))/(cos(dec)*cos(lat)));
   
   //atan2(-cos(el)*sin(az)/cos(dec), sin(el)*cos(lat)+cos(el)*cos(az)*sin(lat) );
}

// DEC = asind(sind(El).*sind(lat)+cosd(El).*cosd(lat).*cosd(Az));
// LHA = atan2(-sind(Az).*cosd(El)./cosd(DEC), ...
//     (sind(El)-sind(DEC).*sind(lat))./(cosd(DEC).*cosd(lat))).*(180/pi);



void azel_to_hadec(mx::Vectord &ha, mx::Vectord &dec, const mx::Vectord &az, const mx::Vectord &el, double lat)
{
   double tha, tdec;
   
   for(size_t i=0; i< az.length(0); i++)
   {
      azel_to_hadec(tha, tdec, az(i), el(i), lat);
      ha(i) = tha;
      dec(i) = tdec;
   }
}



double get_ParAng_deg(double lat, double dec, double ha)
{
   return RTOD(atan2(cos(DTOR(lat))*sin(DTOR(ha)), sin(DTOR(lat))*cos(DTOR(dec)) - cos(DTOR(lat))*sin(DTOR(dec))*cos(DTOR(ha))));
   
   
   //RTOD( atan(-sin(d2r*had), cos(d2r*dec)*tan(d2r*latitude)-sin(d2r*dec)*cos(d2r*had))
}

double get_ParAng_rad(double lat, double dec, double ha)
{
   return atan2(cos(lat)*sin(ha), sin(lat)*cos(dec) - cos(lat)*sin(dec)*cos(ha));
}


int latlon_to_ECI(double &x, double &y, double &z, double lat, double lon, double alt, double lst)
{

   z = (RAD_EARTH + alt) * sin(lat);
   double R = (RAD_EARTH + alt) * cos(lat);

   x = R*cos(lst);
   y = R*sin(lst);
   
   return 0;
}



int latlon_to_ECI(mx::Vectord &x, mx::Vectord &y, double &z, double lat, double lon, double alt, const mx::Vectord &lst)
{
   
   z = (RAD_EARTH + alt) * sin(lat);
   double R = (RAD_EARTH + alt) * cos(lat);

   for(size_t i=0; i< lst.length(0); i++)
   {
      x(i) = R*cos(lst(i));
      y(i) = R*sin(lst(i));
   }
   
   return 0;
}


int ECI_to_TCH( mx::Vectord & az,
                mx::Vectord & el,
                mx::Vectord & r,
                const double lat,
                const mx::Vectord & lst,
                const mx::Vectord & obs_x,
                const mx::Vectord & obs_y,
                const double obs_z,
                const mx::Vectord & tgt_x,
                const mx::Vectord & tgt_y,
                const mx::Vectord & tgt_z )
{
   double cos_lat = cos(lat);
   double sin_lat = sin(lat);
   
   double dx, dy, dz, rs, re, rZ;
   
   for( size_t i=0; i< lst.length(0); i++)
   {
      dx = tgt_x(i) - obs_x(i);
      dy = tgt_y(i) - obs_y(i);
      dz = tgt_z(i) - obs_z;
      
      rs = sin_lat*cos(lst(i))*dx + sin_lat*sin(lst(i))*dy - cos_lat*dz;
      re = -sin(lst(i))*dx + cos(lst(i))*dy;
      rZ = cos_lat*cos(lst(i))*dx + cos_lat*sin(lst(i))*dy + sin_lat*dz;
      
      r(i) = sqrt(rs*rs + re*re + rZ*rZ);
      az(i) = atan2(-re, rs);
      el(i) = asin(rZ/r(i));
   }
   
   return 0;
}
   
   





