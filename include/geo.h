

#ifndef __geo_h__
#define __geo_h__

#include <sofa.h>

/** \addtogroup geo
  * @{
  */

///Convert from degrees to radians
#define DTOR(q) ((q)*DPI/180.)

///Convert from radians to degrees
#define RTOD(q) ((q)*180./DPI)

///Calculate the difference between two angles, correctly across 0/360.
/** Calculates \f$ dq = q2- q1 \f$, but accounts for crossing 0/360.  This implies
  * that \f$ dq \le 180 \f$.
  * 
  * \param q1 angle to subtract from q2, in degrees.
  * \param q2 angle to subtract q1 from, in degrees.
  *
  * \returns the difference of q2 and q1
  *
  * \tparam arithT is the type in which to do arithmetic
  * 
  */
template<int degrad = 0, typename arithT>
arithT angleDiff(arithT q1, arithT q2)
{ 
   arithT dq = q2-q1;

   double half;
   if(!degrad) half = 180.0;
   if(degrad) half = DPI;
      
   if (abs(dq) > half)
   {
      if(dq < 0)
      {
         dq = dq + 2.*half;
      }
      else
      {
         dq = dq - 2.*half;
      }
   }
   
   return dq;
}


/// @}


#endif //__geo_h__
