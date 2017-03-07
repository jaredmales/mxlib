/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __wooferTweeterFilter_hpp__
#define __wooferTweeterFilter_hpp__

#include "wooferTweeterDM.hpp"

namespace mx
{
namespace AO 
{
namespace sim 
{
   

template< typename wooferFilterT, typename tweeterFilterT>
class wooferTweeterFilter
{

public:
   
   typedef typename wooferFilterT::floatT floatT;
     
   typedef std::complex<floatT> complexT;
   
   ///The wavefront data type
   typedef wavefront<floatT> wavefrontT;
   
   ///The command type
   typedef wooferTweeterCommand<floatT> commandT;
   
   
   wooferTweeterFilter();
   
   wooferFilterT woofer;
   tweeterFilterT tweeter;
   
   int wooferNModes;
   int tweeterNModes;
   
public:

   template<typename dmT>
   void initialize(dmT & dm);

   void filterCommands(commandT & filtAmps, commandT & rawAmps);

   int _lowOrders;
   
   void initMeasurements(commandT & filtAmps, commandT & rawAmps)
   {
      filtAmps.wooferVect.resize(1, wooferNModes);
      filtAmps.wooferVect.setZero();
      filtAmps.tweeterVect.resize(1, tweeterNModes);
      filtAmps.tweeterVect.setZero();
      
      rawAmps.wooferVect.resize(1, wooferNModes);
      rawAmps.wooferVect.setZero();
      rawAmps.tweeterVect.resize(1, tweeterNModes);
      rawAmps.tweeterVect.setZero();
   }
   
};


template< typename wooferFilterT, typename tweeterFilterT>
wooferTweeterFilter<wooferFilterT, tweeterFilterT>::wooferTweeterFilter()
{
}




template< typename wooferFilterT, typename tweeterFilterT>
template<typename dmT>
void wooferTweeterFilter<wooferFilterT, tweeterFilterT>::initialize(dmT & dm)
{
   woofer.initialize(dm.woofer);
   wooferNModes = dm.woofer.nModes();
   
   tweeter.initialize(dm.tweeter);
   tweeterNModes = dm.tweeter.nModes();
}



template< typename wooferFilterT, typename tweeterFilterT>
void wooferTweeterFilter<wooferFilterT, tweeterFilterT>::filterCommands(commandT & filtAmps, commandT & rawAmps)
{
   static int called = 0;
   
   if(called  == 51)
   {
      woofer.setGains(0.1);
   }
   
   if(called == 101)
   {
      for(int i=78; i< 1680; ++i)
      tweeter.setLeak(i, 0.01);
   }
   
   woofer.filterCommands(filtAmps.wooferVect, rawAmps.wooferVect);
   tweeter.filterCommands(filtAmps.tweeterVect, rawAmps.tweeterVect);
   
   ++called;
}


} //namespace sim 
} //namespace AO
} //namespace mx 

#endif //__wooferTweeterFilter_hpp__
