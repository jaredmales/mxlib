/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __leakyIntegrator_hpp__
#define __leakyIntegrator_hpp__

namespace mx
{
namespace AO 
{
namespace sim 
{
   
template<class _realT>
struct wfMeasurement
{
   typedef _realT realT;
   
   realT iterNo;
   
   typedef Eigen::Array< _realT, Eigen::Dynamic, Eigen::Dynamic> commandT;
   
   commandT measurement;
};

template<class _realT>
class leakyIntegrator
{

public:
   
   ///The real data type
   typedef _realT realT;
        
   ///The wavefront data type
   typedef wavefront<realT> wavefrontT;
   
   ///The pupil image type
   //typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> commandT;
   typedef wfMeasurement<realT> commandT;
   
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;
   
   leakyIntegrator();

//protected:

   imageT _gains;   
   imageT _leaks;

   imageT _commands;

   bool _openLoop; ///If true, then commands are not integrated.
   
public:

   template<typename dmT>
   void initialize(dmT & dm);

   realT gain(int i)
   {
      return _gains(i);
   }
   
   void setGain(int i, realT g);
   void setGains(realT g);

   void setGains(const std::string & ogainf);
   
   void setLeak(int i, realT l);
   void setLeaks(realT l);


   void filterCommands(int iterNo,
                       commandT & filtAmps, 
                       commandT & rawAmps);

   int _closingDelay;
   int _lowOrders;
   realT _closingGainFact;
   
   int _nModes;
   
   void initMeasurements(commandT & filtAmps, commandT & rawAmps)
   {
      filtAmps.measurement.resize(1, _nModes);
      filtAmps.measurement.setZero();
      
      rawAmps.measurement.resize(1, _nModes);
      rawAmps.measurement.setZero();
      
   }
};


template<class _realT>
leakyIntegrator<_realT>::leakyIntegrator()
{
   _lowOrders = 0;
   
   _openLoop = false;
}



template<class _realT>
template<typename dmT>
void leakyIntegrator<_realT>::initialize(dmT & dm)
{
   _nModes = dm.nModes();
   
   _gains.resize(1,dm.nModes());
   _leaks.resize(1, dm.nModes());
   _commands.resize(1,dm.nModes());
   
   
   _gains.setZero();
   
   _leaks.setZero();
   
   _commands.setZero();
   

}

template<class _realT>
void leakyIntegrator<_realT>::setGain(int i, _realT g)
{
   _gains(0,i) = g;
}

template<class _realT>
void leakyIntegrator<_realT>::setGains(_realT g)
{
   for(int i=0;i<_gains.cols(); ++i) _gains(0,i) = g;
}

template<class _realT>
void leakyIntegrator<_realT>::setGains(const std::string & ogainf)
{
   std::ifstream fin;
   fin.open(ogainf);
   
   if(!fin.good())
   {
      std::cerr << "leakyIntegrator: gain file " << ogainf << " not found.\n";
      exit(-1);
   }
   
   std::string tmpstr;
   for(int i=0;i<_gains.cols();++i)
   {
      fin >> tmpstr;
      fin >> tmpstr;
      fin >> tmpstr;
      fin >> tmpstr;
      
      setGain(i, mx::convertFromString<float>(tmpstr));

   }  
   
}

template<class _realT>
void leakyIntegrator<_realT>::setLeak(int i, _realT l)
{
   _leaks(0,i) = l;
}

template<class _realT>
void leakyIntegrator<_realT>::setLeaks(_realT l)
{
   for(int i=0;i<_leaks.cols(); ++i) _leaks(0,i) = l;
}

template<class realT>
void leakyIntegrator<realT>::filterCommands( int iterNo,
                                               commandT & filtAmps, 
                                               commandT & rawAmps )
{
   
   
   filtAmps.iterNo = rawAmps.iterNo;
   
   if(_openLoop) 
   {
      filtAmps.measurement.setZero();
      return;
   }
   
   int topN = rawAmps.measurement.cols();
   
   
   
   realT gainF = 1.0; 
   if( iterNo < _closingDelay) gainF = ((realT) iterNo)/_closingDelay * _closingGainFact;
      
   if(_lowOrders > 0)
   {
      topN = _lowOrders;   
   }
   
   
   for(int i=0; i< topN; ++i)
   {
      if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0; 

      _commands(0,i) = (1.0-_leaks(0,i))*_commands(0,i) + gainF*_gains(0,i)*rawAmps.measurement(0,i);
      
      filtAmps.measurement(0,i) = _commands(0,i);
   }
   
   for(int i=topN; i< rawAmps.measurement.cols(); ++i)
   {
      filtAmps.measurement(0,i) = 0;
   }
   
   
   
   
}


} //namespace sim 
} //namespace AO
} //namespace mx 

#endif //__leakyIntegrator_hpp__
