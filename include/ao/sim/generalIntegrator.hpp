/** \file 
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief 
  * \ingroup mxAO_sim_files
  * 
  */

#ifndef __generalIntegrator_hpp__
#define __generalIntegrator_hpp__

namespace mx
{
namespace AO 
{
namespace sim 
{
   

template<typename _realT>
class generalIntegrator
{

public:
   
   typedef _realT realT;
     
   typedef std::complex<realT> complexT;
   
   ///The wavefront data type
   typedef wavefront<realT> wavefrontT;
   
   ///The pupil image type
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> commandT;
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;
   
   generalIntegrator();

protected:

   imageT _a;
   int _currA;
   
   imageT _b;
   int _currB;
   
   imageT _gains;   
   
   imageT _commandsIn;
   imageT _commandsOut;
   
public:

   template<typename dmT>
   void initialize(dmT & dm);
   
   void setASize(int n);
   
   void setA(int i, const imageT & a);
   
   void setBSize(int n);
   
   void setB(int i, const imageT & b);
   
   realT gain(int i);
   
   void setGain(int i, realT g);
   
   void setGains(realT g);

   void setGains(const std::string & ogainf);
   
   
   
   void filterCommands( int iterNo,
                        commandT & filtAmps, 
                        commandT & rawAmps);

   int _closingDelay;
   int _lowOrders;
   realT _closingGainFact;
   
   int _nModes;
   
   void initMeasurements(commandT & filtAmps, commandT & rawAmps)
   {
      filtAmps.resize(1, _nModes);
      filtAmps.setZero();
      
      rawAmps.resize(1, _nModes);
      rawAmps.setZero();
      
   }
};


template<typename realT>
generalIntegrator<realT>::generalIntegrator()
{
   _nModes = 0;
   
   _closingDelay = 0;
   _lowOrders = 0;
   _closingGainFact = 1.0;
}



template<typename realT>
template<typename dmT>
void generalIntegrator<realT>::initialize(dmT & dm)
{
   _nModes = dm.nModes();
   
   //If _a has been sized, resize it
   if(_a.cols() > 0)
   {
      int n = _a.cols();
      _a.resize(_nModes, n);
      _a.setZero();
      _currA = 0;
      
      _commandsOut.resize( _nModes, n);
      _commandsOut.setZero();
   }
   
   //If _b has been sized, resize it
   if(_b.cols() > 0)
   {
      int n = _b.cols();
      _b.resize(_nModes, n);
      _b.setZero();
      _currB = 0;
      
      _commandsIn.resize( _nModes, n);
      _commandsIn.setZero();
   }

   _gains.resize(1,_nModes);
   
   _gains.setZero();

}

template<typename realT>
void generalIntegrator<realT>::setASize(int n)
{
   //Resize with _nModes if set
   if( _nModes > 0)
   {
      _a.resize(_nModes, n);
      _commandsOut.resize(_nModes, n);
   }
   else
   {
      _a.resize(1,n);
      _commandsOut.resize(1, n);
   }
   
   _a.setZero();
   _currA = 0;
   _commandsOut.setZero();
      
}
   
template<typename realT>
void generalIntegrator<realT>::setA(int i, const imageT & a)
{
   _a.row(i) = a;
}
   
template<typename realT>
void generalIntegrator<realT>::setBSize(int n)
{
   //Resize with _nModes if set
   if( _nModes > 0)
   {
      _b.resize(_nModes, n);
      _commandsIn.resize(_nModes,n);
   }
   else
   {
      _b.resize(1,n);
      _commandsIn.resize(1,n);
   }
   
   _b.setZero();
   _currB = 0;
   _commandsIn.setZero();
}

  
template<typename realT>
void generalIntegrator<realT>::setB(int i, const imageT & b)
{
   _b.row(i) = b;
}



template<class realT>
realT generalIntegrator<realT>::gain(int i)
{
   return _gains(i);
}

   
template<typename realT>
void generalIntegrator<realT>::setGain(int i, realT g)
{
   _gains(0,i) = g;
}

template<typename realT>
void generalIntegrator<realT>::setGains(realT g)
{
   for(int i=0;i<_gains.cols(); ++i) _gains(0,i) = g;
}

template<typename realT>
void generalIntegrator<realT>::setGains(const std::string & ogainf)
{
   std::ifstream fin;
   fin.open(ogainf);
   
   if(!fin.good())
   {
      std::cerr << "generalIntegrator: gain file " << ogainf << " not found.\n";
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


template<typename realT>
void generalIntegrator<realT>::filterCommands( int iterNo,
                                                 commandT & filtAmps, 
                                                 commandT & rawAmps )
{

   realT aTot, bTot;
   
   
   if( iterNo < _closingDelay )
   {
      BREAD_CRUMB;
      
      for(int i=0; i< _nModes; ++i)
      {
         if( std::isnan( rawAmps(0,i) ) || !std::isfinite(rawAmps(0,i))) rawAmps(0,i) = 0.0; 
         
      
         _commandsIn(i, _currB) = rawAmps(0,i);
      
         aTot = _commandsOut(i, _currA);
         //std::cerr << _currA << " " << aTot << "\n";
         
         
         int cA = _currA + 1;
         if(cA >= _a.cols()) cA = 0;
         
         realT gf = ((realT) iterNo)/_closingDelay * (_closingGainFact);
         
         _commandsOut(i, cA) = aTot + gf*_gains(i) * rawAmps(0,i);
         
         //std::cerr << cA << " " << _commandsOut(i, cA) << "\n";
      
         if( i <= _lowOrders || _lowOrders <= 0)
         {
            filtAmps(0,i) = _commandsOut(i, cA);
         }
         else
         {
            filtAmps(0,i) = 0;
         }
      }  
      
      ++_currB;
      if(_currB >= _b.cols()) _currB = 0;

         
      ++_currA;
      if(_currA >= _a.cols()) _currA = 0;
         
      return;
   }
   

   
   
   for(int i=0; i< _nModes; ++i)
   {
      if(_gains(i) == 0)
      {
         filtAmps(0,i) = 0;
         continue;
      }
      
      
      if( std::isnan( rawAmps(0,i) ) || !std::isfinite(rawAmps(0,i))) rawAmps(0,i) = 0.0; 
      
      aTot = 0;
      for(int j = 0; j < _a.cols(); ++j)
      {
         int k = _currA - j;
         if(k < 0) k += _a.cols();      
         aTot += _a(i,j) * _commandsOut(i,k);
      }
      
      _commandsIn(i, _currB) = rawAmps(0,i);
      
      bTot = 0;
      for(int j = 0; j < _b.cols(); ++j)
      {
         int k = _currB - j;
         if(k < 0) k += _b.cols();
         
         bTot += _b(i,j) * _commandsIn(i,k);
      }
      

      int cA= _currA + 1;
      if( cA >= _a.cols()) cA = 0;
      
      _commandsOut(i, cA) = aTot + _gains(i) * bTot;
      
      filtAmps(0,i) = _commandsOut(i, cA);
   }

   ++_currB;
   if(_currB >= _b.cols()) _currB = 0;
      
   ++_currA;
   if(_currA >= _a.cols()) _currA = 0;
      
}


} //namespace sim 
} //namespace AO
} //namespace mx 

#endif //__generalIntegrator_hpp__
