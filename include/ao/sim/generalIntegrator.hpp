/** \file generalIntegrator.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines a general integrator controller class for AO control.
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
struct wfMeasurement
{
   typedef _realT realT;
   
   realT iterNo;
   
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> commandT;
   
   commandT measurement;
};

///Implements a general integrator controller.
/**
  * \tparam _realT is the floating point type for all calculations.
  */ 
template<typename _realT>
class generalIntegrator
{

public:
   
   ///The real data type
   typedef _realT realT;
     
   ///The real data type
   typedef std::complex<realT> complexT;
   
   ///The wavefront data type
//   typedef wavefront<realT> wavefrontT;
   
   ///The command type
   typedef wfMeasurement<realT> commandT;
   
   ///The image type, used here as a general storage array
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;
   
   ///Default c'tor.
   generalIntegrator();

   imageT _closingGains;   ///< Column-vector of gains used for loop closing as pure integrator
   int _closingRamp; ///< Number of steps over which to ramp gains linearly during loop closing.
   
protected:

   int _nModes; ///< The number of modes being filtered.

   bool _openLoop; ///< If true, then commands are not integrated.
   
   int _closingDelay; ///< If > 0, then the gains are ramped linearly up to this value. Default = 0.
   
   int _lowOrders; ///< If > 0, then this sets the maximum mode number which is filtered. All remaining modes are set to 0.  Default = 0.
   
   
      
   imageT _a;
   int _currA;
   
   imageT _b;
   int _currB;
   
   imageT _gains;   ///< Column-vector of gains
   
   imageT _commandsIn;
   imageT _commandsOut;
   
public:

   ///Allocate and initialize all state.
   /** 
     * \returns 0 on success, a negative integer otherwise.
     */ 
   int initialize(int nModes /**< [in] the number of modes to be filtered */);
   
   ///Get the number of modes
   /** nModes is only set by calling initialize.
     *
     * \returns the current value of _nModes
     */ 
   int nModes();
   
   ///Set the _openLoop flag.
   /** If _openLoop is true, then commands are not filtered.
     *
     * \returns 0 on success, a negative integer on error. 
     */
   int openLoop(bool ol /**< [in] the new value of _openLoop */ );
   
   ///Get the value of the _openLoop flag.
   /** 
     * \returns the current value of _openLoop.
     */
   bool openLoop();
   
   ///Set _closingDelay.
   /** If _closingDelay  > 0, then the gains are ramped linearly up to this value.
     *
     * \returns 0 on success, a negative integer on error. 
     */
   int closingDelay( int cd /**< [in] The new value of _closingDelay */);

   ///Get the value of the _closingDelay.
   /** 
     * \returns the current value of _closingDelay.
     */   
   int closingDelay();

   ///Set _lowOrders.
   /** If _lowOrders > 0, then this sets the maximum mode number which is filtered. All remaining modes are set to 0. 
     *
     * \returns 0 on success, a negative integer on error. 
     */   
   int lowOrders( int lo /**< [in] The new value of _lowOrders */);

   ///Get the value of the _lowOrders.
   /** 
     * \returns the current value of _lowOrders.
     */   
   int lowOrders();
   
   ///Set the size of the IIR vector (the a coefficients)
   /** This allocates _a to be _nModes X n in size.
     * 
     * \returns 0 on success, negative number on error.
     */ 
   int setASize(int n /**< [in] the number of IIR coefficients */ );
   
   ///Set the IIR coefficients for one mode.
   /** 
     * \returns 0 on success, negative number on error.
     */
   int setA( int i,  ///< [in] the mode number
             const imageT & a ///< [in] the IIR coefficients for this mode
           );
   
   ///Set the size of the FIR vector (the b coefficients)
   /** This allocates _b to be _nModes X n in size.
     *
     * \returns 0 on success, negative number on error.
     */ 
   int setBSize(int n /**<  [in] the number of FIR coefficients */ );
   
   ///Set the FIR coefficients for one mode.
   /** 
     * \returns 0 on success, negative number on error.
     */
   int setB( int i,  ///< [in] the mode number
             const imageT & b ///< [in] the IIR coefficients for this mode
           );
   
   ///Get the gain for a single mode.
   /** 
     * \returns the gain value if mode exists.  
     */
   realT gain( int i /**< The mode number*/);
   
   ///Set the gain for a single mode.
   /** 
     * \returns 0 on success, negative number on error.
     */
   int gain( int i,  ///< The mode number 
              realT g ///< The new gain value
            );
   
   ///Set the gain for all modes to a single value
   /** 
     * \returns 0 on success, negative number on error.
     */
   int gains(realT g /**< The new gain value*/ );
   
   ///Set the gains for all modes, using a vector to specify each gain
   /** The vector must be exactly as long as _nModes. 
     * 
     * \returns 0 on success, negative number on error.
     */
   int gains( const std::vector<realT> & gains /**< [in] vector of gains.*/);
   
   int closingGains( const std::vector<realT> & gains /**< [in] vector of gains.*/);
   
   int closingGains(realT g);
   
   ///Set the gains for all modes, using a file to specify each gain
   /** The file format is a simple ASCII single column, with 1 gain per line.
     * Must be exactly as long as _nModes. 
     * 
     * \returns 0 on success, negative number on error.
     */
   int gains( const std::string & ogainf /**<  [in] the name of the file, full path */ );
   
   
   ///Allocate the provided command structures
   /** Used by the calling system to allocate the commands being passed between components.
     *
     * \returns 0 on success, negative number on error.
     */ 
   int initMeasurements( commandT & filtAmps, ///< The structure to contain the filtered commands
                         commandT & rawAmps ///< The structure to contain the raw commands
                       );
   
   int filterCommands( commandT & filtAmps, 
                       commandT & rawAmps,
                       int iterNo
                     );

   
   

};


template<typename realT>
generalIntegrator<realT>::generalIntegrator()
{
   _nModes = 0;
   
   _openLoop = false;
   
   _closingDelay = 0;
   
   _closingRamp = 0;
   
   _lowOrders = 0;
   
}



template<typename realT>
int generalIntegrator<realT>::initialize(int nModes)
{
   _nModes = nModes;
   
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

   _closingGains.resize(1,_nModes);
   
   _closingGains.setZero();
   
   _gains.resize(1,_nModes);
   
   _gains.setZero();

   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::nModes()
{
   return _nModes;
}

template<typename realT>
int generalIntegrator<realT>::openLoop( bool ol )
{
   _openLoop = ol;
   return 0;
}

template<typename realT>
bool generalIntegrator<realT>::openLoop()
{
   return _openLoop;
}

template<typename realT>
int generalIntegrator<realT>::closingDelay( int cd )
{
   _closingDelay = cd;
   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::closingDelay()
{
   return _closingDelay;
}

template<typename realT>
int generalIntegrator<realT>::lowOrders( int lo )
{
   _lowOrders = lo;
   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::lowOrders()
{
   return _lowOrders;
}

template<typename realT>
int generalIntegrator<realT>::setASize(int n)
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
   
   return 0;
}
   
template<typename realT>
int generalIntegrator<realT>::setA(int i, const imageT & a)
{
   _a.row(i) = a;
   
   return 0;
}
   
template<typename realT>
int generalIntegrator<realT>::setBSize(int n)
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
   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::setB(int i, const imageT & b)
{
   _b.row(i) = b;
   
   return 0;
}


template<class realT>
realT generalIntegrator<realT>::gain(int i)
{
   if(i < 0 || i >= _nModes)
   {
      mxError("generalIntegrator::gain", MXE_INVALIDARG, "mode index out of range");
      return 0; ///\retval 0 if the mode doesn't exist.
   }
   
   return _gains(i);
}
   
template<typename realT>
int generalIntegrator<realT>::gain(int i, realT g)
{
   if(i < 0 || i >= _nModes)
   {
      mxError("generalIntegrator::gain", MXE_INVALIDARG, "mode index out of range");
      return 0; ///\retval 0 if the mode doesn't exist.
   }
   
   _gains(0,i) = g;
   
   return 0;
}


template<typename realT>
int generalIntegrator<realT>::gains(realT g)
{
   for(int i=0;i<_gains.cols(); ++i) 
   {
      _gains(0,i) = g;
   }
   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::gains( const std::vector<realT> & gains )
{
   
   if( gains.size() != _gains.cols())
   {
      mxError("generalIntegrator::gains", MXE_SIZEERR, "input gain vector not same size as number of modes");
      return -1; ///\retval -1 on vector size mismatch
   }
   
   for(int i=0;i<_gains.cols(); ++i) 
   {
      _gains(0,i) = gains[i];
   }
   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::closingGains( const std::vector<realT> & gains )
{
   
   if( gains.size() != _closingGains.cols())
   {
      mxError("generalIntegrator::closingGains", MXE_SIZEERR, "input gain vector not same size as number of modes");
      return -1; ///\retval -1 on vector size mismatch
   }
   
   for(int i=0;i<_closingGains.cols(); ++i) 
   {
      _closingGains(0,i) = gains[i];
   }
   
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::closingGains(realT g)
{
   for(int i=0;i<_closingGains.cols(); ++i) 
   {
      _closingGains(0,i) = g;
   }
   
   return 0;
}


template<typename realT>
int generalIntegrator<realT>::gains(const std::string & ogainf)
{
   std::ifstream fin;
   fin.open(ogainf);
   
   if(!fin.good())
   {
      mxError("generalIntegrator::gains", MXE_FILEOERR, "could not open gan file");
      return -1; /// \retval -1 if file open fails
   }
   
   std::string tmpstr;
   realT g;
   for(int i=0;i<_gains.cols();++i)
   {
      fin >> tmpstr;
      g = mx::convertFromString<realT>(tmpstr);
      gain(i, g);

   }
   
   return 0;
}

template<class realT>
int generalIntegrator<realT>::initMeasurements(commandT & filtAmps, commandT & rawAmps)
{
   filtAmps.measurement.resize(1, _nModes);
   filtAmps.measurement.setZero();
      
   rawAmps.measurement.resize(1, _nModes);
   rawAmps.measurement.setZero();
    
   return 0;
}

template<typename realT>
int generalIntegrator<realT>::filterCommands( commandT & filtAmps, 
                                              commandT & rawAmps,
                                              int iterNo )
{
   filtAmps.iterNo = rawAmps.iterNo;
   
   if(_openLoop) 
   {
      filtAmps.measurement.setZero();
      return 0;
   }
   
   
   realT aTot, bTot;

   if( iterNo < _closingDelay)
   {
      BREAD_CRUMB;
      
      for(int i=0; i< _nModes; ++i)
      {
         if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0; 
         
      
         _commandsIn(i, _currB) = rawAmps.measurement(0,i);
      
         aTot = _commandsOut(i, _currA);
         //std::cerr << _currA << " " << aTot << "\n";
         
         
         int cA = _currA + 1;
         if(cA >= _a.cols()) cA = 0;
         
         
         realT gf = 1.0;
         
         if( iterNo < _closingRamp) gf = ((realT) iterNo)/_closingRamp;
         
         
         _commandsOut(i, cA) = aTot + gf*_closingGains(i) * rawAmps.measurement(0,i);
         
         //std::cerr << cA << " " << _commandsOut(i, cA) << "\n";
      
         if( i <= _lowOrders || _lowOrders <= 0)
         {
            filtAmps.measurement(0,i) = _commandsOut(i, cA);
         }
         else
         {
            filtAmps.measurement(0,i) = 0;
         }
      }  
      
      ++_currB;
      if(_currB >= _b.cols()) _currB = 0;

         
      ++_currA;
      if(_currA >= _a.cols()) _currA = 0;
         
      return 0;
   }
   

   
   
   for(int i=0; i< _nModes; ++i)
   {
      if(_gains(i) == 0)
      {
         filtAmps.measurement(0,i) = 0;
         continue;
      }
      
      
      if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0; 
      
      aTot = 0;
      for(int j = 0; j < _a.cols(); ++j)
      {
         int k = _currA - j;
         if(k < 0) k += _a.cols();      
         aTot += _a(i,j) * _commandsOut(i,k);
      }
      
      _commandsIn(i, _currB) = rawAmps.measurement(0,i);
      
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
      
      filtAmps.measurement(0,i) = _commandsOut(i, cA);
   }

   ++_currB;
   if(_currB >= _b.cols()) _currB = 0;
      
   ++_currA;
   if(_currA >= _a.cols()) _currA = 0;
      
   return 0;
}


} //namespace sim 
} //namespace AO
} //namespace mx 

#endif //__generalIntegrator_hpp__
