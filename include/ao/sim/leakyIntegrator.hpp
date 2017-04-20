/** \file leakyIntegrator.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines the leaky integrator controller class for AO control.
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
   
template<typename _realT>
struct wfMeasurement
{
   typedef _realT realT;
   
   realT iterNo;
   
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> commandT;
   
   commandT measurement;
};

///Implements the leaky integrator controller.
/**
  * \tparam _realT is the floating point type for all calculations.
  */  
template<typename _realT>
class leakyIntegrator
{

public:
   
   ///The real data type
   typedef _realT realT;
        
   ///The wavefront data type
   typedef wavefront<realT> wavefrontT;
   
   ///The command type
   typedef wfMeasurement<realT> commandT;
   
   ///The image type, used here as a general storage array
   typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> imageT;
   
   ///Default c'tor.
   leakyIntegrator();

protected:

   int _nModes; ///< The number of modes being filtered.

   bool _openLoop; ///< If true, then commands are not integrated.
   
   int _closingDelay; ///< If > 0, then the gains are ramped linearly up to this value. Default = 0.
   
   int _lowOrders; ///< If > 0, then this sets the maximum mode number which is filtered. All remaining modes are set to 0.  Default = 0.
   
   
   imageT _gains;   ///< Column-vector of gains
   
   imageT _leaks; ///< Column-vector of leaks

   imageT _commands; ///< Column-vector past commands

   
   

public:

   int _lowOrderDelay;
   
   ///Allocate and initialize all state.
   /** 
     * \returns 0 on success, a negative integer otherwise.
     */ 
   int initialize(int nModes /**< [in] the number of modes to be filtered */ );

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
   int openLoop(bool ol /**< The new value of _openLoop */ );
   
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
   int closingDelay( int cd /**< The new value of _closingDelay */);

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
   int lowOrders( int lo /**< The new value of _lowOrders */);

   ///Get the value of the _lowOrders.
   /** 
     * \returns the current value of _lowOrders.
     */   
   int lowOrders();
   
   ///Get the gain for a single mode.
   /** \returns the gain value if mode exists.  
     */
   realT gain( int i /**< [in] the mode number*/);
   
   ///Set the gain for a single mode.
   /** 
     * \returns 0 on success, negative number on error.
     */
   int gain( int i,  ///< [in] the mode number 
             realT g ///< [in] the new gain value
           );
   
   ///Set the gain for all modes to a single value
   /** 
     * \returns 0 on success, negative number on error.
     */
   int gains(realT g /**< [in] the new gain value*/ );

   ///Set the gains for all modes, using a vector to specify each gain
   /** The vector must be exactly as long as _nModes. 
     * 
     * \returns 0 on success, negative number on error.
     */
   int gains(const std::vector<realT> & vgains /**< [in] vector of gains.*/);
   
   ///Set the gains for all modes, using a file to specify each gain
   /** The file format is a simple ASCII single column, with 1 gain per line.
     * Must be exactly as long as _nModes. 
     * 
     * \returns 0 on success, negative number on error.
     */
   int gains(const std::string & ogainf /**< [in] the name of the file, full path */);
   
   
   ///Set the leak for a single mode.
   /**
     * \returns 0 on success, negative number on error.
     */
   int leak( int i, ///< [in] the mode number
             realT l ///< [in] the new leak value for this mode
           );
   
   ///Set a leak for all modes.
   /**
     * \returns 0 on success, negative number on error.
     */
   int leaks( realT l /**< [in] the new leak value to set for all modes */ );


   ///Allocate the provided command structures
   /** Used by the calling system to allocate the commands being passed between components.
     *
     * \returns 0 on success, negative number on error.
     */ 
   int initMeasurements( commandT & filtAmps, ///< The structure to contain the filtered commands 
                         commandT & rawAmps   ///< The structure to contain the raw commands
                       );

      
   ///Apply the leaky integrator
   /**
     * \returns 0 on success, negative number on error.
     */ 
   int filterCommands( commandT & filtAmps, ///< [out] the filtered commands
                       commandT & rawAmps, ///< [in] the raw commands
                       int iterNo ///< [in] The current iteration number
                     );

   

};


template<typename realT>
leakyIntegrator<realT>::leakyIntegrator()
{
   _nModes = 0;

   _openLoop = false;
   
   _closingDelay = 0;
   
   _lowOrders = 0;
   _lowOrderDelay = 0;

}



template<typename realT>
int leakyIntegrator<realT>::initialize(int nModes)
{
   _nModes = nModes;
   
   _gains.resize(1,nModes);
   _leaks.resize(1, nModes);
   _commands.resize(1, nModes);
   
   _gains.setZero();
   
   _leaks.setZero();
   
   _commands.setZero();
   
   return 0;
   
}

template<typename realT>
int leakyIntegrator<realT>::nModes()
{
   return _nModes;
}

template<typename realT>
int leakyIntegrator<realT>::openLoop( bool ol )
{
   _openLoop = ol;
   return 0;
}

template<typename realT>
bool leakyIntegrator<realT>::openLoop()
{
   return _openLoop;
}

template<typename realT>
int leakyIntegrator<realT>::closingDelay( int cd )
{
   _closingDelay = cd;
   
   return 0;
}

template<typename realT>
int leakyIntegrator<realT>::closingDelay()
{
   return _closingDelay;
}

template<typename realT>
int leakyIntegrator<realT>::lowOrders( int lo )
{
   _lowOrders = lo;
   
   return 0;
}

template<typename realT>
int leakyIntegrator<realT>::lowOrders()
{
   return _lowOrders;
}

template<typename realT>
realT leakyIntegrator<realT>::gain(int i)
{
   if(i < 0 || i >= _nModes)
   {
      mxError("leakyIntegrator::gain", MXE_INVALIDARG, "mode index out of range");
      return 0; ///\retval 0 if the mode doesn't exist.
   }
   
   return _gains(i);
}

   
template<typename realT>
int leakyIntegrator<realT>::gain(int i, realT g)
{
   if(i < 0 || i >= _nModes)
   {
      mxError("leakyIntegrator::gain", MXE_INVALIDARG, "mode index out of range");
      return 0; ///\retval 0 if the mode doesn't exist.
   }
   
   _gains(0,i) = g;
   
   return 0;
}

template<typename realT>
int leakyIntegrator<realT>::gains(realT g)
{
   for(int i=0;i<_gains.cols(); ++i) 
   {
      _gains(0,i) = g;
   }
   
   return 0;
}

template<typename realT>
int leakyIntegrator<realT>::gains( const std::vector<realT> & gains )
{
   
   if( gains.size() != _gains.cols())
   {
      mxError("leakyIntegrator::gain", MXE_SIZEERR, "input gain vector not same size as number of modes");
      return -1; ///\retval -1 on vector size mismatch
   }
   
   for(int i=0;i<_gains.cols(); ++i) 
   {
      _gains(0,i) = gains[i];
   }
   
   return 0;
}
   
template<typename realT>
int leakyIntegrator<realT>::gains( const std::string & ogainf )
{
   std::ifstream fin;
   fin.open(ogainf);
   
   if(!fin.good())
   {
      mxError("leakyIntegrator::gains", MXE_FILEOERR, "could not open gan file");
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

template<typename realT>
int leakyIntegrator<realT>::leak(int i, realT l)
{
   _leaks(0,i) = l;
}

template<typename realT>
int leakyIntegrator<realT>::leaks(realT l)
{
   for(int i=0;i<_leaks.cols(); ++i) _leaks(0,i) = l;
   
   return 0;
}

template<class realT>
int leakyIntegrator<realT>::initMeasurements(commandT & filtAmps, commandT & rawAmps)
{
   filtAmps.measurement.resize(1, _nModes);
   filtAmps.measurement.setZero();
      
   rawAmps.measurement.resize(1, _nModes);
   rawAmps.measurement.setZero();
    
   return 0;
}

template<class realT>
int leakyIntegrator<realT>::filterCommands( commandT & filtAmps, 
                                            commandT & rawAmps,
                                            int iterNo )
{
   filtAmps.iterNo = rawAmps.iterNo;
   
   if(_openLoop) 
   {
      filtAmps.measurement.setZero();
      return 0;
   }
   
   int topN = rawAmps.measurement.cols();
   
   
   
   realT gainF = 1.0; 
   realT HOgainF = 0.0;
   
   //leaks(0.01);
   
   if( iterNo < _closingDelay) gainF = ((realT) iterNo)/_closingDelay;
      
   //if( iterNo < 0.55*_closingDelay) leaks(0.1);
   
   if(_lowOrders > 0)
   {
      if( iterNo >= _closingDelay+_lowOrderDelay)
      {
         _lowOrders = 0;
      }
      else
      {
         topN = _lowOrders;   
         
         if( iterNo >= _closingDelay )
         {
            HOgainF = ( (realT) (iterNo - _closingDelay) )/( _lowOrderDelay );
            
         }
         
      }
   }
   
   
   for(int i=0; i< topN; ++i)
   {
      if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0; 

      _commands(0,i) = (1.0-_leaks(0,i))*_commands(0,i) + gainF*_gains(0,i)*rawAmps.measurement(0,i);
      
      filtAmps.measurement(0,i) = _commands(0,i);
   }
   
   for(int i=topN; i< rawAmps.measurement.cols(); ++i)
   {      
      if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0; 

      _commands(0,i) = (1.0-_leaks(0,i))*_commands(0,i) + HOgainF*_gains(0,i)*rawAmps.measurement(0,i);
      
      filtAmps.measurement(0,i) = _commands(0,i);
      
      //filtAmps.measurement(0,i) = 0;
   }
   
   
   return 0;
   
}



} //namespace sim 
} //namespace AO
} //namespace mx 

#endif //__leakyIntegrator_hpp__
