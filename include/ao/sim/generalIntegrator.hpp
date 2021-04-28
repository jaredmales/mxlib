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

   //typedef Eigen::Array< realT, Eigen::Dynamic, Eigen::Dynamic> commandT;
   typedef std::vector<realT> commandT;
   
   commandT measurement;
};

///Implements a general integrator controller.
/** \todo document the math for the G.I.
  * \todo document the filter coefficients, etc.
  * 
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


protected:

   int m_nModes {0};        ///< The number of modes being filtered.

   bool m_openLoop {false}; ///< If true, then commands are not integrated. Default is false.

   int m_openLoopDelay {0}; ///< If > 0, then the loop is open for this time in time steps. Default is 0.
   
   int m_closingRamp {0};   ///< If > 0, then gains are famped linearly up to m_closingGains over this interval in time steps.  Default is 0. This is relative m_openLoopDelay.

   int m_closingDelay {0};  ///< If > 0, then the simple integrator,with m_closingGains is used up to this timestep.   This is relative m_openLoopDelay. Default = 0.

   imageT m_closingGains;   ///< Column-vector of gains used for loop closing as simple integrator

   int m_lowOrders {0};    ///< If > 0, then this sets the maximum mode number which is filtered. All remaining modes are set to 0.  Default = 0.

   imageT m_a;
   int m_currA;

   imageT m_b;
   int m_currB;

   imageT m_gains;   ///< Column-vector of gains

   imageT m_commandsIn;
   imageT m_commandsOut;

public:

   ///Allocate and initialize all state.
   /**
     * \returns 0 on success, a negative integer otherwise.
     */
   int initialize(int nModes /**< [in] the number of modes to be filtered */);

   ///Get the number of modes
   /** nModes is only set by calling initialize.
     *
     * \returns the current value of m_nModes
     */
   int nModes();

   ///Set the m_openLoop flag.
   /** If m_openLoop is true, then commands are not filtered.
     *
     * \returns 0 on success, a negative integer on error.
     */
   int openLoop(bool ol /**< [in] the new value of m_openLoop */ );

   ///Get the value of the m_openLoop flag.
   /**
     * \returns the current value of m_openLoop.
     */
   bool openLoop();

   ///Set the open loop delay.
   /** If m_openLoopDelay  > 0, then the loop is open for this period in time-steps
     * Requires m_closingDelay > 0.
     * 
     * \returns 0 on success, a negative integer on error.
     */
   int openLoopDelay( int cr /**< [in] The new value of m_openLoopDelay */);

   ///Get the value of the m_openLoopDelay.
   /**
     * \returns the current value of m_openLoopDelay.
     */
   int openLoopDelay();
   
   ///Set the closing ramp.
   /** If m_closingRamp  > 0, then the gains are ramped linearly up to m_closingGains over this timer interval in timesteps.
     * Requires m_closingDelay > 0.
     * 
     * \returns 0 on success, a negative integer on error.
     */
   int closingRamp( int cr /**< [in] The new value of m_closingRamp */);

   ///Get the value of the m_closingRamp.
   /**
     * \returns the current value of m_closingRamp.
     */
   int closingRamp();
   
   ///Set the closing delay.
   /** If m_closingDelay  > 0, then the simple integragor is used until this timestep.
     *
     * \returns 0 on success, a negative integer on error.
     */
   int closingDelay( int cd /**< [in] The new value of m_closingDelay */);

   ///Get the value of the m_closingDelay.
   /**
     * \returns the current value of m_closingDelay.
     */
   int closingDelay();

   /// Set the simple integrator gains to use during closing.
   /**
     * \returns 0 on success
     * \returns < 0 on error 
     */
   int closingGains( const std::vector<realT> & gains /**< [in] vector of gains.*/);

   /// Set the simple integrator gain to use for all modes during closing.
   /**
     * \returns 0 on success
     * \returns < 0 on error 
     */
   int closingGains(realT g);

   
   ///Set m_lowOrders.
   /** If m_lowOrders > 0, then this sets the maximum mode number which is filtered. All remaining modes are set to 0.
     *
     * \returns 0 on success, a negative integer on error.
     */
   int lowOrders( int lo /**< [in] The new value of m_lowOrders */);

   ///Get the value of the m_lowOrders.
   /**
     * \returns the current value of m_lowOrders.
     */
   int lowOrders();

   ///Set the size of the IIR vector (the a coefficients)
   /** This allocates m_a to be m_nModes X n in size.
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
   /** This allocates m_b to be m_nModes X n in size.
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
   /** The vector must be exactly as long as m_nModes.
     *
     * \returns 0 on success, negative number on error.
     */
   int gains( const std::vector<realT> & gains /**< [in] vector of gains.*/);


   ///Set the gains for all modes, using a file to specify each gain
   /** The file format is a simple ASCII single column, with 1 gain per line.
     * Must be exactly as long as m_nModes.
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
}



template<typename realT>
int generalIntegrator<realT>::initialize(int nModes)
{
   m_nModes = nModes;

   //If m_a has been sized, resize it
   if(m_a.cols() > 0)
   {
      int n = m_a.cols();
      m_a.resize(m_nModes, n);
      m_a.setZero();
      m_currA = 0;

      m_commandsOut.resize( m_nModes, n);
      m_commandsOut.setZero();
   }

   //If m_b has been sized, resize it
   if(m_b.cols() > 0)
   {
      int n = m_b.cols();
      m_b.resize(m_nModes, n);
      m_b.setZero();
      m_currB = 0;

      m_commandsIn.resize( m_nModes, n);
      m_commandsIn.setZero();
   }

   m_closingGains.resize(1,m_nModes);

   m_closingGains.setZero();

   m_gains.resize(1,m_nModes);

   m_gains.setZero();


   return 0;
}

template<typename realT>
int generalIntegrator<realT>::nModes()
{
   return m_nModes;
}

template<typename realT>
int generalIntegrator<realT>::openLoop( bool ol )
{
   m_openLoop = ol;
   return 0;
}

template<typename realT>
bool generalIntegrator<realT>::openLoop()
{
   return m_openLoop;
}

template<typename realT>
int generalIntegrator<realT>::openLoopDelay( int old )
{
   m_openLoopDelay = old;

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::openLoopDelay()
{
   return m_openLoopDelay;
}

template<typename realT>
int generalIntegrator<realT>::closingRamp( int cr )
{
   m_closingRamp = cr;

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::closingRamp()
{
   return m_closingRamp;
}

template<typename realT>
int generalIntegrator<realT>::closingDelay( int cd )
{
   m_closingDelay = cd;

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::closingDelay()
{
   return m_closingDelay;
}

template<typename realT>
int generalIntegrator<realT>::closingGains( const std::vector<realT> & gains )
{

   if( gains.size() != (size_t) m_closingGains.cols())
   {
      mxError("generalIntegrator::closingGains", MXE_SIZEERR, "input gain vector not same size as number of modes");
      return -1; ///\retval -1 on vector size mismatch
   }

   for(int i=0;i<m_closingGains.cols(); ++i)
   {
      m_closingGains(0,i) = gains[i];
   }

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::closingGains(realT g)
{
   for(int i=0;i<m_closingGains.cols(); ++i)
   {
      m_closingGains(0,i) = g;
   }

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::lowOrders( int lo )
{
   m_lowOrders = lo;

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::lowOrders()
{
   return m_lowOrders;
}

template<typename realT>
int generalIntegrator<realT>::setASize(int n)
{
   //Resize with m_nModes if set
   if( m_nModes > 0)
   {
      m_a.resize(m_nModes, n);
      m_commandsOut.resize(m_nModes, n);
   }
   else
   {
      m_a.resize(1,n);
      m_commandsOut.resize(1, n);
   }

   m_a.setZero();
   m_currA = 0;
   m_commandsOut.setZero();

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::setA(int i, const imageT & a)
{
   m_a.row(i) = a;

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::setBSize(int n)
{
   //Resize with m_nModes if set
   if( m_nModes > 0)
   {
      m_b.resize(m_nModes, n);
      m_commandsIn.resize(m_nModes,n);
   }
   else
   {
      m_b.resize(1,n);
      m_commandsIn.resize(1,n);
   }

   m_b.setZero();
   m_currB = 0;
   m_commandsIn.setZero();

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::setB(int i, const imageT & b)
{
   m_b.row(i) = b;

   return 0;
}


template<class realT>
realT generalIntegrator<realT>::gain(int i)
{
   if(i < 0 || i >= m_nModes)
   {
      mxError("generalIntegrator::gain", MXE_INVALIDARG, "mode index out of range");
      return 0; ///\retval 0 if the mode doesn't exist.
   }

   return m_gains(i);
}

template<typename realT>
int generalIntegrator<realT>::gain(int i, realT g)
{
   if(i < 0 || i >= m_nModes)
   {
      mxError("generalIntegrator::gain", MXE_INVALIDARG, "mode index out of range");
      return 0; ///\retval 0 if the mode doesn't exist.
   }

   m_gains(0,i) = g;

   return 0;
}


template<typename realT>
int generalIntegrator<realT>::gains(realT g)
{
   for(int i=0;i<m_gains.cols(); ++i)
   {
      m_gains(0,i) = g;
   }

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::gains( const std::vector<realT> & gains )
{

   if( gains.size() != (size_t) m_gains.cols())
   {
      mxError("generalIntegrator::gains", MXE_SIZEERR, "input gain vector not same size as number of modes");
      return -1; ///\retval -1 on vector size mismatch
   }

   for(int i=0;i<m_gains.cols(); ++i)
   {
      m_gains(0,i) = gains[i];
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
   for(int i=0;i<m_gains.cols();++i)
   {
      fin >> tmpstr;
      g = ioutils::convertFromString<realT>(tmpstr);
      gain(i, g);

   }

   return 0;
}

template<class realT>
int generalIntegrator<realT>::initMeasurements(commandT & filtAmps, commandT & rawAmps)
{
   filtAmps.measurement.resize(m_nModes);
   for(size_t n=0; n<m_nModes; ++n) filtAmps.measurement[n] = 0;
   //filtAmps.measurement.setZero();

   rawAmps.measurement.resize(m_nModes);
   for(size_t n=0; n<m_nModes; ++n) rawAmps.measurement[n] = 0;
   //rawAmps.measurement.resize(1, m_nModes);
   //rawAmps.measurement.setZero();

   return 0;
}

template<typename realT>
int generalIntegrator<realT>::filterCommands( commandT & filtAmps,
                                              commandT & rawAmps,
                                              int iterNo )
{
   filtAmps.iterNo = rawAmps.iterNo;

   if(m_openLoop || iterNo < m_openLoopDelay)
   {
      for(size_t n=0; n<m_nModes; ++n) filtAmps.measurement[n] = 0;
      //filtAmps.measurement.setZero();
      return 0;
   }


   realT aTot, bTot;

   if( iterNo - m_openLoopDelay < m_closingDelay)
   {
      BREAD_CRUMB;

      for(int i=0; i< m_nModes; ++i)
      {
         //if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0;
         if( std::isnan( rawAmps.measurement[i] ) || !std::isfinite(rawAmps.measurement[i])) rawAmps.measurement[i] = 0.0;

         //m_commandsIn(i, m_currB) = rawAmps.measurement(0,i);
         m_commandsIn(i, m_currB) = rawAmps.measurement[i];

         aTot = m_commandsOut(i, m_currA);
         //std::cerr << m_currA << " " << aTot << "\n";


         int cA = m_currA + 1;
         if(cA >= m_a.cols()) cA = 0;


         realT gf = 1.0;

         if( iterNo - m_openLoopDelay < m_closingRamp) gf = ((realT) iterNo - m_openLoopDelay)/m_closingRamp;


         //m_commandsOut(i, cA) = aTot + gf*m_closingGains(i) * rawAmps.measurement(0,i);
         m_commandsOut(i, cA) = aTot + gf*m_closingGains(i) * rawAmps.measurement[i];

         //std::cerr << cA << " " << m_commandsOut(i, cA) << "\n";

         if( i <= m_lowOrders || m_lowOrders <= 0)
         {
            filtAmps.measurement[i] = m_commandsOut(i, cA);
         }
         else
         {
            filtAmps.measurement[i] = 0;
         }
      }

      ++m_currB;
      if(m_currB >= m_b.cols()) m_currB = 0;


      ++m_currA;
      if(m_currA >= m_a.cols()) m_currA = 0;

      return 0;
   }




   for(int i=0; i< m_nModes; ++i)
   {
      if(m_gains(i) == 0)
      {
         //filtAmps.measurement(0,i) = 0;
         filtAmps.measurement[i] = 0;
         continue;
      }


      //if( std::isnan( rawAmps.measurement(0,i) ) || !std::isfinite(rawAmps.measurement(0,i))) rawAmps.measurement(0,i) = 0.0;
      if( std::isnan( rawAmps.measurement[i] ) || !std::isfinite(rawAmps.measurement[i])) rawAmps.measurement[i] = 0.0;

      aTot = 0;
      for(int j = 0; j < m_a.cols(); ++j)
      {
         int k = m_currA - j;
         if(k < 0) k += m_a.cols();
         aTot += m_a(i,j) * m_commandsOut(i,k);
      }

      m_commandsIn(i, m_currB) = rawAmps.measurement[i];

      bTot = 0;
      for(int j = 0; j < m_b.cols(); ++j)
      {
         int k = m_currB - j;
         if(k < 0) k += m_b.cols();

         bTot += m_b(i,j) * m_commandsIn(i,k);
      }


      int cA= m_currA + 1;
      if( cA >= m_a.cols()) cA = 0;

      m_commandsOut(i, cA) = aTot + m_gains(i) * bTot;

      filtAmps.measurement[i] = m_commandsOut(i, cA);
   }

   ++m_currB;
   if(m_currB >= m_b.cols()) m_currB = 0;

   ++m_currA;
   if(m_currA >= m_a.cols()) m_currA = 0;

   return 0;
}


} //namespace sim
} //namespace AO
} //namespace mx

#endif //__generalIntegrator_hpp__
