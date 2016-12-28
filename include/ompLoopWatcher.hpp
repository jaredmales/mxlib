/** \file ompLoopWatcher.hpp
  * \brief Track iterations in an OMP parallelized looop.
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup utils_files
  *
  */

#ifndef __ompLoopWatcher_hpp__
#define __ompLoopWatcher_hpp__

#include <iostream>

namespace mx
{

///A class to track the number of iterations in an OMP parallelized loop.
/** Uses omp critical directives to deconflict updates by different loops.  Example:
  * \code
    ompLoopWatcher<> watcher(1000, std::cout); //Uses defaults
    
    #pragma omp parallel for 
    int( i=0; i<1000; ++i)
    {
       wathcer.incrementAndOutputStatus();
       //Do loop work
       ...
    }
    \endcode
  * This will result in the following output
   \verbatim
   1 / 1000 (0.1%)
   2 / 1000 (0.2%)
   3 / 1000 (0.3%)
   \endverbatim
  * and so on. 
  *
  * \note This can reduce performance due to critical points it creates if used in fast loops (e.g. do not use inside the inner most loop!).
  *
  * The behavior of the output is controlled through template parameters.  A different output-type can be specified, which needs to accept size_t, and optionally
  * float, and character input using the << operator.
  *
  * \tparam outputT a type which accepts size_t, float, and character input via the  << operator (default is std::ostream).
  * \tparam _printPretty flag to control whether the output is nicely formatted, if false then just the numbers are sent to the output with no spaces or delimiters (default is true).
  * \tparam _printLoops flag to control whether the number of loops is sent to the output each time (default is true).
  * \tparam _printPercent flag to control whether the percentage complete is calculated (default is true). 
  * \tparam _printNLine flag to control whether the newline '\\n' is sent ot output (default is true).
  * 
  * \ingroup mtutils
  */
template<class _outputT=std::ostream, bool _printPretty=true, bool _printLoops=true, bool _printPercent=true, bool _printNLine=true>
class ompLoopWatcher
{
public:
   typedef _outputT outputT;
   
protected:
   size_t _nLoops; ///< The total number of loops
   size_t _counter; ///< The current counter

   
   outputT * _output; ///< Pointer to the instance of type outputT

   ///Increment the counter
   void _increment()
   {
      ++_counter;
   }

   ///Perform the output
   void _outputStatus()
   {
      if(_printPretty)
      {
         (*_output) << _counter;
         if(_printLoops) (*_output) << " / " << _nLoops;
         if(_printPercent) (*_output) << " (" << 100.0*((float) _counter) / _nLoops << "%)";
         if(_printNLine) (*_output) << '\n';
      }
      else
      {
         (*_output) << _counter;
         if(_printLoops) (*_output) << _nLoops;
         if(_printPercent) (*_output) << 100.0*((float) _counter) / _nLoops;
         if(_printNLine) (*_output) << '\n';
      }
   }
   
private:
   //Default C'tor is private since you always have to give the number of loops.
   ompLoopWatcher()
   {
   }

   
public:
      
   ///Constructor
   /** Registers the output and sets the number of loops.
     *
     * \param nLoops is the total number of loops.
     * \param output is the instance of type outputT to which the output will be sent.
     */ 
   ompLoopWatcher( int nLoops, 
                   outputT & output )
   {
      _output = &output;
      _nLoops = nLoops;
      _counter = 0;
   }

   ///Increment the counter.
   /** Call this once per loop.  It contains an omp critical directive.
     */
   void increment()
   {
      #pragma omp critical
      _increment();
   }

   ///Output current status.
   /** Call this whenever you want a status update.  It contains an omp critical directive.
     */
   void outputStatus()
   {
      #pragma omp critical
      _outputStatus();
   }

   ///Increment and output status.
   /** Call this to increment and then give a status update.  Has only one omp critical directive for the two steps.
     */
   void incrementAndOutputStatus()
   {
      #pragma omp critical
      {
         _increment();
         _outputStatus();
      }
   }

};






} //namespace mx

#endif //__ompLoopWatcher_hpp__

