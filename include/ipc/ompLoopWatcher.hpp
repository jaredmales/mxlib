/** \file ompLoopWatcher.hpp
  * \brief Track iterations in an OMP parallelized looop.
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup utils_files
  *
  */

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef ompLoopWatcher_hpp
#define ompLoopWatcher_hpp

#include <iostream>
#include "../sys/timeUtils.hpp"

namespace mx
{
namespace ipc
{
   
///A class to track the number of iterations in an OMP parallelized loop.
/** Uses omp critical directives to deconflict updates by different loops.  Example:
  * \code
    ompLoopWatcher<> watcher(1000, std::cout); //Uses defaults
    
    #pragma omp parallel for 
    int( i=0; i<1000; ++i)
    {
       watcher.incrementAndOutputStatus();
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
  * \tparam _printNLine flag to control whether the newline '\\n' is sent ot output (default is false).  If false, '\\r' is written at end of output.
  * \tparam _time flag to control whether time is tracked.  Default is true.
  * \ingroup mtutils
  */
template<class _outputT=std::ostream, bool _printPretty=true, bool _printLoops=true, bool _printPercent=true, bool _printNLine=false, bool _time = true>
class ompLoopWatcher
{
public:
   typedef _outputT outputT;
   
protected:
   size_t _nLoops; ///< The total number of loops
   size_t _counter; ///< The current counter

   double t0;
   double t1;
   
   outputT * _output; ///< Pointer to the instance of type outputT

   ///Increment the counter
   void _increment()
   {
      ++_counter;
      if(_time) t1 = sys::get_curr_time();
   }

   ///Perform the output
   void _outputStatus()
   {
      if(_printPretty)
      {
         (*_output) << _counter;
         if(_printLoops) (*_output) << " / " << _nLoops;
         if(_printPercent) (*_output) << " (" << 100.0*((float) _counter) / _nLoops << "%)";
         if(_time)
         {
            (*_output) << " " << (t1-t0)/_counter << " s/loop ";
            (*_output) << " ~" << (_nLoops - _counter)*(t1-t0)/_counter << " s left";
         }
         
         if(_printNLine) (*_output) << '\n';
         if(!_printNLine)   (*_output) << "           \r";
         
         (*_output) << std::flush;
      }
      else
      {
         (*_output) << _counter;
         if(_printLoops) (*_output) << _nLoops;
         if(_printPercent) (*_output) << 100.0*((float) _counter) / _nLoops;
         if(_time) (*_output) << " " << t0 << " " << t1;
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
      
      if(_time) t0 = sys::get_curr_time();
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




} //namespace ipc

} //namespace mx

#endif //ompLoopWatcher_hpp

