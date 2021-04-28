/** \file gnuPlot.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declaration and definition of an interface to the gnuplot program
  * \ingroup plotting_files
  * 
*/

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef math_plot_gnuPlot_hpp
#define math_plot_gnuPlot_hpp

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <fcntl.h>
#include <poll.h>

#include <stdlib.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <cmath>

#include "../constants.hpp"
#include "../../mxError.hpp"
#include "../../sys/timeUtils.hpp"

#ifndef MX_GP_FNAME_SZ
///The size of the string for managing the stderr fifo
/** \ingroup plotting
  */
#define MX_GP_FNAME_SZ (128)
#endif
   
#ifndef MX_GP_TEMP_SZ
///The size of the string for managing temporary file names
/** \ingroup plotting
  */
#define MX_GP_TEMP_SZ (128)
#endif

#ifndef MX_GP_FC_TIME
/** \def MX_GP_FC_TIME
  * Time, in microseconds, to wait for gnuplot startup and file creation to complete.
  * Opening is retried after each timeout of this length for MX_GP_FC_RETRIES attempts.
  * \ingroup plotting
  */ 
#define MX_GP_FC_TIME (10000)
#endif
   
#ifndef MX_GP_FC_RETRIES
/** \def MX_GP_FC_RETRIES 
  * Number of times to retry opening the gnuplot stderr file.
  * \ingroup plotting
  */
#define MX_GP_FC_RETRIES (10)
#endif

namespace mx
{
namespace math 
{

   
//Get the gnuplot binary format string for a type
//Specializations below.
template<typename dataT>
std::string gpBinaryFormat()
{
   static_assert(std::is_fundamental<dataT>::value || !std::is_fundamental<dataT>::value, "No gnuplot format specifier available for this type.");
   return "";
}

/// An interactive c++ interface to gnuplot
/** Spawns a gnuplot sesssion and communicates with it.
  * \ingroup plotting 
  * 
  * \todo Use mxError for error reporting.
  * 
  * An example of using gnuPlot to plot data from a file: 
  * \code
  * gnuPlot gp; //This automatically connects, and is now ready to plot.
  * 
  * gp.plot("my_data_file.dat", "u 1:2 w l"); //Plots the data in the file
  * \endcode
  * 
  * Arrays can be plotted directly:
  * \code
  * std::vector<float> x;
  * std::vector<double> y;
  * 
  * //...code to populate x and y
  * 
  * gp.plot(x, y); //Plots y vs x.  
  * \endcode
  * Note that the vector data types are different -- this is handled automatically.
  * 
  * To get a response from the gnuplot session:
  * \code
  * std::cout << gp.getResponse("show terminal") << "\n"; //Sends the command to gnuplot, gets the response.
  * \endcode
  * Error checking and reporting is not straightforward since gnuplot does not return a result if there is no error.  If there is an error,
  * it can take several hundred milliseconds for a response to be available.  So for a typical 
  * plot command, say, one does not want to try to read the stderr ouput just in case as this will just time out with a delay.
  *
  * For cases where a response is expected, the following example shows how to check for errors.
  * \code
  * errno = 0;
  * std::string response = gp.getResponse("show terminal"); //Waits up to 0.5 seconds, or a user specified timeout. 
  * 
  * if(response == "") //This indicates some error occurred, or a timeout.
  * {
  *   if(gp.gpError())
  *   {
  *      std::cerr << "gnuplot returned error:\n" << gp.gpErrorMsg() << "\n";
  *   }
  *   else if(errno)
  *   {
  *      perror("error getting response: ");
  *   }
  *   else
  *   {
  *      std::cerr << "timed out\n";
  *   }
  * }
  * else
  * {
  *   std::cout << response << "\n";
  * }
  * 
  * \endcode
  */
class gnuPlot
{
protected:   

   int _connected {0};
   
   ///Set to true if the response indicates a gnuplot error
   bool _gpError {false};
   std::string _gpErrorMsg;
   
   ///File stream for the gnuplot interface
   FILE * _pipeH {0};
   
   ///Where to create gnuplot stderr fifo
   /** Default is /dev/shm/
     */
   std::string _errLocation;
   
   ///File name of the gnuplot stderr fifo
   std::string _errFName;
   
   ///File descriptor for the gnuplot stderr fifo
   int _errFD {0};
   
   ///Location of temporary files 
   /** Default is /dev/shm/
     */
   std::string _tempLocation;
   
   ///Vector of all temporary file names opened, used for removal on destruction.
   std::vector<std::string> _tempFiles;

   ///Flag to control whether temporary files are deleted on destruction.  Default is true (files deleted).
   bool _deleteTemp {true};

   bool _plotted {false};
   
public:
   
   gnuPlot();
   
   ~gnuPlot();

   ///Connect to gnuplot
   /** Spawns a gnuplot session using popen with stderr redirected to a temporary file.  The temporary file
     * is opened for reading.
     *
     * \retval 0 on success
     * \retval -1 on error
     */ 
   int connect();
   
   ///Return the value of the gpError flag.
   /** This flag is set if gnuplot returns an error message (see checkResponse()).  The error message
     * can be accessed with gpErrorMsg(). 
     *
     * \retval 0 if no error
     * \retval -1 if an error has occurred. 
     */
   bool gpError();
   
   ///Return the gnuplot error message
   /** The error message is extracted by checkResponse() when the response from gnuplot begins with "gnuplot>".
     */
   std::string gpErrorMsg();
   
   ///Send a command to gnuplot
   /** The newline is appended to the command, and then it is sent to gnuplot.
     *
     *
     * \retval 0 on success
     * \retval -1 on error
     */ 
   int command( const std::string & com, ///< [in] the command string
                bool flush = true        ///< [in] [optional] if true (default), then the output stream is flushed once the command is written
              );
     
   ///Check for a response from gnuplot.
   /** It typically takes some time for a response from gnuplot to become readable,
     * whether an error occurs or not. It is often 0.3 seconds or more.  If there is no error, there will be no response so waiting
     * for a response to check for errors after every command can be very time consuming.  Thus it is a choice whether to check 
     * for errors from gnuplot after every command, and what timeout to use. 
     * 
     * \todo investigate having a second thread monitor for responses.
     * 
     * gnuplot terminates all outputs with \c \\n\\n, so this reads up to these two characters.  It then strips any leading and trailing whitespace. 
     * If the response begins with \c gnuplot>, then the response is an error.  In this case, checkResponse returns -1 and the gpError flag is set.  If gpError is not set 
     * but the return value is -1, then some other error occurred (check \c errno).
     * 
     * \retval 0 on timeout or successful read
     * \retval -1 on error, or if the gnuplot response indicates an error.
     */
   int checkResponse( std::string & response, ///< [out] contains the response from gnuplot, but is empty on timeout or error.
                      double timeout = 0      ///< [in] [optional] the length of time, in seconds, to wait for a response. Default is 0, but a minimum of 0.5 if a response is expected.
                    );
   
   ///Get a response from gnuplot for a given command.
   /** This should only be used if a response is expected.
     *
     * \returns "" (empty string) on timeout or error (\c errno and gpError() should be checked), and the response from gnuplot otherwise.
     */ 
   std::string getResponse( const std::string & com, ///< [in] the command string
                            double timeout=0.5       ///< [in] the length of time, in seconds, to wait for a response. Default is 0.5, which is the minimum that should be used.
                          );

   
   /// Issue the \b replot command
   /**
     * \retval 0 on success
     * \retval -1 on error
     */
   int replot();
   
   /// Set the y axis to log scale
   /** Sends the command:
     * \verbatim
       set log y
      \endverbatim
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   int logy();
   
   /// Set the x axis to log scale
   /** Sends the command:
     * \verbatim
       set log x
      \endverbatim
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   int logx();
   
   /// Set the x and y axes to log scale
   /** Sends the commands:
     * \verbatim
       set log x
       set log y
      \endverbatim
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   int logxy();
   
   /// Unset the y axis from log scale
   /** Sends the command:
     * \verbatim
       unset log y
      \endverbatim
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   int ulogy();
   
   /// Unset the x axis from log scale
   /** Sends the command:
     * \verbatim
       unset log x
      \endverbatim
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   int ulogx();
   
   /// Unset the x and y axes from log scale
   /** Sends the command:
     * \verbatim
       unset log x
       unset log y
      \endverbatim
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   int ulogxy();
   
   /// Issue a plot command with a n separate lines specified by gplot_spec structures.
   /** Forms the plot command as follows
     * \verbatim
        plot <plot_specs[0]>, /
             <plot_specs[1]>, /
             ---
             <plot_specs[n-1]>
      \endverbatim
     */
   //int plot( std::vector<gplot_spec> plot_specs ); 
   
   /// Plot from a file
   /** Forms the gnuplot plot command as follows:
     * \verbatim
       plot "<fname>" <modifiers>
       \endverbatim
     * The modifiers string can contain any modifiers such as \a using, \a title, etc.
     * 
     * 
     * \retval 0 on success
     * \retval -1 on error
     */
   int plot( const std::string & fname,        ///< [in] the name (with full path) of the file containing data to plot
             const std::string & modifiers ="" ///< [in] [optional] contains any modifiers to the plot command.
           );
   
   /// Plot data from an array
   /** Copies the data in the array to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataT" u 1 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataT>
   int plot( const dataT * y,                   ///< [in] a pointer to an array of data
             size_t N,                          ///< [in] the length of the array
             const std::string & modifiers ="", ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
             const std::string & title =""      ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
           );      
   
   /// Plot data from a vector
   /** Copies the data in the vector to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataT" u 1 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     * 
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataT>
   int plot( const std::vector<dataT> & y,     ///< [in] the vector containing the data
             const std::string & modifiers="", ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
             const std::string & title = ""    ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
           ); 
   
   
   /// Plot y vs. x data from arrays
   /** Copies the data in the arrays to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataTx%dataTy" u 1:2 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     * 
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataTx, typename dataTy>
   int plot( const dataTx * x,                 ///< [in] a pointer to an array of data for the independent variable
             const dataTy * y,                 ///< [in] a pointer to an array of data for the dependent variable
             size_t N,                         ///< [in] the length of the arrays
             const std::string & modifiers="", ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
             const std::string & title = ""    ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
           );
   
   /// Plot y vs. x data from vectors
   /** Copies the data in the vectors to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataTx%dataTy" u 1:2 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataTx, typename dataTy>
   int plot( const std::vector<dataTx> & x,       ///< [in] a vector of data for the independent variable
             const std::vector<dataTy> & y,       ///< [in] a vector of data for the dependent variable
             const std::string & modifiers = "",  ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
             const std::string & title = ""       ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
           );
   
   /// Plot a single point
   /** Copies the position as a length-1 vector to a temporary binary file, and then plots it.
     * 
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataTx, typename dataTy>
   int point( dataTx x,                           ///< [in] the independent axis (x) coordinate of the point
              dataTy y,                           ///< [in] the dependent axis (y) coordinate of the point
              const std::string & modifiers = "", ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
              const std::string & title = ""      ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
            );
   
   /// Draw a circle on the plot
   /** Creates a circle with specified center and radius, and plots it.
     */
   int circle( double xcen,                        ///< [in] the x-coordinate of the circle center
               double ycen,                        ///< [in] the y-coordinate of the circle center
               double radius,                      ///< [in] the circle radius
               const std::string & modifiers = "", ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
               const std::string & title = "",     ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
               int npoints =  10                   ///< [in] [optional] specifies the number of points in each half circle.  Default 10 is usually sufficient.
             );
   
   /// Draw a circle on the plot around the origin
   /** Creates a circle with radius and center (0,0), and plots it.
     */
   int circle( double radius,                      ///< [in] the circle radius
               const std::string & modifiers = "", ///< [in] [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
               const std::string & title = "",     ///< [in] [optional] contains the title of the data set, default is an empty string and no key on the plot
               int npoints =  10                   ///< [in] [optional] specifies the number of points in each half circle.  Default 10 is usually sufficient.
             );
   
protected:
   
   /// Implementation of 1-D binary plotting.
   int plotImpl( const void * y, 
                 size_t Nbytes,
                 const std::string & binary,
                 const std::string & modifiers, 
                 const std::string & title
               );
   
   /// Implementation of 2-d binary plotting.
   int plotImpl( const void * x, 
                 const void * y, 
                 size_t Npts,
                 size_t sizex,
                 size_t sizey,
                 const std::string & binaryx,
                 const std::string & binaryy,
                 const std::string & modifiers, 
                 const std::string & title
               );
   
   ///Open a temporary binary file, and provide the filename. 
   FILE * openTempFile(char * fname);
   

private:
   ///Initialize an instance, only used at construction
   void init();
   
   
};



template<typename dataT>
int gnuPlot::plot( const dataT * y, 
                   size_t N,  
                   const std::string & modifiers, 
                   const std::string & title
                 )
{
   return plotImpl(y, N*sizeof(dataT), gpBinaryFormat<dataT>(), modifiers, title);
}
  
template<typename dataT>
int gnuPlot::plot( const std::vector<dataT> & y, 
                   const std::string & modifiers, 
                   const std::string & title
                 )
{
   return plot<dataT>( y.data(), y.size(), modifiers, title);
}



template<typename dataTx, typename dataTy>
int gnuPlot::plot( const dataTx * x, 
                   const dataTy * y, 
                   size_t N, 
                   const std::string & modifiers, 
                   const std::string & title
                 )
{
   return plotImpl(x, y, N, sizeof(dataTx), sizeof(dataTy), gpBinaryFormat<dataTx>(), gpBinaryFormat<dataTy>(), modifiers, title);
}

template<typename dataTx, typename dataTy>
int gnuPlot::plot( const std::vector<dataTx> & x, 
                   const std::vector<dataTy> & y, 
                   const std::string & modifiers, 
                   const std::string & title
                 )
{
   return plot( x.data(), y.data(), x.size(), modifiers, title);
}

template<typename dataTx, typename dataTy>
int gnuPlot::point( dataTx x, 
                    dataTy y, 
                    const std::string & modifiers, 
                    const std::string & title
                  )
{
   
   return plot(&x, &y, 1, modifiers, title);
}

template<>
std::string gpBinaryFormat<char>();

template<>
std::string gpBinaryFormat<unsigned char>();

template<>
std::string gpBinaryFormat<short>();

template<>
std::string gpBinaryFormat<unsigned short>();

template<>
std::string gpBinaryFormat<int>();

template<>
std::string gpBinaryFormat<unsigned int>();

template<>
std::string gpBinaryFormat<long>();

template<>
std::string gpBinaryFormat<unsigned long>();

template<>
std::string gpBinaryFormat<float>();

template<>
std::string gpBinaryFormat<double>();

}//namespace math
}//namespace mx

   
#endif //math_plot_gnuPlot_hpp

