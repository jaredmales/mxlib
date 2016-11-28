/** \file gnuPlot.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declaration and definition of an interface to the gnuplot program
  * \ingroup plotting
  * 
*/

#ifndef __gnuPlot_hpp__
#define __gnuPlot_hpp__

//#define _XOPEN_SOURCE 700
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

//For pi<>()
#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include <mx/mxError.hpp>

/** \addtogroup plotting
  * @{
  */
   
namespace mx
{
   
#ifndef MX_GP_FNAME_SZ
///The size of the string for managing the stderr fifo
#define MX_GP_FNAME_SZ (128)
#endif
   
#ifndef MX_GP_TEMP_SZ
///The size of the string for managing temporary file names
#define MX_GP_TEMP_SZ (128)
#endif

#ifndef MX_GP_FC_TIME
/** \def MX_GP_FC_TIME
  * Time, in microseconds, to wait for gnuplot startup and file creation to complete.
  * Opening is retried after each timeout of this length for MX_GP_FC_RETRIES attempts.
  */ 
#define MX_GP_FC_TIME (10000)
#endif
   
#ifndef MX_GP_FC_RETRIES
/** \def MX_GP_FC_RETRIES 
  * Number of times to retry opening the gnuplot stderr file.
  */
#define MX_GP_FC_RETRIES (10)
#endif
   
//Get the gnuplot binary format string for a type
//Specializations below.
template<typename dataT>
std::string gpBinaryFormat()
{
   std::cerr << "No gnuplot format specifier available for this type.\n";
   return "";
}

//Convert a timespec into a double of secons since epoch
double gp_get_curr_time()
{
   struct timespec tsp;
   clock_gettime(CLOCK_REALTIME, &tsp);
   
   return ((double)tsp.tv_sec) + ((double)tsp.tv_nsec)/1e9;
}


/// A c++ interface to gnuplot
/** Spawns a gnuplot sesssion and communicates with it.
  * \ingroup plotting 
  * An example of using gnuPlot to plot data from a file: 
  * \code
  * gnuPlot gp; //This automatically connects, and is now ready to plot.
  * 
  * gp.plot("my_data_file.dat", "u 1:2 w l"); //Plots the data in the file
  * 
  * std::cout << gp.getResponse("show terminal") << "\n"; //Sends the command to gnuplot, gets the response.
  * \endcode
  * 
  * Error checking and reporting is not straightforward since gnuplot does not return a result if there is no error.  If there is an error,
  * it can take several hundred milliseconds for a response to be available.  So for a typical 
  *  plot command, say, one does not want to try to read the stderr ouput just in case as this will just time out with a delay.
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

   int _connected;
   
   ///Set to true if the response indicates a gnuplot error
   bool _gpError;
   std::string _gpErrorMsg;
   
   ///File stream for the gnuplot interface
   FILE * _pipeH;
   
   ///Where to create gnuplot stderr fifo
   /** Default is /dev/shm/
     */
   std::string _errLocation;
   
   ///File name of the gnuplot stderr fifo
   std::string _errFName;
   
   ///File descriptor for the gnuplot stderr fifo
   int _errFD;
   
   ///Location of temporary files 
   /** Default is /dev/shm/
     */
   std::string _tempLocation;
   
   ///Vector of all temporary file names opened, used for removal on destruction.
   std::vector<std::string> _tempFiles;

   ///Flag to control whether temporary files are deleted on destruction.  Default is true (files deleted).
   bool _deleteTemp;

   bool _plotted;
   
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
     * \param[in] com is the command string
     * \param[in] flush [optional] if true, then the output stream is flushed once the command is written.
     *
     * \retval 0 on success
     * \retval -1 on error
     */ 
   int command(const std::string & com, bool flush = true);
     
   ///Check for a response from gnuplot.
   /** It typically takes some time for a response from gnuplot to become readable,
     * whether an error occurs or not. It is often 0.3 seconds or more.  If there is no error, there will be no response so waiting
     * for a response to check for errors after every command can be very time consuming.  Thus it is a choice whether to check 
     * for errors from gnuplot after every command, and what timeout to use. 
     * 
     * gnuplot terminates all outputs with \c \\n\\n, so this reads up to these two characters.  It then strips any leading and trailing whitespace. 
     * 
     * If the response begins with \c gnuplot>, then the response is an error.  In this case, checkResponse returns -1 and the gpError flag is set.  If gpError is not set 
     * but the return value is -1, then some other error occurred (check \c errno).
     * 
     * \param[out] response contains the response from gnuplot, but is empty on timeout or error.
     * \param[in] timeout [optional] is the length of time, in seconds, to wait for a response. Default is 0, but a minimum of 0.5 if a response is expected.
     * 
     * \retval 0 on timeout or successful read
     * \retval -1 on error, or if the gnuplot response indicates an error.
     */
   int checkResponse(std::string & response, double timeout = 0);
   
   ///Get a response from gnuplot for a given command.
   /** This should only be used if a response is expected.
     *
     * \param[in] com is the command string
     * \param[in] timeout [optional] is the length of time, in seconds, to wait for a response. Default is 0.5, which is the minimum that should be used.
     * 
     * \returns "" (empty string) on timeout or error (\c errno and gpError() should be checked), and the response from gnuplot otherwise.
     */ 
   std::string getResponse( const std::string & com, double timeout=0.5);

   
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
     * The modifiers string can contain any modifiers such as \a using, \a title \endverbatim, etc.
     * 
     * \param[in] fname is the name of the file containing data to plot
     * \param[in] modifiers [optional] contains any modifiers to the plot command.
     * 
     * \retval 0 on success
     * \retval -1 on error
     */
   int plot(const std::string & fname, const std::string & modifiers ="");
   
   /// Plot data from an array
   /** Copies the data in the array to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataT" u 1 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     * 
     * \param[in] y is a pointer to an array of data
     * \param[in] N is the length of the array
     * \param[in] modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param[in] title [optional] contains the title of the data set, default is an empty string and no key on the plot
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataT>
   int plot( const dataT * y, size_t N,  const std::string & modifiers ="", const std::string & title ="");      
   
   /// Plot data from a vector
   /** Copies the data in the vector to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataT" u 1 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     * 
     * \param[in] y is the vector containing the data
     * \param[in] modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param[in] title [optional] contains the title of the data set, default is an empty string and no key on the plot
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataT>
   int plot( const std::vector<dataT> & y, const std::string & modifiers="", const std::string & title = "");
   
   
   /// Plot y vs. x data from arrays
   /** Copies the data in the arrays to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataTx%dataTy" u 1:2 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     * 
     * \param[in] x is a pointer to an array of data for the independent variable
     * \param[in] y is a pointer to an array of data for the dependent variable
     * \param[in] N is the length of the arrays
     * \param[in] modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param[in] title [optional] contains the title of the data set, default is an empty string and no key on the plot
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataTx, typename dataTy>
   int plot( const dataTx * x, const dataTy * y, size_t N, const std::string & modifiers="", const std::string & title = "");
   
   /// Plot y vs. x data from vectors
   /** Copies the data in the vectors to a temporary binary file, and then forms the gnuplot plot command as follows:
     * \verbatim 
       plot "temp-file-name" binary format="%dataTx%dataTy" u 1:2 t "title" <modifiers>
      \endverbatim
     * The modifiers string \b must \b NOT contain the \b binary, \b format, \b using, or the \b title  modifiers, but can contain any other modifiers.  Title is
     * specified so that the name of the temporary file name is not printed on the plot. 
     * 
     * \param[in] x is a vector of data for the independent variable
     * \param[in] y is a vector of data for the dependent variable
     * \param[in] modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param[in] title [optional] contains the title of the data set, default is an empty string and no key on the plot
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataTx, typename dataTy>
   int plot( const std::vector<dataTx> & x, const std::vector<dataTy> & y, const std::string & modifiers = "", const std::string & title = "");
   
   template<typename dataT>
   int point( dataT x, dataT y, const std::string & modifiers = "", const std::string & title = "");
   
   template<typename dataT>
   int circle( dataT xcen, dataT ycen, dataT radius, const std::string & modifiers = "", const std::string & title = "", dataT npoints =  10);
   
   template<typename dataT>
   int circle( dataT radius, const std::string & modifiers = "", const std::string & title = "", dataT npoints =  10);
   
protected:
   
      
   ///Open a temporary binary file, and provide the filename. 
   FILE * openTempFile(char * fname);
   

private:
   ///Initialize an instance, only used at construction
   void init();
   
   
};


inline
gnuPlot::gnuPlot()
{
   init();
}

inline
void gnuPlot::init()
{
   _connected = 0;

   _pipeH = 0;

   _errLocation = "/dev/shm";
   _errFD = 0;

   _tempLocation = "/dev/shm";
   
   _plotted = false;
}

inline
gnuPlot::~gnuPlot()
{
   if(_pipeH)
   {
      pclose(_pipeH);
   }
         
   if(_errFD > 0)
   {
      close(_errFD);
   }
   
   if(_errFName.length() > 0)
   {
      remove(_errFName.c_str());
   }
   
   for(int i=0; i < _tempFiles.size(); ++i)
   {
      remove(_tempFiles[i].c_str());
   }
      
}

inline
int gnuPlot::connect()
{

   if(_connected) return 0;

   srand( (unsigned long)  this);  //We don't need random, just different from all other instances in this process
                                   
   
   char tmpStr[MX_GP_FNAME_SZ];

   snprintf(tmpStr, MX_GP_FNAME_SZ, "/gperr%d_%d", getpid(), rand());

   _errFName = _errLocation + tmpStr;
   
   
   std::string comm = "gnuplot -persist 2>";
   comm += _errFName;
   
   errno = 0;
   _pipeH = popen(comm.c_str(), "w");

   if(_pipeH == NULL)
   {
      if(errno)
      {
         mxPError("gnuPlot", errno, "Error starting the gnuplot program with popen.");
      }
      else
      {
         mxError("gnuPlot", MXE_PROCERR, "Error starting the gnuplot program with popen.");
      }
      
      return -1;
   }
   
   
   
   usleep(MX_GP_FC_TIME);//allow file creation to finish before opening stderr file
   
   errno = 0;
   _errFD = open(_errFName.c_str(), O_RDONLY);
   
   int n =0;
   while(_errFD <= 0 && n < MX_GP_FC_RETRIES)
   {
      //First we try again after sleeping again.
      usleep(MX_GP_FC_TIME);//allow file creation to finish before opening stderr file
      _errFD = open(_errFName.c_str(), O_RDONLY);
      ++n;
   }

   if(_errFD <= 0)
   {  
      mxPError("gnuPlot", errno, "gnuPlot failed to open stderr file: ");
      return -1;
   }
      
   _connected = 1;
}

inline
bool gnuPlot::gpError()
{
   return _gpError;
}

inline
std::string gnuPlot::gpErrorMsg()
{
   return _gpErrorMsg;
}
   
inline
int gnuPlot::command(const std::string & com, bool flush)
{
   
   if(!_connected) connect();
   
   fprintf(_pipeH, "%s", (com + "\n").c_str());
   if(flush) fflush(_pipeH);
   
   return 0;
     

}

inline
int gnuPlot::checkResponse(std::string & response, double timeout)
{
   _gpError = false;
   if(_errFD)
   {
      char errstr[1024];
      int rv = 0;
      
      double t0 = gp_get_curr_time();

      errno = 0;
      rv = read(_errFD, errstr, 1024);
      
      while(rv == 0 && gp_get_curr_time() - t0 <= timeout)
      {
         usleep(10); //Give up the thread
         rv = read(_errFD, errstr, 1024);
      }
      
      if(rv < 0)
      {
         mxPError("gnuPlot", errno, "Occurred while reading gnuplot stderr");
         
         response = "";
         return -1;
      }
      
      if(rv == 0)
      {
         mxError("gnuPlot", MXE_TIMEOUT, "Timed out while reading from gnuplot stderr");
         response = "";
         return 0;
      }
      
      errstr[rv] = '\0';
      
      response = errstr;
      
      bool done = false;
      if(response.length() > 1)
      {
         if(response[response.length()-1] == '\n' && response[response.length()-2] == '\n') done = true;
      }
      
      
      while( !done && gp_get_curr_time() - t0 < timeout)
      {
         usleep(10); //Give up the thread
         rv = read(_errFD, errstr, 1024);
         
         if(rv > 0)
         {
            errstr[rv] = '\0';
            response += errstr;
       
            if(response.length() > 1)
            {
               if(response[response.length()-1] == '\n' && response[response.length()-2] == '\n') done = true;
            }
         }
      }

      int first =0;
      while( isspace(response[first]) && first < response.length()-1) ++first;
      response.erase(0, first);
      
      
      int last = response.length()-1;
      while( response[last] == '\n' && last > 0) --last;
      if(response[last] != '\n') ++last;
      
      response.erase(last);
      
      
      if(response.length() > 8)
      {
         if(response.substr(0,8) == "gnuplot>")
         {
            _gpError = true;
            _gpErrorMsg = response;
            
            mxError("gnuPlot", MXE_GNUPLOTERR, "gnuplot says:\n" + response);
            
            response = "";
            return -1;
         }
      }
      
      
      return 0;
   }
   
   //If errFD is not open for some reason
   std::cerr << "gnuplot stderr is not open\n";
   
   
   
   return 1;
   
   
}

inline
std::string gnuPlot::getResponse( const std::string & com, double timeout)
{
   std::string response;
      
   if(command(com, true) < 0) return response;
   
   checkResponse(response, timeout);
   
   return response;
}
   
inline
int gnuPlot::replot()
{
   return command("replot");
}

inline
int gnuPlot::logy()
{
   return command("set logscale y");
}

inline
int gnuPlot::logx()
{
   return command("set logscale x");
}

inline
int gnuPlot::logxy()
{
   int rv;
   rv = logx();
   if(rv !=0) return rv;
   
   return logy();
}

inline
int gnuPlot::ulogy()
{
   return command("unset logscale y");
}

inline
int gnuPlot::ulogx()
{
   return command("unset logscale x");
}

inline
int gnuPlot::ulogxy()
{
   int rv;
   rv = ulogx();
   if(rv !=0) return rv;
   
   return ulogy();
}

int gnuPlot::plot(const std::string & fname, const std::string & modifiers)
{
   std::string com;
   
   if(!_plotted)
   {
      com = "plot ";
   }
   else
   {
      com = "replot ";
   }
   
   com += "\"" + fname;
   
   com += "\" ";
   
   com += modifiers;
 
   _plotted = true;
   
   return command(com, true);
}

template<typename dataT>
int gnuPlot::plot( const dataT * y, size_t N,  const std::string & modifiers, const std::string & title)
{
   FILE * fout;
   char temp[MX_GP_TEMP_SZ];
   
   fout = openTempFile(temp);
   
   if(fout == 0) return -1;
   
   int rv = fwrite(y, sizeof(dataT), N, fout);
   if(rv != N)
   {
      std::cerr << "Error writing to temporary file\n";
      perror("gnuPlot: ");
      return -1;
   }
   fflush(fout);
   
   std::string com;
   
   if(!_plotted)
   {
      com = "plot ";
   }
   else
   {
      com = "replot ";
   }
   com += "\"";
   com += temp;
   com += "\" binary format=\"";
   com +=  gpBinaryFormat<dataT>() + "\" u 1 t \"" + title + "\" " + modifiers;
   
   _plotted = true;
   
   return command(com, true);
}
  
template<typename dataT>
int gnuPlot::plot( const std::vector<dataT> & y, const std::string & modifiers, const std::string & title)
{
   return plot<dataT>( y.data(), y.size(), modifiers, title);
}
   
template<typename dataTx, typename dataTy>
int gnuPlot::plot( const dataTx * x, 
                   const dataTy * y, 
                   size_t N, 
                   const std::string & modifiers, 
                   const std::string & title)
{
   FILE * fout;
   char temp[MX_GP_TEMP_SZ];

   fout = openTempFile(temp);
   
   if(fout == 0) return -1;

   
   for(int i=0; i< N; ++i)
   {
      int rv = fwrite(&x[i], sizeof(dataTx), 1, fout);
      rv += fwrite(&y[i], sizeof(dataTy), 1, fout);
      
      if(rv != 2)
      {
         std::cerr << "Error writing to temporary file\n";
         perror("gnuPlot: ");
         return -1;
      }
   }
   fflush(fout);
   
   std::string com;
   
   if(!_plotted)
   {
      com = "plot ";
   }
   else
   {
      com = "replot ";
   }
   com += "\"";
   com += temp;
   com += "\" binary format=\"";
   com +=  gpBinaryFormat<dataTx>() + gpBinaryFormat<dataTy>() + "\" u 1:2 t \"" + title + "\" " + modifiers;
         
   _plotted = true;
   return command(com, true);
}

template<typename realT>
int gnuPlot::point( realT x, realT y, const std::string & modifiers, const std::string & title)
{
   
   return plot(&x, &y, 1, modifiers, title);
}

template<typename realT>
int gnuPlot::circle( realT xcen, realT ycen, realT radius, const std::string & modifiers, const std::string & title, realT npoints)
{
   int halfCirc = (0.5 * 2.0 * pi<realT>() * radius * npoints + 0.5);
   
   std::vector<realT> xpts(2*halfCirc), ypts(2*halfCirc);
   
   realT x, y;
   
   for(int i=0; i < halfCirc; ++i)
   {
      x = radius * cos(( (realT) i) / ( (realT) halfCirc) * pi<realT>());
      
      y = sqrt( radius*radius - x*x);
      
      xpts[i] = xcen + x;
      ypts[i] = ycen + y;
      
      xpts[halfCirc + i] = xcen - x;
      ypts[halfCirc + i] = ycen - y;
   }
   
   return plot(xpts, ypts, modifiers, title);
}
   
template<typename dataT>
int gnuPlot::circle( dataT radius, const std::string & modifiers, const std::string & title, dataT npoints)
{
   return circle<dataT>(0.0, 0.0, modifiers, title, npoints);
}
   
   
template<typename dataTx, typename dataTy>
int gnuPlot::plot( const std::vector<dataTx> & x, const std::vector<dataTy> & y, const std::string & modifiers, const std::string & title)
{
   return plot( x.data(), y.data(), x.size(), modifiers, title);
}

FILE * gnuPlot::openTempFile(char * temp)
{
   FILE * fout;
   
   snprintf(temp, MX_GP_TEMP_SZ, "%s/gpplot_%d_XXXXXX", _tempLocation.c_str(), getpid());
   mkstemp(temp);
   
   fout = fopen(temp, "wb");

   if(fout == NULL)
   {
      std::cerr << "Could not open tempoary file (" << temp << ") for writing\n";
      perror("gnuPlot: ");
      return 0;
   }
   
   _tempFiles.push_back(temp);
   
   return fout;
} //FILE * gnuPlot::openTempFile(char * temp)

template<>
std::string gpBinaryFormat<char>()
{
   return "%char";
}

template<>
std::string gpBinaryFormat<unsigned char>()
{
   return "%uchar";
}

template<>
std::string gpBinaryFormat<short>()
{
   return "%short";
}

template<>
std::string gpBinaryFormat<unsigned short>()
{
   return "%ushort";
}

template<>
std::string gpBinaryFormat<int>()
{
   return "%int";
}

template<>
std::string gpBinaryFormat<unsigned int>()
{
   return "%uint";
}

template<>
std::string gpBinaryFormat<long>()
{
   return "%long";
}

template<>
std::string gpBinaryFormat<unsigned long>()
{
   return "%ulong";
}

template<>
std::string gpBinaryFormat<float>()
{
   return "%float";
}

template<>
std::string gpBinaryFormat<double>()
{
   return "%double";
}


}//namespace mx

/** @}
  */ //addtogroup plotting
   
#endif //__gnuPlot_hpp__

