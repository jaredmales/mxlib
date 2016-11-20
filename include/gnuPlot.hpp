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

#include <stdlib.h>

#include <string>
#include <iostream>
#include <vector>

#include <cmath>

/** \addtogroup plotting
  * @{
  */
   
namespace mx
{
   
#define MX_GP_FNAME_SZ 64
#define MX_GP_TEMP_SZ 128

template<typename dataT>
std::string gpBinaryFormat()
{
   return "";
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

/// A c++ interface to gnuplot
/** Spawns a gnuplot sesssion and communicates with it.
  */
class gnuPlot
{
protected:   
   
   int _initialized;
   
   FILE * _pipeH;
   
   char _errFName[MX_GP_FNAME_SZ];
   
   int _errFD;
   
   std::vector<std::string> _tempFiles;
   
public:
   
   gnuPlot();
   
   ~gnuPlot();
   
   void init();
   
   int command(const std::string & com, bool flush = true);
   
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
   
   /// Plot from a file
   /** Forms the gnuplot plot command as follows:
     * \verbatim
       plot "<fname>" <modifiers>
       \endverbatim
     * The modifiers string can contain any modifiers such as \a using, \a title \endverbatim, etc.
     * 
     * \param fname [in] is the name of the file containing data to plot
     * \param modifiers [in] [optional] contains any modifiers to the plot command.
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
     * \param y is a pointer to an array of data
     * \param N is the length of the array
     * \param modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param title [optional] contains the title of the data set, default is an empty string and no key on the plot
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
     * \param y is the vector containing the data
     * \param modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param title [optional] contains the title of the data set, default is an empty string and no key on the plot
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
     * \param x is a pointer to an array of data for the independent variable
     * \param y is a pointer to an array of data for the dependent variable
     * \param N is the length of the arrays
     * \param modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param title [optional] contains the title of the data set, default is an empty string and no key on the plot
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
     * \param x is a vector of data for the independent variable
     * \param y is a vector of data for the dependent variable
     * \param modifiers [optional] contains any modifiers to the plot command other than \b binary, \b format, \b using, and \b title.
     * \param title [optional] contains the title of the data set, default is an empty string and no key on the plot
     *
     * \retval 0 on success
     * \retval -1 on error
     */
   template<typename dataTx, typename dataTy>
   int plot( const std::vector<dataTx> & x, const std::vector<dataTy> & y, const std::string & modifiers = "", const std::string & title = "");
   
   
   
   
};


inline
gnuPlot::gnuPlot()
{
   _pipeH = 0;
   
   _errFName[0] = 0;
   _errFD = 0;

   _initialized = 0;
   
   init();
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
   
   if(_errFName[0] > 0)
   {
      remove(_errFName);
   }
   
   for(int i=0; i < _tempFiles.size(); ++i)
   {
      remove(_tempFiles[i].c_str());
   }
      
}

inline
void gnuPlot::init()
{

   if(_initialized) return;
   
   snprintf(_errFName, MX_GP_FNAME_SZ, "/tmp/gperr%d", getpid());
   
   if (mkfifo(_errFName, 0600)) 
   {
      if (errno != EEXIST) 
      {
         fprintf(stderr, "could not create gnuplot stderr fifo\n");
      }
   }
   else
   {
      _errFD = open(_errFName, O_RDONLY  | O_NONBLOCK); 
   }
   
   char comm[2*MX_GP_FNAME_SZ];
   snprintf(comm, 2*MX_GP_FNAME_SZ, "gnuplot -persist 2>%s", _errFName);
   
   _pipeH = popen(comm, "w");
   
   _initialized = 1;
}

inline
int gnuPlot::command(const std::string & com, bool flush)
{
   
   if(!_initialized) init();
   
   fprintf(_pipeH, "%s", (com + "\n").c_str());
   if(flush) fflush(_pipeH);
   usleep(1e6/10); //give time for response
         
   if(_errFD)
   {
      char errstr[1024];
      int rv = read(_errFD, errstr, 1024);
   
      if(rv < 0 && !( errno == EAGAIN || errno == EWOULDBLOCK))
      {
         perror("Could not read gnuplot stderr");
         //std::cerr << "Could not read gnuplot stderr\n";
         return 1;
      }
      
      errstr[rv] = '\0';
      
      if(rv > 0) 
      {
         std::cerr << "gnuplot error\n";
         std::cerr << errstr << "\n";
         
         return -1;
      }
      return 0;
   }
   
   //If errFD is not open for some reason
   std::cerr << "gnuplot stderr is not open\n";
   return 1;
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
   std::string com = "plot \"" + fname;
   
   com += "\" ";
   
   com += modifiers;
   
   return command(com, true);
}

template<typename dataT>
int gnuPlot::plot( const dataT * y, size_t N,  const std::string & modifiers, const std::string & title)
{
   FILE * fout;
   char temp[MX_GP_TEMP_SZ];
   
   snprintf(temp, MX_GP_TEMP_SZ, "/dev/shm/gpplot_%d_XXXXXX", getpid());
   mkstemp(temp);
   
   fout = fopen(temp, "wb");

   if(fout == NULL)
   {
      std::cerr << "Could not open /dev/shm file for writing\n";
      perror("gnuPlot: ");
      return -1;
   }
   
   _tempFiles.push_back(temp);
   
   int rv = fwrite(y, sizeof(dataT), N, fout);
   if(rv != N)
   {
      std::cerr << "Error writing to /dev/shm\n";
      perror("gnuPlot: ");
      return -1;
   }
   fflush(fout);
   
   std::string com = "plot \"";
   com += temp;
   com += "\" binary format=\"";
   com +=  gpBinaryFormat<dataT>() + "\" u 1 t \"" + title + "\" " + modifiers;
         
   return command(com, true);
}
  
template<typename dataT>
int gnuPlot::plot( const std::vector<dataT> & y, const std::string & modifiers, const std::string & title)
{
   return plot<dataT>( y.data(), y.size(), modifiers, title);
}
   
template<typename dataTx, typename dataTy>
int gnuPlot::plot( const dataTx * x, const dataTy * y, size_t N, const std::string & modifiers, const std::string & title)
{
   FILE * fout;
   char temp[MX_GP_TEMP_SZ];
   
   snprintf(temp, MX_GP_TEMP_SZ, "/dev/shm/gpplot_%d_XXXXXX", getpid());
   mkstemp(temp);
   
   fout = fopen(temp, "wb");

   if(fout == NULL)
   {
      std::cerr << "Could not open /dev/shm file for writing\n";
      perror("gnuPlot: ");
      return -1;
   }
   
   _tempFiles.push_back(temp);
   
   for(int i=0; i< N; ++i)
   {
      int rv = fwrite(&x[i], sizeof(dataTx), 1, fout);
      rv += fwrite(&y[i], sizeof(dataTy), 1, fout);
      
      if(rv != 2)
      {
         std::cerr << "Error writing to /dev/shm\n";
         perror("gnuPlot: ");
         return -1;
      }
   }
   fflush(fout);
   
   std::string com = "plot \"";
   com += temp;
   com += "\" binary format=\"";
   com +=  gpBinaryFormat<dataTx>() + gpBinaryFormat<dataTy>() + "\" u 1:2 t \"" + title + "\" " + modifiers;
         
   return command(com, true);
   
   
}
 
template<typename dataTx, typename dataTy>
int gnuPlot::plot( const std::vector<dataTx> & x, const std::vector<dataTy> & y, const std::string & modifiers, const std::string & title)
{
   return plot( x.data(), y.data(), x.size(), modifiers, title);
}

}//namespace mx

/** @}
  */ //addtogroup plotting
   
#endif //__gnuPlot_hpp__

