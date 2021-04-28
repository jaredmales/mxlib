/** \file gnuPlot.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Implementation of an interface to the gnuplot program
  * \ingroup plotting_files
  * 
*/

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "math/plot/gnuPlot.hpp"
   
namespace mx
{
namespace math
{

gnuPlot::gnuPlot()
{
   init();
}

void gnuPlot::init()
{
   _errLocation = "/dev/shm";
   _tempLocation = "/dev/shm";

}

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
   
   for(size_t i=0; i < _tempFiles.size(); ++i)
   {
      remove(_tempFiles[i].c_str());
   }
      
}

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
   
   return 0;
}

bool gnuPlot::gpError()
{
   return _gpError;
}

std::string gnuPlot::gpErrorMsg()
{
   return _gpErrorMsg;
}
   
int gnuPlot::command( const std::string & com, 
                      bool flush
                    )
{
   
   if(!_connected) connect();
   
   fprintf(_pipeH, "%s", (com + "\n").c_str());
   if(flush) fflush(_pipeH);
   
   return 0;
     

}

int gnuPlot::checkResponse( std::string & response, 
                            double timeout
                          )
{
   _gpError = false;
   if(_errFD)
   {
      char errstr[1024];
      int rv = 0;
      
      double t0 = sys::get_curr_time();

      errno = 0;
      rv = read(_errFD, errstr, 1024);
      
      while(rv == 0 && sys::get_curr_time() - t0 <= timeout)
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
      
      
      while( !done && sys::get_curr_time() - t0 < timeout)
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

      size_t first =0;
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

std::string gnuPlot::getResponse( const std::string & com, 
                                  double timeout
                                )
{
   std::string response;
      
   if(command(com, true) < 0) return response;
   
   checkResponse(response, timeout);
   
   return response;
}
   
int gnuPlot::replot()
{
   return command("replot");
}

int gnuPlot::logy()
{
   return command("set logscale y");
}

int gnuPlot::logx()
{
   return command("set logscale x");
}

int gnuPlot::logxy()
{
   int rv;
   rv = logx();
   if(rv !=0) return rv;
   
   return logy();
}

int gnuPlot::ulogy()
{
   return command("unset logscale y");
}

int gnuPlot::ulogx()
{
   return command("unset logscale x");
}

int gnuPlot::ulogxy()
{
   int rv;
   rv = ulogx();
   if(rv !=0) return rv;
   
   return ulogy();
}

int gnuPlot::plot( const std::string & fname, 
                   const std::string & modifiers
                 )
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

int gnuPlot::circle( double xcen, 
                     double ycen, 
                     double radius, 
                     const std::string & modifiers, 
                     const std::string & title, 
                     int npoints
                   )
{
   int halfCirc = (0.5 * 2.0 * math::pi<double>() * radius * npoints + 0.5);
   
   std::vector<double> xpts(2*halfCirc), ypts(2*halfCirc);
   
   
   for(int i=0; i < halfCirc; ++i)
   {
      double x = radius * cos(( (double) i) / ( (double) halfCirc) * math::pi<double>());
      
      double y = sqrt( radius*radius - x*x);
      
      xpts[i] = xcen + x;
      ypts[i] = ycen + y;
      
      xpts[halfCirc + i] = xcen - x;
      ypts[halfCirc + i] = ycen - y;
   }
   
   return plot(xpts, ypts, modifiers, title);
}
   
int gnuPlot::circle( double radius, 
                     const std::string & modifiers, 
                     const std::string & title, 
                     int npoints
                   )
{
   return circle(0.0, 0.0, radius, modifiers, title, npoints);
}

int gnuPlot::plotImpl( const void * y, 
                       size_t Nbytes,
                       const std::string & binary,
                       const std::string & modifiers, 
                       const std::string & title
                     )
{
   FILE * fout;
   char temp[MX_GP_TEMP_SZ];
   
   fout = openTempFile(temp);
   
   if(fout == 0) return -1;
   
   int rv = fwrite(y, 1, Nbytes, fout);
   if(rv != Nbytes)
   {
      std::cerr << "Error writing to temporary file\n";
      perror("gnuPlot: ");
      return -1;
   }
   fflush(fout);
   fclose(fout);
   
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
   com +=  binary + "\" u 1 t \"" + title + "\" " + modifiers;
   
   _plotted = true;
   
   return command(com, true);
}

int gnuPlot::plotImpl( const void * x, 
                       const void * y, 
                       size_t Npts,
                       size_t sizex,
                       size_t sizey,
                       const std::string & binaryx,
                       const std::string & binaryy,
                       const std::string & modifiers, 
                       const std::string & title
                     )
{
   FILE * fout;
   char temp[MX_GP_TEMP_SZ];

   fout = openTempFile(temp);
   
   if(fout == 0) return -1;

   
   for(int i=0; i< Npts; ++i)
   {
      int rv = fwrite( (char *) x + i*sizex, sizex, 1, fout);
      rv += fwrite( (char *) y + i*sizey, sizey, 1, fout);
      
      if(rv != 2)
      {
         std::cerr << "Error writing to temporary file\n";
         perror("gnuPlot: ");
         return -1;
      }
   }
   fflush(fout);
   fclose(fout);
   
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
   com +=  binaryx + binaryy + "\" u 1:2 t \"" + title + "\" " + modifiers;
         
   _plotted = true;
   return command(com, true);
}

FILE * gnuPlot::openTempFile(char * temp)
{
   FILE * fout;
   
   snprintf(temp, MX_GP_TEMP_SZ, "%s/gpplot_%d_XXXXXX", _tempLocation.c_str(), getpid());
   int rv = mkstemp(temp);
   if(rv < 0) return 0;
      
   close(rv);
   
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

} //namespace math
} //namespace mx

   

