/** \file ds9Interface.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Implementation of the mxlib c++ DS9 interface
  * \ingroup image_processing_files
  *
*/

//***********************************************************************//
// Copyright 2021 Jared R. Males (jaredmales@gmail.com)
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


#include "improc/ds9Interface.hpp"

namespace mx
{
namespace improc
{

#ifndef DS9INTERFACE_SPAWN_SLEEP
/// The time to sleep after spawning, in msecs, while checking for the new instance to be ready.
/**
  * \ingroup image_processing
  * \ingroup plotting
  */
#define DS9INTERFACE_SPAWN_SLEEP (100)
#endif

#ifndef DS9INTERFACE_SPAWN_TIMEOUT
/// The maximum time to wait after spawning, in msecs, before giving up.
/**
  * \ingroup image_processing
  * \ingroup plotting
  */
#define DS9INTERFACE_SPAWN_TIMEOUT (2000)
#endif

#ifndef DS9INTERFACE_CMD_MAX_LENGTH
/// The maximum length of a ds9 command
/**
  * \ingroup image_processing
  * \ingroup plotting
  */
#define DS9INTERFACE_CMD_MAX_LENGTH (512)
#endif

ds9Interface::ds9Interface()
{
   initialize();
}

ds9Interface::ds9Interface(const std::string & nn)
{
   initialize();
   title(nn);
}

ds9Interface::ds9Interface(const char * nn)
{
   initialize();
   title(nn);
}



#ifndef DS9INTERFACE_NO_EIGEN

template ds9Interface::ds9Interface<eigenImage<signed char>>( const eigenImage<signed char> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<unsigned char>>( const eigenImage<unsigned char> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<short>>( const eigenImage<short> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<unsigned short>>( const eigenImage<unsigned short> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<int>>( const eigenImage<int> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<unsigned int>>( const eigenImage<unsigned int> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<long>>( const eigenImage<long> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<unsigned long>>( const eigenImage<unsigned long> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<float>>( const eigenImage<float> & array, int frame);
template ds9Interface::ds9Interface<eigenImage<double>>( const eigenImage<double> & array, int frame);

#endif //DS9INTERFACE_NO_EIGEN

ds9Interface::~ds9Interface()
{
   shutdown();
   XPAClose(xpa);
}

void ds9Interface::initialize()
{
   xpa = XPAOpen(NULL);
}

std::string ds9Interface::title()
{
   return m_title;
}

int ds9Interface::title(const std::string & nn)
{
   m_title = nn;
   m_connected = false;
   
   return 0;
}

int ds9Interface::connect()
{
   
   if(xpa) 
   {
      shutdown();
      xpa = NULL;
   }
   
   if(xpa == NULL)
   {
      xpa = XPAOpen(NULL);
   }
   
   int  n = 1;
   char *names[1];
   names[0] = NULL;

   std::string tmpl;
   if( m_title.find(':', 0) == std::string::npos)
   {
      if(m_title == "") m_title = "ds9";

      tmpl = "DS9:";
      tmpl += m_title;
   }
   else
   {
      tmpl = m_title;
   }

   char paramlist[] = "gsi";

   int rv = XPAAccess(xpa, const_cast<char *>(tmpl.c_str()), paramlist, NULL, names, NULL, n);

   if(rv == 0)
   {
      if( spawn() != 0) return -1;
      if(names[0]) 
      {
         free(names[0]);
         names[0] = NULL;
      }
      
      int slept = 0;
      while(rv == 0 && slept < DS9INTERFACE_SPAWN_TIMEOUT)
      {
         rv = XPAAccess(xpa, const_cast<char *>(tmpl.c_str()), paramlist, NULL, names, NULL, n);

         usleep(DS9INTERFACE_SPAWN_SLEEP*1000);
         slept += DS9INTERFACE_SPAWN_SLEEP;
      }

      if(rv == 0)
      {
         std::cerr << "ds9Interface: failed to connect after attempting to spawn.  Timed out.\n";
         if(names[0]) free(names[0]);
         return -1;
      }
   }

   char * st = strchr(names[0], ' ');

   m_ipAndPort = st + 1;

   free(names[0]);

   m_connected = true;

   return 0;
}

int ds9Interface::spawn()
{
   if(m_title == "" || m_title == "*") m_title = "ds9";

   int pid = fork();
   if (pid==0)
   {
      /* Create a new SID for the child process so it is detached*/
      pid_t sid = setsid();
     
      if (sid < 0) 
      {
         std::cerr << "ds9Interface: failed to detach.\n";
         perror("ds9Interface");
      }
      
      errno = 0;
      execlp("ds9", "ds9", "-title", m_title.c_str(), (char *) 0);

      std::cerr << "ds9Interface: spawning failed, execlp returned.\n";
      perror("ds9Interface");

      return -1;
   }

   return 0;
}

int ds9Interface::XPASet( const char * cmd )
{
   if(!m_connected) if(connect() < 0) return -1;

   int rv = ::XPASet(xpa, const_cast<char *>(m_ipAndPort.c_str()), const_cast<char *>(cmd), NULL, NULL, 0, NULL, NULL, 1);

   std::cerr << cmd << "\n";
   if(rv != 1)
   {
      std::cerr << "ds9Interface::XPASet: did not send cmd properly.\n";
      return -1;
   }
   
   return 0;
}

int ds9Interface::addsegment( size_t frame )
{
   size_t curr_n;

   if(frame == 0) return -1;
   
   if((size_t)(frame-1) < m_segs.size()) return 0;

   curr_n = m_segs.size();

   m_segs.resize(frame);

   for(size_t i = curr_n; i< m_segs.size(); ++i)
   {
      m_segs[i].initialize();
      m_segs[i].setKey(0, IPC_PRIVATE);
   }

   return 0;

}

int ds9Interface::addframe( size_t frame )
{
   char cmd[DS9INTERFACE_CMD_MAX_LENGTH];

   addsegment( frame );

   snprintf(cmd, DS9INTERFACE_CMD_MAX_LENGTH, "frame %zu", frame);

   std::cerr << cmd << "\n";
   
   if(!m_connected) if(connect() < 0) return -1;

   int rv = XPASet(cmd);

   if(rv != 0)
   {
      std::cerr << "ds9Interface: could not add frame.\n";
      m_connected = false;
      return -1;
   }

   return 0;
}

int ds9Interface::togglePreserveRegions( bool onoff)
{
   
   for(size_t frame=1; frame< m_segs.size()+1; ++frame)
   {
      int rv = togglePreserveRegions(frame, onoff);
      if(rv < 0) return -1;
   }

   return 0;
}

int ds9Interface::togglePreserveRegions( size_t frame,
                                         bool onoff
                                       )
{
   int rv;

   char cmd[DS9INTERFACE_CMD_MAX_LENGTH];
   
   snprintf(cmd, DS9INTERFACE_CMD_MAX_LENGTH, "frame %zu", frame);
   rv = XPASet(cmd);
 
   if(rv < 0)
   {
      std::cerr << "ds9Interface::preserveRegions: error sending frame." << "\n";
      return -1;
   }

   if(onoff == true)
   {
      rv = XPASet("preserve regions yes");
      m_regionsPreserved = true;
   }
   else
   {
      rv = XPASet("preserve regions no");
      m_regionsPreserved = false;
   }

   if(rv < 0)
   {
      std::cerr << "ds9Interface::preserveRegions: error sending preserve regions." << "\n";
      return -1;
   }

   return 0;
}

int ds9Interface::togglePreservePan(bool onoff)
{
   int rv;

   if(onoff == true)
   {
      rv = XPASet("preserve pan yes");
      m_panPreserved = true;
   }
   else
   {
      rv = XPASet("preserve pan no");
      m_panPreserved = false;
   }

   if(rv < 0)
   {
      std::cerr << "ds9Interface::preserveRegions: error sending preserve regions." << "\n";
      return -1;
   }

   return 0;
}

int ds9Interface::display( const void * im,
                           int bitpix,
                           size_t pixsz,
                           size_t dim1,
                           size_t dim2,
                           size_t dim3,
                           int frame
                          )
{
   size_t tot_size;
   char cmd[DS9INTERFACE_CMD_MAX_LENGTH];

   if(frame < 1)
   {
      std::cerr <<  "ds9Interface: frame must >= 1\n" << "\n";
      return -1;
   }

   if(!m_connected) if(connect() < 0) return -1;

   if(addframe(frame) < 0) 
   {
      m_connected = false;
      return -1;
   }
   //Calculate total size
   tot_size= pixsz;
   tot_size*=dim1;
   tot_size*=dim2;
   tot_size*=dim3;
   
   bool realloc = false;

   //Re-allocate shared memory if necessary
   if(tot_size > m_segs[frame-1].size)
   {
      if( m_segs[frame-1].size > 0 )
      {
         m_segs[frame-1].detach();
      }
      m_segs[frame-1].create(tot_size);
      
      realloc = true;
   }
   else
   {
      if( dim1 != m_segs[frame-1].dim1 || dim2 != m_segs[frame-1].dim2 || dim3 != m_segs[frame-1].dim3 || bitpix != m_segs[frame-1].bitpix)
      {
         realloc = true; //force a new shm command
      }
   }

   memcpy( m_segs[frame-1].addr, im, tot_size );

   m_segs[frame-1].dim1 = dim1;
   m_segs[frame-1].dim2 = dim2;
   m_segs[frame-1].dim3 = dim3;
   m_segs[frame-1].bitpix = bitpix;
   
   if(realloc)
   {
      //Handle single image so that the cube dialog doesn't open up if dim3=1
      if(dim3 == 1)
      {
         snprintf(cmd, DS9INTERFACE_CMD_MAX_LENGTH, "shm array shmid %i [xdim=%zu,ydim=%zu,bitpix=%i]",
                                         m_segs[frame-1].shmemid,
                                        dim1, dim2, bitpix);
      }
      else
      {
         snprintf(cmd, DS9INTERFACE_CMD_MAX_LENGTH, "shm array shmid %i [xdim=%zu,ydim=%zu,zdim=%zu,bitpix=%i]",
                                         m_segs[frame-1].shmemid,
                                        dim1, dim2, dim3, bitpix);
      }
   }
   else
   {
      snprintf(cmd, DS9INTERFACE_CMD_MAX_LENGTH, "update");
   }

   int rv = XPASet(cmd);

   if(rv != 0)
   {
      std::cerr << "ds9Interface: sending shm array command to ds9 failed.\n";
      m_connected = false;
      return -1;
   }

   if( m_regionsPreserved != m_preserveRegions ) togglePreserveRegions(m_preserveRegions);
   if( m_panPreserved != m_preservePan ) togglePreservePan(m_preservePan);


   return 0;

}

template int ds9Interface::display<signed char>(const signed char * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<unsigned char>(const unsigned char * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<short>(const short * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<unsigned short>(const unsigned short * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<int>(const int * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<unsigned int>(const unsigned int * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<long>(const long * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<unsigned long>(const unsigned long * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<float>(const float * im, size_t dim1, size_t dim2, size_t dim3, int frame);
template int ds9Interface::display<double>(const double * im, size_t dim1, size_t dim2, size_t dim3, int frame);


#ifndef DS9INTERFACE_NO_EIGEN

template int ds9Interface::display<eigenImage<signed char>>( const eigenImage<signed char> & im, int frame);
template int ds9Interface::display<eigenImage<unsigned char>>( const eigenImage<unsigned char> & im, int frame);
template int ds9Interface::display<eigenImage<short>>( const eigenImage<short> & im, int frame);
template int ds9Interface::display<eigenImage<unsigned short>>( const eigenImage<unsigned short> & im, int frame);
template int ds9Interface::display<eigenImage<int>>( const eigenImage<int> & im, int frame);
template int ds9Interface::display<eigenImage<unsigned int>>( const eigenImage<unsigned int> & im, int frame);
template int ds9Interface::display<eigenImage<long>>( const eigenImage<long> & im, int frame);
template int ds9Interface::display<eigenImage<unsigned long>>( const eigenImage<unsigned long> & im, int frame);
template int ds9Interface::display<eigenImage<float>>( const eigenImage<float> & im, int frame);
template int ds9Interface::display<eigenImage<double>>( const eigenImage<double> & im, int frame);

template int ds9Interface::operator()<eigenImage<signed char>>( const eigenImage<signed char> & im, int frame);
template int ds9Interface::operator()<eigenImage<unsigned char>>( const eigenImage<unsigned char> & im, int frame);
template int ds9Interface::operator()<eigenImage<short>>( const eigenImage<short> & im, int frame);
template int ds9Interface::operator()<eigenImage<unsigned short>>( const eigenImage<unsigned short> & im, int frame);
template int ds9Interface::operator()<eigenImage<int>>( const eigenImage<int> & im, int frame);
template int ds9Interface::operator()<eigenImage<unsigned int>>( const eigenImage<unsigned int> & im, int frame);
template int ds9Interface::operator()<eigenImage<long>>( const eigenImage<long> & im, int frame);
template int ds9Interface::operator()<eigenImage<unsigned long>>( const eigenImage<unsigned long> & im, int frame);
template int ds9Interface::operator()<eigenImage<float>>( const eigenImage<float> & im, int frame);
template int ds9Interface::operator()<eigenImage<double>>( const eigenImage<double> & im, int frame);

#endif //DS9INTERFACE_NO_EIGEN

int ds9Interface::loadRegion( size_t frame,
                              const std::string & fname
                            )
{
   std::string cmd;
   
   cmd = "frame " + std::to_string(frame);
   int rv = XPASet(cmd.c_str());
 
   if(rv < 0)
   {
      std::cerr << "ds9Interface::loadRegion: error sending frame." << "\n";
      return -1;
   }
   
   cmd = "regions load " + fname;
   
   rv = XPASet(cmd.c_str());
   
   if(rv < 0)
   {
      std::cerr << "ds9Interface::loadRegion: error loading region." << "\n";
      return -1;
   }
   
   return 0;
}

int ds9Interface::loadRegion( const std::string & fname )
{
   std::string cmd;
      
   cmd = "regions load all " + fname;
   
   int rv = XPASet(cmd.c_str());
   
   if(rv < 0)
   {
      std::cerr << "ds9Interface::loadRegion: error loading region." << "\n";
      return -1;
   }
   
   return 0;
}
   
int ds9Interface::shutdown()
{
   size_t i;

   for(i=0; i < m_segs.size(); i++) m_segs[i].detach();

   m_segs.clear();

   return 0;
}

} //namespace improc
} //namespace mx
