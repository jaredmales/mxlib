/** \file ds9Interface.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib c++ DS9 interface
  * \ingroup image_processing_files
  *
*/

//***********************************************************************//
// Copyright 2015, 2016, 2017, 2018 Jared R. Males (jaredmales@gmail.com)
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


#ifndef improc_ds9Interface_hpp
#define improc_ds9Interface_hpp

#include <cstring>
#include <iostream>
#include <vector>

#include <unistd.h>

#include <xpa.h>


#include "../ipc/sharedMemSegment.hpp"
#include "fitsUtils.hpp"

#ifndef DS9INTERFACE_NO_EIGEN
#include "eigenImage.hpp"
#endif

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


class ds9Segment : public ipc::sharedMemSegment
{
public:
   size_t dim1 {0};
   size_t dim2 {0};
   size_t dim3 {0};
   int bitpix {0};
   
};

/// An interface to the ds9 image viewer.
/** Handles spawning the ds9 window, and manages shared memory segments,
  * one for each frame.
  *
  *
  * \ingroup image_processing
  * \ingroup plotting
  */
class ds9Interface
{

protected:

   /// The XPA structure to hold connections.
   XPA xpa {NULL};

   ///The title of the window, which sets the XPA access point
   std::string m_title;

   ///The ip:port string to uniquely identify a single XPA access point
   /** We use this so that we don't send to multiple "DS9:ds9" windows, and
     * and always send to the same one.
     */
   std::string m_ipAndPort;

   ///Whether or not the connect() procedure has completed successfully.
   bool m_connected {false};

   ///A vector of shared memory segments, one per frame
   //std::vector<ipc::sharedMemSegment> m_segs;
   std::vector<ds9Segment> m_segs;

   bool m_preserveRegions{true};
   bool m_regionsPreserved {false};

   bool m_preservePan{true};
   bool m_panPreserved {false};


public:

   ///Default c'tor
   ds9Interface();

   ///Constructor which initializes the access point title.
   explicit ds9Interface(const std::string & nn);

   ///Constructor which initializes the access point title.
   explicit ds9Interface(const char * nn);
   
   #ifndef DS9INTERFACE_NO_EIGEN
   /// Constructor which, after initialization, proceeds to display an Eigen-like array in ds9.
   /**
     * see \ref display<typename arrayT>(arrayT & array, int frame=1) for more info.
     */
   template<typename arrayT>
   ds9Interface( const arrayT & array, ///< [in] An array containing an image or cube of images
                 int frame = 1 ///<  [in] [optional] the number of the new frame to initialize.  \note frame must be >= 1.
               );
   #endif //DS9INTERFACE_NO_EIGEN

   ///Destructor.
   /** Calls shutdown, and closes the XPA connections.
     */
   ~ds9Interface();

   ///Initialize the ds9Interface structure
   /** Established a persistent XPA connection.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   void initialize();

   ///Get the title of the ds9Interface
   /**
     * \returns the current value of m_title
     */
   std::string title();

   ///Set the title of the ds9Interface
   /** The title is used as the name of the XPA access point
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int title(const std::string & nn /**< [in] the title to set. */);

   ///Establish the existence of the desired DS9 XPA access point, spawning a new instance if needed.
   /** This isn't really a "connection", the main point is to get the unique ip:port name of the
     * DS9 instance we will be communicating with.
     *
     * \retval 0 on success.
     * \retval -1 on error.
     */
   int connect();

   int XPASet( const char * cmd );

protected:
   ///Spawn (open) the ds9 image viewer
   /** This forks and execs.  An error is returned if exec fails.
     *
     * \retval 0 on sucess
     * \retval -1 on an error.
     */
   int spawn();

   ///Add a segment corresponding to a particular frame in ds9
   /** Nothing is done if the frame already exists.  Note that this does not open a new frame in ds9.
     *
     * \param frame is the number of the new frame to initialize.  Note frame must be >= 1.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int addsegment(size_t frame /**< [in] the number of the new frame to initialize.  \note frame must be >= 1. */);

public:
   ///Open a frame in ds9
   /** Nothing is done if the frame already exists.  First calls \ref addsegment.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int addframe( size_t frame /**< [in] the number of the new frame to initialize.  \note frame must be >= 1. */);

   int togglePreserveRegions( size_t frame,
                              bool onoff
                            );
   
   int togglePreserveRegions(bool onoff);

   int togglePreservePan(bool onoff);

   ///Display an image in ds9.
   /** A new ds9 instance is opened if necessary, and a new sharedmemory segment is added if necessary.
     * The image is described by a pointer and its 2 or 3 dimensions.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     *
     */
   int display( const void *im, ///< [in] the address of the image
                int bitpix,
                size_t pixsz,
                size_t dim1,     ///< [in] the first dimension of the image (in pixels)
                size_t dim2,     ///< [in] the second dimension of the image (in pixels)
                size_t dim3,     ///< [in] the third dimension of the image (in pixels), set to 1 if not a cube.
                int frame = 1    ///< [in] [optional] the number of the new frame to initialize.  \note frame must be >= 1.
              );

   ///Display an image in ds9.
   /** A new ds9 instance is opened if necessary, and a new sharedmemory segment is added if necessary.
     * The image is described by a pointer and its 2 or 3 dimensions.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     *
     */
   template<typename dataT>
   int display( const dataT *im, ///< [in] the address of the image
                size_t dim1,     ///< [in] the first dimension of the image (in pixels)
                size_t dim2,     ///< [in] the second dimension of the image (in pixels)
                size_t dim3,     ///< [in] the third dimension of the image (in pixels), set to 1 if not a cube.
                int frame = 1    ///< [in] [optional] the number of the new frame to initialize.  \note frame must be >= 1.
              );

   #ifndef DS9INTERFACE_NO_EIGEN

   /// Display an Eigen-like array in ds9.
   /** Uses the rows(), cols(), and possibly planes(), methods of arrayT.
     *
     * see \ref display<typename dataT>(const dataT *im, size_t dim1, size_t dim2, size_t dim3, int frame=1) for more.
     *
     * \overload
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   template<typename arrayT>
   int display( const arrayT & array, ///< [in] An array containing an image or cube of images
                int frame = 1 ///<  [in] [optional] the number of the new frame to initialize.  \note frame must be >= 1.
              );

   /// Display an Eigen-like array in ds9.
   /**
     * see \ref display<typename arrayT>(arrayT & array, int frame=1) for more.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   template<typename arrayT>
   int operator()( const arrayT & array, ///< [in] An array containing an image or cube of images
                   int frame = 1 ///<  [in] [optional] the number of the new frame to initialize.  \note frame must be >= 1.
                 );
   #endif //DS9INTERFACE_NO_EIGEN

   
   int loadRegion( size_t frame,
                   const std::string & fname
                 );
   
   int loadRegion( const std::string & fname
                 );
   
   
   ///Shutdown the ds9 interface
   /** Detaches from the shared memory segments.
     *
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int shutdown();

};

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
template<typename arrayT>
ds9Interface::ds9Interface( const arrayT & array,
                            int frame
                          )
{
   initialize();
   display(array, frame);
}
#endif //DS9INTERFACE_NO_EIGEN

ds9Interface::~ds9Interface()
{
   shutdown();
   XPAClose(xpa);
}

inline
void ds9Interface::initialize()
{
   xpa = XPAOpen(NULL);
}

inline
std::string ds9Interface::title()
{
   return m_title;
}

inline
int ds9Interface::title(const std::string & nn)
{
   m_title = nn;
   m_connected = false;
   
   return 0;
}

inline
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

inline
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

inline
int ds9Interface::XPASet( const char * cmd )
{
   if(!m_connected) if(connect() < 0) return -1;

   int rv = ::XPASet(xpa, const_cast<char *>(m_ipAndPort.c_str()), const_cast<char *>(cmd), NULL, NULL, 0, NULL, NULL, 1);

   if(rv != 1)
   {
      std::cerr << "ds9Interface::XPASet: did not send cmd properly.\n";
      return -1;
   }
   
   return 0;
}

inline
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

inline
int ds9Interface::addframe( size_t frame )
{
   char cmd[DS9INTERFACE_CMD_MAX_LENGTH];

   addsegment( frame );

   snprintf(cmd, DS9INTERFACE_CMD_MAX_LENGTH, "frame %zu", frame);

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

inline
int ds9Interface::togglePreserveRegions( bool onoff)
{
   
   for(size_t frame=1; frame< m_segs.size()+1; ++frame)
   {
      int rv = togglePreserveRegions(frame, onoff);
      if(rv < 0) return -1;
   }

   return 0;
}

inline
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

inline
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

template<typename dataT>
int ds9Interface::display( const dataT * im,
                           size_t dim1,
                           size_t dim2,
                           size_t dim3,
                           int frame
                          )
{
   return display(im, getFitsBITPIX<dataT>(), sizeof(dataT), dim1, dim2, dim3, frame);
}

#ifndef DS9INTERFACE_NO_EIGEN
template<typename arrayT>
int ds9Interface::display( const arrayT & array,
                           int frame
                         )
{
   eigenArrPlanes<arrayT> planes;

   return display(array.data(), array.rows(), array.cols(), planes(array), frame);
}

template<typename arrayT>
int ds9Interface::operator()( const arrayT & array,
                              int frame
                            )
{
   return display(array, frame);
}

#endif //DS9INTERFACE_NO_EIGEN

inline
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

inline
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
   
inline
int ds9Interface::shutdown()
{
   size_t i;

   for(i=0; i < m_segs.size(); i++) m_segs[i].detach();

   m_segs.clear();

   return 0;
}

} //namespace improc
} //namespace mx


#endif //improc_ds9Interface_hpp
