/** \file ds9Interface.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib c++ ds9 interface
  * \ingroup image_processing_files
  * 
*/

#ifndef __ds9Interface_hpp__
#define __ds9Interface_hpp__

#include <signal.h>
#include <sys/wait.h>

#include "fitsUtils.hpp"
#include "eigenImage.hpp"

//#include "../eigenUtils.hpp"

#include "../IPC.h"


namespace mx
{
namespace improc 
{
   
///The maximum length of a ds9 command
/** 
  * \ingroup image_processing
  * \ingroup plotting
  */ 
#define DS9_CMD_MAX_LENGTH 512
   
   
/// An interface to the ds9 image viewer.
/** Handles spawning the ds9 window, and manages shared memory segments, 
  * one for each frame.  
  * 
  * This makes use of system() calls.  In each case the result is checked and if SIGINT or SIGQUIT was caught
  * by the ds9 child process, it will be re-raised in the calling process.
  * 
  * \ingroup image_processing
  * \ingroup plotting
  */   
class ds9Interface
{
   
protected:
   
   ///The title of the window, which sets the XPA access point 
   std::string _title;
   
   ///The port to use (normally 0)
   int _port;
   
   ///An array of shared memory segments, one per frame
   sharedmem_segment * segs;
   
   ///The number of segments in \ref segs
   size_t nsegs;
   
public:
   
   ///Default c'tor
   ds9Interface();
   
   /// Constructor which, after initialization, proceeds to display an Eigen-like array in ds9.
   /** 
     * see \ref display<typename arrayT>(arrayT & array, int frame=1) for more.
     */
   template<typename arrayT>
   ds9Interface( const arrayT & array, ///< [in] An array containing an image or cube of images
                 int frame = 1 ///<  [in] [optional] the number of the new frame to initialize.  \note frame must be >= 1.
               );
   
   
   ///Destructor
   ~ds9Interface();
   
   ///Initialize the ds9Interface structure
   /** The title is set to "ds9", the port is set to 0, segs is set to NULL, and nsegs to 0.
     * 
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int init();

   
   ///Get the title of the ds9Interface
   /**
     * \returns the current value of _title
     */ 
   std::string title();
   
   ///Set the title of the ds9Interface 
   /** The title is used as the name of the XPA access point
     * 
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int title( const std::string & new_title /**< [in] the title to set. */);

   ///Get the port of the ds9Interface
   /**
     * \returns the current value of _port
     */ 
   int port();
   
   ///Set the port of the ds9Interface 
   /** 
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int port( const int & new_port /**< [in] the port to set. */);
   
   ///Spawn (open) the ds9 image viewer
   /** This makes use of system() calls.  In each case the result is checked and if SIGINT or SIGQUIT was caught
     * by the ds9 child process, it will be re-raised in the calling process.
     * 
     * \retval 0 on sucess
     * \retval -1 on an error, which may be caused by a response timeout.
     */
   int spawn();

   ///Add a segment corresponding to a particular frame in ds9
   /** Nothing is done if the frame already exists.  Note that this does not open a new frame in ds9.
     * 
     * This makes use of system() calls.  In each case the result is checked and if SIGINT or SIGQUIT was caught
     * by the ds9 child process, it will be re-raised in the calling process.
     * 
     * \param frame is the number of the new frame to initialize.  Note frame must be >= 1.
     * 
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int addsegment(int frame /**< [in] the number of the new frame to initialize.  \note frame must be >= 1. */);

   ///Open a frame in ds9
   /** Nothing is done if the frame already exists.  First calls \ref addsegment.
     * 
     * This makes use of system() calls.  In each case the result is checked and if SIGINT or SIGQUIT was caught
     * by the ds9 child process, it will be re-raised in the calling process.
     * 
     * \param frame is 
     * 
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int addframe( int frame /**< [in] the number of the new frame to initialize.  \note frame must be >= 1. */);

   ///Display an image in ds9.
   /** A new ds9 instance is opened if necessary, and a new sharedmemory segment is added if necessary.
     * The image is described by a pointer and its 2 or 3 dimensions.
     * 
     * This makes use of system() calls.  In each case the result is checked and if SIGINT or SIGQUIT was caught
     * by the ds9 child process, it will be re-raised in the calling process.
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
   
   ///Shutdown the ds9 interface
   /** Mainly detaches from the shared memory segments
     *  
     * \retval 0 on sucess
     * \retval -1 on an error
     */
   int shutdown();

};

inline
ds9Interface::ds9Interface()
{
   init();
}

template<typename arrayT>
ds9Interface::ds9Interface( const arrayT & array,
                            int frame
                          )
{
   init();
   display(array, frame);
}

inline
ds9Interface::~ds9Interface()
{
   shutdown();
}

inline
int ds9Interface::init()
{
   title("ds9");
   
   _port = 0;

   nsegs = 0;
   segs = 0;

   return 0;
}

inline
std::string ds9Interface::title()
{
   return _title;
}

inline
int ds9Interface::title( const std::string & new_title)
{
   _title = new_title;

   return 0;
}

inline
int ds9Interface::port()
{
   return _port;
}

inline
int ds9Interface::port( const int & new_port)
{
   _port = new_port;

   return 0;
}


inline
int ds9Interface::spawn()
{
   int i;
   char resp[32];
   char cmd[DS9_CMD_MAX_LENGTH];
   int ret;
   
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaaccess %s", _title.c_str());
   
   resp[0] = 0;
   if( command_response(cmd, resp, 128) ) return -1;
   
   //Don't respawn if it already exists.
   if(strcmp(resp, "yes\n") == 0) return 0;
   
   ret = system("xpans &");
   if (WIFSIGNALED(ret))
   {
      if( WTERMSIG(ret) == SIGINT) raise(SIGINT);
      if( WTERMSIG(ret) == SIGQUIT) raise(SIGQUIT);
   }
   
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "ds9 -title %s -port %i &", _title.c_str(), _port);
   
   ret = system(cmd);
   if (WIFSIGNALED(ret))
   {
      if( WTERMSIG(ret) == SIGINT) raise(SIGINT);
      if( WTERMSIG(ret) == SIGQUIT) raise(SIGQUIT);
   }
   
   //Now wait for ds9 to respond or timeout.
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaaccess %s",  _title.c_str());
   
   for(i=0;i<10000;i++)
   {
      resp[0] = 0;
      if( command_response(cmd, resp, 128) ) return -1;
      if(strcmp(resp, "yes\n") == 0) return 0;
      usleep(100);
   }
      
   return -1; //a timeout
}

inline
int ds9Interface::addsegment( int frame )
{
   int i;
   size_t curr_n;
   
   if(frame-1 < nsegs) return 0;
   
   curr_n = nsegs;
   
   nsegs = frame;
   
   segs = (sharedmem_segment *) realloc( segs, sizeof(sharedmem_segment) * nsegs);
      
   for(i = curr_n; i< nsegs; i++)
   {
      sharedmem_segment_initialize( &segs[i] );
      sharedmem_segment_set_key( &segs[i], 0, IPC_PRIVATE);
   }
   
   return 0;
   
}

inline
int ds9Interface::addframe( int frame )
{
   char cmd[DS9_CMD_MAX_LENGTH];
   int ret;
   
   addsegment( frame );

   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s frame %i", _title.c_str(), frame);
   ret = system(cmd);
   if (WIFSIGNALED(ret))
   {
      if( WTERMSIG(ret) == SIGINT) raise(SIGINT);
      if( WTERMSIG(ret) == SIGQUIT) raise(SIGQUIT);
   }
   
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
   size_t i, tot_size;
   char cmd[DS9_CMD_MAX_LENGTH];
   size_t pixsz;
   int bitpix;
   int ret;
   
   if(frame < 1)
   {
      fprintf(stderr, "frame must >= 1\n");
      return -1;
   }
   
   spawn();
   
   bitpix = getFitsBITPIX<dataT>();
   
   pixsz = sizeof(dataT);
      
   //If needed add a shared memory segment for this frame
   if(frame-1 >= nsegs)
   {
      addsegment( frame );
   }

   //Calculate total size
   tot_size= pixsz;
   tot_size*=dim1;
   tot_size*=dim2;
   tot_size*=dim3;

   //Re-allocate shared memory if necessary
   if(tot_size > segs[frame-1].size)
   {
      if( segs[frame-1].size > 0 )
      {
         sharedmem_segment_detach( &segs[frame-1]);
      }
      sharedmem_segment_create( &segs[frame-1], tot_size);
   }
   
   memcpy( segs[frame-1].addr, im, tot_size );
   
   
   snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s frame %i", _title.c_str(), frame);
   ret = system(cmd);
   if (WIFSIGNALED(ret))
   {
      if( WTERMSIG(ret) == SIGINT) raise(SIGINT);
      if( WTERMSIG(ret) == SIGQUIT) raise(SIGQUIT);
   }
   
   
   //Handle single image so that the cube dialog doesn't open up if dim3=1
   if(dim3 == 1)
   {
      snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s shm array shmid %i [xdim=%zu,ydim=%zu,bitpix=%i] &", 
                                      _title.c_str(), segs[frame-1].shmemid,
                                     dim1, dim2, bitpix);
   }
   else
   {
      snprintf(cmd, DS9_CMD_MAX_LENGTH, "xpaset -p %s shm array shmid %i [xdim=%zu,ydim=%zu,zdim=%zu,bitpix=%i] &", 
                                      _title.c_str(), segs[frame-1].shmemid,
                                     dim1, dim2, dim3, bitpix);
   }
   
   ret = system(cmd);
   if (WIFSIGNALED(ret))
   {
      if( WTERMSIG(ret) == SIGINT) raise(SIGINT);
      if( WTERMSIG(ret) == SIGQUIT) raise(SIGQUIT);
   }
   
   return 0;
   
}

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
   
inline
int ds9Interface::shutdown()
{
   size_t i;
   
   for(i=0; i < nsegs; i++) sharedmem_segment_detach( &segs[i] );
   
   free( segs );
   
   nsegs = 0;
   
   return 0;
}


} //namespace improc 
} //namespace mx

#endif //__ds9Interface_hpp__

