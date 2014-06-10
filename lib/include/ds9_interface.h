/** \file ds9_interface.h
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declarations for the mxlib c ds9 interface
  * \ingroup image_processing
  * 
*/

#ifndef __ds9_interface_h__
#define __ds9_interface_h__

#include "IPC.h"


#ifdef __cplusplus
extern "C"
{
#endif

/** \addtogroup image_processing
  * @{
  */

/** \ingroup plotting
  * @{
  */

///The maximum length of the ds9 window title
#define DS9_TITLE_LENGTH 128   
   
///The maximum length of a ds9 command
#define DS9_CMD_MAX_LENGTH 512
   
///Structure to hold the details of an interface to the ds9 image viewer
/** Handles spawning the ds9 window, and manages shared memory segments, 
  * one for each frame.  
  * 
  * \sa Functions for working with ds9_interface include \ref ds9_interface_init,
  * \ref ds9_interface_set_title, \ref ds9_interface_spawn, \ref ds9_interface_addsegment,
  * \ref ds9_interface_addframe, \ref ds9_interface_display_raw, \ref ds9_interface_shutdown
  */   
typedef struct
{
   ///The title of the window, which sets the XPA access point 
   char title[DS9_TITLE_LENGTH];
   
   ///The port to use (normally 0)
   int port;
   
   ///An array of shared memory segments, one per frame
   sharedmem_segment * segs;
   
   ///The number of segments in \ref segs
   size_t nsegs;
   
} ds9_interface;
   
///Initialize the ds9_interface structure
/** The title is set to "ds9", the port is set to 0, segs is set to NULL, and nsegs to 0.
  *
  * \param ds9i is the interface to initialize
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  */
int ds9_interface_init(ds9_interface * ds9i);

///Set the title of an ds9_interface structure
/** The title is used as the name of the XPA access point
  *
  * \param ds9i is the interface which will have its title set
  * \param new_title is the title to set.  Maximum length is \ref DS9_TITLE_LENGTH.
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  */
int ds9_interface_set_title(ds9_interface * ds9i, const char * new_title);

///Spawn (open) the ds9 image viewer
/**
  * \param ds9i is the interface to spawn 
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  */
int ds9_interface_spawn(ds9_interface * ds9i);

///Add a segment corresponding to a particular frame in ds9
/** Nothing is done if the frame already exists.  Note that this does not open a new frame in ds9.
  * 
  * \param ds9i is the interface to add a segment too
  * \param frame is the number of the new frame to initialize.  Note frame must be >= 1.
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  */
int ds9_interface_addsegment(ds9_interface *ds9i, int frame);

///Open a frame in ds9
/** Nothing is done if the frame already exists.  First calls \ref ds9_interface_addsegment.
  * 
  * \param ds9i is the interface to add a frame too
  * \param frame is the number of the new frame to initialize.  Note frame must be >= 1.
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  */
int ds9_interface_addframe(ds9_interface *ds9i, int frame);

///Display an image in ds9
/** A new ds9 instance is opened if necessary, and a new sharedmemory segment is added if necessary.
  * The image is described by its address, its 3 dimensions, and the parameters of its data type.
  * 
  * \param ds9i is the interface to use for display
  * \param frame is the number of the new frame to initialize.  Note frame must be >= 1.
  * \param im is the address of the image
  * \param dim1 is the first dimension of the image (in pixels)
  * \param dim2 is the second dimension of the image (in pixels)
  * \param dim3 is the thirst dimension of the image (in pixels), set to 1 if not a cube.
  * \param bitpix corresponds to the datatype, using the FITS specification.
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  *
  */
int ds9_interface_display_raw(ds9_interface *ds9i, int frame, const void *im, size_t dim1, size_t dim2, size_t dim3, int bitpix);

///Shutdown the ds9 interface
/** Mainly detaches from the shared memory segments
  *  
  * \param ds9i is the interface to shutdown
  * 
  * \retval 0 on sucess
  * \retval -1 on an error
  */
int ds9_interface_shutdown(ds9_interface *ds9i);

///Display an image using a static ds9_interface
/** For convenience to skip the declaration and initialization steps.
  * \note you must call \ref ds9_display_shutdown when you are done with ds9.
  * 
  * \param frame is the number of the new frame to initialize.  Note frame must be >= 1.
  * \param im is the address of the image
  * \param dim1 is the first dimension of the image (in pixels)
  * \param dim2 is the second dimension of the image (in pixels)
  * \param dim3 is the thirst dimension of the image (in pixels), set to 1 if not a cube.
  * \param bitpix corresponds to the datatype, using the FITS specification.
  */
int ds9_display(int frame, const void *im, size_t dim1, size_t dim2, size_t dim3, int bitpix);

///Shutdown the static ds9_interface used by ds9_display.
int ds9_display_shutdown();

///@}
///@}

#ifdef __cplusplus
} //extern "C"
#endif
   
#endif //__ds9_interface_h__

