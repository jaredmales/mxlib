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
#include "../ioutils/fits/fitsUtils.hpp"

#ifndef DS9INTERFACE_NO_EIGEN
#include "eigenImage.hpp"
#endif

namespace mx
{
namespace improc
{

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

#ifndef DS9INTERFACE_NO_EIGEN
template<typename arrayT>
ds9Interface::ds9Interface( const arrayT & array,
                            int frame
                          )
{
   initialize();
   display(array, frame);
}

extern template ds9Interface::ds9Interface<eigenImage<signed char>>( const eigenImage<signed char> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<unsigned char>>( const eigenImage<unsigned char> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<short>>( const eigenImage<short> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<unsigned short>>( const eigenImage<unsigned short> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<int>>( const eigenImage<int> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<unsigned int>>( const eigenImage<unsigned int> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<long>>( const eigenImage<long> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<unsigned long>>( const eigenImage<unsigned long> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<float>>( const eigenImage<float> & array, int frame);
extern template ds9Interface::ds9Interface<eigenImage<double>>( const eigenImage<double> & array, int frame);

#endif //DS9INTERFACE_NO_EIGEN


template<typename dataT>
int ds9Interface::display( const dataT * im,
                           size_t dim1,
                           size_t dim2,
                           size_t dim3,
                           int frame
                          )
{
   return display(im, fits::getFitsBITPIX<dataT>(), sizeof(dataT), dim1, dim2, dim3, frame);
}

extern template int ds9Interface::display<signed char>(const signed char * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<unsigned char>(const unsigned char * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<short>(const short * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<unsigned short>(const unsigned short * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<int>(const int * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<unsigned int>(const unsigned int * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<long>(const long * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<unsigned long>(const unsigned long * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<float>(const float * im, size_t dim1, size_t dim2, size_t dim3, int frame);
extern template int ds9Interface::display<double>(const double * im, size_t dim1, size_t dim2, size_t dim3, int frame);

#ifndef DS9INTERFACE_NO_EIGEN
template<typename arrayT>
int ds9Interface::display( const arrayT & array,
                           int frame
                         )
{
   eigenArrPlanes<arrayT> planes;

   return display(array.data(), array.rows(), array.cols(), planes(array), frame);
}

extern template int ds9Interface::display<eigenImage<signed char>>( const eigenImage<signed char> & im, int frame);
extern template int ds9Interface::display<eigenImage<unsigned char>>( const eigenImage<unsigned char> & im, int frame);
extern template int ds9Interface::display<eigenImage<short>>( const eigenImage<short> & im, int frame);
extern template int ds9Interface::display<eigenImage<unsigned short>>( const eigenImage<unsigned short> & im, int frame);
extern template int ds9Interface::display<eigenImage<int>>( const eigenImage<int> & im, int frame);
extern template int ds9Interface::display<eigenImage<unsigned int>>( const eigenImage<unsigned int> & im, int frame);
extern template int ds9Interface::display<eigenImage<long>>( const eigenImage<long> & im, int frame);
extern template int ds9Interface::display<eigenImage<unsigned long>>( const eigenImage<unsigned long> & im, int frame);
extern template int ds9Interface::display<eigenImage<float>>( const eigenImage<float> & im, int frame);
extern template int ds9Interface::display<eigenImage<double>>( const eigenImage<double> & im, int frame);

template<typename arrayT>
int ds9Interface::operator()( const arrayT & array,
                              int frame
                            )
{
   return display(array, frame);
}

extern template int ds9Interface::operator()<eigenImage<signed char>>( const eigenImage<signed char> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<unsigned char>>( const eigenImage<unsigned char> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<short>>( const eigenImage<short> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<unsigned short>>( const eigenImage<unsigned short> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<int>>( const eigenImage<int> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<unsigned int>>( const eigenImage<unsigned int> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<long>>( const eigenImage<long> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<unsigned long>>( const eigenImage<unsigned long> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<float>>( const eigenImage<float> & im, int frame);
extern template int ds9Interface::operator()<eigenImage<double>>( const eigenImage<double> & im, int frame);

#endif //DS9INTERFACE_NO_EIGEN


} //namespace improc
} //namespace mx


#endif //improc_ds9Interface_hpp
