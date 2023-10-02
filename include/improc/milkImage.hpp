/** \file milkImage.hpp
  * \brief Interface to MILK::ImageStreamIO shared memory streams
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
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

#ifndef improc_milkImage_hpp
#define improc_milkImage_hpp

#ifdef MXLIB_MILK

#ifdef MXLIB_CUDA
#define HAVE_CUDA
#endif

#include <ImageStreamIO/ImageStreamIO.h>

#include "../mxException.hpp"
#include "eigenImage.hpp"

namespace mx
{
namespace improc 
{

template<typename typeT>
struct ImageStructTypeCode;

template<>
struct ImageStructTypeCode<uint8_t>
{
    constexpr static int TypeCode = _DATATYPE_UINT8;
};

template<>
struct ImageStructTypeCode<int8_t>
{
    constexpr static int TypeCode = _DATATYPE_INT8;
};

template<>
struct ImageStructTypeCode<char>
{
    constexpr static int TypeCode = _DATATYPE_INT8;
};

template<>
struct ImageStructTypeCode<uint16_t>
{
    constexpr static int TypeCode = _DATATYPE_UINT16;
};

template<>
struct ImageStructTypeCode<int16_t>
{
    constexpr static int TypeCode = _DATATYPE_INT16;
};

template<>
struct ImageStructTypeCode<uint32_t>
{
    constexpr static int TypeCode = _DATATYPE_UINT32;
};

template<>
struct ImageStructTypeCode<int32_t>
{
    constexpr static int TypeCode = _DATATYPE_INT32;
};

template<>
struct ImageStructTypeCode<uint64_t>
{
    constexpr static int TypeCode = _DATATYPE_UINT64;
};
template<>
struct ImageStructTypeCode<int64_t>
{
    constexpr static int TypeCode = _DATATYPE_INT64;
};
template<>
struct ImageStructTypeCode<float>
{
    constexpr static int TypeCode = _DATATYPE_FLOAT;
};

template<>
struct ImageStructTypeCode<double>
{
    constexpr static int TypeCode = _DATATYPE_DOUBLE;
};

template<>
struct ImageStructTypeCode<std::complex<float>>
{
    constexpr static int TypeCode = _DATATYPE_COMPLEX_FLOAT;
};

template<>
struct ImageStructTypeCode<std::complex<double>>
{
    constexpr static int TypeCode = _DATATYPE_COMPLEX_DOUBLE;
};

/// Class to interface with an ImageStreamIO image in shared memory
/**
  *
  * Use with Eigen::Map (aliased as mx::improc::eigenMap)
  * \code
  * using namespace mx::improc;
  * milkImage<float> mim("image"); //connects to image.im.shm
  * eigenMap<float> im(mim); //the conversion operator passes a reference to the internal map
  * im.setRandom(); // im is now a map pointed at the image data.  Eigen functions are now available.
  * im(64,64) *= 2000;
  * im /= 0.2;
  * \endcode
  * Once you have changed something via the Eigen::Map you want to notify others connected to the stream 
  * via
  * \code
  * mim.post();
  * \endcode
  * 
  * \ingroup eigen_image_processing
  */  
template<typename _dataT>
class milkImage
{
public:

    typedef _dataT dataT; ///< The data type 

protected:

    std::string m_name; ///< The image name, from name.im.shm (the .im.shm should not be given). 

    IMAGE * m_image {nullptr}; ///< Pointer to the ImageStreamIO IMAGE structure.

    dataT * m_raw {nullptr}; //The raw pointer address is stored here for checks

    uint64_t m_size_0 {0}; ///< The size[0] of the image when last opened.

    uint64_t m_size_1 {0}; ///< The size[1] of the image when last opened.

public:

    /// Default c'tor
    milkImage();

    /// Constructor which opens the specified image
    milkImage( const std::string & imname /**< [in] The image name, from name.im.shm (the .im.shm should not be given).*/ );

    /// Open and connect to an image, allocating the eigenMap.
    /**
      * \throws std::invalid_argument if the image type_code does not match dataT.
      */  
    void open( const std::string & imname /**< [in] The image name, from name.im.shm (the .im.shm should not be given).*/);

    /// Checks if the image is still the same as when connected.
    /** Checks on pointer value, size[], and data_type.
      *
      * \returns true if no changes
      * \returns false if something changed.  all maps are now invalid. 
      */
    bool unchanged();

    /// Reopens the image. 
    /** Same as
      * \code
      * close();
      * open(m_name);
      * \endcode
      */ 
    void reopen();

    // Close the image
    void close();

    /// Conversion operator returns an eigenMap
    /** Use this like
      * \code
      * milkImage<float> mim("imname");
      * eigenMap<float> im(mim); //you now have an Eigen interface to the image
      * \endcode
      * but with caution:
      * - there is no way to know if the image has changed
      * 
      * Note: we assume this does a move
      * 
      * \returns an Eigen::Map<Array,-1,-1> object
      * 
      * \throws mx::err::mxException if the image is not opened
      * 
      */
    operator eigenMap<dataT>();

    /// Update the metadata and post all semaphores
    /**
      * \throws mx::err::mxException if the image is not opened
      */ 
    void post();

};

template<typename dataT>
milkImage<dataT>::milkImage()
{
}

template<typename dataT>
milkImage<dataT>::milkImage( const std::string & imname )
{
    open(imname);
}

template<typename dataT>
void milkImage<dataT>::open( const std::string & imname )
{
    if(m_image)
    {
        ImageStreamIO_closeIm(m_image);
        delete m_image;
        m_image = nullptr;
    }
    
    m_image = new IMAGE;
    
    errno_t rv = ImageStreamIO_openIm(m_image, imname.c_str());

    if(rv != IMAGESTREAMIO_SUCCESS) 
    {
        delete m_image;
        m_image = nullptr;
        throw err::mxException("", 0, "", 0, "", 0, "ImageStreamIO_openIm returned an error");     
    }

    if(ImageStructTypeCode<dataT>::TypeCode != m_image->md->datatype)
    {
        delete m_image;
        m_image = nullptr;
        throw std::invalid_argument("shmim datatype does not match template type");
    }

    m_name = imname;
    m_raw = (dataT*) m_image->array.raw;
    m_size_0 = m_image->md->size[0];
    m_size_1 = m_image->md->size[1];
}

template<typename dataT>
milkImage<dataT>::operator eigenMap<dataT>()
{
    if(m_image == nullptr)
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open");
    }

    if(!unchanged())
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image has changed");
    }

    return eigenMap<dataT>((dataT*)m_image->array.raw, m_image->md->size[0], m_image->md->size[1]);
}

template<typename dataT>
void milkImage<dataT>::reopen()
{
    close();
    open(m_name);
}

template<typename dataT>
void milkImage<dataT>::close()
{
    if(m_image == nullptr) 
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open");
    }

    errno_t rv = ImageStreamIO_closeIm(m_image);
    delete m_image;
    m_image = nullptr;

    if(rv != IMAGESTREAMIO_SUCCESS) throw err::mxException("", 0, "", 0, "", 0, "ImageStreamIO_closeIm returned an error");
}

template<typename dataT>
bool milkImage<dataT>::unchanged()
{
    if(m_image == nullptr)
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open");
    }

    if( m_raw != (dataT*) m_image->array.raw || 
            m_size_0 != m_image->md->size[0] ||
                m_size_1 != m_image->md->size[1] ||
                    ImageStructTypeCode<dataT>::TypeCode != m_image->md->datatype )
    {
        return false;
    }

    return true;
}

template<typename dataT>
void milkImage<dataT>::post()
{
    if(m_image == nullptr)
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open");
    }

    errno_t rv = ImageStreamIO_UpdateIm(m_image);

    if(rv != IMAGESTREAMIO_SUCCESS) throw err::mxException("", 0, "", 0, "", 0, "ImageStreamIO_closeIm returned an error");
}

}
}

#endif //MXLIB_MILK
#endif //improc_milkImage_hpp

