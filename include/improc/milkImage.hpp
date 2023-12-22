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

    eigenMap<dataT> * m_map {nullptr}; //An Eigen::Map of the array
    
    uint64_t m_size_0 {0}; ///< The size[0] of the image when last opened.

    uint64_t m_size_1 {0}; ///< The size[1] of the image when last opened.

public:

    /// Default c'tor
    milkImage();

    /// Constructor which opens the specified image
    milkImage( const std::string & imname /**< [in] The image name, from name.im.shm (the .im.shm should not be given).*/ );

    /// D'tor
    ~milkImage();

    /// Open and connect to an image, allocating the eigenMap.
    /**
      * \throws std::invalid_argument if the image type_code does not match dataT.
      */  
    void open( const std::string & imname /**< [in] The image name, from name.im.shm (the .im.shm should not be given).*/);

    /// Create and connect to an image, allocating the eigenMap.
    /**
      * \throws std::invalid_argument if the image type_code does not match dataT.
      */  
    void create( const std::string & imname, ///< [in] The image name, for name.im.shm (the .im.shm should not be given).
                 uint32_t sz0,               ///< [in] the x size of the image
                 uint32_t sz1                ///< [in] the y size of the image
               );

    /// Get the width of the image
    /**
      * \returns the current value of m_size_0 
      */
    uint32_t rows();

    /// Get the height of the image
    /**
      * \returns the current value of m_size_1
      */
    uint32_t cols();

    /// Get the size of a dimension of the image
    /**
      * \returns the current value of m_size_0 or m_size_1 depending on n
      */
    uint32_t size( unsigned n /**< [in] the dimension to get the size of*/);

    /// Checks if the image is connected and is still the same format as when connected.
    /** Checks on pointer value, size[], and data_type.
      *
      * \returns true if connected and no changes
      * \returns false if not connected or something changed.  All maps are now invalid. 
      */
    bool valid();

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

    /// Get an eigenMap 
    /** Use with caution:
      * - there is no way to know if the image has changed
      * - you should check \ref valid() before using
      * 
      * 
      * \returns an Eigen::Map<Array,-1,-1> reference
      * 
      * \throws mx::err::mxException if the image is not opened
      * 
      */
    eigenMap<dataT> & operator()();

    /// Conversion operator returns an eigenMap
    /** Use this like
      * \code
      * milkImage<float> mim("imname");
      * eigenMap<float> im(mim); //you now have an Eigen interface to the image
      * \endcode
      * but with caution:
      * - there is no way to know if the image has changed
      * - you should check \ref valid() before using
      * 
      * 
      * \returns an Eigen::Map<Array,-1,-1> object
      * 
      * \throws mx::err::mxException if the image is not opened
      * 
      */
    operator eigenMap<dataT>();

    

    /// Copy data from an Eigen Array type to the shared memory stream
    /** Sets the write flag, copies using the Eigen assigment to map, unsets the write flag, then posts.
      * 
      * \throws mxException on an error
      * 
      */
    template<typename eigenT> 
    milkImage & operator=(const eigenT & im /**< [in] the eigen array to copy to the stream*/);

    /// Set the write flag
    /** The write flag is set to indicate whether or not the the data is being changed.
      * The write flag will be set to false by \ref post().
      *  
      * \throws mx::err::mxException if the image is not opened
      */
    void setWrite( bool wrflag = true /**< [in] [optional] the desired value of the write flag.  Default is true.*/);

    /// Update the metadata and post all semaphores
    /** 
      * \todo need to set wtime, have a version with atime
      * 
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
milkImage<dataT>::~milkImage()
{
    close();
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
    
    if(m_map)
    {
        delete m_map;
        m_map = nullptr;
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

    m_map = new eigenMap<dataT>(m_raw, m_size_0, m_size_1);
}

template<typename dataT>
void milkImage<dataT>::create( const std::string & imname,
                               uint32_t sz0, 
                               uint32_t sz1  
                             )
{
    if(m_image)
    {
        ImageStreamIO_closeIm(m_image);
        delete m_image;
        m_image = nullptr;
    }
    
    if(m_map)
    {
        delete m_map;
        m_map = nullptr;
    }

    m_image = new IMAGE;
    
    uint32_t imsize[3];
    imsize[0] = sz0; 
    imsize[1] = sz1;
    imsize[2] = 1;

    errno_t rv =  ImageStreamIO_createIm_gpu(m_image, imname.c_str(), 3, imsize, ImageStructTypeCode<dataT>::TypeCode, 
                                                   -1, 1, IMAGE_NB_SEMAPHORE, 0, CIRCULAR_BUFFER | ZAXIS_TEMPORAL, 0);

    if(rv != IMAGESTREAMIO_SUCCESS) 
    {
        delete m_image;
        m_image = nullptr;
        throw err::mxException("", 0, "", 0, "", 0, "ImageStreamIO_createIm_gpu returned an error");     
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

    m_map = new eigenMap<dataT>(m_raw, m_size_0, m_size_1);          

}

template<typename dataT>
eigenMap<dataT> & milkImage<dataT>::operator()()
{
    if(m_map == nullptr)
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open (map is null)");
    }

    return *m_map;
}

template<typename dataT>
milkImage<dataT>::operator eigenMap<dataT>()
{
    return operator()();
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
    if(m_map != nullptr)
    {
        delete m_map;
        m_map = nullptr;
    }

    if(m_image == nullptr) 
    {
        return;
    }

    errno_t rv = ImageStreamIO_closeIm(m_image);
    delete m_image;
    m_image = nullptr;

    if(rv != IMAGESTREAMIO_SUCCESS) throw err::mxException("", 0, "", 0, "", 0, "ImageStreamIO_closeIm returned an error");
}

template<typename dataT>
uint32_t milkImage<dataT>::rows()
{
    return m_size_0;
}

template<typename dataT>
uint32_t milkImage<dataT>::cols()
{
    return m_size_1;
}

template<typename dataT>
uint32_t milkImage<dataT>::size( unsigned n )
{
    if(n == 0)
    {
        return m_size_0;
    }
    else if(n == 1)
    {
        return m_size_1;
    }

    return 0;
}

template<typename dataT>
bool milkImage<dataT>::valid()
{
    if(m_image == nullptr)
    {
        return false;
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
template<typename eigenT> 
milkImage<dataT> & milkImage<dataT>::operator=(const eigenT & im)
{
    eigenMap<dataT> map((dataT*)m_image->array.raw, m_image->md->size[0], m_image->md->size[1]);
        
    setWrite(true);

    map = im;
        
    setWrite(false);
        
    post();

    return *this;
}

template<typename dataT>
void milkImage<dataT>::setWrite(bool wrflag)
{
    if(m_image == nullptr)
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open");
    }

    m_image->md->write = wrflag;
}

template<typename dataT>
void milkImage<dataT>::post()
{
    if(m_image == nullptr)
    {
        throw err::mxException("", 0, "", 0, "", 0, "Image is not open");
    }

    errno_t rv = ImageStreamIO_UpdateIm(m_image);
    m_image->md->atime=m_image->md->writetime;

    if(rv != IMAGESTREAMIO_SUCCESS) 
    {
        throw err::mxException("", 0, "", 0, "", 0, "ImageStreamIO_UpdateIm returned an error");
    }
}

}
}

#endif //MXLIB_MILK
#endif //improc_milkImage_hpp

