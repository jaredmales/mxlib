/** \file milkImage.hpp
  * \brief Tools for using the eigen library for image processing
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


template<typename _dataT>
class milkImage
{
public:

    typedef _dataT dataT;

protected:

    std::string m_name; ///< The image name, from name.im.shm (the .im.shm is silent). 

    IMAGE * m_image {nullptr};

    eigenImageMap<dataT> * m_map {nullptr};

public:

    milkImage();

    milkImage( const std::string & imname );

    void open( const std::string & imname );

    void reopen();

    void close();

    eigenImageMap<dataT> & map(); 

    eigenImageMap<dataT> * operator->();

    operator eigenImageMap<dataT>&();

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
    
    ImageStreamIO_openIm(m_image, imname.c_str());
     
    if(ImageStructTypeCode<dataT>::TypeCode != m_image->md->datatype)
    {
        delete m_image;
        m_image = nullptr;
        throw std::invalid_argument("shmim datatype does not match template type");
    }

    m_name = imname;

    if(m_map)
    {
        delete m_map;
        m_map = nullptr;
    }

    m_map = new eigenImageMap<dataT>((dataT*)m_image->array.raw, m_image->md->size[0], m_image->md->size[1]);
}

template<typename dataT>
eigenImageMap<dataT> & milkImage<dataT>::map() 
{
    return *m_map;
}

template<typename dataT>
eigenImageMap<dataT> * milkImage<dataT>::operator->() 
{
    return m_map;
}

template<typename dataT>
milkImage<dataT>::operator eigenImageMap<dataT>&()
{
    return *m_map;
}

template<typename dataT>
void milkImage<dataT>::post()
{
    if(!m_image) return;
    ImageStreamIO_UpdateIm(m_image);
}

}
}

#endif //MXLIB_MILK
#endif //improc_milkImage_hpp

