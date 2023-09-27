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

template<typename _dataT>
class milkImage
{
public:

    typedef _dataT dataT;

    IMAGE * m_image {nullptr};

    eigenImageMap<dataT> * m_map {nullptr};

    void open( const std::string & imname )
    {
        if(m_image)
        {
            ImageStreamIO_closeIm(m_image);
            delete m_image;
            m_image = nullptr;
        }

        m_image = new IMAGE;

        ImageStreamIO_openIm(m_image, imname.c_str());

        if(m_map)
        {
            delete m_map;
            m_map = nullptr;
        }

        m_map = new eigenImageMap<dataT>((dataT*)m_image->array.raw, m_image->md->size[0], m_image->md->size[1]);
    }

    eigenImageMap<dataT> & image() {return *m_map;}
    
    eigenImageMap<dataT> * operator->() {return m_map;}

    operator eigenImageMap<dataT>&() {return *m_map;}



    void post()
    {
        if(!m_image) return;

        ImageStreamIO_UpdateIm(m_image);
    }
};


}
}

#endif //MXLIB_MILK
#endif //improc_milkImage_hpp

