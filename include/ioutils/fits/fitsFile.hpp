/** \file fitsFile.hpp
  * \brief Declares and defines a class to work with a FITS file
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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

#ifndef ioutils_fits_fitsFile_hpp
#define ioutils_fits_fitsFile_hpp


#include "../../mxError.hpp"

#include "../../improc/eigenImage.hpp"

#include "fitsUtils.hpp"
#include "fitsHeader.hpp"

namespace mx
{
namespace fits 
{
   


/// Class to manage interactions with a FITS file
/** This class wraps the functionality of cfitsio.
  *
  * \tparam dataT is the datatype to use for in-memory storage of the image.  This does not have to match the data type stored on disk for reading, but will be the type used for writing.
  *
  * \ingroup fits_processing
  */
template<typename dataT> class fitsFile
{

protected:

   ///The path to the file
   std::string m_fileName;

   ///The cfitsio data structure
   fitsfile * m_fptr {nullptr};

   ///The dimensions of the image (1D, 2D, 3D etc)
   int m_naxis;

   ///The size of each dimension
   long * m_naxes {nullptr};

   ///Flag indicating whether the file is open or not
   bool m_isOpen;

   ///The value to replace null values with
   dataT m_nulval;

   ///Records whether any null values were replaced
   int m_anynul;

   ///Flag to control whether the comment string is read.
   int m_noComment;


   ///The starting x-pixel to read from
   long m_x0;

   ///The starting y-pixel to read from
   long m_y0;

   ///The number of x-pixels to read
   long m_xpix;

   ///The number of y-pixels to read
   long m_ypix;

   ///One time initialization common to all constructors
   void construct();

public:

   ///Default constructor
   fitsFile();

   ///Constructor with m_fileName, and option to open.
   fitsFile( const std::string & fname, ///<File name to set on construction
             bool doopen = true ///< If true, then the file is opened (the default).
           );

   ///Destructor
   ~fitsFile();

   ///Get the current value of m_fileName
   /**
     * \returns the current value of m_fileName.
     */
   std::string fileName();

   ///Set the file path, and optionally open the file.
   /**
     * \returns 0 on success
     * \returns -1 on success
     */
   int fileName( const std::string & fname,  ///< The new file name.
                 bool doopen = true ///< If true, then the file is opened (the default).
               );

   ///Get the current value of m_naxis
   /**
     * \returns the current value of m_naxis
     */
   int naxis();

   ///Get the current value of m_naxes for the specified dimension
   /**
     * \returns the current value of m_naxes for the specified dimension
     */
   long naxes( int dim /**< [in] the dimension */);



   ///Open the file
   /** File name needs to already have been set.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int open();

   ///Open the file, first setting the file path.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */
   int open( const std::string & fname /**< The name of the file to open. */);

   ///Close the file.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */
   int close();

   ///Get the number of dimensions (i.e. m_naxis)
   int getDimensions();

   ///Get the total size
   long getSize();

   ///Get the size of a specific dimension
   long getSize(size_t axis);

   /** \name Reading Basic Arrays
     * These methods read FITS data into basic or raw arrays specified by a pointer.
     * @{
     */

   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int read( dataT * data /**< [out] an allocated arrray large enough to hold the entire image */ );

   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int read( dataT * data, ///< [out] an allocated arrray large enough to hold the entire image
             fitsHeader &head ///< [out] a fitsHeader object which is passed to \ref readHeader
           );

   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int read( dataT * data, ///< [out] is an allocated arrray large enough to hold the entire image
             const std::string &fname ///< [in] is the file path, which is passed to \ref fileName
           );

   ///Read the contents of the FITS file into an array and read the header.
   /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int read( dataT * data,             ///< [out] an allocated arrray large enough to hold the entire image
             fitsHeader &head,         ///< [out] a fitsHeader object which is passed to \ref readHeader
             const std::string &fname ///< [in] the file path, which is passed to \ref fileName
           );

   ///Read data from a vector list of files into an image cube
   int read( dataT * im, ///< [out] An allocated array large enough to hold all the images
             const std::vector<std::string> & flist ///< [in] The list of files to read.
           );

   ///Read data from a vector of files into an image cube with individual headers
   int read( dataT * im, ///< [out] An allocated array large enough to hold all the images
             std::vector<fitsHeader> &heads, ///< [in/out] The vector of fits headers, allocated to contain one per image.
             const std::vector<std::string> & flist  ///< [in] The list of files to read.
           );


   ///@}

   /** \name Reading Eigen Arrays
     * These methods read FITS data into array types with an Eigen-like interface.
     * @{
     */

   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     *
     * \tparam arrT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename arrT>
   int read( arrT & data /**< [out] is the array, which will be resized as necessary using its resize(int, int) member */);

   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     *
     * \tparam arrT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename arrT>
   int read( arrT & data,      ///< [out] is the array, which will be resized as necessary using its resize(int, int) member
             fitsHeader &head  ///< [out] is a fitsHeader object which is passed to \ref readHeader
           );

   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     *
     * \tparam arrT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename arrT>
   int read( arrT & data,              ///< [out] is the array, which will be resized as necessary using its resize(int, int) member
             const std::string &fname  ///< [in] is the file path, which is passed to \ref fileName
           );

   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     *
     * \tparam arrT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename arrT>
   int read( arrT & data,  ///< [out] the array, which will be resized as necessary using its resize(int, int) member
             fitsHeader &head, ///< [out] a fitsHeader object which is passed to \ref readHeader
             const std::string &fname ///< [in] the file path, which is passed to \ref fileName
           );

   ///Read data from a vector list of files into an image cube
   /** The type cubeT can be any type with the following members defined:
     * - resize(int, int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef cubeT::Scalar (is the data type, does not have to be dataT)
     *
     * \note The images must all be the same size, or all be as large or larger than the subset specified.
     *
     * \tparam cubeT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename cubeT>
   int read( cubeT & cube,    ///< [out] A cube which will be resized using its resize(int, int, int) member.
             const std::vector<std::string> & flist, ///< [in] The list of files to read.
             std::vector<fitsHeader> * heads = 0 ///< [out] [optional] A vector of fits headers, allocated to contain one per image.
           );

   ///Read data from a vector of files into an image cube with individual headers
   /** The type cubeT can be any type with the following members defined:
     * - resize(int, int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef cubeT::Scalar (is the data type, does not have to be dataT)
     *
     * \tparam cubeT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename cubeT>
   int read( cubeT & cube, ///< [out] A cube which will be resized using its resize(int, int, int) member.
             std::vector<fitsHeader> &heads, ///< [out] The vector of fits headers, allocated to contain one per image.
             const std::vector<std::string> & flist  ///< [in] The list of files to read.
           );


   ///@}


   /** \name Reading Headers
     * @{
     */

   ///Read the header from the fits file.
   /** If head is not empty, then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int readHeader(fitsHeader &head /**< [out] a fitsHeader object */);

   ///Read the header from the fits file.
   /** If head is not empty, then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int readHeader( fitsHeader &head,  ///< [out] a fitsHeader object
                   const std::string &fname ///< [in] the file path, which is passed to \ref fileName
                 );

   /// Read the headers from a list of FITS files.
   /** In each case, if the header is not empty then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int readHeader( std::vector<fitsHeader> & heads, /// A vector of fitsHeader objects to read into.
                   const std::vector<std::string> & flist ///< [in] A list of files, each of which is passed to \ref fileName
                 );

   ///@}

   /** \name Writing Basic Arrays
     * These methods write basic arrays specified by a pointer to FITS.
     * @{
     */

   ///Write the contents of a raw array to the FITS file.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */
   int write( const dataT * im, ///< [in] is the array
              int d1,  ///< [in] is the first dimension
              int d2,  ///< [in] is the second dimension
              int d3,  ///< [in] is the third dimenesion (minimum value is 1)
              fitsHeader * head ///< [in] a pointer to the header.  Set to 0 if not used.
            );

   ///Write the contents of a raw array to the FITS file.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */
   int write( const dataT * im, ///< [in] is the array
              int d1,  ///< [in] is the first dimension
              int d2,  ///< [in] is the second dimension
              int d3  ///< [in] is the third dimenesion (minimum value is 1)
            );

   ///Write the contents of a raw array to the FITS file.
   /**
     * Note: the type of the array must match dataT
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   int write( const dataT * im, ///< [in] is the array
              int d1,  ///< [in] is the first dimension
              int d2,  ///< [in] is the second dimension
              int d3,  ///< [in] is the third dimenesion (minimum value is 1)
              fitsHeader & head ///< [in] is the header
            );

   ///Write the contents of a raw array to the FITS file.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */
   int write( const std::string & fname, ///< [in] is the name of the file.
              const dataT * im, ///< [in] is the array
              int d1,  ///< [in] is the first dimension
              int d2,  ///< [in] is the second dimension
              int d3   ///< [in] is the third dimenesion (minimum value is 1)
            );

   ///Write the contents of a raw array to the FITS file.
   /**
     * \returns 0 on success
     * \returns -1 on error
     */
   int write( const std::string  & fname, ///< [in] is the name of the file.
              const dataT * im, ///< [in] is the array
              int d1,  ///< [in] is the first dimension
              int d2,  ///< [in] is the second dimension
              int d3,  ///< [in] is the third dimenesion (minimum value is 1)
              fitsHeader & head ///< [in] is the header
            );

   ///@}

   /** \name Writing Eigen Arrays
     * These methods write array types with an Eigen-like interface.
     * @{
     */

   ///Write the contents of an Eigen-type array to a FITS file.
   /** The type arrT can be any type with the following members defined:
     * - data() (returns a pointer to the underlying array)
     * - rows() (returrns the number of rows)
     * - cols() (returns the number of columns)
     * - may have planes() defined
     *
     * Note: as with all write methods, the Scalar type of the array must match dataT
     *
     * \tparam arrT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename arrT>
   int write( const std::string & fname, ///< [in] is the name of the file.
              const arrT & im ///< [in] is the array
            );

   ///Write the contents of an Eigen-type array to a FITS file.
   /** The type arrT can be any type with the following members defined:
     * - data() (returns a pointer to the underlying array)
     * - rows() (returrns the number of rows)
     * - cols() (returns the number of columns)
     * - may have planes() defined.
     *
     * Note: as with all write methods, the Scalar type of the array must match dataT
     *
     * \tparam arrT is the type of array, see requirements above.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
   template<typename arrT>
   int write( const std::string & fname,  ///< [in] is the file path, which is passed to \ref fileName
              const arrT & im, ///< [in] is the array
              fitsHeader & head ///< [in] is a fitsHeader object which is passed to \ref readHeader
            );

   //int writeHeader( fitsHeader &head );

   ///@}

   /** \name Reading Subsets
     * It is often desirable to read only a subset of an image or images into memory.  These methods
     * allow you to specify this.
     * @{
     */

   ///Set to read all the pixels in the file
   void setReadSize();

   ///Set to read only a subset of the pixels in the file
   /**
     */
   void setReadSize( long x0,   ///< is the starting x-pixel to read
                     long y0,   ///< is the starting y-pixel to read
                     long xpix, ///< is the number of x-pixels to read
                     long ypix  ///< is the number of y-pixels to read
                   );

   ///@}

}; // fitsFile


template<typename dataT>
void fitsFile<dataT>::construct()
{   
   m_isOpen = 0;

   m_nulval = 0;
   m_anynul = 0;

   m_noComment = 0;

   setReadSize();

}

template<typename dataT>
fitsFile<dataT>::fitsFile()
{
   construct();
}

template<typename dataT>
fitsFile<dataT>::fitsFile(const std::string &fname, bool doopen)
{
   construct();
   fileName(fname, doopen);
}

template<typename dataT>
fitsFile<dataT>::~fitsFile()
{
   if(m_isOpen) close();

   if(m_naxes) delete[] m_naxes;
}

template<typename dataT>
std::string fitsFile<dataT>::fileName()
{
   return m_fileName;
}

template<typename dataT>
int fitsFile<dataT>::fileName(const std::string &fname, bool doopen)
{
   if(m_isOpen)
   {
      close();
   }

   m_fileName = fname;

   if(doopen)
   {
      return open();
   }

   return 0;

}

template<typename dataT>
int fitsFile<dataT>::naxis()
{
   return m_naxis;
}

template<typename dataT>
long fitsFile<dataT>::naxes( int dim)
{
   if(m_naxes == nullptr) return -1;
   
   return m_naxes[dim];
}

template<typename dataT>
int fitsFile<dataT>::open()
{
   int fstatus = 0;

   fits_open_file(&m_fptr, m_fileName.c_str(), READONLY, &fstatus);

   if (fstatus)
   {
      std::string explan = "Error opening file";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILEOERR, explan);

      return -1;

   }

   fits_get_img_dim(m_fptr, &m_naxis, &fstatus);
   if (fstatus)
   {
      std::string explan = "Error getting number of axes in file";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILERERR, explan);

      return -1;

   }

   if(m_naxes) delete[] m_naxes;
   m_naxes = new long[m_naxis];

   fits_get_img_size(m_fptr, m_naxis, m_naxes, &fstatus);
   if (fstatus)
   {
      std::string explan = "Error getting dimensions in file";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILERERR, explan);

      return -1;
   }

   m_isOpen = true; //Only set this after opening is complete.

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::open(const std::string & fname)
{
   return fileName(fname, true);
}

template<typename dataT>
int fitsFile<dataT>::close()
{
   //int stat;
   int fstatus = 0;

   if(!m_isOpen)
   {
      return 0; //No error.
   }

   //stat =
   fits_close_file(m_fptr, &fstatus);

   if (fstatus)
   {
      std::string explan = "Error closing file";
      fitsErrText( explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILECERR, explan);

      return -1;
   }

   m_isOpen = 0;
   fstatus = 0;

   if(m_naxes) delete[] m_naxes;

   m_naxes = nullptr;

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::getDimensions()
{
   if(!m_isOpen or !m_naxes) return -1;

   return m_naxis;
}

template<typename dataT>
long fitsFile<dataT>::getSize()
{
   if(!m_isOpen or !m_naxes) return -1;

   long sz = 1;

   if(m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis ==2)
   {
      return m_xpix*m_ypix;
   }
   else
   {
      for(int i=0;i<m_naxis;++i) sz*= m_naxes[i];
   }

   return sz;
}

template<typename dataT>
long fitsFile<dataT>::getSize(size_t axis)
{
   if(!m_isOpen or !m_naxes) return -1;

   if(m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis ==2)
   {
      if(axis == 0) return m_xpix;
      return m_ypix;
   }
   else
   {
      return m_naxes[axis];
   }
}

/************************************************************/
/***                      Basic Arrays                    ***/
/************************************************************/

template<typename dataT>
int fitsFile<dataT>::read(dataT * data)
{
   int fstatus = 0;

   if(!m_isOpen)
   {
      if( open() < 0) return -1;
   }

//   long long nelements = 1;
   long *fpix = new long[m_naxis];
   long *lpix = new long[m_naxis];
   long *inc = new long[m_naxis];

   if(m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis == 2)
   {
      fpix[0] = m_x0+1;
      lpix[0] = fpix[0] + m_xpix-1;
      fpix[1] = m_y0+1;
      lpix[1] = fpix[1] + m_ypix-1;

      inc[0] = 1;
      inc[1] = 1;
   }
   else
   {
      for(int i=0;i<m_naxis; i++)
      {
         fpix[i] = 1;
         lpix[i] = m_naxes[i];
         inc[i] = 1;
         //nelements *= m_naxes[i];
      }
   }

   //fits_read_pix(m_fptr, fitsType<dataT>(), fpix, nelements, (void *) &m_nulval,
                                     //(void *) data, &m_anynul, &fstatus);
   fits_read_subset(m_fptr, fitsType<dataT>(), fpix, lpix, inc, (void *) &m_nulval,
                                     (void *) data, &m_anynul, &fstatus);

   if (fstatus && fstatus != 107)
   {
      std::string explan = "Error reading data from file";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILERERR, explan);

      delete[] fpix;
      delete[] lpix;
      delete[] inc;

      return -1;
   }

   delete[] fpix;
   delete[] lpix;
   delete[] inc;

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::read( dataT * data,
                           fitsHeader &head
                         )
{
   if( read(data) < 0 ) return -1;
   if( readHeader(head) < 0) return -1;

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::read( dataT * data,
                           const std::string &fname
                         )
{
   if( fileName(fname) < 0 ) return -1;
   if( read(data) < 0 ) return -1;
}

template<typename dataT>
int fitsFile<dataT>::read( dataT * data,
                           fitsHeader &head,
                           const std::string &fname
                         )
{
   if( fileName(fname) < 0 ) return -1;
   if( read(data) < 0 ) return -1;
   if( readHeader(head) < 0 ) return -1;
}

template<typename dataT>
int fitsFile<dataT>::read( dataT * im,
                           const std::vector<std::string> & flist
                         )
{
   if(flist.size() == 0)
   {
      mxError("fitsFile", MXE_PARAMNOTSET, "Empty file list");
      return - 1;
   }

   long sz0 =0, sz1=0;

   for(int i=0;i<flist.size(); ++i)
   {
      if( fileName(flist[i], 1) < 0 ) return -1;

      if( read(im + i*sz0*sz1) < 0 ) return -1;

      sz0 = getSize(0);
      sz1 = getSize(1);
   }

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::read( dataT * im,
                           std::vector<fitsHeader> &heads,
                           const std::vector<std::string> & flist
                         )
{
   if(flist.size() == 0)
   {
      mxError("fitsFile", MXE_PARAMNOTSET, "Empty file list");
      return -1;
   }

   long sz0 =0, sz1=0;

   for(size_t i=0;i<flist.size(); ++i)
   {
      if( fileName(flist[i], 1) < 0 ) return -1;

      if( read(im + i*sz0*sz1) < 0 ) return -1;

      if( readHeader(heads[i]) < 0 ) return  -1;

      sz0 = getSize(0);
      sz1 = getSize(1);

   }

   return 0;
}


/************************************************************/
/***                      Eigen Arrays                    ***/
/************************************************************/

template<typename arrT, bool isCube=improc::is_eigenCube<arrT>::value>
struct eigenArrResize
{
   //If it's a cube, always pass zsz
   void resize(arrT & arr, int xsz, int ysz, int zsz)
   {
      arr.resize(xsz, ysz, zsz);
   }
};

template<typename arrT>
struct eigenArrResize<arrT, false>
{
   //If it's not a cube, never pass zsz
   void resize(arrT & arr, int xsz, int ysz, int zsz)
   {
      static_cast<void>(zsz);
      
      arr.resize(xsz, ysz);
   }
};

template<typename dataT>
template<typename arrT>
int fitsFile<dataT>::read(arrT & im)
{
   int fstatus = 0;

   if(!m_isOpen)
   {
      if( open() < 0 ) return -1;
   }

   //long long nelements = 1;
   long *fpix = new long[m_naxis];
   long *lpix = new long[m_naxis];
   long *inc = new long[m_naxis];
   eigenArrResize<arrT> arrresz;

   if(m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis == 2)
   {
      fpix[0] = m_x0+1;
      lpix[0] = fpix[0] + m_xpix-1;

      fpix[1] = m_y0+1;
      lpix[1] = fpix[1] + m_ypix-1;
      //nelements = m_xpix*m_ypix;

      inc[0] = 1;
      inc[1] = 1;
      arrresz.resize(im, m_xpix, m_ypix,1);
   }
   else
   {
      for(int i=0;i<m_naxis; i++)
      {
         fpix[i] = 1;
         lpix[i] = m_naxes[i];
         inc[i] = 1;
         //nelements *= m_naxes[i];
      }

      if(m_naxis > 2)
      {
         arrresz.resize(im, m_naxes[0], m_naxes[1], m_naxes[2]);
      }
      else if(m_naxis > 1)
      {
         arrresz.resize(im, m_naxes[0], m_naxes[1],1);
      }
      else
      {
         arrresz.resize(im, m_naxes[0], 1,1);
      }
   }

   if( fits_read_subset(m_fptr, fitsType<typename arrT::Scalar>(), fpix, lpix, inc, (void *) &m_nulval,
                                     (void *) im.data(), &m_anynul, &fstatus) < 0 ) return -1;

   if (fstatus && fstatus != 107)
   {
      std::string explan = "Error reading data from file";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILERERR, explan);

      delete[] fpix;
      delete[] lpix;
      delete[] inc;

      return -1;
   }
   delete[] fpix;
   delete[] lpix;
   delete[] inc;

   return 0;

}

template<typename dataT>
template<typename arrT>
int fitsFile<dataT>::read(arrT & data, fitsHeader &head)
{
   if( read(data) < 0 ) return -1;
   if( readHeader(head) < 0 ) return -1;
   return 0;
}

template<typename dataT>
template<typename arrT>
int fitsFile<dataT>::read( arrT & data,
                           const std::string &fname
                         )
{
   if( fileName(fname) < 0 ) return -1;
   if( read(data) < 0 ) return -1;
   return 0;
}

template<typename dataT>
template<typename arrT>
int fitsFile<dataT>::read( arrT & data,
                           fitsHeader &head,
                           const std::string &fname
                         )
{
   if( fileName(fname) < 0 ) return -1;
   if( read(data) < 0 ) return -1;
   if( readHeader(head) < 0 ) return  -1;

   return 0;
}


template<typename dataT>
template<typename cubeT>
int fitsFile<dataT>::read( cubeT & cube,
                           const std::vector<std::string> & flist,
                           std::vector<fitsHeader> * heads
                         )
{
   int fstatus = 0;

   if( flist.size() == 0 )
   {
      mxError("fitsFile", MXE_PARAMNOTSET, "Empty file list");
      return -1;
   }


   //Open the first file to get the dimensions.
   if( fileName(flist[0], 1) < 0) return -1;

   long *fpix = new long[m_naxis];
   long *lpix = new long[m_naxis];
   long *inc = new long[m_naxis];

   //Check if we're reading subsets
   if(m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis == 2)
   {
      fpix[0] = m_x0+1;
      lpix[0] = fpix[0] + m_xpix-1;

      fpix[1] = m_y0+1;
      lpix[1] = fpix[1] + m_ypix-1;
      //nelements = m_xpix*m_ypix;

      inc[0] = 1;
      inc[1] = 1;
      cube.resize(m_xpix, m_ypix, flist.size());
   }
   else
   {
      for(int i=0;i<m_naxis; i++)
      {
         fpix[i] = 1;
         lpix[i] = m_naxes[i];
         inc[i] = 1;
      }

      cube.resize( m_naxes[0], m_naxes[1], flist.size());
   }


   //Now read first image.
   fits_read_subset(m_fptr, fitsType<typename cubeT::Scalar>(), fpix, lpix, inc, (void *) &m_nulval,
                                     (void *) cube.image(0).data(), &m_anynul, &fstatus);


   if (fstatus && fstatus != 107)
   {

      std::string explan = "Error reading data from file";
      fitsErrText(explan, m_fileName, fstatus);

      mxError("cfitsio", MXE_FILERERR, explan);

      delete[] fpix;
      delete[] lpix;
      delete[] inc;

      return -1;
   }

   if(heads)
   {
      if( readHeader( (*heads)[0]) < 0) return -1;
   }


   //Now read in the rest.
   for(int i=1; i< flist.size(); ++i)
   {
      fileName(flist[i], 1);

      //Now read image.
      fits_read_subset(m_fptr, fitsType<typename cubeT::Scalar>(), fpix, lpix, inc, (void *) &m_nulval,
                                     (void *) cube.image(i).data(), &m_anynul, &fstatus);

      if (fstatus && fstatus != 107)
      {
         std::string explan = "Error reading data from file";
         fitsErrText(explan, m_fileName, fstatus);
         mxError("cfitsio", MXE_FILERERR, explan);

         delete[] fpix;
         delete[] lpix;
         delete[] inc;

         return -1;
      }

      if(heads)
      {
         if( readHeader( (*heads)[i]) < 0) return -1;
      }
   }


   delete[] fpix;
   delete[] lpix;
   delete[] inc;

   return 0;
}


template<typename dataT>
template<typename cubeT>
int fitsFile<dataT>::read( cubeT & cube,
                           std::vector<fitsHeader> & heads,
                           const std::vector<std::string> & flist
                         )
{
   return read(cube, flist, &heads);
}

template<typename dataT>
int fitsFile<dataT>::readHeader(fitsHeader & head)
{
   int fstatus = 0;

   char keyword[FLEN_KEYWORD];
   char value[FLEN_VALUE];
   char * comment;

   //The keys to look for if head is already populated
   std::list<fitsHeader::headerIterator> head_keys;
   std::list<fitsHeader::headerIterator>::iterator head_keys_it;
//   int num_head_keys;

   bool head_keys_only = false;
   if(head.size() > 0)
   {
      head_keys_only = true;
      fitsHeader::headerIterator headIt = head.begin();
      while(headIt != head.end())
      {
         head_keys.push_back(headIt);
         ++headIt;
      }
//      num_head_keys = head.size();
   }

   //If m_noComment is set, then we don't read in the comment
   if(m_noComment)
   {
      comment = 0;
   }
   else
   {
      comment = new char[FLEN_COMMENT];
   }

   int keysexist;
   int morekeys;

   if(!m_isOpen)
   {
      open();
   }

   //This gets the number of header keys to read
   fits_get_hdrspace(m_fptr, &keysexist, &morekeys, &fstatus);

   if (fstatus)
   {
      std::string explan = "Error reading header from file";
      fitsErrText( explan, m_fileName, fstatus);
      mxError( "fitsFile", MXE_FILERERR, explan);

      if(comment) delete[] comment;

      return -1;
   }

   for(int i=0; i<keysexist; i++)
   {
      fits_read_keyn(m_fptr, i+1, keyword, value, comment, &fstatus);

      if (fstatus)
      {
         std::string explan = "Error reading header from file";
         fitsErrText( explan, m_fileName, fstatus);
         mxError( "fitsFile", MXE_FILERERR, explan);

         if(comment) delete[] comment;

         return -1;
      }

      if(!head_keys_only)
      {
         if(strcmp(keyword, "COMMENT") == 0)
         {
            head.append<fitsCommentType>(keyword,  fitsCommentType(value), comment);
         }
         else if(strcmp(keyword, "HISTORY") == 0)
         {
            head.append<fitsHistoryType>(keyword,  fitsHistoryType(value), comment);
         }
         else
         {
            //Otherwise we append it as an unknown type
            head.append(keyword, value, comment);
         }
      }
      else
      {
         head_keys_it = head_keys.begin();
         while(head_keys_it != head_keys.end())
         {
            if( (*(*head_keys_it)).keyword() == keyword)
            {
               head[keyword].value(value);
               if(comment) head[keyword].comment(comment);

               head_keys.erase(head_keys_it);

               break;
            }
            ++head_keys_it;
         }

         //Quit if we're done.
         if(head_keys.empty()) break;
      }
   }

   if(comment) delete[] comment;

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::readHeader( fitsHeader &head,
                                 const std::string &fname
                               )
{
   if( fileName(fname) < 0 ) return -1;
   if( readHeader(head) < 0 ) return -1;
}

template<typename dataT>
int fitsFile<dataT>::readHeader( std::vector<fitsHeader> & heads,
                                 const std::vector<std::string> & flist
                               )
{
   for(int i=0; i< flist.size(); ++i)
   {
      fileName(flist[i], 1);

      if( readHeader( heads[i]) < 0) return -1;
   }

   return 0;
}


template<typename dataT>
int fitsFile<dataT>::write( const dataT * im,
                            int d1,
                            int d2,
                            int d3,
                            fitsHeader * head
                          )
{
   int fstatus = 0;

   if(m_isOpen) close();
   if(m_naxes) delete[] m_naxes;

   fstatus = 0;
   m_naxis = 1;
   if(d2 > 0)
   {
      if(d3 > 1) m_naxis = 3;
      else m_naxis = 2;
   }

   m_naxes = new long[m_naxis];

   m_naxes[0] = d1;
   if(m_naxis > 1) m_naxes[1] = d2;
   if(m_naxis > 2) m_naxes[2] = d3;

   std::string forceFileName = "!"+m_fileName;

   fits_create_file(&m_fptr, forceFileName.c_str(), &fstatus);
   if (fstatus)
   {
      std::string explan = "Error creating file";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILEOERR, explan);

      return -1;
   }
   m_isOpen = true;

   fits_create_img( m_fptr, fitsBITPIX<dataT>(), m_naxis, m_naxes, &fstatus);
   if (fstatus)
   {
      std::string explan = "Error creating image";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILEWERR, explan);

      return -1;
   }

   long fpixel[3];
   fpixel[0] = 1;
   fpixel[1] = 1;
   fpixel[2] = 1;

   LONGLONG nelements = 1;

   for(int i=0;i<m_naxis;++i) nelements *= m_naxes[i];

   fits_write_pix( m_fptr,  fitsType<dataT>(), fpixel, nelements, (void *) im, &fstatus);
   if (fstatus)
   {
      std::string explan = "Error writing data";
      fitsErrText(explan, m_fileName, fstatus);
      mxError("fitsFile", MXE_FILEWERR, explan);

      return -1;

   }

   if(head != 0)
   {
      fitsHeader::headerIterator it;

      for(it = head->begin(); it != head->end(); ++it)
      {
         int wrv = it->write(m_fptr);
         if(wrv != 0)
         {
            std::string explan = "Error writing keyword";
            fitsErrText(explan, m_fileName, wrv);
            mxError("fitsFile", MXE_FILEWERR, explan);
         }

      }
   }

   if( close() < 0) return -1;

   return 0;
}

template<typename dataT>
int fitsFile<dataT>::write( const dataT * im,
                            int d1,
                            int d2,
                            int d3
                          )
{
   return write(im, d1, d2, d3, (fitsHeader *) 0);
}

template<typename dataT>
int fitsFile<dataT>::write( const dataT * im,
                            int d1,
                            int d2,
                            int d3,
                            fitsHeader & head)
{
   return write(im, d1, d2, d3, &head);
}


template<typename dataT>
int fitsFile<dataT>::write( const std::string & fname,
                            const dataT * im,
                            int d1,
                            int d2,
                            int d3
                          )
{
   if( fileName(fname, false) < 0) return -1;
   return write(im, d1, d2, d3, (fitsHeader *) 0);
}

template<typename dataT>
int fitsFile<dataT>::write( const std::string & fname,
                            const dataT * im,
                            int d1,
                            int d2,
                            int d3,
                            fitsHeader & head)
{
   if( fileName(fname, false) < 0) return -1;
   return write(im, d1, d2, d3, &head);
}



template<typename dataT>
template<typename arrT>
int fitsFile<dataT>::write( const std::string & fname,
                            const arrT & im
                          )
{
   improc::eigenArrPlanes<arrT> planes;

   return write(fname, im.data(), im.rows(), im.cols(), planes(im));

}

template<typename dataT>
template<typename arrT>
int fitsFile<dataT>::write( const std::string & fname,
                            const arrT & im,
                            fitsHeader & head)
{
   improc::eigenArrPlanes<arrT> planes;

   return write(fname, im.data(), im.rows(), im.cols(), planes(im), head);

}

template<typename dataT>
void fitsFile<dataT>::setReadSize()
{
   m_x0 = -1;
   m_y0 = -1;
   m_xpix = -1;
   m_ypix = -1;
}

template<typename dataT>
void fitsFile<dataT>::setReadSize(long x0, long y0, long xpix, long ypix)
{
   m_x0 = x0;
   m_y0 = y0;
   m_xpix = xpix;
   m_ypix = ypix;
}


/** \ingroup fits_processing_typedefs
  * @{
  */

/// A \ref fitsFile to work in signed characters
typedef fitsFile<char> fitsFilec;

/// A \ref fitsFile to work in unsigned characters
typedef fitsFile<unsigned char> fitsFileuc;

/// A \ref fitsFile to work in signed short integers
typedef fitsFile<short> fitsFiles;

/// A \ref fitsFile to work in unsigned short integers
typedef fitsFile<unsigned short> fitsFileus;

/// A \ref fitsFile to work in signed integers
typedef fitsFile<int> fitsFilei;

/// A \ref fitsFile to work in unsigned integers
typedef fitsFile<unsigned int> fitsFileui;

/// A \ref fitsFile to work in signed long integers
typedef fitsFile<long> fitsFilel;

/// A \ref fitsFile to work in single precision floats
typedef fitsFile<float> fitsFilef;

/// A \ref fitsFile to work in double precision
typedef fitsFile<double> fitsFiled;

///@}

} //namespace fits
} //namespace mx

#endif //ioutils_fits_fitsFile_hpp
