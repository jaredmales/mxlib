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

#include <memory>

#include "../../mxlib.hpp"

#include "../../improc/eigenImage.hpp"

#include "fitsUtils.hpp"
#include "fitsHeader.hpp"

namespace mx
{
namespace fits
{

namespace fitsFileDetail
{

/** \cond */
/// Resettable CFITSIO operation table used by fitsFile and its deterministic failure tests.
struct fitsFileCfitsioOps
{
    int ( *openFile )( fitsfile **, const char *, int, int * ); ///< Open an existing FITS file.
    decltype( &ffgidm ) getImageDim;                            ///< Read the number of image axes.
    decltype( &ffgisz ) getImageSize;                           ///< Read image-axis sizes.
    decltype( &ffclos ) closeFile;                              ///< Close a FITS file.
    decltype( &ffgsv ) readSubset;                              ///< Read a FITS pixel subset.
    decltype( &ffghsp ) getHeaderSpace;                         ///< Read the number of FITS header cards.
    decltype( &ffgkyn ) readKey;                                ///< Read one FITS header card.
    decltype( &ffinit ) createFile;                             ///< Create a FITS file.
    decltype( &ffcrim ) createImage;                            ///< Create the primary FITS image.
    decltype( &ffppx ) writePixels;                             ///< Write FITS pixels.
    long *( *allocateLongs )( size_t );                         ///< Allocate a pixel-coordinate array.
};

/// Access the process-wide CFITSIO operation table.
fitsFileCfitsioOps &fitsFileCfitsioOpsInstance();

/// Restore every CFITSIO operation to its production implementation.
void resetFitsFileCfitsioOps();
/** \endcond */

} // namespace fitsFileDetail

/// Class to manage interactions with a FITS file
/** This class wraps the functionality of cfitsio.
 *
 * \tparam dataT is the datatype to use for in-memory storage of the image.
 * This does not have to match the data type
 * stored on disk for reading, but will be the type used for writing.
 *
 * \ingroup fits_processing
 */
template <typename dataT, class verboseT = verbose::d>
class fitsFile
{

    friend class fitsFile_test;

  public:
    typedef typename fitsHeader<verboseT>::headerIteratorT headerIteratorT;

  protected:
    /// The path to the file
    std::string m_fileName;

    /// The cfitsio data structure
    fitsfile *m_fptr{ nullptr };

    /// The dimensions of the image (1D, 2D, or 3D)
    int m_naxis{ 0 };

    /// The size of each dimension
    long m_naxes[3];

    /// Flag indicating whether the file is open or not
    bool m_isOpen{ false };

    /// The value to replace null values with
    dataT m_nulval{ 0 };

    /// Records whether any null values were replaced
    int m_anynul{ 0 };

    /// Flag to control whether the comment string is read.
    int m_noComment{ 0 };

    /// The starting x-pixel to read from
    long m_x0{ -1 };

    /// The starting y-pixel to read from
    long m_y0{ -1 };

    /// The number of x-pixels to read
    long m_xpix{ -1 };

    /// The number of y-pixels to read
    long m_ypix{ -1 };

    /// The starting frame to read from a cube
    long m_z0{ -1 };

    /// The number of frames to read from a cube
    long m_zframes{ -1 };

    /// One time initialization common to all constructors
    void construct();

  public:
    /// Default constructor
    fitsFile();

    /// Default constructor with error code
    fitsFile( error_t &errc /**< [out] error_t code indicating success or error */ );

    /// Constructor with file name and error code
    /** The file is not opened.
     *
     */
    fitsFile( const std::string &fname, ///< [in] File name to set on construction
              error_t &errc             /**< [out] error_t code indicating success or error */
    );

    /// Constructor with file name, and option to open.
    fitsFile( const std::string &fname, ///< [in] File name to set on construction
              bool doopen = true        ///< [in] [optional] If true, then the file is opened (the default).
    );

    /// Constructor with file name, option to open, and error code
    fitsFile( const std::string &fname, ///< File name to set on construction
              bool doopen,              ///< If true, then the file is opened (the default).
              error_t &errc             /**< [out] error_t code indicating success or error */
    );

    /// Destructor
    ~fitsFile();

    /// Get the current value of m_fileName
    /**
     * \returns the current value of m_fileName.
     */
    std::string fileName();

    /// Set the file path, and optionally open the file.
    /**
     * \returns mx::error_t::noerror on success
     * \returns mx::error_t codes from \ref close or \ref open
     */
    error_t fileName( const std::string &fname, ///< The new file name.
                      bool doopen = true        ///< If true, then the file is opened (the default).
    );

    /// Get the current value of m_naxis
    /**
     * \returns the current value of m_naxis
     *
     */
    int naxis();

    /// Get the current value of m_naxes for the specified dimension
    /**
     * \returns the current value of m_naxes for the specified dimension. -1 if no such dimension
     *
     */
    long naxes( int dim /**< [in] the dimension */ );

    /// Open the file and gets its dimensions.
    /** File name needs to already have been set.
     * If the file has already been opened, this returns immediately with no re-open.
     *
     * \returns mx::error_t::noerror on success
     * \returns mx::invalidconfig if m_filename is empty.
     * \returns mx::error_t::bad_alloc on an allocation error
     * \returns mx::error_t::std_exception on other errors during allocations
     * \returns mx::error_t::exception on unexpected errors during allocations
     * \returns mx::error_t::allocerr if allocation fails without an exception
     * \returns mx::error_t::fits_* codes from cfitsio functions
     *
     */
    error_t open();

    /// Open the file, first setting the file path.
    /**
     * \returns mx::error_t::noerror on success
     * \returns mx::invalidconfig if m_filename is empty.
     * \returns mx::error_t::bad_alloc on an allocation error
     * \returns mx::error_t::std_exception on other errors during allocations
     * \returns mx::error_t::exception on unexpected errors during allocations
     * \returns mx::error_t::allocerr if allocation fails without an exception
     * \returns mx::error_t::fits_* codes from cfitsio functions
     */
    error_t open( const std::string &fname /**< The name of the file to open. */ );

    /// Close the file.
    /**
     * \returns mx::error_t::noerror on success
     * \returns mx::error_t::fits_* codes from cfitsio functions
     */
    error_t close();

    /// Get the number of dimensions (i.e. m_naxis)
    int getDimensions( error_t &errc );

    /// Get the total size
    long getSize( error_t &errc );

    /// Get the size of a specific dimension
    long getSize( size_t axis, error_t &errc );

    /** \name Reading Basic Arrays
     * These methods read FITS data into basic or raw arrays specified by a pointer.
     * @{
     */

  protected:
    struct pixarrT
    {
        long *fpix{ nullptr }; ///< Populated with the lower left pixel to read.
        long *lpix{ nullptr }; ///< Populated with the upper right pixel to read.
        long *inc{ nullptr };  ///< The increment.

        error_t allocate( int naxis )
        {
            if( naxis <= 0 )
            {
                return internal::mxlib_error_report<verboseT>( error_t::paramnotset, "naxis" );
            }

            if( fpix )
            {
                delete[] fpix;
                fpix = nullptr;
            }

            if( lpix )
            {
                delete[] lpix;
                lpix = nullptr;
            }

            if( inc )
            {
                delete[] inc;
                inc = nullptr;
            }

            try
            {
                fpix = fitsFileDetail::fitsFileCfitsioOpsInstance().allocateLongs( naxis );
                lpix = fitsFileDetail::fitsFileCfitsioOpsInstance().allocateLongs( naxis );
                inc = fitsFileDetail::fitsFileCfitsioOpsInstance().allocateLongs( naxis );
            }
            catch( const std::bad_alloc &e )
            {
                internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                        std::string( "allocating pixel read arrays: " ) + e.what() );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
                return error_t::std_bad_alloc;
#else
                throw;
#endif
            }
            catch( const std::exception &e )
            {
                internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                        std::string( "allocating pixel read arrays: " ) + e.what() );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
                return error_t::exception;
#else
                throw;
#endif
            }
            catch( ... )
            {
                internal::mxlib_error_report<verboseT>( error_t::exception, "allocating pixel read arrays" );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
                return error_t::exception;
#else
                throw;
#endif
            }

            // also check for null allocations (in case compiled with exceptions off)
            if( fpix == nullptr )
            {
                return internal::mxlib_error_report<verboseT>( error_t::allocerr, "fpix is nullptr" );
            }

            if( lpix == nullptr )
            {
                return internal::mxlib_error_report<verboseT>( error_t::allocerr, "lpix is nullptr" );
            }

            if( inc == nullptr )
            {
                return internal::mxlib_error_report<verboseT>( error_t::allocerr, "inc is nullptr" );
            }

            return error_t::noerror;
        }

        ~pixarrT()
        {
            if( fpix )
            {
                delete[] fpix;
            }

            if( lpix )
            {
                delete[] lpix;
            }

            if( inc )
            {
                delete[] inc;
            }
        }
    };

    /// Fill in the read-size arrays for reading a subset (always used)
    /**
     */
    error_t calcPixarrs( pixarrT &pixarr /**< [out] Populated with the allocated read-size arrays*/ );

  public:
    /// Read the contents of the FITS file into an array.
    /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t read( dataT *data /**< [out] an allocated arrray large enough to hold the entire image */ );

    /// Read the contents of the FITS file into an array.
    /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t read( dataT *data,               ///< [out] an allocated arrray large enough to hold the entire image
                  fitsHeader<verboseT> &head ///< [out] a fitsHeader object which is passed to \ref readHeader
    );

    /// Read the contents of the FITS file into an array.
    /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t read( dataT *data,             ///< [out] is an allocated arrray large enough to hold the entire image
                  const std::string &fname ///< [in] is the file path, which is passed to \ref fileName
    );

    /// Read the contents of the FITS file into an array and read the header.
    /** The array pointed to by data must have been allocated.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t read( dataT *data,                ///< [out] an allocated arrray large enough to hold the entire image
                  fitsHeader<verboseT> &head, ///< [out] a fitsHeader object which is passed to \ref readHeader
                  const std::string &fname    ///< [in] the file path, which is passed to \ref fileName
    );

    /// Read data from a vector list of files into an image cube
    error_t read( dataT *im, ///< [out] An allocated array large enough to hold all the images
                  const std::vector<std::string> &flist ///< [in] The list of files to read.
    );

    /// Read data from a vector of files into an image cube with individual headers
    error_t read( dataT *im, /**< [out] An allocated array large enough to hold all the images */
                  std::vector<fitsHeader<verboseT>> &heads, /**< [in/out] The vector of fits headers, allocated
                                                                to contain one per image. */
                  const std::vector<std::string> &flist     /**< [in] The list of files to read. */
    );

    ///@}

    /** \name Reading Eigen Arrays
     * These methods read FITS data into array types with an Eigen-like interface.
     * @{
     */

    /// Read the contents of the FITS file into an Eigen array type (not a simple pointer).
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
    template <typename arrT>
    error_t read( arrT &data /**< [out] is the array, which will be resized as
        necessary using its resize(int, int) member */ );

    /// Read the contents of the FITS file into an Eigen array type (not a simple pointer).
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
    template <typename arrT>
    error_t read( arrT &data,                /**< [out] is the array, which will be resized as necessary
                                                        using its resize(int, int) member*/
                  fitsHeader<verboseT> &head ///< [out] is a fitsHeader object which is passed to \ref readHeader
    );

    /// Read the contents of the FITS file into an Eigen array type (not a simple pointer).
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
    template <typename arrT>
    error_t read( arrT &data,              /**< [out] is the array, which will be resized as
                                                      necessary using its resize(int, int) member*/
                  const std::string &fname ///< [in] is the file path, which is passed to \ref fileName
    );

    /// Read the contents of the FITS file into an Eigen array type (not a simple pointer).
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
    template <typename arrT>
    error_t read( arrT &data,                 /**< [out] the array, which will be resized as
                                                         necessary using its resize(int, int) member*/
                  fitsHeader<verboseT> &head, ///< [out] a fitsHeader object which is passed to \ref readHeader
                  const std::string &fname    ///< [in] the file path, which is passed to \ref fileName
    );

    /// Read data from a vector list of files into an image cube
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
    template <typename cubeT>
    error_t read( cubeT &cube, ///< [out] A cube which will be resized using its resize(int, int, int) member.
                  const std::vector<std::string> &flist, ///< [in] The list of files to read.
                  std::vector<fitsHeader<verboseT>> *heads =
                      0 ///< [out] [optional] A vector of fits headers, allocated to contain one per image.
    );

    /// Read data from a vector of files into an image cube with individual headers
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
    template <typename cubeT>
    error_t read( cubeT &cube,                              /**< [out] A cube which will be resized
                                                                       using its resize(int, int, int) member.*/
                  std::vector<fitsHeader<verboseT>> &heads, /**< [out] The vector of fits headers, allocated to
                                                             contain one per image.*/
                  const std::vector<std::string> &flist     ///< [in] The list of files to read.
    );

    ///@}

    /** \name Reading Headers
     * @{
     */

    /// Read the header from the fits file.
    /** If head is not empty, then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t readHeader( fitsHeader<verboseT> &head /**< [out] a fitsHeader object */ );

    /// Read the header from the fits file.
    /** If head is not empty, then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t readHeader( fitsHeader<verboseT> &head, ///< [out] a fitsHeader object
                        const std::string &fname    ///< [in] the file path, which is passed to \ref fileName
    );

    /// Read the headers from a list of FITS files.
    /** In each case, if the header is not empty then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t readHeader(
        std::vector<fitsHeader<verboseT>> &heads, /// A vector of fitsHeader objects to read into.
        const std::vector<std::string> &flist     ///< [in] A list of files, each of which is passed to \ref fileName
    );

    ///@}

    /** \name Writing Basic Arrays
     * These methods write basic arrays specified by a pointer to FITS.
     * @{
     */

    /// Write the contents of a raw array to the FITS file.
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t write( const dataT *im,           ///< [in] is the array
                   int d1,                    ///< [in] is the first dimension
                   int d2,                    ///< [in] is the second dimension
                   int d3,                    ///< [in] is the third dimenesion (minimum value is 1)
                   fitsHeader<verboseT> *head ///< [in] a pointer to the header.  Set to 0 if not used.
    );

    /// Write the contents of a raw array to the FITS file.
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t write( const dataT *im, ///< [in] is the array
                   int d1,          ///< [in] is the first dimension
                   int d2,          ///< [in] is the second dimension
                   int d3           ///< [in] is the third dimenesion (minimum value is 1)
    );

    /// Write the contents of a raw array to the FITS file.
    /**
     * Note: the type of the array must match dataT
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t write( const dataT *im,           ///< [in] is the array
                   int d1,                    ///< [in] is the first dimension
                   int d2,                    ///< [in] is the second dimension
                   int d3,                    ///< [in] is the third dimenesion (minimum value is 1)
                   fitsHeader<verboseT> &head ///< [in] is the header
    );

    /// Write the contents of a raw array to the FITS file.
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t write( const std::string &fname, ///< [in] is the name of the file.
                   const dataT *im,          ///< [in] is the array
                   int d1,                   ///< [in] is the first dimension
                   int d2,                   ///< [in] is the second dimension
                   int d3                    ///< [in] is the third dimenesion (minimum value is 1)
    );

    /// Write the contents of a raw array to the FITS file.
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    error_t write( const std::string &fname,  ///< [in] is the name of the file.
                   const dataT *im,           ///< [in] is the array
                   int d1,                    ///< [in] is the first dimension
                   int d2,                    ///< [in] is the second dimension
                   int d3,                    ///< [in] is the third dimenesion (minimum value is 1)
                   fitsHeader<verboseT> &head ///< [in] is the header
    );

    ///@}

    /** \name Writing Eigen Arrays
     * These methods write array types with an Eigen-like interface.
     * @{
     */

    /// Write the contents of an Eigen-type array to a FITS file.
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
    template <typename arrT>
    error_t write( const std::string &fname, ///< [in] is the name of the file.
                   const arrT &im            ///< [in] is the array
    );

    /// Write the contents of an Eigen-type array to a FITS file.
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
    template <typename arrT>
    error_t write( const std::string &fname,  ///< [in] is the file path, which is passed to \ref fileName
                   const arrT &im,            ///< [in] is the array
                   fitsHeader<verboseT> &head ///< [in] is a fitsHeader object which is passed to \ref readHeader
    );

    // int writeHeader( fitsHeader<verboseT> &head );

    ///@}

    /** \name Reading Subsets
     * It is often desirable to read only a subset of an image or images into memory.  These methods
     * allow you to specify this.
     * @{
     */

    /// Set to read all the pixels in the file
    void setReadSize();

    /// Set to read only a subset of the pixels in the file
    /**
     */
    void setReadSize( long x0,   ///< is the starting x-pixel to read
                      long y0,   ///< is the starting y-pixel to read
                      long xpix, ///< is the number of x-pixels to read
                      long ypix  ///< is the number of y-pixels to read
    );

    /// Set to read all frames from a cube.
    /**
     */
    void setCubeReadSize();

    /// Set the number of frames to read from a cube.
    /**
     *
     */
    void setCubeReadSize( long z0,     ///< is the starting frame to read
                          long zframes ///< is the number of frames to read
    );

    ///@}

}; // fitsFile

template <typename dataT, class verboseT>
fitsFile<dataT, verboseT>::fitsFile()
{
}

template <typename dataT, class verboseT>
fitsFile<dataT, verboseT>::fitsFile( error_t &errc )
{
    errc = error_t::noerror;
}

template <typename dataT, class verboseT>
fitsFile<dataT, verboseT>::fitsFile( const std::string &fname, error_t &errc )
{
    // no errors are actually possible
    errc = internal::mxlib_error_report<verboseT>( fileName( fname, false ) ); // nothing printed if noerror
}

template <typename dataT, class verboseT>
fitsFile<dataT, verboseT>::fitsFile( const std::string &fname, bool doopen )
{
    // If an error happens on open(), then m_open will be false and this will persist
    // so no need to throw an exception
    internal::mxlib_error_report<verboseT>( fileName( fname, doopen ) ); // nothing printed if noerror
}

template <typename dataT, class verboseT>
fitsFile<dataT, verboseT>::fitsFile( const std::string &fname, bool doopen, error_t &errc )
{
    // If an error happens on open(), then m_open will be false and this will persist
    // so no need to throw an exception
    errc = internal::mxlib_error_report<verboseT>( fileName( fname, doopen ) ); // nothing printed if noerror
}

template <typename dataT, class verboseT>
fitsFile<dataT, verboseT>::~fitsFile()
{
    if( m_isOpen )
    {
        internal::mxlib_error_report<verboseT>( close() );
    }
}

template <typename dataT, class verboseT>
std::string fitsFile<dataT, verboseT>::fileName()
{
    return m_fileName;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::fileName( const std::string &fname, bool doopen )
{
    if( m_isOpen )
    {
        mxlib_error_check( close() );
    }

    m_fileName = fname;

    if( doopen )
    {
        mxlib_error_check( open() );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
int fitsFile<dataT, verboseT>::naxis()
{
    return m_naxis;
}

template <typename dataT, class verboseT>
long fitsFile<dataT, verboseT>::naxes( int dim )
{
    if( dim >= m_naxis || dim > 2 )
    {
        return -1;
    }

    return m_naxes[dim];
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::open()
{
    if( m_isOpen ) // no-op
    {
        return error_t::noerror;
    }

    if( m_fileName == "" )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidconfig, "File name is not set" );
    }

    int fstatus = 0;

    fitsFileDetail::fitsFileCfitsioOpsInstance().openFile( &m_fptr, m_fileName.c_str(), READONLY, &fstatus );

    if( fstatus )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ), "Opening file " + m_fileName );
    }

    fstatus = 0;
    fitsFileDetail::fitsFileCfitsioOpsInstance().getImageDim( m_fptr, &m_naxis, &fstatus );
    if( fstatus )
    {
        int closeStatus = 0;
        ffclos( m_fptr, &closeStatus );
        m_fptr = nullptr;
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Getting number of axes in file " + m_fileName );
    }

    fstatus = 0;
    fitsFileDetail::fitsFileCfitsioOpsInstance().getImageSize( m_fptr, m_naxis, m_naxes, &fstatus );
    if( fstatus )
    {
        int closeStatus = 0;
        ffclos( m_fptr, &closeStatus );
        m_fptr = nullptr;
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Getting dimensions in file " + m_fileName );
    }

    m_isOpen = true; // Only set this after opening is complete.

    return error_t::noerror;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::open( const std::string &fname )
{
    mxlib_error_return( fileName( fname, true ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::close()
{
    if( !m_isOpen )
    {
        return error_t::noerror; // No error.
    }

    int fstatus = 0;
    fitsFileDetail::fitsFileCfitsioOpsInstance().closeFile( m_fptr, &fstatus );

    if( fstatus )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ), "Closing file " + m_fileName );
    }

    m_fptr = nullptr;
    m_isOpen = false;

    return error_t::noerror;
}

template <typename dataT, class verboseT>
int fitsFile<dataT, verboseT>::getDimensions( error_t &errc )
{
    if( !m_isOpen )
    {
        errc = error_t::invalidconfig;
        return -1;
    }

    return m_naxis;
}

template <typename dataT, class verboseT>
long fitsFile<dataT, verboseT>::getSize( error_t &errc )
{
    if( !m_isOpen )
    {
        errc = error_t::invalidconfig;
        return -1;
    }

    long sz = 1;

    errc = error_t::noerror;

    if( m_naxis == 1 )
    {
        if( m_x0 > -1 && m_xpix > -1 )
        {
            return m_xpix;
        }
        return m_naxes[0];
    }
    else if( m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis == 2 )
    {
        return m_xpix * m_ypix;
    }
    else
    {
        for( int i = 0; i < m_naxis && i < 3; ++i )
        {
            sz *= m_naxes[i];
        }
    }

    return sz;
}

template <typename dataT, class verboseT>
long fitsFile<dataT, verboseT>::getSize( size_t axis, error_t &errc )
{
    if( !m_isOpen )
    {
        errc = error_t::invalidconfig;
        return -1;
    }

    if( axis >= m_naxis || axis > 2 )
    {
        errc = error_t::invalidarg;
        return -1;
    }

    errc = error_t::noerror;

    if( m_naxis == 1 && m_x0 > -1 && m_xpix > -1 )
    {
        return m_xpix;
    }
    else if( m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis == 2 )
    {
        if( axis == 0 )
        {
            return m_xpix;
        }
        return m_ypix;
    }
    else
    {
        return m_naxes[axis];
    }
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::calcPixarrs( pixarrT &pixarr )
{
    error_t errc = pixarr.allocate( m_naxis );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    if( m_naxis == 1 )
    {
        if( m_x0 > -1 && m_xpix > -1 )
        {
            pixarr.fpix[0] = m_x0 + 1;
            pixarr.lpix[0] = pixarr.fpix[0] + m_xpix - 1;
        }
        else
        {
            pixarr.fpix[0] = 1;
            pixarr.lpix[0] = m_naxes[0];
        }
        pixarr.inc[0] = 1;
    }
    else if( m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 && m_naxis == 2 )
    {
        pixarr.fpix[0] = m_x0 + 1;
        pixarr.lpix[0] = pixarr.fpix[0] + m_xpix - 1;
        pixarr.fpix[1] = m_y0 + 1;
        pixarr.lpix[1] = pixarr.fpix[1] + m_ypix - 1;

        pixarr.inc[0] = 1;
        pixarr.inc[1] = 1;
    }
    else
    {
        if( m_x0 < 0 && m_y0 < 0 && m_xpix < 0 && m_ypix < 0 && m_z0 < 0 && m_zframes < 0 )
        {
            for( int i = 0; i < m_naxis && i < 3; i++ )
            {
                pixarr.fpix[i] = 1;
                pixarr.lpix[i] = m_naxes[i];
                pixarr.inc[i] = 1;
            }
        }
        else
        {
            if( m_x0 > -1 && m_y0 > -1 && m_xpix > -1 && m_ypix > -1 )
            {
                pixarr.fpix[0] = m_x0 + 1;
                pixarr.lpix[0] = pixarr.fpix[0] + m_xpix - 1;
                pixarr.fpix[1] = m_y0 + 1;
                pixarr.lpix[1] = pixarr.fpix[1] + m_ypix - 1;

                pixarr.inc[0] = 1;
                pixarr.inc[1] = 1;
            }
            else
            {
                pixarr.fpix[0] = 1;
                pixarr.lpix[0] = m_naxes[0];
                pixarr.fpix[1] = 1;
                pixarr.lpix[1] = m_naxes[1];

                pixarr.inc[0] = 1;
                pixarr.inc[1] = 1;
            }

            if( m_z0 > -1 && m_zframes > -1 )
            {
                pixarr.fpix[2] = m_z0 + 1;
                pixarr.lpix[2] = pixarr.fpix[2] + m_zframes - 1;
                pixarr.inc[2] = 1;
            }
            else
            {
                pixarr.fpix[2] = 1;
                pixarr.lpix[2] = m_naxes[2];
                pixarr.inc[2] = 1;
            }
        }
    }

    return error_t::noerror;
}

/************************************************************/
/***                      Basic Arrays                    ***/
/************************************************************/

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::read( dataT *data )
{
    if( !m_isOpen )
    {
        mxlib_error_check( open() );
    }

    pixarrT pixarrs;

    mxlib_error_check( calcPixarrs( pixarrs ) );

    ///\todo test if there is a speed difference for full reads for fits_read_pix/subset

    int fstatus = 0;

    fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset( m_fptr,
                                                             fitsType<dataT>(),
                                                             pixarrs.fpix,
                                                             pixarrs.lpix,
                                                             pixarrs.inc,
                                                             (void *)&m_nulval,
                                                             (void *)data,
                                                             &m_anynul,
                                                             &fstatus );

    if( fstatus && fstatus != END_OF_FILE )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Reading data from " + m_fileName );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::read( dataT *data, fitsHeader<verboseT> &head )
{
    mxlib_error_check( read( data ) );

    mxlib_error_return( readHeader( head ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::read( dataT *data, const std::string &fname )
{
    mxlib_error_check( fileName( fname ) );

    mxlib_error_return( read( data ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::read( dataT *data, fitsHeader<verboseT> &head, const std::string &fname )
{
    mxlib_error_check( fileName( fname ) );

    mxlib_error_check( read( data ) );

    mxlib_error_return( readHeader( head ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::read( dataT *im, const std::vector<std::string> &flist )
{
    if( flist.size() == 0 )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "Empty file list" );
    }

    long sz0 = 0, sz1 = 0;

    for( int i = 0; i < flist.size(); ++i )
    {
        mxlib_error_check( fileName( flist[i], 1 ) );

        mxlib_error_check( read( im + i * sz0 * sz1 ) );

        error_t errc;
        sz0 = getSize( 0, errc );
        mxlib_error_check( errc );

        sz1 = getSize( 1, errc );
        mxlib_error_check( errc );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::read( dataT *im,
                                         std::vector<fitsHeader<verboseT>> &heads,
                                         const std::vector<std::string> &flist )
{
    if( flist.size() == 0 )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "Empty file list" );
    }

    if( heads.size() != flist.size() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                       "Header vector is not the same size as the file list" );
    }

    long sz0 = 0, sz1 = 0;

    for( size_t i = 0; i < flist.size(); ++i )
    {

        mxlib_error_check( fileName( flist[i], 1 ) );

        mxlib_error_check( read( im + i * sz0 * sz1 ) );

        mxlib_error_check( readHeader( heads[i] ) );

        error_t errc;
        sz0 = getSize( 0, errc );
        mxlib_error_check( errc );

        sz1 = getSize( 1, errc );
        mxlib_error_check( errc );
    }

    return error_t::noerror;
}

/************************************************************/
/***                      Eigen Arrays                    ***/
/************************************************************/

template <typename arrT, class verboseT, bool isCube = improc::is_eigenCube<arrT>::value>
struct eigenArrResize
{
    // If it's a cube, always pass zsz
    error_t resize( arrT &arr, int xsz, int ysz, int zsz )
    {
        try
        {
            arr.resize( xsz, ysz, zsz );
        }
        catch( const std::bad_alloc &e )
        {
            internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                    std::string( "resizing array: " ) + e.what() );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
            return error_t::std_bad_alloc;
#else
            throw;
#endif
        }
        catch( const std::exception &e )
        {
            internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                    std::string( "resizing array: " ) + e.what() );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
            return error_t::std_exception;
#else
            throw;
#endif
        }
        catch( ... )
        {
            internal::mxlib_error_report<verboseT>( error_t::exception, "resizing array" );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
            return error_t::exception;
#else
            throw;
#endif
        }

        return error_t::noerror;
    }
};

template <typename arrT, class verboseT>
struct eigenArrResize<arrT, verboseT, false>
{
    // If it's not a cube, never pass zsz
    error_t resize( arrT &arr, int xsz, int ysz, [[maybe_unused]] int zsz )
    {
        try
        {
            arr.resize( xsz, ysz );
        }
        catch( const std::bad_alloc &e )
        {
            internal::mxlib_error_report<verboseT>( error_t::std_bad_alloc,
                                                    std::string( "resizing array: " ) + e.what() );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
            return error_t::std_bad_alloc;
#else
            throw;
#endif
        }
        catch( const std::exception &e )
        {
            internal::mxlib_error_report<verboseT>( error_t::std_exception,
                                                    std::string( "resizing array: " ) + e.what() );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
            return error_t::std_exception;
#else
            throw;
#endif
        }
        catch( ... )
        {
            internal::mxlib_error_report<verboseT>( error_t::exception, "resizing array" );
#ifdef MXLIB_TRAP_ALLOC_ERRORS
            return error_t::exception;
#else
            throw;
#endif
        }

        return error_t::noerror;
    }
};

template <typename dataT, class verboseT>
template <typename arrT>
error_t fitsFile<dataT, verboseT>::read( arrT &im )
{
    ///\todo this can probably be made part of one read function (or call read(data *)) with a call to resize with
    /// SFINAE
    int fstatus = 0;

    if( !m_isOpen )
    {
        mxlib_error_check( open() );
    }

    pixarrT pixarrs;
    mxlib_error_check( calcPixarrs( pixarrs ) );

    eigenArrResize<arrT, verboseT> arrresz;
    if( m_naxis > 2 )
    {
        mxlib_error_check( arrresz.resize( im,
                                           pixarrs.lpix[0] - pixarrs.fpix[0] + 1,
                                           pixarrs.lpix[1] - pixarrs.fpix[1] + 1,
                                           pixarrs.lpix[2] - pixarrs.fpix[2] + 1 ) );
    }
    else if( m_naxis > 1 )
    {
        mxlib_error_check(
            arrresz.resize( im, pixarrs.lpix[0] - pixarrs.fpix[0] + 1, pixarrs.lpix[1] - pixarrs.fpix[1] + 1, 1 ) );
    }
    else
    {
        mxlib_error_check( arrresz.resize( im, pixarrs.lpix[0] - pixarrs.fpix[0] + 1, 1, 1 ) );
    }

    fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset( m_fptr,
                                                             fitsType<typename arrT::Scalar>(),
                                                             pixarrs.fpix,
                                                             pixarrs.lpix,
                                                             pixarrs.inc,
                                                             (void *)&m_nulval,
                                                             (void *)im.data(),
                                                             &m_anynul,
                                                             &fstatus );

    if( fstatus && fstatus != END_OF_FILE )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Reading data from " + m_fileName );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
template <typename arrT>
error_t fitsFile<dataT, verboseT>::read( arrT &data, fitsHeader<verboseT> &head )
{
    error_t errc;
    errc = read( data );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    errc = readHeader( head );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
template <typename arrT>
error_t fitsFile<dataT, verboseT>::read( arrT &data, const std::string &fname )
{
    error_t errc;
    errc = fileName( fname );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    errc = read( data );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }
    return error_t::noerror;
}

template <typename dataT, class verboseT>
template <typename arrT>
error_t fitsFile<dataT, verboseT>::read( arrT &data, fitsHeader<verboseT> &head, const std::string &fname )
{
    error_t errc;
    errc = fileName( fname );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    errc = read( data );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    errc = readHeader( head );
    if( errc != error_t::noerror )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
template <typename cubeT>
error_t fitsFile<dataT, verboseT>::read( cubeT &cube,
                                         const std::vector<std::string> &flist,
                                         std::vector<fitsHeader<verboseT>> *heads )
{
    error_t errc;
    int fstatus = 0;

    if( flist.size() == 0 )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "Empty file list" );
    }

    if( heads && heads->size() != flist.size() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                       "Header vector is not the same size as the file list" );
    }

    // Open the first file to get the dimensions.
    errc = fileName( flist[0], 1 );
    if( !!errc )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    pixarrT pixarrs;
    errc = calcPixarrs( pixarrs );

    if( !!errc )
    {
        return internal::mxlib_error_report<verboseT>( errc );
    }

    cube.resize( pixarrs.lpix[0] - pixarrs.fpix[0] + 1, pixarrs.lpix[1] - pixarrs.fpix[1] + 1, flist.size() );

    // Now read first image.
    fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset( m_fptr,
                                                             fitsType<typename cubeT::Scalar>(),
                                                             pixarrs.fpix,
                                                             pixarrs.lpix,
                                                             pixarrs.inc,
                                                             (void *)&m_nulval,
                                                             (void *)cube.image( 0 ).data(),
                                                             &m_anynul,
                                                             &fstatus );

    if( fstatus && fstatus != END_OF_FILE )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Reading data from " + m_fileName );
    }

    if( heads )
    {
        errc = readHeader( ( *heads )[0] );
        if( errc != error_t::noerror )
        {
            return internal::mxlib_error_report<verboseT>( errc );
        }
    }

    // Now read in the rest.
    for( int i = 1; i < flist.size(); ++i )
    {
        errc = fileName( flist[i], 1 );
        if( errc != error_t::noerror )
        {
            return internal::mxlib_error_report<verboseT>( errc );
        }

        // Now read image.
        fitsFileDetail::fitsFileCfitsioOpsInstance().readSubset( m_fptr,
                                                                 fitsType<typename cubeT::Scalar>(),
                                                                 pixarrs.fpix,
                                                                 pixarrs.lpix,
                                                                 pixarrs.inc,
                                                                 (void *)&m_nulval,
                                                                 (void *)cube.image( i ).data(),
                                                                 &m_anynul,
                                                                 &fstatus );

        if( fstatus && fstatus != END_OF_FILE )
        {
            return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                           "Reading data from " + m_fileName );
        }

        if( heads )
        {
            errc = readHeader( ( *heads )[i] );
            if( errc != error_t::noerror )
            {
                return internal::mxlib_error_report<verboseT>( errc );
            }
        }
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
template <typename cubeT>
error_t fitsFile<dataT, verboseT>::read( cubeT &cube,
                                         std::vector<fitsHeader<verboseT>> &heads,
                                         const std::vector<std::string> &flist )
{
    mxlib_error_return( read( cube, flist, &heads ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::readHeader( fitsHeader<verboseT> &head )
{
    int fstatus = 0;

    char keyword[FLEN_KEYWORD];
    char value[FLEN_VALUE];
    std::unique_ptr<char[]> commentStorage;
    char *comment = nullptr;

    // The keys to look for if head is already populated
    typename std::list<headerIteratorT> head_keys;
    typename std::list<headerIteratorT>::iterator head_keys_it;
    //   int num_head_keys;

    bool head_keys_only = false;
    if( head.size() > 0 )
    {
        head_keys_only = true;
        headerIteratorT headIt = head.begin();
        while( headIt != head.end() )
        {
            head_keys.push_back( headIt );
            ++headIt;
        }
        //      num_head_keys = head.size();
    }

    // If m_noComment is set, then we don't read in the comment
    if( m_noComment )
    {
        comment = nullptr;
    }
    else
    {
        commentStorage = std::make_unique<char[]>( FLEN_COMMENT );
        comment = commentStorage.get();
    }

    int keysexist;
    int morekeys;

    if( !m_isOpen )
    {
        mxlib_error_check( open() );
    }

    // This gets the number of header keys to read
    fitsFileDetail::fitsFileCfitsioOpsInstance().getHeaderSpace( m_fptr, &keysexist, &morekeys, &fstatus );

    if( fstatus )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Reading header from " + m_fileName );
    }

    for( int i = 0; i < keysexist; i++ )
    {
        fitsFileDetail::fitsFileCfitsioOpsInstance().readKey( m_fptr, i + 1, keyword, value, comment, &fstatus );

        if( fstatus )
        {
            return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                           "Reading header from " + m_fileName );
        }

        if( !head_keys_only )
        {
            if( strcmp( keyword, "COMMENT" ) == 0 )
            {
                head.template append<fitsCommentType>( keyword, fitsCommentType( value ), comment );
            }
            else if( strcmp( keyword, "HISTORY" ) == 0 )
            {
                head.template append<fitsHistoryType>( keyword, fitsHistoryType( value ), comment );
            }
            else
            {
                // Otherwise we append it as an unknown type
                head.append( keyword, value, comment );
            }
        }
        else
        {
            head_keys_it = head_keys.begin();
            while( head_keys_it != head_keys.end() )
            {
                if( ( *( *head_keys_it ) ).keyword() == keyword )
                {
                    head[keyword].value( (const char *)value );
                    if( comment )
                    {
                        head[keyword].comment( comment );
                    }

                    head_keys.erase( head_keys_it );

                    break;
                }
                ++head_keys_it;
            }

            // Quit if we're done.
            if( head_keys.empty() )
            {
                break;
            }
        }
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::readHeader( fitsHeader<verboseT> &head, const std::string &fname )
{
    error_t errc;
    mxlib_error_check( fileName( fname ) );
    mxlib_error_check( readHeader( head ) );
    return error_t::noerror;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::readHeader( std::vector<fitsHeader<verboseT>> &heads,
                                               const std::vector<std::string> &flist )
{
    if( heads.size() != 0 && heads.size() != flist.size() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                       "head vector is not empty and not same size as file list" );
    }

    if( heads.size() == 0 )
    {
        heads.resize( flist.size() );
    }

    error_t errc;
    for( int i = 0; i < flist.size(); ++i )
    {
        mxlib_error_check( fileName( flist[i], 1 ) );

        mxlib_error_check( readHeader( heads[i] ) );
    }

    return error_t::noerror;
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::write( const dataT *im, int d1, int d2, int d3, fitsHeader<verboseT> *head )
{
    int fstatus = 0;

    if( m_isOpen )
    {
        mxlib_error_check( close() );
    }

    m_naxis = 1;
    if( d2 > 0 )
    {
        if( d3 > 1 )
        {
            m_naxis = 3;
        }
        else
        {
            m_naxis = 2;
        }
    }

    m_naxes[0] = d1;
    if( m_naxis > 1 )
    {
        m_naxes[1] = d2;
    }
    if( m_naxis > 2 )
    {
        m_naxes[2] = d3;
    }

    std::string forceFileName = "!" + m_fileName;

    fitsFileDetail::fitsFileCfitsioOpsInstance().createFile( &m_fptr, forceFileName.c_str(), &fstatus );
    if( fstatus )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ), "Creating " + m_fileName );
    }
    m_isOpen = true;

    fstatus = 0;
    fitsFileDetail::fitsFileCfitsioOpsInstance().createImage( m_fptr, fitsBITPIX<dataT>(), m_naxis, m_naxes, &fstatus );
    if( fstatus )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ),
                                                       "Creating image in" + m_fileName );
    }

    long fpixel[3];
    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 1;

    LONGLONG nelements = 1;

    for( int i = 0; i < m_naxis && i < 3; ++i )
    {
        nelements *= m_naxes[i];
    }

    fstatus = 0;
    fitsFileDetail::fitsFileCfitsioOpsInstance()
        .writePixels( m_fptr, fitsType<dataT>(), fpixel, nelements, (void *)im, &fstatus );
    if( fstatus )
    {
        return internal::mxlib_error_report<verboseT>( fits_status2error_t( fstatus ), "Writing data " + m_fileName );
    }

    if( head != 0 )
    {
        headerIteratorT it;

        for( it = head->begin(); it != head->end(); ++it )
        {
            error_t errc = it->write( m_fptr );
            if( errc != error_t::noerror )
            {
                return internal::mxlib_error_report<verboseT>( errc, "Writing keyword " + m_fileName );
            }
        }
    }

    mxlib_error_return( close() )
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::write( const dataT *im, int d1, int d2, int d3 )
{
    mxlib_error_return( write( im, d1, d2, d3, (fitsHeader<verboseT> *)0 ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::write( const dataT *im, int d1, int d2, int d3, fitsHeader<verboseT> &head )
{
    mxlib_error_return( write( im, d1, d2, d3, &head ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::write( const std::string &fname, const dataT *im, int d1, int d2, int d3 )
{
    mxlib_error_check( fileName( fname, false ) );

    mxlib_error_return( write( im, d1, d2, d3, (fitsHeader<verboseT> *)0 ) );
}

template <typename dataT, class verboseT>
error_t fitsFile<dataT, verboseT>::write(
    const std::string &fname, const dataT *im, int d1, int d2, int d3, fitsHeader<verboseT> &head )
{
    mxlib_error_check( fileName( fname, false ) );

    mxlib_error_return( write( im, d1, d2, d3, &head ) );
}

template <typename dataT, class verboseT>
template <typename arrT>
error_t fitsFile<dataT, verboseT>::write( const std::string &fname, const arrT &im )
{
    improc::eigenArrPlanes<arrT> planes;

    mxlib_error_return( write( fname, im.data(), im.rows(), im.cols(), planes( im ) ) );
}

template <typename dataT, class verboseT>
template <typename arrT>
error_t fitsFile<dataT, verboseT>::write( const std::string &fname, const arrT &im, fitsHeader<verboseT> &head )
{
    improc::eigenArrPlanes<arrT> planes;

    return write( fname, im.data(), im.rows(), im.cols(), planes( im ), head );
}

template <typename dataT, class verboseT>
void fitsFile<dataT, verboseT>::setReadSize()
{
    m_x0 = -1;
    m_y0 = -1;
    m_xpix = -1;
    m_ypix = -1;
}

template <typename dataT, class verboseT>
void fitsFile<dataT, verboseT>::setReadSize( long x0, long y0, long xpix, long ypix )
{
    m_x0 = x0;
    m_y0 = y0;
    m_xpix = xpix;
    m_ypix = ypix;
}

template <typename dataT, class verboseT>
void fitsFile<dataT, verboseT>::setCubeReadSize()
{
    m_z0 = -1;
    m_zframes = -1;
}

template <typename dataT, class verboseT>
void fitsFile<dataT, verboseT>::setCubeReadSize( long z0, long zframes )
{
    m_z0 = z0;
    m_zframes = zframes;
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

extern template class fitsFile<short, verbose::d>;
extern template class fitsFile<unsigned short, verbose::d>;
extern template class fitsFile<int, verbose::d>;
extern template class fitsFile<unsigned int, verbose::d>;
extern template class fitsFile<float, verbose::d>;
extern template class fitsFile<double, verbose::d>;

} // namespace fits
} // namespace mx

#endif // ioutils_fits_fitsFile_hpp
