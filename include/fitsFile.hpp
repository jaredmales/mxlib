/** \file fitsFile.hpp
  * \brief Declares and defines a class to work with a FITS file
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __fitsFile_hpp__
#define __fitsFile_hpp__


#include "fitsUtils.hpp"

#include "fitsHeader.hpp"
#include "mxException.hpp"
#include "eigenCube.hpp"

namespace mx
{

/** \addtogroup fits_processing
  * @{
  */


/// Class to manage interactions with a FITS file
/** This class conveniently wraps the functionality of cfitsio.
  * 
  * \tparam dataT is the datatype to use for in-memory storage of the image.  This does not have to match the data type stored on disk.
  * 
  */
template<typename dataT> class fitsFile
{

protected:
   
   ///The path to the file
   std::string fileName;
   
   ///The cfitsio data structure
   fitsfile * fptr;
   
   ///The status code returned by the cfitsio library
   int fstatus;
   
   ///The dimensions of the image (1D, 2D, 3D etc)
   int naxis;
   
   ///The size of each dimension
   long * naxes;

   ///Flag indicating whether the file is open or not
   bool isOpen;

   ///The value to replace null values with
   dataT nulval;
   
   ///Records whether any null values were replaced
   int anynul;
   
   ///Flag to control whether the comment string is read.
   int noComment; 
   
   ///One time initialization common to all constructors
   void construct();
  
   ///The starting x-pixel to read from
   long _x0;
   
   ///The starting y-pixel to read from
   long _y0;
   
   ///The number of x-pixels to read
   long _xpix;
   
   ///The number of y-pixels to read
   long _ypix;
   
public:
   
   ///Default constructor
   fitsFile();
   
   ///Constructor with fileName, and option to open.
   fitsFile(const std::string & fname, bool doopen = 1);
   
   ///Destructor
   ~fitsFile();
   
   ///Set the file path, and optionally open the file.
   void setFilename(const std::string & fname, bool doopen = 1);
   
   ///Open the file
   void open();
   
   ///Open the file, first setting the file path.
   void open(const std::string & fname);
   
   ///Close the file.
   void close();

   ///Get the number of dimensions (i.e. naxis)
   int getDimensions();
   
   ///Get the total size
   long getSize();
   
   ///Get the size of a specific dimension
   long getSize(size_t axis);
   
   /** \name Reading Basic Arrays
     * @{
     */
   
   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \param data an allocated arrray large enough to hold the entire image 
     * 
     * \throws mxException on error
     * 
     */ 
   void read(dataT * data);

   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \param data is an allocated arrray large enough to hold the entire image 
     * \param head is a fitsHeader object which is passed to \ref readHeader
     * 
     * \throws mxException is thrown on error
     * 
     */ 
   void read(dataT * data, fitsHeader &head);
   
   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \param fname is the file path, which is passed to \ref setFilename
     * \param data is an allocated arrray large enough to hold the entire image 
     * 
     * \throws mxException is thrown on error
     * 
     */ 
   void read(const std::string &fname, dataT * data);
   
   ///Read the contents of the FITS file into an array.
   /** The array pointed to by data must have been allocated.
     *
     * \param fname is the file path, which is passed to \ref setFilename
     * \param data is an allocated arrray large enough to hold the entire image 
     * \param head is a fitsHeader object which is passed to \ref readHeader
     * 
     * \throws mxException is thrown on error
     * 
     */ 
   void read(const std::string &fname, dataT * data, fitsHeader &head);

   
   ///@}
   
   /** \name Reading Eigen Arrays
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
     * \param data is the array, which will be resized as necessary using its resize(int, int) member
     * 
     * \throws mxException on error
     * 
     */ 
   template<typename arrT>
   void read(arrT & data);
   
   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     * 
     * \tparam arrT is the type of array, see requirements above.
     * 
     * \param data is the array, which will be resized as necessary using its resize(int, int) member
     * \param head is a fitsHeader object which is passed to \ref readHeader
     * 
     * \throws mxException on error
     * 
     */ 
   template<typename arrT>
   void read(arrT & data, fitsHeader &head);
   
   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     * 
     * \tparam arrT is the type of array, see requirements above.
     * 
     * \param fname is the file path, which is passed to \ref setFilename
     * \param data is the array, which will be resized as necessary using its resize(int, int) member
     * 
     * \throws mxException on error
     * 
     */ 
   template<typename arrT>
   void read(const std::string &fname, arrT & data);
   
   ///Read the contents of the FITS file into an Eigen array type (not a simple pointer).
   /** The type arrT can be any type with the following members defined:
     * - resize(int, int) (allocates memory)
     * - data() (returns a pointer to the underlying array)
     * - typedef arrT::Scalar (is the data type, does not have to be dataT)
     * 
     * \tparam arrT is the type of array, see requirements above.
     * 
     * \param fname is the file path, which is passed to \ref setFilename
     * \param data is the array, which will be resized as necessary using its resize(int, int) member
     * \param head is a fitsHeader object which is passed to \ref readHeader
     * 
     * \throws mxException on error
     * 
     */ 
   template<typename arrT>
   void read(const std::string &fname, arrT & data, fitsHeader &head);
   
   ///@}
   
   
   void read(const std::vector<std::string> & flist, dataT * im);
   
   void read(const std::vector<std::string> & flist, dataT * im, std::vector<fitsHeader> &heads);
   
   
   ///Read the header from the fits file.
   /** If head is not empty, then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     * 
     * \param head is a fitsHeader object
     * 
     * \throws mxException is thrown on error
     */
   void readHeader(fitsHeader &head);
   
   ///Read the header from the fits file.
   /** If head is not empty, then only the keywords already in head are updated.  Otherwise
     * the complete header is read.
     * 
     * \param fname is the file path, which is passed to \ref setFilename
     * \param head is a fitsHeader object
     * 
     * \throws mxException is thrown on error
     */
   void readHeader(const std::string &fname, fitsHeader &head);
   
   
   /** \name Writing Data Arrays
     * @{
     */
   
   ///Write the contents of a raw array to the FITS file.
   /** 
     * Note: the type of the array must match dataT
     *  
     * \param im is the array
     * \param d1 is the first dimension
     * \param d2 is the second dimension
     * \param d3 is the third dimenesion (minimum 1)
     * \param head is an (optiona) pointer to the header.  Set to 0 if not used.
     * 
     * \throws mxException on error
     * 
     */ 
   void write(dataT * im, int d1, int d2, int d3, fitsHeader * head);
   
   ///Write the contents of a raw array to the FITS file.
   /** 
     * Note: the type of the array must match dataT
     *  
     * \param im is the array
     * \param d1 is the first dimension
     * \param d2 is the second dimension
     * \param d3 is the third dimenesion (minimum 1)
     * 
     * \throws mxException on error
     * 
     */
   void write(dataT * im, int d1, int d2, int d3);
   
   ///Write the contents of a raw array to the FITS file.
   /** 
     * Note: the type of the array must match dataT
     *  
     * \param im is the array
     * \param d1 is the first dimension
     * \param d2 is the second dimension
     * \param d3 is the third dimenesion (minimum 1)
     * \param head is the header. 
     * 
     * \throws mxException on error
     * 
     */
   void write(dataT * im, int d1, int d2, int d3, fitsHeader & head);
   
   ///Write the contents of a raw array to the FITS file.
   /** 
     * Note: the type of the array must match dataT
     *  
     * \param fname is the name of the file
     * \param im is the array
     * \param d1 is the first dimension
     * \param d2 is the second dimension
     * \param d3 is the third dimenesion (minimum 1)
     * 
     * \throws mxException on error
     * 
     */
   void write(std::string fname, dataT * im, int d1, int d2, int d3);
   
   ///Write the contents of a raw array to the FITS file.
   /** 
     * Note: the type of the array must match dataT
     *  
     * \param fname is the name of the file
     * \param im is the array
     * \param d1 is the first dimension
     * \param d2 is the second dimension
     * \param d3 is the third dimenesion (minimum 1)
     * \param head is the header. 
     * 
     * \throws mxException on error
     * 
     */
   void write(std::string fname, dataT * im, int d1, int d2, int d3, fitsHeader & head);

   ///Write the contents of and Eigen-type array to a FITS file.
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
     * \param fname is the file path, which is passed to \ref setFilename
     * \param im is the array
     * 
     * \throws mxException on error
     * 
     */ 
   template<typename arrT>
   void write(std::string fname, arrT & im); 
  
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
     * \param fname is the file path, which is passed to \ref setFilename
     * \param im is the array
     * \param head is a fitsHeader object which is passed to \ref readHeader
     * 
     * \throws mxException on error
     * 
     */ 
   template<typename arrT>
   void write(std::string fname, arrT & im, fitsHeader & head); 
   
   void writeHeader(fitsHeader &head);
   
   ///@}
   
   ///Set to read all the pixels in the file
   void setReadSize();
//    {
//       _x0 = -1;
//       _y0 = -1;
//       _xpix = -1;
//       _ypix = -1;
//    }

   ///Set to read only a subet of the pixels in the file
   /**
     * \param x0  is the starting x-pixel to read
     * \param y0  is the starting y-pixel to read
     * \param xpix is the number of x-pixels to read
     * \param ypix is the number of y-pixels to read
     */ 
   void setReadSize(long x0, long y0, long xpix, long ypix);
//    {
//       _x0 = x0;
//       _y0 = y0;
//       _xpix = xpix;
//       _ypix = ypix;
//       
//    }
   
   
}; // fitsFile

///@}

template<typename dataT>
void fitsFile<dataT>::construct()
{
   naxes = 0;
   isOpen = 0;
   fstatus = 0;
   
   nulval = 0;
   anynul = 0;
   
   noComment = 0;
   
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
   setFilename(fname, doopen);
}

template<typename dataT>
fitsFile<dataT>::~fitsFile()
{
   if(isOpen) close();
   
   if(naxes) delete naxes;
}

template<typename dataT>
void fitsFile<dataT>::setFilename(const std::string &fname, bool doopen)
{
   if(isOpen) 
   {
      close();
   }
   
   fileName = fname;
   
   if(doopen)
   {
      open();
   }

} 

template<typename dataT>
void fitsFile<dataT>::open()
{
   fits_open_file(&fptr, fileName.c_str(), READONLY, &fstatus);
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error opening file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);

      throw e;

   }
   fits_get_img_dim(fptr, &naxis, &fstatus);
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error getting number of axes in file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
      throw e;
      

   }

   if(naxes) delete naxes;
   naxes = new long[naxis];

   fits_get_img_size(fptr, naxis, naxes, &fstatus);
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error getting dimensions in file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
      throw e;
   }

   isOpen = true; //Only set this after opening is complete.
}

template<typename dataT>
void fitsFile<dataT>::open(const std::string & fname)
{
   setFilename(fname);
   open();
}

template<typename dataT>
void fitsFile<dataT>::close()
{
   int stat;
   fstatus = 0;
      
   if(!isOpen) 
   {
      std::cerr << "Tried to close but not open\n";
      return;
   }
   stat = fits_close_file(fptr, &fstatus);
      
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error closing file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
      return;
   }
   
   isOpen = 0;
   fstatus = 0;
      
   if(naxes) delete naxes;

   naxes = 0;
}

template<typename dataT>
int fitsFile<dataT>::getDimensions()
{
   if(!isOpen or !naxes) return -1;
   
   return naxis;
}

template<typename dataT>
long fitsFile<dataT>::getSize()
{
   if(!isOpen or !naxes) return -1;

   long sz = 1;
   
   if(_x0 > -1 && _y0 > -1 && _xpix > -1 && _ypix > -1 && naxis ==2)
   {
      return _xpix*_ypix;
   }
   else
   {
      for(int i=0;i<naxis;++i) sz*= naxes[i];
   }
   
   return sz;
}

template<typename dataT>
long fitsFile<dataT>::getSize(size_t axis)
{
   if(!isOpen or !naxes) return -1;

   if(_x0 > -1 && _y0 > -1 && _xpix > -1 && _ypix > -1 && naxis ==2)
   {
      if(axis == 0) return _xpix;
      return _ypix;
   }
   else
   {
      return naxes[axis];
   }
}

/************************************************************/
/***                      Basic Arrays                    ***/
/************************************************************/

template<typename dataT>
void fitsFile<dataT>::read(dataT * data)
{
   if(!isOpen)
   {
      open();
   
      if(!isOpen)
      {
         std::cerr << "File not open.\n";
         return;
      }
   }
   
   long long nelements = 1;
   long *fpix = new long[naxis];
   long *lpix = new long[naxis];
   long *inc = new long[naxis];
   
   if(_x0 > -1 && _y0 > -1 && _xpix > -1 && _ypix > -1 && naxis == 2)
   {
      fpix[0] = _x0+1;
      lpix[0] = fpix[0] + _xpix-1;
      fpix[1] = _y0+1;
      lpix[1] = fpix[1] + _ypix-1;
   
      inc[0] = 1;
      inc[1] = 1;
   }
   else
   {
      for(int i=0;i<naxis; i++)
      {
         fpix[i] = 1;
         lpix[i] = naxes[i];
         inc[i] = 1;
         //nelements *= naxes[i];
      }
   }
   
   //fits_read_pix(fptr, getFitsType<dataT>(), fpix, nelements, (void *) &nulval, 
                                     //(void *) data, &anynul, &fstatus);
   fits_read_subset(fptr, getFitsType<dataT>(), fpix, lpix, inc, (void *) &nulval, 
                                     (void *) data, &anynul, &fstatus);
   
   if (fstatus && fstatus != 107)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error reading data from file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
      delete fpix;
      delete lpix;
      delete inc;
   
      throw e;
   }

   delete fpix;
   delete lpix;
   delete inc;
}

template<typename dataT>
void fitsFile<dataT>::read(dataT * data, fitsHeader &head)
{
   read(data);
   readHeader(head);
}

template<typename dataT>
void fitsFile<dataT>::read(const std::string &fname, dataT * data)
{
   setFilename(fname);
   read(data);
}

template<typename dataT>
void fitsFile<dataT>::read(const std::string &fname, dataT * data, fitsHeader &head)
{
   setFilename(fname);
   read(data);
   readHeader(head);
}

/************************************************************/
/***                      Eigen Arrays                    ***/
/************************************************************/


template<typename arrT, bool isCube=is_eigenCube<arrT>::value>
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
      arr.resize(xsz, ysz);
   }
};

template<typename dataT>
template<typename arrT>
void fitsFile<dataT>::read(arrT & im)
{
   if(!isOpen)
   {
      open();
   
      if(!isOpen)
      {
         std::cerr << "File not open.\n";
         return;
      }
   }
      
   //long long nelements = 1;
   long *fpix = new long[naxis];
   long *lpix = new long[naxis];
   long *inc = new long[naxis];
   eigenArrResize<arrT> arrresz;
 
      
   if(_x0 > -1 && _y0 > -1 && _xpix > -1 && _ypix > -1 && naxis == 2)
   {      
      fpix[0] = _x0+1;
      lpix[0] = fpix[0] + _xpix-1;
      
      fpix[1] = _y0+1;
      lpix[1] = fpix[1] + _ypix-1;
      //nelements = _xpix*_ypix;
     
      inc[0] = 1;
      inc[1] = 1;
      arrresz.resize(im, _xpix, _ypix,1);
   }
   else
   {
      for(int i=0;i<naxis; i++)
      {
         fpix[i] = 1;
         lpix[i] = naxes[i];
         inc[i] = 1;
         //nelements *= naxes[i];
      }

   
      if(naxis > 2)
      {
         arrresz.resize(im, naxes[0], naxes[1], naxes[2]);
      }
      else
      {
         arrresz.resize(im, naxes[0], naxes[1],1);
      }
   }

   fits_read_subset(fptr, getFitsType<typename arrT::Scalar>(), fpix, lpix, inc, (void *) &nulval, 
                                     (void *) im.data(), &anynul, &fstatus);
 
   if (fstatus && fstatus != 107)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error reading data from file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
      delete fpix;
      
      throw e;
   }
   delete fpix;
      
}

template<typename dataT>
template<typename arrT>
void fitsFile<dataT>::read(arrT & data, fitsHeader &head)
{
   read(data);
   readHeader(head);
}
   
template<typename dataT>
template<typename arrT>
void fitsFile<dataT>::read(const std::string &fname, arrT & data)
{
   setFilename(fname);
   read(data);
}
   
template<typename dataT>
template<typename arrT>
void fitsFile<dataT>::read(const std::string &fname, arrT & data, fitsHeader &head)
{
   setFilename(fname);
   read(data);
   readHeader(head);
}
   
template<typename dataT>
void fitsFile<dataT>::read(const std::vector<std::string> & flist, dataT * im)
{
   if(flist.size() == 0)
   {
      std::cerr << "Empty file list\n";
      return;
   }
   
   long sz0 =0, sz1=0;
   
   for(int i=0;i<flist.size(); ++i)
   {
      setFilename(flist[i], 1);
      
      read(im + i*sz0*sz1);
      
      sz0 = getSize(0);
      sz1 = getSize(1);      
   }
}
   
template<typename dataT>
void fitsFile<dataT>::read(const std::vector<std::string> & flist, dataT * im, std::vector<fitsHeader> &heads)
{
   if(flist.size() == 0)
   {
      std::cerr << "Empty file list\n";
      return;
   }
   
   long sz0 =0, sz1=0;
   
   for(int i=0;i<flist.size(); ++i)
   {
      setFilename(flist[i], 1);
      
      read(im + i*sz0*sz1);
    
      readHeader(heads[i]);
      //pout(heads[i]["ROTOFF"].value);
      
      sz0 = getSize(0);
      sz1 = getSize(1);
      
   }
}


template<typename dataT>
void fitsFile<dataT>::readHeader(fitsHeader & head)
{
   char keyword[FLEN_KEYWORD];
   char value[FLEN_VALUE];
   char * comment;

   //The keys to look for if head is already populated
   std::list<fitsHeader::headerIterator> head_keys;
   std::list<fitsHeader::headerIterator>::iterator head_keys_it;
   int num_head_keys;
   
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
      num_head_keys = head.size();
   }
   
   //If noComment is set, then we don't read in the comment
   if(noComment)
   {
      comment = 0;
   }
   else
   {
      comment = new char[FLEN_COMMENT];
   }
  
   int keysexist;
   int morekeys;
   
   if(!isOpen)
   {
      open();
   }
   
   //This gets the number of header keys to read
   fits_get_hdrspace(fptr, &keysexist, &morekeys, &fstatus);
   
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error reading header from file: ";
      explan +=  fileName + ". ";      
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
      delete comment;
      
      throw e;
   }
   
   for(int i=0; i<keysexist; i++)
   {
      fits_read_keyn(fptr, i+1, keyword, value, comment, &fstatus);
      
      if (fstatus)
      {
         char emnem[31];
         std::string explan;
         fits_get_errstatus(fstatus, emnem);
         explan = "Error reading header from file: ";
         explan +=  fileName + ". ";      
         mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      
         delete comment;
         
         throw e;
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
            //Otherwise we append it as an unknown type, using void *
            head.append(keyword, (void *) value, comment);
         }
      }
      else
      {
         head_keys_it = head_keys.begin();
         while(head_keys_it != head_keys.end())
         {
            if( (*(*head_keys_it)).keyword == keyword)
            {
               head[keyword].value = value;
               if(comment) head[keyword].comment = comment;
               
               head_keys.erase(head_keys_it);
               
               break;
            }
            ++head_keys_it;
         }
         
         //Quit if we're done.
         if(head_keys.size() == 0) break;
      }
   }
   
   
   
   
   delete comment;
}

template<typename dataT>
void fitsFile<dataT>::readHeader(const std::string &fname, fitsHeader &head)
{
   setFilename(fname);
   readHeader(head);
}

template<typename dataT>
void fitsFile<dataT>::write(dataT * im, int d1, int d2, int d3, fitsHeader * head)
{
   if(isOpen) close();
   if(naxes) delete naxes;

   fstatus = 0;
   naxis = 1;
   if(d2 > 0) 
   {
      if(d3 > 1) naxis = 3;
      else naxis = 2;
   }
   
   naxes = new long[naxis];
   
   naxes[0] = d1;
   if(naxis > 1) naxes[1] = d2;
   if(naxis > 2) naxes[2] = d3;
   
   std::string forceFileName = "!"+fileName;
   
   fits_create_file(&fptr, forceFileName.c_str(), &fstatus);
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error creating file: ";
      explan +=  fileName + ". ";
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      throw e;

   }
   isOpen = true;
   
   fits_create_img( fptr, getFitsBITPIX<dataT>(), naxis, naxes, &fstatus);
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error creating image: ";
      explan +=  fileName + ". ";
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      throw e;

   }
   
   long fpixel[3];
   fpixel[0] = 1;
   fpixel[1] = 1;
   fpixel[2] = 1;
   
   LONGLONG nelements = 1;
   
   for(int i=0;i<naxis;++i) nelements *= naxes[i];
      
   fits_write_pix( fptr,  getFitsType<dataT>(), fpixel, nelements, im, &fstatus);
   if (fstatus)
   {
      char emnem[31];
      std::string explan;
      fits_get_errstatus(fstatus, emnem);
      explan = "Error writing data: ";
      explan +=  fileName + ". ";
      mxException e("cfitsio", fstatus, emnem, __FILE__, __LINE__, explan);
      throw e;

   }
   
   if(head != 0)
   {
      fitsHeader::headerIterator it;
      
      for(it = head->begin(); it != head->end(); ++it)
      {
         it->write(fptr);
         
         // fits_update_key(fptr, TSTRING, it->keyword.c_str(), (void *)it->value.c_str(),(char *)it->comment.c_str(),  &fstatus);
      }
   }
   
   close();
}
      
template<typename dataT>
void fitsFile<dataT>::write(dataT * im, int d1, int d2, int d3)
{
   write(im, d1, d2, d3, (fitsHeader *) 0);
}

template<typename dataT>
void fitsFile<dataT>::write(dataT * im, int d1, int d2, int d3, fitsHeader & head)
{
   write(im, d1, d2, d3, &head);
}


template<typename dataT>
void fitsFile<dataT>::write(std::string fname, dataT * im, int d1, int d2, int d3)
{
   setFilename(fname, false);
   write(im, d1, d2, d3, (fitsHeader *) 0);
}

template<typename dataT>
void fitsFile<dataT>::write(std::string fname, dataT * im, int d1, int d2, int d3, fitsHeader & head)
{
   setFilename(fname, false);
   write(im, d1, d2, d3, &head);
}


/************************************************************/
/***                      Eigen Arrays                    ***/
/************************************************************/

//SFINAE to check for 3D eigenCube
template<typename arrT, bool isCube=is_eigenCube<arrT>::value>
struct eigenArrPlanes
{
   //If it's an eigenCube, call planes planes()
   int operator()(arrT & arr)
   {
      return arr.planes();
   }
};

template<typename arrT>
struct eigenArrPlanes<arrT, false>
{
   //If it's not an eigenCube, never call planes()
   int operator()(arrT & arr)
   {
      return 1;
   }
};


template<typename dataT>
template<typename arrT>
void fitsFile<dataT>::write(std::string fname, arrT & im)
{
   eigenArrPlanes<arrT> planes;
   
   write(fname, im.data(), im.rows(), im.cols(), planes(im));
   
}

template<typename dataT>
template<typename arrT>
void fitsFile<dataT>::write(std::string fname, arrT & im, fitsHeader & head)
{
   eigenArrPlanes<arrT> planes;
   
   write(fname, im.data(), im.rows(), im.cols(), planes(im), head);
   
}

template<typename dataT>
void fitsFile<dataT>::setReadSize()
{
   _x0 = -1;
   _y0 = -1;
   _xpix = -1;
   _ypix = -1;
}

template<typename dataT>
void fitsFile<dataT>::setReadSize(long x0, long y0, long xpix, long ypix)
{
   _x0 = x0;
   _y0 = y0;
   _xpix = xpix;
   _ypix = ypix;
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

//@}

} //namespace mx

#endif //__fitsFile__



