/** \file binVector.hpp
  * \author Jared R. Males
  * \brief A utility to read/write vectors of data from/to a binary file.
  * \ingroup ioutils
  */

#ifndef __binVector_hpp__
#define __binVector_hpp__



namespace mx
{
   

/** \defgroup binvector BinVector Files
  * \ingroup ioutils
  * \brief A simple binary file format for storing vectors of data on disk.
  * 
  * The BinVector file format is very simple: the first 8 bytes contain a unit64_t integer which specifies
  * the length L of the data vector.  The second 8 bytes contain a uint64_t integer which records sizeof(dataT), 
  * where dataT is the type. The remaining L*sizeof(dataT) bytes contain the data.
  * No provision is made for recording the data type.  The programmer must know this ahead of time, though 
  * the consistency between the recorded data size and sizeof(dataT) is checked.
  * 
  * The suggested extension for BinVector files is ".binv".
  * 
  * @{
  */

/// Read a BinVector file from disk.
/** 
  * \param [out] vec is a vector which will be resized and populated.
  * \param [in] fname is the name (full-path) of the file.
  * 
  * \note dataT must match what was stored in the file.
  * 
  * \returns 0 on success.
  * \returns -1 if an error occurs.
  */
template<typename dataT>
int readBinVector( std::vector<dataT> & vec, 
                   const std::string & fname )
{
   FILE *fin;
   uint64_t sz, ssz;
   int nrd;

   fin = fopen(fname.c_str(), "r");
   if(fin == 0)
   {
      mxPError("readBinVector", errno, "Error from fopen.");
      return -1;
   }
   
   errno = 0;
   nrd = fread(&sz, sizeof(uint64_t), 1, fin);
   
   if(nrd != 1)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("readBinVector", errno, "Error reading vector size.");
      }
      else
      {
         mxError("readBinVector", MXE_FILERERR, "Error reading vector size, did not read enough bytes.");
      }
      fclose(fin);
      return -1;
   }

   errno = 0;
   nrd = fread(&ssz, sizeof(uint64_t), 1, fin);
   
   if(nrd != 1)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("readBinVector", errno, "Error reading data size.");
      }
      else
      {
         mxError("readBinVector", MXE_FILERERR, "Error reading data size, did not read enough bytes.");
      }
      fclose(fin);
      return -1;
   }
   
   if( ssz != sizeof(dataT) )
   {
      mxError("readBinVector", MXE_SIZEERR, "Mismatch between sizeof(dataT) and data size in file.");
      return -1;
   }
   
   
   vec.resize(sz);
   
   errno = 0;
   nrd = fread(vec.data(), sizeof(dataT), sz, fin);
   
   if(nrd != sz)
   {
      //Have to handle case where EOF reached but no error.
      if(errno != 0)
      {
         mxPError("readBinVector", errno, "Error reading data.");
      }
      else
      {
         mxError("readBinVector", MXE_FILERERR, "Did not read enough data.");
      }
      fclose(fin);
      return -1;
   }
   
   fclose(fin);
   
   return 0;
}

/// Write a BinVector file to disk.
/** 
  * \param [in] fname is the name (full-path) of the file.
  * \param [in] vec is the vector which will be written to disk.
  * 
  * \returns 0 on success.
  * \returns -1 if an error occurs.
  */
template<typename dataT>
int writeBinVector( const std::string & fname,  
                    std::vector<dataT> & vec )
{
     
   FILE *fout;
   int nwr;
   uint64_t sz = vec.size();
   uint64_t ssz = sizeof(dataT);

   fout = fopen(fname.c_str(), "wb");
   if(fout == 0)
   {
      mxPError("writeBinVector", errno, "Error from fopen.");
      return -1;
   }

   nwr = fwrite( &sz, sizeof(uint64_t), 1, fout);
   if(nwr != 1)
   {
      mxPError("writeBinVector", errno, "Error writing vector size.");
      fclose(fout);
      return -1;
   }
   
   nwr = fwrite( &ssz, sizeof(uint64_t), 1, fout);
   if(nwr != 1)
   {
      mxPError("writeBinVector", errno, "Error writing data size.");
      fclose(fout);
      return -1;
   }
   
   nwr = fwrite( vec.data(), sizeof(dataT), sz, fout);
   if(nwr != sz)
   {
      mxPError("writeBinVector", errno, "Error writing data.");
      fclose(fout);
      return -1;
   }
   
   fclose(fout);
   
   return 0;
}

      
/// @}

} //namespace mx

#endif //__binVector_hpp__
