/** \file binVector.hpp
  * \author Jared R. Males
  * \brief A utility to read/write vectors of data from/to a binary file.
  * \ingroup ioutils
  */

#ifndef __binVector_hpp__
#define __binVector_hpp__



namespace mx
{
   

/** \defgroup binvector binVector Files
  * \ingroup ioutils
  * \brief A simple binary file format for storing vectors of data on disk.
  * 
  * The binVector file format is very simple: the first 8 bytes contain a unit64_t integer which specifies
  * the type of the data. The second 8 bytes contain a uint64_t integer which specifies the length L of the data vector.  
  * The remaining L*sizeof(dataT) bytes contain the data.
  * 
  * The suggested extension for BinVector files is ".binv".
  * 
  * @{
  */

///Get the integer type code corresponding to the type.
/**
  * \returns an integer which uniquely identifies the type.
  * 
  * \tparam dataT is the type 
  */ 
template<typename dataT>
uint64_t binVectorTypeCode();

template<>
uint64_t binVectorTypeCode<bool>()
{
   return 0;
}

template<>
uint64_t binVectorTypeCode<signed char>()
{
   return 1;
}

template<>
uint64_t binVectorTypeCode<unsigned char>()
{
   return 2;
}

template<>
uint64_t binVectorTypeCode<char>()
{
   return 3;
}

template<>
uint64_t binVectorTypeCode<wchar_t>()
{
   return 4;
}

template<>
uint64_t binVectorTypeCode<char16_t>()
{
   return 5;
}

template<>
uint64_t binVectorTypeCode<char32_t>()
{
   return 6;
}

template<>
uint64_t binVectorTypeCode<int>()
{
   return 7;
}

template<>
uint64_t binVectorTypeCode<unsigned int>()
{
   return 8;
}

template<>
uint64_t binVectorTypeCode<short int>()
{
   return 9;
}

template<>
uint64_t binVectorTypeCode<short unsigned int>()
{
   return 10;
}

template<>
uint64_t binVectorTypeCode<long int>()
{
   return 11;
}

template<>
uint64_t binVectorTypeCode<long unsigned int>()
{
   return 12;
}

template<>
uint64_t binVectorTypeCode<long long int>()
{
   return 13;
}

template<>
uint64_t binVectorTypeCode<long long unsigned int>()
{
   return 14;
}


template<>
uint64_t binVectorTypeCode<float>()
{
   return 15;
}

template<>
uint64_t binVectorTypeCode<double>()
{
   return 16;
}

template<>
uint64_t binVectorTypeCode<long double>()
{
   return 17;
}


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
   uint64_t typecode, sz;
   int nrd;

   fin = fopen(fname.c_str(), "r");
   if(fin == 0)
   {
      mxPError("readBinVector", errno, "Error from fopen.");
      return -1;
   }

   errno = 0;
   nrd = fread(&typecode, sizeof(uint64_t), 1, fin);
   
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
   
   if( typecode != binVectorTypeCode<dataT> )
   {
      mxError("readBinVector", MXE_SIZEERR, "Mismatch between type dataT and type in file.");
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
   uint64_t typecode = binVectorTypeCode<dataT>(); 
   uint64_t sz = vec.size();
   
   fout = fopen(fname.c_str(), "wb");
   if(fout == 0)
   {
      mxPError("writeBinVector", errno, "Error from fopen.");
      return -1;
   }

   nwr = fwrite( &typecode, sizeof(uint64_t), 1, fout);
   if(nwr != 1)
   {
      mxPError("writeBinVector", errno, "Error writing typecode");
      fclose(fout);
      return -1;
   }

   nwr = fwrite( &sz, sizeof(uint64_t), 1, fout);
   if(nwr != 1)
   {
      mxPError("writeBinVector", errno, "Error writing vector size.");
      fclose(fout);
      return -1;
   }
      
   nwr = fwrite( vec.data(), sizeof(dataT), vec.size(), fout);
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
