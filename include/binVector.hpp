/** \file binVector.hpp
  * \author Jared R. Males
  * \brief A utility to read/write vectors of data from/to a binary file.
  * \ingroup ioutils
  */

#ifndef __binVector_hpp__
#define __binVector_hpp__



namespace mx
{
   

/** \addtogroup ioutils
  * @{
  */

///dataT must match what was stored in the file
template<typename dataT>
void readBinVector( std::vector<dataT> & vec, const std::string & fname)
{
   FILE *fin;
   uint64_t sz;
   int nrd;
   
   fin = fopen(fname.c_str(), "r");

   nrd = fread(&sz, sizeof(uint64_t), 1, fin);

   vec.resize(sz);
   
   nrd = fread(vec.data(), sizeof(dataT), sz, fin);
   
   fclose(fin);
}


template<typename dataT>
void writeBinVector( const std::string & fname,  std::vector<dataT> & vec)
{
     
   FILE *fout;
   fout = fopen(fname.c_str(), "wb");
   uint64_t sz = vec.size();
   fwrite( &sz, sizeof(uint64_t), 1, fout);
   fwrite( vec.data(), sizeof(dataT), sz, fout);
   fclose(fout);
}

      
/// @}

} //namespace mx

#endif //__binVector_hpp__