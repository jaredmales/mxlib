/** \file readBinVector.hpp
  * \author Jared R. Males
  * \brief A utility to read in a vector of data from a binary file.
  * \ingroup ioutils
  */

#ifndef __readBinVector_hpp__
#define __readBinVector_hpp__



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
/// @}

} //namespace mx

#endif //__readColumns_hpp__